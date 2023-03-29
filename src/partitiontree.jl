export 
    partlen,
    PartVals, ConstVal, FunVals,
    TagDict,
    val,
    clear!,
    indicator,
    index_partitions,
    treecolors,
    PartitionTree, lchild!, rchild!, islchild, isrchild

# PartVals are function values over indices in some discrete partition
abstract type PartVals end
Base.copy(pv::T) where T<:PartVals = T(copy(val(pv)))

PartValType = (Number | Vector{<:Number} | SparseVector{<:Number})
PartVals(vals::PartValType) = convert(PartVals, vals)

# represents a constant value over an array of indices
struct ConstVal <: PartVals
    val::Float64
end

# any PartVals constructed from a single number is a ConstVal
Base.convert(::Type{PartVals}, x::Number) = ConstVal(x)
Base.getindex(fv::ConstVal, ::Int) = fv.val
val(fv::ConstVal) = fv.val

# represents some function values as an array that is dense over the root domain of the partition
struct FunVals <: PartVals
    vals::Vector{Float64}
end

# any PartVals constructed from a vector is a FunVals
Base.convert(::Type{PartVals}, v::Vector{<:Number}) = FunVals(v)
Base.getindex(fv::FunVals, ind::Int) = fv.vals[ind]
val(fv::FunVals) = fv.vals

# represents some function values as a sparse array, with nonzeros on one region of the partition
struct SparseVals <: PartVals
    vals::SparseVector{Float64}
end

# any PartVals constructed from a sparse vector is a SparseVals
Base.convert(::Type{PartVals}, v::SparseVector{<:Number}) = SparseVals(v)
Base.getindex(fv::SparseVals, ind::Int) = fv.vals[ind]
val(sv::SparseVals) = sv.vals

Tag = Int
TagDict = Dict{Tag, PartVals}

Base.setindex!(d::TagDict, t::Tag, x::Number) = setindex!(d, t, ConstVal(x))
Base.setindex!(d::TagDict, t::Tag, v::Vector{<:Number}) = setindex!(d, t, FunVals(v))
Base.setindex!(d::TagDict, t::Tag, v::SparseVector{<:Number}) = setindex!(d, t, SparseVals(v))

# represents a subtree of the partitioning of the k-simplices of a graph
Base.@kwdef mutable struct PartitionTree <: Tree
    sinds::Vector{Int}
    parent::(Nothing|PartitionTree) = nothing
    lchild::(Nothing|PartitionTree) = nothing
    rchild::(Nothing|PartitionTree) = nothing
    fvals::(Nothing|TagDict) = nothing
end

PartitionTree(n::Int, fvals::(Nothing|TagDict)=nothing) = PartitionTree(; sinds = 1:n, fvals)
PartitionTree(g::AbstractGraph, ksimpfun::Function = g -> 1:nv(g)) = PartitionTree(length(ksimpfun(g)))

function lchild!(tree::PartitionTree, child::PartitionTree)
    tree.lchild = child
    child.parent = tree
end

function rchild!(tree::PartitionTree, child::PartitionTree)
    tree.rchild = child
    child.parent = tree
end

partlen(tree::PartitionTree) = length(tree.sinds)
# don't define a length -- it would just specialize the version already defined for Tree
# Base.length(tree::PartitionTree) = length(tree.sinds)

isleaf(tree::PartitionTree) = partlen(tree) == 1

Base.copy(tree::PartitionTree; keep_vals = false) = keep_vals ?
    PartitionTree(partlen(tree), copy(tree.fvals)) : 
    PartitionTree(partlen(tree))

function Base.deepcopy(tree::PartitionTree; keep_vals = false)
    new_tree = copy(tree; keep_vals)
    isnothing(tree.lchild) || lchild!(new_tree, deepcopy(tree.lchild; keep_vals))
    isnothing(tree.rchild) || rchild!(new_tree, deepcopy(tree.rchild; keep_vals))
    new_tree
end

clear!(tree::PartitionTree) = dfs(tsn -> (tsn.node.fvals = nothing), tree)

children(tree::PartitionTree) = filter(!isnothing, [tree.lchild, tree.rchild])
islchild(tree::PartitionTree) = !isnothing(tree.parent) && tree.parent.lchild == tree
isrchild(tree::PartitionTree) = !isnothing(tree.parent) && tree.parent.rchild == tree

function treecolors(tree::PartitionTree; maxlevel = Inf)
	colors = ones(partlen(tree))
    function colorleft(tsn)
        if islchild(tsn.node)
            colors[tsn.node.sinds] .-= 2.0^-tsn.level
        end
    end

    bfs(colorleft, tree; maxlevel)
	colors
end

# N = 2 for TagArrays of coefficients, (n x L+1)
# N = 3 for TagArrays of basis vectors, (n x n x L+1) so each vector is a column
# notice that getindex signature is ergonomic for notation, but storage is different for efficiency
struct TagArray{R, N}
    data::Array{R, N}
    tagmap::Dict{NTuple{3, Int}, Int}
    offsets::Vector{Vector{Int}}
    sindmap::Dict{Vector{Int}, NTuple{2, Int}}
end

TagMatrix{R} = TagArray{R, 2}
TagTensor{R} = TagArray{R, 3}

# get a single coefficient
Base.getindex(M::TagMatrix, j::Int, k::Int, l::Int) = M.data[M.tagmap[(j, k, l)], j+1]
# get a single basis vector
Base.getindex(M::TagTensor, j::Int, k::Int, l::Int) = M.data[:, M.tagmap[(j, k, l)], j+1]

# get coefficients corresponding to subregion as a vector
Base.getindex(M::TagMatrix, j::Int, k::Int, ::Colon) = M.data[
    M.offsets[j+1][k+1]+1 : (length(M.offsets[j+1]) == k+1 ? lastindex(M.data, 1) : M.offsets[j+1][k+2]),
    j+1
]
# get basis vectors corresponding to subregion as a matrix
Base.getindex(M::TagTensor, j::Int, k::Int, ::Colon) = M.data[
    :,
    M.offsets[j+1][k+1]+1 : (length(M.offsets[j+1]) == k+1 ? lastindex(M.data, 2) : M.offsets[j+1][k+2]), 
    j+1
]

# get all coefficients for given level as a vector
Base.getindex(M::TagMatrix, j::Int, ::Colon) = M.data[:, j+1]
# get all basis vectors for given level as a matrix
Base.getindex(M::TagTensor, j::Int, ::Colon) = M.data[:, :, j+1]

# get either all coefficients as a matrix, or all basis vectors as a tensor
Base.getindex(M::TagArray, ::Colon, ::Colon) = M.data

function index_partitions(tree::PartitionTree)
    maxlevel = dfs(tsn -> tsn.level, tree) |> maximum
    isbasis = isa(val(first(values(tree.fvals))), AbstractVector)

    res = [SortedDict{Tag, PartValType}[] for _ ∈ 0:maxlevel]
    ss = [Vector{Int}[] for _ ∈ 0:maxlevel]
    function propagate_vals(tsn)
        push!(ss[tsn.level+1], tsn.node.sinds)
        for level ∈ tsn.level : (isempty(children(tsn.node)) ? maxlevel : tsn.level)
            push!(res[level+1], SortedDict(tag => val(vals) for (tag, vals) ∈ tsn.node.fvals))
        end
    end
    dfs(propagate_vals, tree)

    M = if isbasis
        reduce(
            (x,y) -> cat(x, y; dims=3), 
            [reduce(hcat, [reduce(hcat, values(d)) for d ∈ ds]) for ds ∈ res]
        )
    else
        reduce(
            hcat, 
            [reduce(vcat, [reduce(vcat, values(d)) for d ∈ ds]) for ds ∈ res]
        )
    end

    tagmap = Dict{NTuple{3, Int}, Int}()
    sindmap = Dict{Vector{Int}, NTuple{2, Int}}()
    offsets = Vector{Int}[]
    for j = 0:maxlevel
        partlens = length.(res[j+1])
        offset = [0; cumsum(partlens)[1:end-1]]
        push!(offsets, offset)
        for (k, d) ∈ enumerate(res[j+1])
            length(ss[j+1]) ≥ k && (sindmap[ss[j+1][k]] = (j, k-1))
            for (l, tag) ∈ enumerate(keys(d))
                tagmap[(j, k-1, tag)] = offset[k] + l
            end
        end
    end
    
    TagArray(M, tagmap, offsets, sindmap)
end

function indicator(tree::PartitionTree, n::Int)
    indicator = sparsevec(tree.sinds, ones(partlen(tree)), n)
    lc = tree.lchild
    if isnothing(lc)
        Vector(indicator)
    else
        Vector(indicator - 2 .* sparsevec(lc.sinds, ones(partlen(lc)), n))
    end
end
