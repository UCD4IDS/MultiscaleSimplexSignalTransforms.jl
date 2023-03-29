export 
    Representation, SubRepresentation, Basis, PartitionInput, PartitionMethod, EigenMethod,
    SCPartition,
    root,
    partition!, hierarchical_partition!,
    part_by_sign, part_by_kmeans

# how are the adjacencies and degrees represented as a Laplacian-like matrix?
# leaving room for, e.g., method using 0-lap on graph formed by k-simplices
@enumx Representation KLaplacian Distance

# what is the technique used to construct representation matrices
# on the partition sub-complexes?
# the SubmatrixPartition just computes L once and then uses submatrices to get 
#   approximate, but more stable Fiedler vectors
# the FullPartition recomputes L explicitly for each partition level
@enumx SubRepresentation Submatrix Full

# how to compute the eigenvector basis?
# leaving room for, e.g., searching for bottom non-zero values
@enumx Basis Entire Fiedler

# how do we augment the input basis before passing it through
# the partition method?
# Identity: pass the basis through as is
# DCOrientation: apply an orientation such that the bottom
#   eigenvector becomes a DC component
@enumx PartitionInput Identity DCOrientation

# given a representation matrix, how is the partition performed?
@enumx PartitionMethod FiedlerSign Kmeans

# which eigenvalue is sought by the Fiedler basis method?
# first: simply the first one
# positive: an exponential search for the first nonzero eig
@enumx EigenMethod First Positive

abstract type SCPartition end
# interface:
root(::SCPartition) = error("needs to be implemented")
n(::SCPartition) = error("needs to be implemented")
_partition!(::SCPartition, ::PartitionTree) = error("needs to be implemented")

SCPartition(subr::SubRepresentation.T, region::R; kwargs...) where R<:Region = @match subr begin
    SubRepresentation.Submatrix => SubmatrixPartition(region; kwargs...)
    SubRepresentation.Full      => FullPartition(region; kwargs...)
end

Base.@kwdef struct PartitionOutput{M<:AbstractMatrix}
    basis::(M|Nothing)=nothing
    parts::NTuple{2, Vector{Int}}
    vals::(Vector{Float64}|Nothing)=nothing
end

PartitionOutput(parts::NTuple{2, Vector{Int}}) = PartitionOutput{Matrix{Float64}}(; parts)

function part_components(ccs)
    lpart, rpart = Int[], Int[]
    for cc in sort!(ccs, by=cc->length(cc), rev=true)
        append!(length(lpart) < length(rpart) ? lpart : rpart, cc)
    end
    lpart, rpart
end

function partition!(
    part::P,
    node::PartitionTree = root(part)
) where P<:SCPartition
    if isleaf(node)
		return ([0.0], [1.0])
    end

    if partlen(node) == 2
        lchild!(node, PartitionTree(; sinds=[node.sinds[1]]))
        rchild!(node, PartitionTree(; sinds=[node.sinds[2]]))
        return ([0.0,0.0], [1.0 0.0; 0.0 1.0])
    end

    # SCPartition implementation-specific partitioning should
    # return a basis -- this is forwarded to hierarchical_partition!
    (; basis, parts, vals) =  _partition!(part, node)
    (lpart, rpart) = parts

    any(isempty.((lpart, rpart))) && error("partition was trivial")

    lchild!(node, PartitionTree(; sinds=node.sinds[lpart]))
    rchild!(node, PartitionTree(; sinds=node.sinds[rpart]))
    
    # return the basis so that it can be used by maps which
    # iterate along a hierarchical partition, like the HGLET
    (vals, basis)
end

indicator(part::P) where P<:SCPartition = indicator(root(part), partlen(root(part)))

iterate_tree(iter_fn, part_fn, part, n) = function(tsn)
    vals, basis = part_fn(part, tsn.node)
    !isnothing(iter_fn) && iter_fn(tsn.node; n, basis)
    (tsn.node.sinds, vals)
end

recurse_tree(recurse_fn) = function recurse_tree_inner(node)
    lout, rout = isleaf(node) ?
        (nothing, nothing) :
        recurse_tree_inner.(children(node))
    recurse_fn(node, lout, rout)
end

resolve_missing_values(eigvals) = function(node, _, _)
    if !isleaf(node) && isnothing(node.fvals)
        nind = findfirst(v -> v[1] == node.sinds, eigvals)
        lind = findfirst(v -> v[1] == node.lchild.sinds, eigvals)
        rind = findfirst(v -> v[1] == node.rchild.sinds, eigvals)
        stackvals = [eigvals[lind][2]; eigvals[rind][2]]
        p = sortperm(stackvals)
        node.fvals = TagDict(
            i-1 => c
            for (i, c) = enumerate([
                td[tdi]
                for td ∈ (node.lchild.fvals, node.rchild.fvals)
                for tdi ∈ sort(collect(keys(td)))
            ][p])
        )
        eigvals[nind] = (eigvals[nind][1], stackvals[p])
    end
    nothing
end

function hierarchical_partition!(
    part::P;
    maxlevel::Real=Inf,
    iterate_with_partition::(Nothing|Function)=nothing,
    recurse_after_partition::(Nothing|Function)=nothing
) where P<:SCPartition
    iter_f = !isnothing(iterate_with_partition)
    recurse_f = !isnothing(recurse_after_partition)
    no_f = !iter_f && !recurse_f

    partfun = iterate_tree(iterate_with_partition, partition!, part, n(part))
        # iterate_with_partition(tsn.node; n, basis=partition!(part, tsn.node))
    itercond = _->true #tsn -> partlen(tsn.node) > 1
    eigvals = dfs(partfun, root(part); itercond, maxlevel)
    
    recurse_f && recurse_tree(recurse_after_partition)(root(part))
    !no_f && recurse_tree(resolve_missing_values(eigvals))(root(part))
end

function resolve_zeros_randomly!(signs::Vector; tol=1e-14, check_nonuniform=false)
    # break ties randomly
    zs = abs.(signs) .< tol
    nzs = sum(zs)
    nzs == 0 && return

    signs[zs] = rand((-1,1), nzs)

    # don't allow uniform output
    if check_nonuniform && ((signs .> 0) |> x -> (all(x) || all(.!x)))
        signs[rand(findall(zs))] *= -1
    end
end

function dc_orientation(basis::M, resolve_zeros=resolve_zeros_randomly!) where M<:AbstractMatrix
    signs = sign.(basis[:,1])
    resolve_zeros(signs)
    basis .* signs
end

function part_by_sign(signs::Vector; resolve_zeros=resolve_zeros_randomly!, fix_first=true)
    resolve_zeros(signs; check_nonuniform=true)
    if fix_first
        signs .*= sign(signs[1])
    end
    findall.([>(0), <(0)], Ref(signs)) |> Tuple
end

function part_by_kmeans(X::Array)
    labels = kmeans(X', 2).assignments
    findall.(
        [==(labels[1]), !=(labels[1])],
        Ref(labels)
    ) |> Tuple
end

reversed_eigen(M) = eigen(M) |> ((Λ, X), ) -> (reverse(Λ), X[:,end:-1:1])

function fiedler_eigs(representation::Representation.T, basis::Basis.T, eigenmethod::EigenMethod.T; tol = 1e-4, maxiter=30000)
    whicheig, eignum = @match representation begin
        Representation.KLaplacian => (:SM, 2)
        Representation.Distance => (:LR, 1)
    end

    entirebasis = basis == Basis.Entire
    keeptrying = eigenmethod == EigenMethod.Positive
    keeptrying && entirebasis && error("EigenMethod.Positive and Basis.Entire are incompatible")
    base_eignum = eignum

    function poseigs(M, eignum=eignum; keeptrying=keeptrying, entirebasis=entirebasis)
        takeall = entirebasis || eignum > size(M,1)/4
        vals, vecs = takeall ?
            eigen(Matrix(M)) :
            eigs(M+I; nev=eignum, which=whicheig, maxiter)
        !takeall && (vals .-= 1)

        if keeptrying && !takeall && maximum(vals[1:end-1]) < tol
            return poseigs(M, min(2*eignum, size(M, 1)); keeptrying)
        end

        eignum = base_eignum
        lim = keeptrying ? findfirst(≥(tol), vals)+1 : eignum
        ret_range = entirebasis ? ((lim-eignum+1):lastindex(vals)) : ((lim-eignum+1):lim)

        (vals[ret_range], vecs[:,ret_range], eignum) |> vv -> (
            whicheig==:SM ? vv : (reverse(vv[1]), vv[2][:,end:-1:1], vv[3])
        )
    end
end
