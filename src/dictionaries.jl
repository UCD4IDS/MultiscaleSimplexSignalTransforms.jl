export
    MultiscaleSimplicialTransform, kGHWT, kHGLET,
    SimplicialSignal, BasisSignal, CoefSignal,
    analyze, basis!, transform!,
    best_basis, alt_best_basis

abstract type SimplicialSignal end
val(signal::SimplicialSignal) = val(signal.val)

# used to ask an MSST to construct a basis of dense vectors
struct BasisSignal{P<:PartVals} <: SimplicialSignal
    val::P
    n::Int
end

BasisSignal(v::PartValType, n::Int) = BasisSignal(PartVals(v), n)
Base.getindex(signal::BasisSignal, ind::Int) = signal.val[ind]
haar_repr(signal::BasisSignal, ind::Int) = delta(ind, signal.n, signal[ind])

# used to ask an MSST to construct coefficients 
struct CoefSignal{P<:PartVals} <: SimplicialSignal
    val::P
end

CoefSignal(v::PartValType) = CoefSignal(PartVals(v))
Base.getindex(signal::CoefSignal, ind::Int) = signal.val[ind]
haar_repr(signal::CoefSignal, ind::Int) = signal[ind]

abstract type MultiscaleSimplicialTransform end
MSST = MultiscaleSimplicialTransform

part(::MSST) = error("must implement `part`")
inds(::MSST) = error("must implement `inds`")
root(transform::MSST) = root(part(transform))
clear!(transform::MSST) = clear!(root(transform))

Base.getindex(transform::MSST, jkl...) = getindex(inds(transform), jkl...)
Base.getindex(transform::MSST, II::NTuple{3,Int}) = getindex(transform, II...)
Base.show(io::IO, ::MIME"text/plain", transform::MSST) = print(
    io,
    repr("text/plain", typeof(transform); context=:compact => true)
)

# intended to be part of the API, but not currently used, because
# - basis is only ever constructed when initializing the first MSST instance for a given complex
# - basis is automatically constructed on this initialization, for both kGHWT and kHGLET
basis!(
    T::Type{<:MSST}, part::SCPartition, signal::PartValType=1.0
) = transform!(T, part, BasisSignal(signal, n(part)))

analyze(
    transform::MSST, signal::PartValType
) = typeof(transform)(deepcopy(part(transform)), CoefSignal(signal); dict=transform)

function best_basis(transform::MSST, signal::PartValType; p=0.5, cost=v -> norm(v, p))
    coefs = analyze(transform, signal)

    function compare_bases(node, lout, rout)
        jk = transform.inds.sindmap[node.sinds]
        # new_coefs = coefs[jk..., :]
        new_coefs = val.(values(node.fvals))

        isleaf(node) && return (; locs=[jk], coefs=new_coefs)

        join_locs = [lout.locs; rout.locs]
        join_coefs = [lout.coefs; rout.coefs]
        if length(new_coefs) != length(join_coefs)
            @info ("jk:", jk)
            @info ("locs:", lout, rout)
            @info new_coefs
            @info ("coefs:", length(new_coefs), length(join_coefs))
            @info node.sinds
            @info coefs[jk..., :]
            error("bah")
        end

        cost(new_coefs) < cost(join_coefs) ?
        (; locs=[jk], coefs=new_coefs) :
        (; locs=join_locs, coefs=join_coefs)
    end

    recurse_tree(compare_bases)(root(coefs))
end

function alt_best_basis(transform::MSST, signal::PartValType; p=0.5, cost=v -> norm(v, p))
    coefs = analyze(transform, signal)

    function compare_bases(node, lout, rout)
        new_locs = Dict(node.sinds => node.fvals)

        isleaf(node) && return new_locs

        new_coefs = val.(values(node.fvals))

        join_locs = merge(lout, rout)
        join_coefs = @chain join_locs begin
            values
            @. values
            @. collect
            reduce(vcat, _)
            @. val
        end
        # @info ( new_locs, new_coefs, join_locs, join_coefs)

        cost(new_coefs) < cost(join_coefs) ? new_locs : join_locs
    end

    recurse_tree(compare_bases)(root(coefs))
end


struct kGHWT{P<:SCPartition} <: MSST
    part::P
    inds::TagArray
    function kGHWT(
        part::P, signal::SimplicialSignal=BasisSignal(1.0, n(part)); dict::(kGHWT | Nothing)=nothing
    ) where {P}
        transform!(isnothing(dict) ? kGHWT : dict, part, signal)
        new{P}(part, index_partitions(root(part)))
    end
end

kGHWT{P}(args...; kwargs...) where {P} = kGHWT(args...; kwargs...)

kGHWT(
    region::R, subr::SubRepresentation.T=SubRepresentation.Submatrix; kwargs...
) where {R<:Region} = kGHWT(SCPartition(subr, region; kwargs...))

part(transform::kGHWT) = transform.part
inds(transform::kGHWT) = transform.inds

function kghwt!(node::PartitionTree, signal::SimplicialSignal)
    node.fvals = TagDict()

    if isleaf(node)
        node.fvals[0] = haar_repr(signal, node.sinds[1])
        return
    end

    nl = partlen(node.lchild)
    nr = partlen(node.rchild)

    fl = node.lchild.fvals
    fr = node.rchild.fvals

    node.fvals[0] = (sqrt(nl) .* val(fl[0]) + sqrt(nr) .* val(fr[0])) ./ sqrt(nl + nr)
    node.fvals[1] = (sqrt(nr) .* val(fl[0]) - sqrt(nl) .* val(fr[0])) ./ sqrt(nl + nr)

    for ℓ ∈ delete!(keys(fl) ∪ keys(fr), 0)
        if haskey(fl, ℓ)
            if haskey(fr, ℓ)
                node.fvals[2ℓ] = (val(fl[ℓ]) .+ val(fr[ℓ])) ./ sqrt(2)
                node.fvals[2ℓ+1] = (val(fl[ℓ]) .- val(fr[ℓ])) ./ sqrt(2)
            else
                node.fvals[2ℓ] = val(fl[ℓ])
            end
        else
            node.fvals[2ℓ] = val(fr[ℓ])
        end
    end

    nothing
end

# forming a basis requires performing the hierarchical partition
transform!(::Type{kGHWT}, part::SCPartition, signal::BasisSignal) = hierarchical_partition!(
    part; recurse_after_partition=(node, _, _) -> kghwt!(node, signal)
)

# if just analyzing, only perform the tree recursion
transform!(::Type{kGHWT}, part::SCPartition, signal::CoefSignal) = recurse_tree(
    (node, _, _) -> kghwt!(node, signal)
)(root(part))

# if we receive an existing dictionary for kGHWT, ignore it
transform!(::kGHWT, args...) = transform!(kGHWT, args...)


struct kHGLET{P<:SCPartition} <: MSST
    part::P
    inds::TagArray
    function kHGLET(
        part::P, signal::SimplicialSignal=BasisSignal(1.0, n(part)); dict::(kHGLET | Nothing)=nothing
    ) where {P}
        transform!(isnothing(dict) ? kHGLET : dict, part, signal)
        new{P}(part, index_partitions(root(part)))
    end
end

kHGLET{P}(args...; kwargs...) where {P} = kHGLET(args...; kwargs...)

kHGLET(
    region::R, subr::SubRepresentation.T=SubRepresentation.Submatrix; kwargs...
) where {R<:Region} = kHGLET(SCPartition(subr, region; basis=Basis.Entire, kwargs...))

part(transform::kHGLET) = transform.part
inds(transform::kHGLET) = transform.inds

function khglet!(
    node::PartitionTree; basis::(M | Nothing), n::Int
) where {M<:AbstractArray}
    isnothing(basis) && return
    node.fvals = TagDict(
        i - 1 => Vector(sparsevec(node.sinds, c, n))
        for (i, c) = enumerate(eachcol(basis))
    )
end

# currently forming a basis with kHGLET ignores any input signal
transform!(::Type{kHGLET}, part::SCPartition, ::BasisSignal) = hierarchical_partition!(
    part; iterate_with_partition=khglet!
)

# analyzing requires taking inner products against each basis function
transform!(::Type{kHGLET}, part::SCPartition, signal::CoefSignal) =
    error("Producing kHGLET coefficients requires an existing dictionary")

# ignore existing dictionary when forming a basis
transform!(::kHGLET, part::SCPartition, signal::BasisSignal) = transform!(kHGLET, part, signal)

# analyzing requires taking inner products against each basis function
function transform!(dict::kHGLET, part::SCPartition, signal::CoefSignal)
    dict_xs = dfs(tsn -> tsn.node.fvals, root(dict))
    for (i, tsn) ∈ enumerate(SearchTree(root(part), DFS))
        tsn.node.fvals = TagDict()
        for (tag, x) ∈ dict_xs[i]
            tsn.node.fvals[tag] = sum(val(x) .* val(signal))
        end
    end
end
