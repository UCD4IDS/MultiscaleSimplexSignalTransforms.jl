export
    KRegion, ZeroRegion,
    n_simp, n_bd, n_hull,
    boundaries, hulls,
    getweights,
    absolute_boundary_matrix

abstract type WeightMap{K} end

struct Weights{K} <: WeightMap{K}
    weights::OrderedDict{NodeTuple{K},Float64}
end

@forward Weights.weights Base.get, Base.getindex, Base.setindex!, Base.length, Base.haskey
setfn!(wts::Weights{K}, weightfn, simplex::NodeTuple{K}) where {K} = (wts[simplex] = weightfn(simplex))
setfn!(wts::Weights, weightfn, simplex::Nodes) = setfn!(wts, weightfn, Tuple(simplex))
simplices(wts::Weights) = collect(keys(wts.weights))
getweights(wts::Weights) = collect(values(wts.weights))

struct ConstantWeight{K} <: WeightMap{K}
    members::OrderedSet{NodeTuple{K}}
    wt::Float64
end

@forward ConstantWeight.members Base.length
Base.haskey(wts::ConstantWeight{K}, nt::NodeTuple{K}) where {K} = nt ∈ wts.members
Base.get(wts::ConstantWeight{K}, nt::NodeTuple{K}, default::Float64) where {K} = nt ∈ wts.members ? wts.wt : default
Base.getindex(wts::ConstantWeight{K}, nt::NodeTuple{K}) where {K} = get(wts.members, nt, wts.wt)
setfn!(wts::ConstantWeight{K}, ::Nothing, simplex::NodeTuple{K}) where {K} = (push!(wts.members, simplex); wts.wt)
setfn!(wts::ConstantWeight, ::Nothing, simplex::Nodes) = setfn!(wts, nothing, Tuple(simplex))
simplices(wts::ConstantWeight) = collect(wts.members)
getweights(wts::ConstantWeight) = ones(Float64, length(wts.members))

Base.getindex(t::NTuple{N,Int64}, ind::InvertedIndices.InvertedIndex) where {N} = collect(t)[ind]

abstract type Region end

Base.show(io::IO, ::MIME"text/plain", region::R) where {R<:Region} = print(
    io,
    repr("text/plain", typeof(region); context=:compact => true)
)

struct ZeroRegion <: Region
    weights::SimpleWeightedGraph
    diagmodifier::Float64 # the weak adj weight applies to all nodes, so it just modifies the diagonal
end

k(::ZeroRegion) = 0
n_simp(region::ZeroRegion) = nv(region.weights)
n_bd(::ZeroRegion) = 0
n_hull(region::ZeroRegion) = ne(region.weights)

complex(region::ZeroRegion) = cliquecomplex(region.weights, 1)

simplices(region::ZeroRegion) = [(i,) for i in Graphs.vertices(region.weights)]
boundaries(::ZeroRegion) = []
hulls(region::ZeroRegion) = map(e -> Tuple(e)[1:2], edges(region.weights))

function ZeroRegion(
    g::AbstractGraph; hullweightfn=nothing, weakadjweight::Float64=0.0, subregion_inds=nothing
)
    if !isnothing(subregion_inds)
        g, _ = induced_subgraph(g, subregion_inds)
    end

    es = Tuple.(edges(g))

    # if !isnothing(subregion_inds)
    #     filter!(e -> all(e .∈ subregion_inds), es)
    #     g = SimpleGraph(Edge.(es))
    # end

    ZeroRegion(
        isnothing(hullweightfn) ?
        SimpleWeightedGraph(g) :
        SimpleWeightedGraph(getindex.(es, 1), getindex.(es, 2), hullweightfn.(es)),
        weakadjweight
    )
end

ZeroRegion(tree::SimplexTree; kwargs...) = ZeroRegion(
    SimpleGraph(Edge.(Tuple.(simplices(tree, 1))));
    kwargs...
)

ZeroRegion(M::AbstractMatrix; kwargs...) = ZeroRegion(SimpleWeightedGraph(M); kwargs...)

struct AdjacentKSimplices{KU}
    hull::NodeTuple{KU}
    leaves::NTuple{2,Node}
    natconsistent::Bool

    function AdjacentKSimplices{KU}(hinge::NodeTuple{KL}, leaves::NTuple{2,Node}) where {KL,KU}
        ns = [leaves..., hinge...]
        p = sortperm(ns)
        ip = invperm(p)
        new{KU}(NodeTuple{KU}(ns[p]), leaves, (ip[1] + ip[2]) % 2 == 0)
    end
end

HingeAdjacency{K} = @NamedTuple{leaves::Dict{Node,Int}, pairs::Vector{AdjacentKSimplices{K}}}

struct KRegion{KL,K,KU} <: Region
    tree::SimplexTree
    hingemap::DefaultDict{NodeTuple{KL},HingeAdjacency{KU}}
    bdweights::Weights{KL}
    simplexweights::Weights{K}
    hullweights::WeightMap{KU}
end

k(::KRegion{KL}) where {KL} = KL
n_simp(region::KRegion) = length(region.simplexweights)
n_bd(region::KRegion) = length(region.bdweights)
n_hull(region::KRegion) = length(region.hullweights)

complex(region::KRegion) = region.tree

simplices(region::KRegion) = simplices(region.simplexweights)
boundaries(region::KRegion) = simplices(region.bdweights)
hulls(region::KRegion) = simplices(region.hullweights)

function KRegion(
    tree::SimplexTree, k::Int;
    hullweightfn=nothing,
    weakadjweight::Float64=0.0,
    subregion_inds::(Nothing | Vector{Int})=nothing
)
    k == 0 && return ZeroRegion(tree; hullweightfn, weakadjweight, subregion_inds)

    KL, K, KU = k:k+2 # number of nodes in the k-1 boundary simplex, k-simplex, k+1 hull simplex
    simps = collect(SimplexIterator(tree; maxlevel=k + 2)) # avoid iterating through higher-than-needed simplices

    # set a consistent indexing for the simplices
    bds, ksimps, hulls = filter.(((s -> length(s) == ki) for ki ∈ k:k+2), Ref(simps))
    if !isnothing(subregion_inds)
        ksimps = ksimps[subregion_inds]
        hulls = filter(
            hl -> all(face ∈ ksimps for face in faces(hl)),
            hulls
        )
    end

    hingemap = DefaultDict{NodeTuple{KL},HingeAdjacency{KU}}(
        () -> HingeAdjacency{KU}((Dict(), AdjacentKSimplices{KU}[]))
    )
    bdweights = Weights{KL}(OrderedDict(Tuple(bd) => 0.0 for bd in bds))
    simplexweights = Weights{K}(OrderedDict(Tuple(ksimp) => weakadjweight for ksimp in ksimps))
    hullweights = isnothing(hullweightfn) ?
                  ConstantWeight{KU}(OrderedSet(Tuple.(hulls)), 1.0) :
                  Weights{KU}(OrderedDict(Tuple(hull) => 0.0 for hull in hulls))

    # first build up adjacency structure
    for (simplex_ind, simplex) ∈ enumerate(ksimps)
        for (hinge, leaf) ∈ faces(simplex)
            for adjleaf ∈ keys(hingemap[hinge].leaves)
                as = AdjacentKSimplices{KU}(hinge, (leaf, adjleaf))
                push!(hingemap[hinge].pairs, as)
            end
            hingemap[hinge].leaves[leaf] = simplex_ind
        end
    end

    # then assign weights via hulls
    for hull ∈ hulls
        # handle hull weights
        wt = setfn!(hullweights, hullweightfn, hull)

        # accumulate to simplex degrees ...
        for (simplex, _) ∈ faces(hull)
            simplexweights[simplex] += wt

            # and to boundary degrees   
            for (hinge, _) ∈ faces(simplex)
                bdweights[hinge] += wt / 2 # each hull contributes twice to each boundary
            end
        end
    end

    # finally accumulate weakadjweight to boundary degrees
    for (hinge, adj) ∈ hingemap
        bdweights[hinge] += weakadjweight * length(adj.leaves)
    end

    KRegion(tree, hingemap, bdweights, simplexweights, hullweights)
end

KRegion(g::AbstractGraph, k::Int; kwargs...) = KRegion(cliquecomplex(g, k + 1), k; kwargs...)

function absolute_boundary_matrix(region::KRegion{KL}, k::Int) where {KL}
    @match k begin
        $(KL - 1) => sparse( # lower boundary matrix
            unzip([
                (i, j)
                for (i, b) ∈ enumerate(boundaries(region))
                for j ∈ values(region.hingemap[b].leaves)
            ])...,
            1.0, n_bd(region), n_simp(region)
        )
        $(KL) => error("upper boundary matrix not implemented yet")
        _ => error("for a k-region, boundary matrix k must be k-1 for lower or k for upper")
    end
end
