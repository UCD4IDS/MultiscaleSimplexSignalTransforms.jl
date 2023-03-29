export 
    Node, Nodes, NodeTuple,
    SimplexTree, SimplexIterator,
    simplices,
    cliquecomplex,
    relabeled

"""
SimplexTree represents a simplicial complex via a tree data structure
"""
mutable struct SimplexTree <: Tree
    parent::(Nothing|SimplexTree)
    node::(Nothing|Int)
    children::SortedDict{Int, SimplexTree}
end

children(tree::SimplexTree) = collect(values(tree.children))

SimplexTree() = SimplexTree(nothing, nothing, Dict())
SimplexTree(parent::SimplexTree, node::Int) = SimplexTree(parent, node, Dict())

child!(tree::SimplexTree, n::Int) = (
    haskey(tree.children, n),
    get!(() -> SimplexTree(tree, n), tree.children, n)
)
child(tree::SimplexTree, n::Int) = get(tree.children, n, nothing)

# starts a SimplexTree as a disjoint collection of n nodes (zero-simplices)
function SimplexTree(v::AbstractVector{Int})
    tree = SimplexTree()
    for i ∈ v
        child!(tree, i)
    end
    tree
end

SimplexTree(n::Int) = SimplexTree(1:n)

Node = Int
Nodes = Vector{Node}
NodeTuple{K} = NTuple{K, Node}
# NodeContainer{K} = Union{Nodes, NodeTuple{K}}

struct SimplexIterator
    st::SearchTree{SimplexTree}
    parentmap::Vector{Nodes}
    SimplexIterator(tree::SimplexTree; kwargs...) = new(SearchTree(tree, DFS; kwargs...), [])
end

@forward SimplexIterator.st Base.length
Base.filter(f, si::SimplexIterator) = [ s for s ∈ si if f(s) ]

function Base.iterate(si::SimplexIterator, state::TreeSearchState{SimplexTree}=TreeSearchState(si.st.tree))
    newiter = iterate(si.st, state)
    isnothing(newiter) && return nothing
    (subtree, level), newstate = newiter
    
    simp = @match level begin
        0 => Node[]
        1 => [subtree.node]
        _ => push!(copy(si.parentmap[level-1]), subtree.node)
    end

    if level > 0
        if level > length(si.parentmap)
            push!(si.parentmap, simp)
        else
            si.parentmap[level] = simp
        end
    end

    (simp, newstate)
end

simplices(tree::SimplexTree) = collect(SimplexIterator(tree))

simplices(tree::SimplexTree, k::Int)::Vector{Nodes} = filter(
    s -> length(s)==k+1,
    SimplexIterator(tree; maxlevel=k+1)
)

function relabeled(tree::SimplexTree, switched::NTuple{2, Int}...)
    reindex = DefaultDict(
		k -> k,
		Dict(
			p1 => p2
			for p ∈ switched
			for (p1, p2) ∈ (p, reverse(p))
		);
		passkey=true
	)

    ret = SimplexTree()
    # reinserting _every_ simplex is super inefficient but currently getting the job done
    for simp ∈ simplices(tree)
        insert!(ret, getindex.(Ref(reindex), simp))
    end
    ret
end

function get_sorted(tree::SimplexTree, simplex::Nodes)
    t = tree
    for ind ∈ simplex
         t = child(t, ind)
         isnothing(t) && break
    end
    t
end

has_sorted(tree::SimplexTree, simplex::Nodes) = !isnothing(get_sorted(tree, simplex))

function insert_sorted!(tree::SimplexTree, simplex::Nodes; level::Int=1)
    rest = simplex
    ret = []
    for _=1:length(simplex)
        ind, rest = rest[1], rest[2:end]
        hasind, t = child!(tree, ind)
        push!(ret, (t, hasind, level))
        append!(ret, insert_sorted!(t, rest; level=level+1))
    end
    ret
end

Base.get(tree::SimplexTree, simplex::Nodes) = get_sorted(tree, sort(simplex))
has(tree::SimplexTree, simplex::Nodes) = has_sorted(tree, sort(simplex))
Base.insert!(tree::SimplexTree, simplex::Nodes) = insert_sorted!(tree, sort(simplex))
Base.insert!(tree::SimplexTree, simplex::NodeTuple{K}) where K = insert!(tree, collect(simplex))
Base.insert!(tree::SimplexTree, edge::Edge) = insert!(tree, NodeTuple{2}(edge))

struct FaceIterator{C}
    simplex::C
end

function Base.iterate(fi::FaceIterator, state::Int=1)
    state > length(fi.simplex) && return nothing
    ((Tuple(fi.simplex[InvertedIndices.Not(state)]), fi.simplex[state]), state+1)
end

faces(simplex::NodeTuple{K}) where K = FaceIterator{NodeTuple{K}}(simplex)
faces(simplex::Nodes) = FaceIterator{Nodes}(simplex)

function Graphs.adjacency_matrix(tree::SimplexTree)
    nvs = length(simplices(tree, 0))
    es = simplices(tree, 1)
    I, J = [ getindex.(es, i) for i ∈ 1:2 ]
    M = sparse(I, J, ones(length(I)), nvs, nvs)
    M+M' # slapdash symmetrization
end

function cliquecomplex(A::AbstractMatrix, kmax::R=2) where R<:Real
    t = SimplexTree(size(A, 1))
    vr = vietorisrips(A, kmax)
    for simplices ∈ getfield.(vr[2:end], :simplices)
        insert!.(Ref(t), simplices)
    end
    t
end

cliquecomplex(g::AbstractGraph, kmax::R=2) where R<:Real = cliquecomplex(adjacency_matrix(g), kmax)
cliquecomplex(tree::SimplexTree, kmax::R=2) where R<:Real = cliquecomplex(adjacency_matrix(tree), kmax)
