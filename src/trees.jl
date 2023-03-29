export SearchTree,
    dfs, bfs,
    DFS, ReversedDFS, BFS, ReversedBFS,
    GTree, children

abstract type Tree end
# interface
# Base.length(::Tree)::Int = error("implement!")
children(::Tree)::Vector{<:Tree} = error("implement!")

Base.show(io::IO, ::MIME"text/plain", tree::Tree) = print(
    io,
    repr("text/plain", typeof(tree); context = :compact=>true)
)

@enum SearchMethod DFS BFS ReversedDFS ReversedBFS

"""
for DFS methods, should pop! (from the end)
for BFS methods, should popfirst! (from the front)
"""
shouldpop(st::SearchMethod) = st ∈ [DFS, ReversedDFS]

"""
for DFS, default is reversed (so ReversedDFS should not reverse)
for BFS, default is normal (so ReversedBFS should reverse)
"""
shouldreverse(st::SearchMethod) = st ∈ [DFS, ReversedBFS]

struct TreeSearchNode{T<:Tree}
    node::T
    level::Int
end

# convenience functions to auto-deconstruct a tuple from a TreeSearchNode
Base.iterate(tsn::TreeSearchNode) = (tsn.node, false)
Base.iterate(tsn::TreeSearchNode, done::Bool) = done ? nothing : (tsn.level, true)

# keep prev node in order to wait until after results of search operations to continue
# iterating into children. This allows building the tree on the fly with the search
mutable struct TreeSearchState{T<:Tree}
    next::Vector{TreeSearchNode{T}}
    prev::(Nothing|TreeSearchNode{T})
end

TreeSearchState(t::T) where T<:Tree = TreeSearchState{T}([TreeSearchNode(t, 0)], nothing)

next!(tss::TreeSearchState; popside::Bool) = popside ? pop!(tss.next) : popfirst!(tss.next)
Base.isempty(tss::TreeSearchState) = isempty(tss.next)
Base.push!(tss::TreeSearchState{T}, item::TreeSearchNode{T}) where T<:Tree = 
    push!(tss.next, item)

struct SearchTree{T<:Tree, R<:Real}
    tree::T
    method::SearchMethod
    maxlevel::R
    itercond::Function
    SearchTree(
        tree::T, method::SearchMethod;
        maxlevel::R=Inf, itercond::Function=(_)->true
    ) where {T<:Tree, R<:Real} = new{T, R}(tree, method, maxlevel, itercond)
end

function Base.iterate(
    tree::SearchTree{T},
    state::TreeSearchState{T}=TreeSearchState(tree.tree)
) where T<:Tree
    if !isnothing(state.prev)
        cs = children(state.prev.node)
        shouldreverse(tree.method) && reverse!(cs)

        for c in cs
            cnode = TreeSearchNode(c, state.prev.level+1)
            tree.itercond(cnode) && push!(state, cnode)
        end
    end

    popside = shouldpop(tree.method)

    while !isempty(state)
        tsn = next!(state; popside)
        tsn.level > tree.maxlevel && continue
        state.prev = tsn
        return (tsn, state)
    end

    nothing
end

Base.eltype(::T) where T<:Tree = TreeSearchNode{T}
Base.eltype(::SearchTree{T}) where T<:Tree = TreeSearchNode{T}

# in order to use a Tree in generator expressions, need to compute all descendants separately
# at least once. I am going for a lazy/expensive solution instead of caching/tracking
function Base.length(st::SearchTree)
    count = 0
    for _ ∈ st
        count += 1
    end
    count
end

Base.length(tree::T) where T<:Tree = length(SearchTree(tree, ReversedDFS))

function treemap(
    f::Function, tree::T;
    method::SearchMethod, maxlevel::Real, itercond::Function
) where T<:Tree
    res = []
    for tsn ∈ SearchTree(tree, method; maxlevel, itercond)
        push!(res, f(tsn))
    end
    res
end

# dfs and itercond functions are applied to the iterator: a TreeSearchNode with node and level fields
dfs(
    f::Function, tree::T;
    righttoleft::Bool=false, maxlevel::Real=Inf, itercond::Function =(_)->true
) where T<:Tree = treemap(f, tree; method = righttoleft ? ReversedDFS : DFS, maxlevel, itercond)

# bfs and itercond functions are applied to the iterator: a TreeSearchNode with node and level fields
bfs(
    f::Function, tree::T;
    righttoleft::Bool=false, maxlevel::Real=Inf, itercond::Function =(_)->true
) where T<:Tree = treemap(f, tree; method = righttoleft ? ReversedBFS : BFS, maxlevel, itercond)

struct GTree{G<:AbstractGraph} <: Tree
    g::G
    visited::Vector{Bool}
    node::Int
end

function GTree(g::AbstractGraph)
    @assert ne(g) == nv(g)-1 "This graph doesn't look like a tree"
    @assert length(connected_components(g)) == 1 "This graph looks disconnected"
    GTree(g, [true; fill(false, nv(g)-1)], 1)
end

function children(tree::GTree{G}) where G
    res = GTree{G}[]
    for node ∈ neighbors(tree.g, tree.node)
        tree.visited[node] && continue
        tree.visited[node] = true
        push!(res, GTree(tree.g, tree.visited, node))
    end
    res
end

# re-order the nodes of a tree, in order to represent a rooted tree
# the first node is the root, and the rest have only earlier nodes as parents
# returns the adjacency matrix of the rooted tree
function rooted(g::AbstractGraph, usedfs=true)
    searchfn = usedfs ? dfs : bfs
    perm = searchfn(tsn -> tsn.node.node, GTree(g)) .|> Int
    permute(adjacency_matrix(g), perm, perm)
end
