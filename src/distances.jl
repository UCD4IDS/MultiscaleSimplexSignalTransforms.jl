export
    Sandwich,
    distance_kernel,
    dijkstra,
    censored_distmat, distmat

struct Sandwich{T<:Real,X<:AbstractMatrix{T},V<:AbstractVector{T}} <: AbstractMatrix{T}
    M::X
    p::V
    u::V
    v::V
    s::Float64
    function Sandwich(M, v=ones(size(M, 1)))
        p = v ./ norm(v)
        u = M * p
        v = M' * p
        new{eltype(M),typeof(M),typeof(p)}(M, p, u, v, u' * p)
    end
end

Base.show(io::IO, ::MIME"text/plain", A::Sandwich) = print(io, join([
        join(string.(size(A)), '×'),
        repr("text/plain", typeof(A); context=:compact => true)
    ], ' '))

function LinearAlgebra.mul!(y, A::Sandwich, x)
    px = A.p' * x
    vx = A.v' * x
    y = A.M * x + (A.s * px - vx) .* A.p - px .* A.u
end

Base.size(A::Sandwich) = size(A.M)
Base.getindex(A::Sandwich, i, j) = A.M[i, j] + (A.s * A.p[j] - A.v[j]) * A.p[i] - A.p[j] * A.u[i]
LinearAlgebra.issymmetric(A::Sandwich) = issymmetric(A.M)
SparseArrays.nnz(A::Sandwich) = nnz(A.M)
SparseArrays.nnz(A::Sandwich{T,S}) where {T,S<:Symmetric} = nnz(parent(A.M))
Arpack.eigs(A::Sandwich; kwargs...) = issymmetric(A) ?
                                      eigs(Symmetric(A); kwargs...) :
                                      error("for some reason asymmetric eigs output incorrect values, so this is disabled")

# why are these not already in SparseArrays??
SparseArrays.nnz(M::Symmetric{T,S}) where {T,S<:AbstractSparseArray} = nnz(parent(M))
SparseArrays.nonzeros(M::Symmetric{T,S}) where {T,S<:AbstractSparseArray} = nonzeros(parent(M))

function dijkstra(
    g::AbstractGraph,
    src::Int,
    distmx::AbstractMatrix=Graphs.weights(g);
    dmax=Inf
)
    dmax = Float64(dmax)
    dists = sparsevec([src], [0.0], nv(g))

    H = PriorityQueue{Int,Float64}()
    H[src] = 0.0

    while !isempty(H)
        u = dequeue!(H)
        d = dists[u]
        for v in outneighbors(g, u)
            alt = min(d + distmx[u, v], dmax)
            if alt < dmax && (v ∉ findnz(dists)[1] || alt < dists[v])
                dists[v] = alt
                H[v] = alt
            end
        end
    end

    dists
end

distmat(g::AbstractGraph; dmax=Inf) = reduce(
    hcat,
    (dijkstra(g, i; dmax) for i ∈ 1:nv(g))
)'

distmat(W::AbstractMatrix; dmax=Inf) = distmat(SimpleWeightedGraph(W); dmax)

# NOTE: censored_distmat already adds dmax*I !
function censored_distmat(g::AbstractGraph, dmax::Real)
    dmax == Inf && error("cannot censor with dmax=Inf, use distmat instead")
    nvg = nv(g)

    reduce(
        hcat,
        (
            findnz(dists) |> nzds -> sparsevec(nzds[1], dmax .- nzds[2], nvg)
            for dists ∈ (dijkstra(g, i; dmax) for i ∈ 1:nvg)
        )
    )'
end

censored_distmat(W::AbstractMatrix, dmax::Real; issymmetric=true) = censored_distmat(
    issymmetric ? SimpleWeightedGraph(W) : SimpleWeightedDiGraph(W),
    dmax
)

distance_kernel(region::Region; kwargs...) = error("can only use Representation.Distance for k=0 for now")

function distance_kernel(
    region::ZeroRegion;
    subregion_inds::(Nothing | Vector{Int})=nothing,
    normalization::Normalization.T=Normalization.Weighted,
    withdistmat=false,
    dmax=Inf,
    issymmetric=true
)
    W = @chain Graphs.weights(region.weights) begin
        ifelse.(
            _ .== 0,
            0.0,
            normalization == Normalization.Combinatorial ? 1.0 : _ .^ -1
        )
        isnothing(subregion_inds) ? _ : _[subregion_inds, subregion_inds]
    end

    if issymmetric
        W = Symmetric(W)
    end

    K = 0.5 * (isinf(dmax) ? -distmat(W) : censored_distmat(W, dmax; issymmetric))

    if issymmetric
        K = Symmetric(K)
    end

    d = Vector(sum(W, dims=2)[:]) .^ 0.5
    D = Diagonal(d)

    norm(d) == 0 && error("uh oh bad degree vector for $(size(region.weights))")

    kernel = @match normalization begin
        $(Normalization.Symmetric) => Sandwich(Symmetric(D * K * D), d)
        _ => Sandwich(K)
    end

    withdistmat ? (K, kernel) : kernel
end

distance_kernel(X; issymmetric=true, kwargs...) = distance_kernel(
    ZeroRegion(
        issymmetric ? SimpleWeightedGraph(X) : SimpleWeightedDiGraph(X);
        issymmetric
    );
    issymmetric, kwargs...
)
