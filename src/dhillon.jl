export BipartiteRegion, dhillon_representation

struct BipartiteRegion{A<:AbstractMatrix} <: Region
    data::A
    BipartiteRegion(data::A) where {A<:AbstractMatrix} = new{A}(data)
end

entire_adjacency_matrix(region::BipartiteRegion) = [
    zeros(size(region.data, 1), size(region.data, 1)) region.data;
    region.data' zeros(size(region.data, 2), size(region.data, 2))
]

function separate_rows_cols(inds, nrows::Int64)
    imax = findfirst(>(nrows), inds)
    isnothing(imax) && return inds, Int64[]
    inds[1:imax-1], inds[imax:end] .- nrows
end

separate_rows_cols(inds, M::AbstractMatrix) = separate_rows_cols(inds, size(M, 1))

k(::BipartiteRegion) = 0
n_simp(region::BipartiteRegion) = sum(size(region.data))
n_bd(::BipartiteRegion) = 0
n_hull(region::BipartiteRegion) = sum(region.data .> 0)

complex(region::BipartiteRegion) = cliquecomplex(entire_adjacency_matrix(region), 1)

simplices(region::BipartiteRegion) = [(i,) for i in 1:n_simp(region)]
boundaries(::BipartiteRegion) = []
hulls(region::BipartiteRegion) = map(
    e -> (e[1], e[2]+size(region.data, 1)),
    Iterators.filter(
        x -> x[3]>0,
        ((II[1], II[2], region.data[II]) for II in CartesianIndices(region.data))
    )
)

function dhillon_representation(
    region::BipartiteRegion;
    subregion_inds::(Nothing | Vector{Int})=nothing,
    normalization::Normalization.T=Normalization.Symmetric
)
    @assert normalization == Normalization.Symmetric "Dhillon method requires symmetric normalization"

    W = @chain region.data begin
        isnothing(subregion_inds) ? _ : _[subregion_inds, subregion_inds]
    end

    d_r = sum(W, dims=2)[:]
    d_c = sum(W, dims=1)[:]
    Diagonal(d_r .^ -0.5) * W * Diagonal(d_c .^ -0.5)
end
