# abstract type AdjMat end

# struct ZeroAdj{M<:AbstractMatrix} <: AdjMat
#     W::M
#     ZeroAdj(W::M) where M<:AbstractMatrix = new{M}(W)
# end

# adjmat(reg::ZeroRegion) = ZeroAdj(Graphs.weights(reg.weights))
# laplacian(adj::ZeroAdj, normalization::Normalization.T) = k_laplacian(
#     ZeroRegion(adj.W)
# )

# struct KAdj{M<:AbstractMatrix} <: AdjMat
#     Wup::M
#     Wdown::M
#     KAdj(W::M,V::M) where M<:AbstractMatrix = new{M}(W, V)
# end

# adjmat(reg::KRegion; lap_kwargs...) = KAdj(.-k_laplacian(
#     reg; 
#     with_diagonal=false,
#     separate=true,
#     lap_kwargs...
# )...)

export
    FullPartition

struct FullPartition{T<:SimplexTree} <: SCPartition
    root::PartitionTree
    complex::T
    regionfn
    representation::Representation.T
    normalization::Normalization.T
    basis::Basis.T
    input::PartitionInput.T
    method::PartitionMethod.T
    eigenmethod::EigenMethod.T
    repr_kwargs::Dict
    function FullPartition(
        complex::T, regionfn;
        representation::Representation.T=Representation.KLaplacian,
        normalization::Normalization.T=Normalization.Symmetric,
        basis::Basis.T=Basis.Fiedler,
        input::PartitionInput.T=PartitionInput.DCOrientation,
        method::PartitionMethod.T=PartitionMethod.FiedlerSign,
        eigenmethod::EigenMethod.T=EigenMethod.First,
        repr_kwargs...
    ) where {T<:SimplexTree}
        new{T}(PartitionTree(
            n_simp(regionfn(; complex, subregion_inds=nothing))),
            complex, regionfn, representation, normalization, basis, input, method, eigenmethod, repr_kwargs
        )
    end
end

FullPartition(
    region::R,
    regionfn = (; complex, subregion_inds) -> KRegion(
        complex, k(region); subregion_inds, weakadjweight=(k(region) == 0 ? 0.0 : 1.0)
    );
    kwargs...
) where {R<:Region} = FullPartition(complex(region), regionfn; kwargs...)

root(part::FullPartition) = part.root
n(part::FullPartition) = partlen(root(part))

function _partition!(part::FullPartition, node::PartitionTree)
    region = part.regionfn(complex = part.complex, subregion_inds = node.sinds)

    M = @match part.representation begin
        Representation.KLaplacian => k_laplacian(region; part.normalization, part.repr_kwargs...)
        Representation.Distance => distance_kernel(region; part.normalization, part.repr_kwargs...)
    end

    # Ll, Lu = k_laplacian(region; part.normalization, separate=true)
    # M = Ll+Lu

    vals, basis, eignum = fiedler_eigs(part.representation, part.basis, part.eigenmethod)(M)

    # ccs = connected_components(Graph(abs.(Ll)+abs.(Lu)))
    # if length(ccs) > 1
    #     # @info "disconnected $(length(ccs))"
	# 	return PartitionOutput(part_components(ccs))
	# end

    modbasis = copy(basis) |> B -> @match part.input begin
        PartitionInput.Identity => B
        PartitionInput.DCOrientation => dc_orientation(B)
    end

    PartitionOutput(;
        vals,
        basis = modbasis,
        parts = @match part.method begin
            PartitionMethod.FiedlerSign => part_by_sign(modbasis[:,eignum])
            PartitionMethod.Kmeans => part_by_kmeans(modbasis[:,1:eignum])
        end
    )
end
