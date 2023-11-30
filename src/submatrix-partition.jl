export SubmatrixPartition

struct SubmatrixPartition{M<:AbstractMatrix} <: SCPartition
    root::PartitionTree
    matrix::M
    representation::Representation.T
    basis::Basis.T
    input::PartitionInput.T
    method::PartitionMethod.T
    eigenmethod::EigenMethod.T
    function SubmatrixPartition(
        region::R;
        representation::Representation.T=Representation.KLaplacian,
        normalization::Normalization.T=Normalization.Symmetric,
        basis::Basis.T=Basis.Fiedler,
        input::PartitionInput.T=PartitionInput.DCOrientation,
        method::PartitionMethod.T=PartitionMethod.FiedlerSign,
        eigenmethod::EigenMethod.T=EigenMethod.First,
        repr_kwargs...
    ) where {R<:Region}
        M = @match representation begin
            $(Representation.KLaplacian) => k_laplacian(region; normalization, repr_kwargs...)
            # Representation.Distance => distance_kernel(region; normalization, repr_kwargs...)
            $(Representation.Distance) => error("can't use Submatrix Subrepresentation on a Distance Representation")
        end
        new{typeof(M)}(PartitionTree(n_simp(region)), M, representation, basis, input, method, eigenmethod)
    end
end

root(part::SubmatrixPartition) = part.root
n(part::SubmatrixPartition) = partlen(part.root)

function _partition!(part::SubmatrixPartition, node::PartitionTree)
    M = part.matrix[node.sinds, node.sinds]
    part.representation == Representation.Distance && error("can't use Submatrix Subrepresentation on a Distance Representation")

    # ccs = connected_components(Graph(M))
    # if length(ccs) > 1
    #     # @info "disconnected $(length(ccs))"
    # 	return PartitionOutput(part_components(ccs))
    # end

    vals, basis, eignum = fiedler_eigs(part.representation, part.basis, part.eigenmethod)(M)

    modbasis = copy(basis) |> B -> @match part.input begin
        $(PartitionInput.Identity) => B
        $(PartitionInput.DCOrientation) => dc_orientation(B)
    end

    PartitionOutput(;
        vals,
        basis=modbasis,
        parts=@match part.method begin
            $(PartitionMethod.FiedlerSign) => part_by_sign(modbasis[:, eignum])
            $(PartitionMethod.Kmeans) => part_by_kmeans(modbasis[:, 1:eignum])
        end
    )
end
