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
            # $(Representation.Distance) => distance_kernel(region; normalization, repr_kwargs...)
            $(Representation.Distance) => error("can't use Submatrix Subrepresentation on a Distance Representation")
            $(Representation.Dhillon) => dhillon_representation(region; normalization, repr_kwargs...)
        end
        new{typeof(M)}(PartitionTree(n_simp(region)), M, representation, basis, input, method, eigenmethod)
    end
end

root(part::SubmatrixPartition) = part.root
n(part::SubmatrixPartition) = partlen(part.root)

function is_trivial_leaf(part::SubmatrixPartition, node::PartitionTree)
    if part.representation == Representation.Dhillon
        II, JJ = separate_rows_cols(node.sinds, part.matrix)
        return length(II) == 1 || length(JJ) == 1
    end

    return isleaf(node)
end

is_trivial_pair(::SubmatrixPartition, node::PartitionTree) = partlen(node) == 2

function _partition!(part::SubmatrixPartition, node::PartitionTree)
    if part.representation == Representation.Dhillon
        i_inds, j_inds = separate_rows_cols(node.sinds, part.matrix)
        M = part.matrix[i_inds, j_inds]
    else
        M = part.matrix[node.sinds, node.sinds]
    end
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

    parts = @match part.method begin
        $(PartitionMethod.FiedlerSign) => part_by_sign(modbasis[:, eignum])
        $(PartitionMethod.Kmeans) => part_by_kmeans(modbasis[:, 1:eignum])
    end

    PartitionOutput(; vals, basis=modbasis, parts)
end
