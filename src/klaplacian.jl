export
    Normalization, WhichLap,
    k_laplacian,
    HullPartition

@enumx Normalization begin
    Combinatorial    # B'B + BB'                    maintains homology, ignores weights
    Decoupled        # B'UB + BVB'                  maintains homology
    Weighted         # WB'U⁻¹BW + BVB'              breaks homology
    Symmetric        # √W B'U⁻¹B √W + √W⁻¹BVB'√W⁻¹  maintains homology
end

@enumx HullPartition BothKeep BothDrop

@enumx WhichLap Lower Upper Both

function k_laplacian(
    region::KRegion;
    oris::BitVector=trues(n_simp(region)),
    normalization::Normalization.T=Normalization.Weighted,
    which::WhichLap.T=WhichLap.Both,
    Clower = 1.0,
    Cupper = 1.0,
    with_diagonal=true,
    separate::Bool=false,
    oriented=true,
    tol=1e-14
)
    @assert !separate || which==WhichLap.Both "separate=true requires which=WhichLap.Both"
    L = @NamedTuple{i::Int, j::Int, lower::Float64, upper::Float64, natcon::Bool}[]
    sizehint!(L, sum(length(adj.pairs) for (_, adj) ∈ region.hingemap))

    Dlower = zeros(n_simp(region))

    for (hinge, adj) ∈ region.hingemap
        wt_lower = @match normalization begin
            Normalization.Combinatorial => 1.0
            Normalization.Decoupled => region.bdweights[hinge]
            _ => region.bdweights[hinge] ^ -1.0
        end
        Dlower[collect(values(adj.leaves))] .+= wt_lower

        for adjsimps ∈ adj.pairs
            wt_upper = @match normalization begin
                Normalization.Combinatorial => haskey(region.hullweights, adjsimps.hull) ? 1.0 : 0.0
                _ => get(region.hullweights, adjsimps.hull, 0.0)
            end
            push!(
                L, 
                (
                    i = adj.leaves[adjsimps.leaves[2]],
                    j = adj.leaves[adjsimps.leaves[1]],
                    lower = wt_lower,
                    upper = wt_upper,
                    natcon = adjsimps.natconsistent
                )
            )
        end
    end

    sparse_constructor(v::Vector{Float64}) = sparse(getfield.(L, :i), getfield.(L, :j), v, n_simp(region), n_simp(region))
    sparse_constructor(v) = isempty(v) ? 
        sparse_constructor(Float64[]) : # encountered problem for "dust"; empty A matrix.
        error("should not be here")

    Alower = sparse_constructor(getfield.(L, :lower)) .* Clower
    Aupper = sparse_constructor(getfield.(L, :upper)) .* Cupper

    S = sparse_constructor(2.0*getfield.(L, :natcon) .- 1)

    dk = getweights(region.simplexweights)
    sandwich(M, p=1; v=dk) = Diagonal(v.^p) |> D -> D*M*D

    # trick to get row or col sum for just the triangular matrix
    Dlower = with_diagonal ? Dlower : spzeros(n_simp(region))
    Dupper = with_diagonal ? (sum(Aupper, dims=2)[:] + sum(Aupper, dims=1)[:]) ./ (k(region)+1) : spzeros(n_simp(region))

    S = oriented ?
        sandwich(S; v=ifelse.(oris, 1, -1)) :
        abs.(S)

    # usually upper and lower Laplacians contribute oppositely; ori_correction accounts
    # for this depending on whether we use the oriented or unoriented Laplacian
    ori_correction = oriented ? -1 : 1
    Llower = droptol!(Diagonal(Dlower) - S .* Alower, tol)
    Lupper = droptol!(Diagonal(Dupper) - ori_correction .* S .* Aupper, tol)

    Llower, Lupper = @chain begin
        (
            @match normalization begin
                Normalization.Combinatorial => (Llower, Lupper)
                Normalization.Decoupled => (Llower, Lupper)
                Normalization.Weighted => (sandwich(Llower), Lupper)
                Normalization.Symmetric => (sandwich(Llower, 0.5), sandwich(Lupper, -0.5))
            end
        )
        Symmetric.(_)
    end

    @match which begin
        WhichLap.Lower => Llower
        WhichLap.Upper => Lupper
        WhichLap.Both => separate ? (Llower, Lupper) : Llower + Lupper
    end
end

function k_laplacian(
    region::ZeroRegion;
    subregion_inds::(Nothing|Vector{Int})=nothing,
    normalization::Normalization.T=Normalization.Weighted
)
    W = @chain Graphs.weights(region.weights) begin
        normalization == Normalization.Combinatorial ? ifelse.(_.==0.0, 0.0, 1.0) : _
        isnothing(subregion_inds) ? _ : _[subregion_inds, subregion_inds]
        Symmetric
    end
    d = sum(W, dims=2)[:]
    sandwich(M, p=1; v=d) = Diagonal(v.^p) |> D -> D*M*D |> Symmetric

    L = Diagonal(d)-W
    @match normalization begin
        Normalization.Symmetric => sandwich(L, -0.5; v = d .+ region.diagmodifier)
        _ => L
    end
end
