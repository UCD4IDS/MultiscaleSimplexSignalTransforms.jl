export
    lap_from_adj

lap_from_adj(M, k) = KRegion(cliquecomplex(M, k+1), k; weakadjweight=1.0) |>
    reg -> (
		k_laplacian(
			reg;
			normalization=Normalization.Combinatorial
		),
		simplices(reg)
	)
