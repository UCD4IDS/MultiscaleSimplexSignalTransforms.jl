function sponge(W, nev, τ₊ = 1, τ₋ = 1)
	# here W should represent a signed graph
	# symmetric, with either positive or negative entry
	usesparse = nev*3 < size(W, 1)

	# obtain the bottom generalized eigenvectors
	signedlaps(W, τ₊, τ₋) |> (
		Ms -> usesparse ?
			lobpcg(Ms..., false, nev) :
			eigen(Ms...)
	) |> (
		eigstruct -> usesparse ?
			(eigstruct.λ, eigstruct.X) :
			(eigstruct.values[1:nev], eigstruct.vectors[1:nev])
	)
end

function signedlaps(W, τ₊ = 1, τ₋ = 1)
	# here W should represent a signed graph
	# symmetric, with either positive or negative entry
	W₊ = (W .> 0) .* W
	D₊ = Diagonal(sum(W₊, dims=2)[:])
	L₊ = D₊ - W₊

	W₋ = -(W .< 0) .* W
	D₋ = Diagonal(sum(W₋, dims=2)[:])
	L₋ = D₋ - W₋

	(L₊ + τ₋*D₋, L₋ + τ₊*D₊)
end

function clusteredges()
	# here W should represent a directed graph
	# asymmetric



end