### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 7421b50c-9270-11ef-1912-614e3be9a306
begin
	using Pkg, Revise
	Pkg.activate(".")
end

# ╔═╡ 4c8ec920-0212-46ae-91c5-57a86f4d3896
using LinearAlgebra, Arpack, MultiscaleSimplexSignalTransforms, Plots, Graphs, CSV, DataFrames, DelimitedFiles

# ╔═╡ 21349d65-3e23-4df1-abed-fd5fde27f803
n = 10

# ╔═╡ 956ff7d0-3267-4dca-8449-b6427694a0e7

	spy(adjacency_matrix(
		stochastic_block_model(49, 0, [50,50,50])
	))


# ╔═╡ 2155ccde-1e2d-489e-bc8d-c920a173c4be
begin
	Mblock = [ zeros(n,n) ones(n,n); ones(n,n) ones(n,n); ones(n,n) zeros(n,n) ]
	heatmap(Mblock)
end

# ╔═╡ 5fec5407-e6b1-4614-8d4f-d691d0cd553f
r = BipartiteRegion(Mblock)

# ╔═╡ e3470d3c-586c-40b8-8d7a-10c224ef2fd1
d = kGHWT(r, representation = Representation.Dhillon)

# ╔═╡ 2639b575-71dc-4e48-b266-0832f8a8d27a
r2 = ZeroRegion(
	[zeros(3n,3n) Mblock; Mblock' zeros(2n,2n)]
)

# ╔═╡ 98d6c815-39af-40a2-914f-0ad24d4ae659
d2 = kGHWT(r2)

# ╔═╡ 0e502d95-c6df-48cf-9030-180dcf5d0255
grouppos(n, k, offset) = [
	n/2*(3k-1)-v+offset
	for i in 0:k-1
	for v in i*n/2 .+ (n*i+1 : n*(i+1))
]

# ╔═╡ b25948bf-2b9f-4ca3-831a-bdc2c76296ab
pos = Pos(
	[zeros(3n); 50*ones(2n)],
	[grouppos(n, 3, 0); grouppos(n, 2, 3n/4)]
)

# ╔═╡ da85495e-38da-4eda-9e3a-6267cb40bda3
(5, 3) |> ((n, k),) -> scatter(zeros(n*k), grouppos(n, k, 0), zcolor = repeat(1:k, inner=n)); scatter!(ones(10), grouppos(5, 2, 15/4))

# ╔═╡ 4496133d-62af-45f1-a0e5-a0c85b005759
splot(r, d2[0,0,1], pos, simpsize=10)

# ╔═╡ 3a0ca20e-4461-46a8-8c35-d979230eaba6
dfs(
	tsn -> MultiscaleSimplexSignalTransforms.separate_rows_cols(tsn.node.sinds, 30),
	d2.part.root
)

# ╔═╡ 63fd2cd3-9a39-483d-864f-cb342e3f9ffc
M = readdlm("../data/small_layer_6.csv", ',', Float64)

# ╔═╡ Cell order:
# ╠═7421b50c-9270-11ef-1912-614e3be9a306
# ╠═4c8ec920-0212-46ae-91c5-57a86f4d3896
# ╠═21349d65-3e23-4df1-abed-fd5fde27f803
# ╠═956ff7d0-3267-4dca-8449-b6427694a0e7
# ╠═2155ccde-1e2d-489e-bc8d-c920a173c4be
# ╠═5fec5407-e6b1-4614-8d4f-d691d0cd553f
# ╠═e3470d3c-586c-40b8-8d7a-10c224ef2fd1
# ╠═2639b575-71dc-4e48-b266-0832f8a8d27a
# ╠═98d6c815-39af-40a2-914f-0ad24d4ae659
# ╠═0e502d95-c6df-48cf-9030-180dcf5d0255
# ╠═b25948bf-2b9f-4ca3-831a-bdc2c76296ab
# ╠═da85495e-38da-4eda-9e3a-6267cb40bda3
# ╠═4496133d-62af-45f1-a0e5-a0c85b005759
# ╠═3a0ca20e-4461-46a8-8c35-d979230eaba6
# ╠═63fd2cd3-9a39-483d-864f-cb342e3f9ffc
