### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ d8cba18c-2f63-11ee-062b-337a52a86175
using Pkg; Pkg.activate(".")

# ╔═╡ 2d1ddd7f-f3dc-4310-bfb1-84a197e32dc9
using LinearAlgebra

# ╔═╡ 04f67a57-3a5d-479b-85d2-f2f3deac35e3
using MultiscaleSimplexSignalTransforms, Graphs, Arpack

# ╔═╡ bead31a0-1589-4517-bea4-943128676de7
n = 100

# ╔═╡ 0a4c2fb0-dbf8-440a-9d9e-a0dd666377f2


# ╔═╡ ef5f18e2-8448-4edf-b0aa-2d9679bf3109
g = path_graph(n)

# ╔═╡ c8a2b78e-5025-4f4e-b7ca-1fa3e0606957


# ╔═╡ 4d1715c9-e7d1-4b8f-8f0c-c43a949eb80d
region = ZeroRegion(g)

# ╔═╡ d26bf1d9-15a9-4486-959d-ca0cc65b9c84
L = k_laplacian(region)

# ╔═╡ 7b8f5f94-3c8b-4b35-babe-d51f871fcefc
pl = f -> splot(
	region, f; 
	xs = collect(Float64.(1:n)), 
	ys = zeros(n), 
	palette=:redblue,
	size=10
)

# ╔═╡ b8791d74-0008-455b-aaef-99118f0571be
pls = M -> splot_many(
	region, collect(eachcol(M))...; 
	xs = collect(Float64.(1:n)), 
	ys = zeros(n), 
	palette=:redblue,
	size=10
)

# ╔═╡ a9102099-b4eb-4842-8c18-bb1761243c7f
Λ, X = eigs(L, which=:SM)

# ╔═╡ 8e134977-e7b8-4adb-8a0b-780b08d6a687
pls(X)

# ╔═╡ Cell order:
# ╠═d8cba18c-2f63-11ee-062b-337a52a86175
# ╠═2d1ddd7f-f3dc-4310-bfb1-84a197e32dc9
# ╠═bead31a0-1589-4517-bea4-943128676de7
# ╠═04f67a57-3a5d-479b-85d2-f2f3deac35e3
# ╠═0a4c2fb0-dbf8-440a-9d9e-a0dd666377f2
# ╠═ef5f18e2-8448-4edf-b0aa-2d9679bf3109
# ╠═c8a2b78e-5025-4f4e-b7ca-1fa3e0606957
# ╠═4d1715c9-e7d1-4b8f-8f0c-c43a949eb80d
# ╠═d26bf1d9-15a9-4486-959d-ca0cc65b9c84
# ╠═7b8f5f94-3c8b-4b35-babe-d51f871fcefc
# ╠═b8791d74-0008-455b-aaef-99118f0571be
# ╠═a9102099-b4eb-4842-8c18-bb1761243c7f
# ╠═8e134977-e7b8-4adb-8a0b-780b08d6a687
