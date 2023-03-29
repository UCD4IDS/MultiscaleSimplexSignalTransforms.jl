### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ b3b28169-06a3-4274-b965-161a293a5362
using Pkg; Pkg.activate(".")

# ╔═╡ 64b484a1-b5c9-4558-987c-a1a20f1ffc9b
using Revise

# ╔═╡ 60a5352d-e090-4073-9f98-792caa9cf200
using MultiscaleSimplexSignalTransforms

# ╔═╡ ae69c53e-cb07-11ec-207c-e7380f19ca22
using Images, StatsBase, Plots, Graphs, GraphRecipes, SparseArrays, LinearAlgebra, Arpack, Clustering

# ╔═╡ b6baf5ec-6cbb-4ebd-bbf9-a18deceb49f0
im = load("../data/cameraman.png") .|> Gray

# ╔═╡ f0c38d49-1657-4abb-97a4-46fa8c63c969
[pixel.val |> reinterpret for pixel in im[:]] |> histogram # countmap

# ╔═╡ c5bc4665-a3fd-409e-a147-61680fad0d65
function pixelmesh(n, m=n; tribackslash=true)
	g = Graphs.grid((n, m))
	tris = Simplex[]
	for i=1:n-1, j=1:m-1
		if tribackslash
			add_edge!(g, i+1+(j-1)*n, i+j*n)   # diagonal edges look like /
			push!(
				tris, 
				Simplex([i+1+(j-1)*n, i+j*n, i+(j-1)*n]), 
				Simplex([i+1+(j-1)*n, i+j*n, i+1+j*n])
			)
		else
			add_edge!(g, i+(j-1)*n, i+1+j*n)   # diagonal edges look like \
			push!(
				tris, 
				Simplex([i+1+(j-1)*n, i+1+j*n, i+(j-1)*n]), 
				Simplex([i+(j-1)*n, i+j*n, i+1+j*n])
			)
		end
	end
	(g, tris)
end

# ╔═╡ b09fface-feb0-459a-9bd2-a3de46c64ed7
(4, 5) |> sz -> graphplot(
	pixelmesh(sz...)[1];
	names=1:prod(sz), 
	x = repeat(1:sz[1], sz[2]), 
	y = repeat(sz[2]:-1:1, 1, sz[1])'[:],
	curves=false
)

# ╔═╡ 22b04ce8-83a5-4a4d-93ee-d7d47658ade7
edgesignal(f_node, edges; oris = [src(e) < dst(e) for e in edges]) = [
	(f_node[dst(e)]-f_node[src(e)]) * (oris[i] ? 1 : -1)
	for (i, e) in enumerate(edges)
]

# ╔═╡ 7faf58a1-7f74-4cfc-83b9-9ee9a638da02
trisignal(f_edge, tris; oris = trues(length(f_edge))) = [
	collect(t.nodes) |> 
		ns ->
			f_edge[Edge(ns[1], ns[3])] -
			f_edge[Edge(ns[1], ns[2])] -
			f_edge[Edge(ns[2], ns[3])]
	for (i, t) in enumerate(tris)
]

# ╔═╡ ded12198-437c-41c3-88b4-04d3c8ab6852
grid_kernel_nz(n, m=n; σ=1.0, r=5) = CartesianIndices((n, m)) |> 
	I -> [
		(
			i, 
			jx + (jy-1)*n,
			I[i] - CartesianIndex(jx, jy) |> Tuple |> norm |> x -> exp(-x^2/σ)
		)
		for i=1:n*m
		for jx = max(1, I[i][1]-r):min(n, I[i][1]+r)
		for jy = floor(Int, sqrt(r^2-(I[i][1]-jx)^2)) |> 
			s -> max(1, I[i][2]-s):min(m, I[i][2]+s)
		if i != jx + (jy-1)*n # eliminate the diagonal
	] #|>
	#T -> sparse((getindex.(T, i) for i=1:3)..., n*m, n*m)

# ╔═╡ 4dd6ed35-814c-4d19-863b-6ca4fc68a06d
function image_kernel(im; σ_prox=1.0, σ_color=1.0, r=5)
	sz = size(im)
	psz = prod(sz)
	I = CartesianIndices(sz)
	Wnz = grid_kernel_nz(sz...; σ=σ_prox, r)
	kernelvals = [
		v * exp(-(im[I[i]]-im[I[j]])^2/σ_color)
		for (i, j, v) in Wnz
	]
	sparse(getindex.(Wnz, 1), getindex.(Wnz, 2), kernelvals, psz, psz)
	
end

# ╔═╡ 60a838e4-cb59-47e9-af10-4bc5c0fae87f
Wim = image_kernel(getfield.(im, :val) .|> float64; σ_prox=4.0, σ_color=0.005)

# ╔═╡ a5df41ff-6099-4610-9d8e-3ba6c49c34d4
K = distance_kernel(ZeroRegion(Wim))

# ╔═╡ 0cf1c29b-5e92-4ea4-ad03-6b093833f48b


# ╔═╡ e3dee408-51d6-489a-a5fa-48d08a98d105
Lim = Diagonal(sum(Wim; dims=2)[:])^-0.5 |>
	D -> I - D * Wim * D |>
	Symmetric

# ╔═╡ f8e4b572-790d-4c5d-99d3-5d5a8545eaf1
k = 3

# ╔═╡ 36a0eae7-f2c4-4523-aeec-45a9e30e533d
E = eigs(Lim+I; nev=10, which=:SM, maxiter=10000)

# ╔═╡ 81ed8d28-64e5-4008-ba8d-235ca0e3bc58
clusters = kmeans(E[2][:,1:k]', k).assignments

# ╔═╡ 5f497f29-9b49-4fda-b7d1-c1a85a4762dc
Bool.([1 1 0]) |>
flipped -> Tuple(
	im .* [
		(clusters[i] == kk) |> b -> flipped[kk] ? !b : b
		for (i, I) in enumerate(CartesianIndices(size(im)))
	]
	for kk=1:k
)

# ╔═╡ 08ef307a-e8d2-4778-9b8a-10ee40c4153b
es = pixelmesh(256)[1] |> edges |> collect

# ╔═╡ e45b15b5-f142-4af4-9957-4fe4a45db437
emap = Dict(e => i for (i, e) in enumerate(es))

# ╔═╡ 5998dfe3-2cd3-4b65-89df-b90c9a27e945
esig = edgesignal(getfield.(im'[:], :val) .|> float64, es)

# ╔═╡ fb9d596f-9252-4aeb-b266-a8053397a5d0
Le = hodgelaplacian(es, [src(e)<dst(e) ? 1 : -1 for e in es])

# ╔═╡ 7549afef-22ea-4d98-9413-c2ae786bf930
We = -hodgelaplacian(es, [src(e)<dst(e) ? 1 : -1 for e in es]; with_diagonal=false)

# ╔═╡ 8048bf11-7785-48de-996a-f91173759cd2
function edge_kernel_nz(n, m=n; σ, r)
	I = CartesianIndices((n, m))
	grid_kernel_nz(n, m; σ, r) |>
		Ts -> [
			()
			for (i, j, v) in Ts
			for e in [()]
	]
end

# ╔═╡ Cell order:
# ╠═b3b28169-06a3-4274-b965-161a293a5362
# ╠═64b484a1-b5c9-4558-987c-a1a20f1ffc9b
# ╠═60a5352d-e090-4073-9f98-792caa9cf200
# ╠═ae69c53e-cb07-11ec-207c-e7380f19ca22
# ╠═b6baf5ec-6cbb-4ebd-bbf9-a18deceb49f0
# ╠═f0c38d49-1657-4abb-97a4-46fa8c63c969
# ╠═c5bc4665-a3fd-409e-a147-61680fad0d65
# ╠═b09fface-feb0-459a-9bd2-a3de46c64ed7
# ╠═22b04ce8-83a5-4a4d-93ee-d7d47658ade7
# ╠═7faf58a1-7f74-4cfc-83b9-9ee9a638da02
# ╠═ded12198-437c-41c3-88b4-04d3c8ab6852
# ╠═4dd6ed35-814c-4d19-863b-6ca4fc68a06d
# ╠═60a838e4-cb59-47e9-af10-4bc5c0fae87f
# ╠═a5df41ff-6099-4610-9d8e-3ba6c49c34d4
# ╠═0cf1c29b-5e92-4ea4-ad03-6b093833f48b
# ╠═e3dee408-51d6-489a-a5fa-48d08a98d105
# ╠═f8e4b572-790d-4c5d-99d3-5d5a8545eaf1
# ╠═36a0eae7-f2c4-4523-aeec-45a9e30e533d
# ╠═81ed8d28-64e5-4008-ba8d-235ca0e3bc58
# ╠═5f497f29-9b49-4fda-b7d1-c1a85a4762dc
# ╠═08ef307a-e8d2-4778-9b8a-10ee40c4153b
# ╠═e45b15b5-f142-4af4-9957-4fe4a45db437
# ╠═5998dfe3-2cd3-4b65-89df-b90c9a27e945
# ╠═fb9d596f-9252-4aeb-b266-a8053397a5d0
# ╠═7549afef-22ea-4d98-9413-c2ae786bf930
# ╠═8048bf11-7785-48de-996a-f91173759cd2
