### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 19746b80-9afb-11ee-38d4-43a96fc94583
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ fda9797e-5a17-4686-a7e9-f42133bcff0b
using Revise

# ╔═╡ 8a2ca6cf-a29c-43fd-916e-fec16b0f7557
using MultiscaleSimplexSignalTransforms

# ╔═╡ 6de9e4cf-366f-4992-8b41-1eec760466a3
using SparseArrays, LinearAlgebra, StatsBase, Arpack, Graphs, Plots, Chain, DataStructures, MAT, GraphRecipes, InvertedIndices, SimpleWeightedGraphs

# ╔═╡ 821aeaad-7b39-48f8-8d60-7909c107554d
function remove_diag!(M)
	for ij ∈ CartesianIndices(M)
		if ij.I[1] == ij.I[2]
			M[ij] = 0
		end
	end
	M
end

# ╔═╡ 2cfc36d2-eca3-436a-a1c3-a69edd63d226
dendrite_W() = @chain begin
	matread("../data/RGC60_100.mat")
	-_["L100"]
	remove_diag!
end

# ╔═╡ 2225db66-2fc5-4b77-8148-1b090592cbfc
begin
	W = dendrite_W()
	g = Graph(W)
	h = DiGraph(UpperTriangular(W))
end

# ╔═╡ 35488aaa-fb62-4066-bd43-6f81fa30bf1f
adjacency_matrix(h)

# ╔═╡ f8fa892e-9209-4785-99b4-f4e0806c288a
d = distmat(h)

# ╔═╡ d8b77bc6-8a61-4c91-8c59-65202522fe09
adjacency_matrix(h)[1:3,1:3]

# ╔═╡ 27ae267a-4c6d-4fad-812c-7dc85ee49238
outneighbors(h, 1)

# ╔═╡ Cell order:
# ╠═19746b80-9afb-11ee-38d4-43a96fc94583
# ╠═fda9797e-5a17-4686-a7e9-f42133bcff0b
# ╠═8a2ca6cf-a29c-43fd-916e-fec16b0f7557
# ╠═6de9e4cf-366f-4992-8b41-1eec760466a3
# ╠═821aeaad-7b39-48f8-8d60-7909c107554d
# ╠═2cfc36d2-eca3-436a-a1c3-a69edd63d226
# ╠═2225db66-2fc5-4b77-8148-1b090592cbfc
# ╠═35488aaa-fb62-4066-bd43-6f81fa30bf1f
# ╠═f8fa892e-9209-4785-99b4-f4e0806c288a
# ╠═d8b77bc6-8a61-4c91-8c59-65202522fe09
# ╠═27ae267a-4c6d-4fad-812c-7dc85ee49238
