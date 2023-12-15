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
using SparseArrays, LinearAlgebra, StatsBase, Arpack, Graphs, Plots, Chain, DataStructures, MAT, GraphRecipes, InvertedIndices, SimpleWeightedGraphs, Distances

# ╔═╡ 821aeaad-7b39-48f8-8d60-7909c107554d
function remove_diag!(M)
	for ij ∈ CartesianIndices(M)
		if ij.I[1] == ij.I[2]
			M[ij] = 0
		end
	end
	M
end

# ╔═╡ 619417e1-fe3e-4e6e-90c4-6d2f240d6f12
dendrite_data = matread("../data/RGC60_100.mat")

# ╔═╡ eed4eca1-dab9-4675-b2a8-480404461153
xs, ys = @chain begin
	dendrite_data["xyz100"]
	(_[:,1], _[:,2])
end

# ╔═╡ b7bd92c5-22a3-4b45-82e8-a4a4a3175bf7
xyz_dists = Distances.pairwise(Euclidean(), dendrite_data["xyz100"], dims=1)

# ╔═╡ 2225db66-2fc5-4b77-8148-1b090592cbfc
begin
	W = -remove_diag!(dendrite_data["L100"])
	g = SimpleWeightedGraph(W .* xyz_dists)
	h = SimpleWeightedDiGraph(UpperTriangular(W) .* xyz_dists)
	dmax = 10000
	# Dh = distmat(h)
	# Dg, Dh = distmat(g), distmat(h)
end

# ╔═╡ 7fee970b-1f97-4470-8e51-f9100575e4e7
Dg = Matrix(distmat(g))

# ╔═╡ e14bfac7-3940-4f31-bcee-d61c32b75693
Kg = remove_diag!(1 ./ Dg)

# ╔═╡ ce9af7c4-8e56-47ef-9968-295d700de5f3
begin
	Dh = remove_diag!(fill(Float64(dmax), (nv(h), nv(h))))
	jj, ii, vv = findnz(distmat(h)')
	for (i,j,v) in zip(ii, jj, vv)
		Dh[i,j] = v
	end
	
	Kh = Matrix(sparse(ii, jj, 1 ./ vv, size(Dh)...))
	(Dh, Kh)
end

# ╔═╡ 152d2cdf-de1f-4ec2-9933-af4210c48d85
begin
	Pn(n) = I - ones(n,n)/n
	P(M::AbstractMatrix) = Pn(size(M, 1)) |> PP -> PP*M*PP
end

# ╔═╡ 8c7bf95a-0f78-474c-995b-e2302cf8e54c
# our eight matrices
# is it directed? (sym/dir)
# is it original distances or their reciprocals? (d/k)
# is the matrix projected against? (orig/proj)
matrices = Dict(
	:sym_d_orig => Dg,
	:sym_d_proj => P(Dg),
	:sym_k_orig => Kg,
	:sym_k_proj => P(Kg),
	:dir_d_orig => Dh,
	:dir_d_proj => P(Dh),
	:dir_k_orig => Kh,
	:dir_k_proj => P(Kh)
)

# ╔═╡ 4fdfee57-3353-4fd2-8953-4d00ce9dea72
decomps = Dict(
	k => svd(M)
	for (k,M) in matrices
)

# ╔═╡ 71fd7e5c-04b2-483a-b53b-fb807e856c2e
pos = Pos(xs, ys)

# ╔═╡ 3a78d48d-4484-4230-84a8-84090b3f85a6
zreg = ZeroRegion(W)

# ╔═╡ bde81468-8b03-4f75-8f2c-ce0a791683ee
splot(zreg, decomps[:sym_d_proj].U[:,1], pos)

# ╔═╡ b2ffd9bf-ff2e-4995-8bea-6270d8c5d8f3
function save_splot(
	n;
	is_directed, is_reciprocal, is_projected, from_end, use_left_sv
)
	dirsymbol = is_directed ? :dir : :sym
	recsymbol = is_reciprocal ? :k : :d
	projsymbol = is_projected ? :proj : :orig
	handle = Symbol(dirsymbol, :_, recsymbol, :_, projsymbol)
	fn = "../output/$handle.png"

	source = (
		use_left_sv ? decomps[handle].U : decomps[handle].Vt'
	)

	signal = from_end ? source[:, end-n+1] : source[:, n]

	# savefig(splot(zreg, signal, pos), fn)
	splot(zreg, signal, pos)
end

# ╔═╡ 508853b8-0e30-409d-ad5b-621be6bcd6e0
save_splot(4; is_directed=false, is_reciprocal=false, is_projected=true, use_left_sv=true, from_end=false)

# ╔═╡ Cell order:
# ╠═19746b80-9afb-11ee-38d4-43a96fc94583
# ╠═fda9797e-5a17-4686-a7e9-f42133bcff0b
# ╠═8a2ca6cf-a29c-43fd-916e-fec16b0f7557
# ╠═6de9e4cf-366f-4992-8b41-1eec760466a3
# ╠═821aeaad-7b39-48f8-8d60-7909c107554d
# ╠═619417e1-fe3e-4e6e-90c4-6d2f240d6f12
# ╠═eed4eca1-dab9-4675-b2a8-480404461153
# ╠═b7bd92c5-22a3-4b45-82e8-a4a4a3175bf7
# ╠═2225db66-2fc5-4b77-8148-1b090592cbfc
# ╠═7fee970b-1f97-4470-8e51-f9100575e4e7
# ╠═e14bfac7-3940-4f31-bcee-d61c32b75693
# ╠═ce9af7c4-8e56-47ef-9968-295d700de5f3
# ╠═152d2cdf-de1f-4ec2-9933-af4210c48d85
# ╠═8c7bf95a-0f78-474c-995b-e2302cf8e54c
# ╠═4fdfee57-3353-4fd2-8953-4d00ce9dea72
# ╠═71fd7e5c-04b2-483a-b53b-fb807e856c2e
# ╠═3a78d48d-4484-4230-84a8-84090b3f85a6
# ╠═bde81468-8b03-4f75-8f2c-ce0a791683ee
# ╠═b2ffd9bf-ff2e-4995-8bea-6270d8c5d8f3
# ╠═508853b8-0e30-409d-ad5b-621be6bcd6e0
