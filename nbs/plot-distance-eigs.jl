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

# ╔═╡ 31246c91-5d30-429f-935f-0ed069c4dbd9
using Unzip

# ╔═╡ 45bdfde0-9677-4cc6-86ec-ee1f52a793da
using Glob

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

# ╔═╡ d8830422-03b9-4358-8865-b86b80a8f76a
randpath = cumsum(rand(100))

# ╔═╡ b27f62ed-851b-46d5-bc51-0603f97b4777
cyclecoor = Dict(
	1 => cos,
	2 => sin,
	3 => ((_) -> 0)
)

# ╔═╡ f11ffc25-b0ed-4dc6-b810-c16a18363462
g_data = Dict(
	:dendrite => (
		W=remove_diag!(-dendrite_data["L100"]),
		coor = dendrite_data["xyz100"]
	),
	:path => (
		W=adjacency_matrix(path_graph(100)),
		coor = [
			i
			for i=1:100, j=1:3
		]
	),
	:wt_path => (
		W=adjacency_matrix(path_graph(100)),
		coor = [
			randpath[i]
			for i=1:100, j=1:3
		]
	),
	:cycle => (
		W = Symmetric(adjacency_matrix(cycle_graph(100))),
		coor = [
			cyclecoor[j](2π*i/100)
			for i=0:99, j=1:3
		]
	)
)

# ╔═╡ d7ab7b8e-8a6b-4778-854d-6258cbec5707
which_g = :dendrite

# ╔═╡ 6810e7a3-d52b-474c-8ea5-2280b93177af
use_weights = false

# ╔═╡ e5d21394-56b7-46ad-a3d9-ef7f649bf323
use_dmax_on_diag = true

# ╔═╡ c04c85f9-727a-4032-99f3-9cc3ecb59be5
dir_mod = true

# ╔═╡ eed4eca1-dab9-4675-b2a8-480404461153
xs, ys = @chain begin
	g_data[which_g].coor
	(_[:,1], _[:,2])
end

# ╔═╡ b7bd92c5-22a3-4b45-82e8-a4a4a3175bf7
xyz_dists = Distances.pairwise(Euclidean(), g_data[which_g].coor, dims=1)

# ╔═╡ 2225db66-2fc5-4b77-8148-1b090592cbfc
begin
	W = g_data[which_g].W
	g = SimpleWeightedGraph(
		use_weights ? W .* xyz_dists : W
	)
	h_mat = use_weights ? UpperTriangular(W) .* xyz_dists : UpperTriangular(W)
	if which_g == :cycle
		h_mat = Matrix(h_mat)
		h_mat[1, end] = 0
		h_mat[end, 1] = 1
		if dir_mod
			for i=2:size(h_mat,1)÷2+1
				h_mat[i,i+1] = 0
				h_mat[i,i-1] = 1
			end
		end
	end

	if which_g == :dendrite
		if dir_mod
			h_mat = Matrix(h_mat)
			h_mat[11:97, 11:97] = h_mat[11:97,11:97]'
		end
	end
	h = SimpleWeightedDiGraph(h_mat)
	dmax = 10000
end

# ╔═╡ a38a91a4-d806-476e-8edc-e77906e411f0
g3d_edges = @chain begin
	W
	Graph
	edges
	[(e.src, e.dst) for e in _]
	[ [(g_data[which_g].coor[i,1], g_data[which_g].coor[i,2], g_data[which_g].coor[i,3]) for i in ei] for ei in _]
	[
		pt
		for pts in _
		for pt in [pts; (NaN, NaN, NaN)]
	]
	unzip
end

# ╔═╡ 6c3e4bf5-9f76-4c77-865f-6bd8a623a0dc
begin
	plot(g3d_edges...)
	scatter!(eachcol(g_data[which_g].coor)..., ms=1)
end

# ╔═╡ 150a5f2a-5b79-4a9f-a09d-6a4da54c7fd4
# begin
# 	g_rec = SimpleWeightedGraph(
# 		use_weights ? W .* remove_diag!(1 ./ xyz_dists) : W
# 	)
# 	h_rec = SimpleWeightedDiGraph(
# 		use_weights ? UpperTriangular(W) .*  remove_diag!(1 ./ xyz_dists) : UpperTriangular(W)
# 	)
# end

# ╔═╡ 840eddb8-0182-4ccf-b1c6-5f035608bba1
# plot_sv(sparsevec(12:95, 1, nv(g)))

# ╔═╡ cb47d3a4-2c21-4a30-a05a-d5ea0a0660f0
# findall(degree(g).==3)

# ╔═╡ 357bfc0a-2583-4c3d-9657-f886f2e6e1ef
# findnz(sparse(W)[11,:])

# ╔═╡ 7fee970b-1f97-4470-8e51-f9100575e4e7
Dg = Matrix(distmat(g))

# ╔═╡ e14bfac7-3940-4f31-bcee-d61c32b75693
# Kg = remove_diag!(1 ./ Matrix(distmat(g_rec))) + (use_dmax_on_diag ? dmax*I : 0*I)
Kg = remove_diag!(1 ./ Dg) + (use_dmax_on_diag ? dmax*I : 0*I)

# ╔═╡ 282e0dc3-30b2-4ce4-9fc0-63e6ff5914b8
# histogram(Kg[Kg .< dmax], yscale=:log10, legend=false)

# ╔═╡ ce9af7c4-8e56-47ef-9968-295d700de5f3
begin
	Dh = remove_diag!(fill(Float64(dmax), (nv(h), nv(h))))
	jj, ii, vv = findnz(distmat(h)')
	for (i,j,v) in zip(ii, jj, vv)
		Dh[i,j] = v
	end

	# jj, ii, vv = findnz(distmat(h_rec)')
	Kh = Matrix(sparse(ii, jj, 1 ./ vv, size(Dh)...)) + (use_dmax_on_diag ? dmax*I : 0*I)
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
	:sym_d_orig => Symmetric(Dg),
	:sym_d_proj => Symmetric(P(Dg)),
	:sym_k_orig => Symmetric(Kg),
	:sym_k_proj => Symmetric(P(Kg)),
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

# ╔═╡ ea14be67-4a57-45aa-a847-0365047bcf23
g_edges = @chain begin
	W
	Graph
	edges
	[(e.src, e.dst) for e in _]
	[ [(xs[i], ys[i]) for i in ei] for ei in _]
	@. Shape
end

# ╔═╡ f84f60ef-4b6d-4b50-a7e9-88b802e613cb
function uscale(v, thresh=1e-12)
	m, M = extrema(v)
	if M-m < thresh
		ones(size(v))
	else
		(v .- m) ./ (M-m)
	end
end

# ╔═╡ f237ffca-ea30-40b2-918a-bb6914b13b8f
cscale(v) = maximum(abs.(extrema(v))) |> M -> (-1,1) .* M
# (v .+ M) ./ 2M

# ╔═╡ 939defff-376a-48be-af29-e162a2d38ba1
function plot_sv(v)
	plot(g_edges; legend=false, axis=([], false), aspectratio=:equal, size=(400, 400))
	scatter!(xs, ys, mz=v.*sign(v[1]), msw=0, ms=4uscale(abs.(v)), c=:viridis, clims=cscale(v))
	# scatter!(xs, ys, mz=v, msw=0, ms=1, c=:seismic)
end

# ╔═╡ 4f408f07-1a71-4658-bc24-71c989988529
Knaive = svd(P(remove_diag!(1 ./xyz_dists)));

# ╔═╡ 4ae6a7c6-c910-40b5-ac07-1dad65180bff
plot_sv(sparsevec(11:97, 1, nv(g)))

# ╔═╡ cd0aff95-619f-4c71-b71a-4e52fff08b4b
1 |> n -> plot(
	# plot([
	# 	plot_sv(decomps[:dir_k_orig].U[:,end-i+1])
	# 	for i in (1:n) .+ 0
	# ]..., layout=(n,1)),
	plot([
		plot_sv(decomps[:dir_k_orig].Vt'[:,end-i+1])
		for i in (1:n) .+ 0
	]..., layout=(n,1)),
	# plot([
	# 	plot_sv(decomps[:dir_k_proj].U[:,end-i+1])
	# 	for i in (1:n) .+ 0
	# ]..., layout=(n,1)),
	plot([
		plot_sv(decomps[:dir_k_proj].Vt'[:,end-i+1])
		for i in (1:n) .+ 0
	]..., layout=(n,1)),
	size=(800,400n),
	layout=(1,2)
)

# ╔═╡ 434fd627-58e2-4fa0-b444-5f4cfd9bc709
# plot(decomps[:dir_d_orig].Vt'[:,2])

# ╔═╡ 6e537699-2004-4388-af79-3c9ef5141a67
# plot(
# 	[
# 		# plot_sv(decomp.Vt'[:,endswith(string(k), "orig") ? 1 : end])
# 		plot_sv(decomp.Vt'[:,1])
# 		for (k, decomp) in decomps
# 		if startswith(string(k), "dir")
# 	]...
# )

# ╔═╡ b2ffd9bf-ff2e-4995-8bea-6270d8c5d8f3
function save_sv_plot(
	n;
	is_directed, is_reciprocal, is_projected, from_end, use_left_sv
)
	dirsymbol = is_directed ? :dir : :sym
	recsymbol = is_reciprocal ? :k : :d
	projsymbol = is_projected ? :proj : :orig
	handle = Symbol(dirsymbol, :_, recsymbol, :_, projsymbol)

	svsymbol = use_left_sv ? :U : :V

	ind = from_end ? size(W, 1)-n+1 : n
	fn = "../output/$(handle)_$(svsymbol)_$ind.png"

	signal = (
		use_left_sv ? decomps[handle].U : decomps[handle].Vt'
	)[:, ind]

	fig = plot_sv(signal)
	savefig(fig, fn)
	# fig
end

# ╔═╡ 508853b8-0e30-409d-ad5b-621be6bcd6e0
# save_splot(1; is_directed=true, is_reciprocal=false, is_projected=false, use_left_sv=true, from_end=true)

# ╔═╡ c15f01cf-5d5d-4385-b42d-9abe94d81db3
function save_all_sv_plots(nfirst, nlast)
	for bools in Iterators.product(fill([true, false], 5)...)
		is_directed, is_reciprocal, is_projected, from_end, use_left_sv = bools
		n = from_end ? nlast : nfirst
		for ni in 1:n
			save_sv_plot(
				ni; is_directed, is_reciprocal, is_projected, from_end, use_left_sv
			)
		end
	end
end

# ╔═╡ d3b0918c-b256-4ff9-bea3-068d2c64c3a7
mpath = rsplit(pathof(MultiscaleSimplexSignalTransforms), '/', limit=3)[1]

# ╔═╡ eee7238c-1daf-4b07-9fc1-85fea61360cb
function get_compressed_plots(name="svplots")
	# first clear the output directory
	@chain begin
		glob("*", "$mpath/output")
		[run(`rm $fi`).exitcode for fi in _]
	end

	# then refill the output directory
	save_all_sv_plots(8, 16)

	# then form the archive
	run(`tar -czvf $mpath/$name.tgz $mpath/output`)
end

# ╔═╡ 21b3f93b-f27b-4785-ad00-c91372fa49cb
# get_compressed_plots("svplots_viridis_switchbranch")

# ╔═╡ Cell order:
# ╠═19746b80-9afb-11ee-38d4-43a96fc94583
# ╠═fda9797e-5a17-4686-a7e9-f42133bcff0b
# ╠═8a2ca6cf-a29c-43fd-916e-fec16b0f7557
# ╠═6de9e4cf-366f-4992-8b41-1eec760466a3
# ╠═821aeaad-7b39-48f8-8d60-7909c107554d
# ╠═619417e1-fe3e-4e6e-90c4-6d2f240d6f12
# ╠═d8830422-03b9-4358-8865-b86b80a8f76a
# ╠═b27f62ed-851b-46d5-bc51-0603f97b4777
# ╠═f11ffc25-b0ed-4dc6-b810-c16a18363462
# ╠═d7ab7b8e-8a6b-4778-854d-6258cbec5707
# ╠═6810e7a3-d52b-474c-8ea5-2280b93177af
# ╠═e5d21394-56b7-46ad-a3d9-ef7f649bf323
# ╠═c04c85f9-727a-4032-99f3-9cc3ecb59be5
# ╠═eed4eca1-dab9-4675-b2a8-480404461153
# ╠═31246c91-5d30-429f-935f-0ed069c4dbd9
# ╠═a38a91a4-d806-476e-8edc-e77906e411f0
# ╠═6c3e4bf5-9f76-4c77-865f-6bd8a623a0dc
# ╠═b7bd92c5-22a3-4b45-82e8-a4a4a3175bf7
# ╠═2225db66-2fc5-4b77-8148-1b090592cbfc
# ╠═150a5f2a-5b79-4a9f-a09d-6a4da54c7fd4
# ╠═840eddb8-0182-4ccf-b1c6-5f035608bba1
# ╠═cb47d3a4-2c21-4a30-a05a-d5ea0a0660f0
# ╠═357bfc0a-2583-4c3d-9657-f886f2e6e1ef
# ╠═7fee970b-1f97-4470-8e51-f9100575e4e7
# ╠═e14bfac7-3940-4f31-bcee-d61c32b75693
# ╠═282e0dc3-30b2-4ce4-9fc0-63e6ff5914b8
# ╠═ce9af7c4-8e56-47ef-9968-295d700de5f3
# ╠═152d2cdf-de1f-4ec2-9933-af4210c48d85
# ╠═8c7bf95a-0f78-474c-995b-e2302cf8e54c
# ╠═4fdfee57-3353-4fd2-8953-4d00ce9dea72
# ╠═ea14be67-4a57-45aa-a847-0365047bcf23
# ╠═f84f60ef-4b6d-4b50-a7e9-88b802e613cb
# ╠═f237ffca-ea30-40b2-918a-bb6914b13b8f
# ╠═939defff-376a-48be-af29-e162a2d38ba1
# ╠═4f408f07-1a71-4658-bc24-71c989988529
# ╠═4ae6a7c6-c910-40b5-ac07-1dad65180bff
# ╠═cd0aff95-619f-4c71-b71a-4e52fff08b4b
# ╠═434fd627-58e2-4fa0-b444-5f4cfd9bc709
# ╠═6e537699-2004-4388-af79-3c9ef5141a67
# ╠═b2ffd9bf-ff2e-4995-8bea-6270d8c5d8f3
# ╠═508853b8-0e30-409d-ad5b-621be6bcd6e0
# ╠═c15f01cf-5d5d-4385-b42d-9abe94d81db3
# ╠═d3b0918c-b256-4ff9-bea3-068d2c64c3a7
# ╠═45bdfde0-9677-4cc6-86ec-ee1f52a793da
# ╠═eee7238c-1daf-4b07-9fc1-85fea61360cb
# ╠═21b3f93b-f27b-4785-ad00-c91372fa49cb
