### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 3ab55f74-30af-11ec-2bfb-9bbf2581fb13
using Revise, Plots

# ╔═╡ 0c11b8b1-1734-4918-b5af-d9e7dedf58e5
using MAT, LinearAlgebra, SparseArrays, Arpack

# ╔═╡ b1e1e11f-ad3f-4995-b9c0-da5c892feff1
plot([1, 2, NaN, 4, 5], line_z = [1, 2, 3, 4, 5], arrows=true, legend=:none)

# ╔═╡ f90e4ed3-7760-46ce-bdb4-5ebfb6fc1463
function edgeplot(x_loc, y_loc, edges, edge_vals)
	xs = []
	ys = []
	zs = []
	for (i, e) in enumerate(edges)
		if edge_vals[i] > 0
			push!(xs, x_loc[e[1]], x_loc[e[2]], NaN)
			push!(ys, y_loc[e[1]], y_loc[e[2]], NaN)
		else
			push!(xs, x_loc[e[2]], x_loc[e[1]], NaN)
			push!(ys, y_loc[e[2]], y_loc[e[1]], NaN)
		end
		# push!(zs, abs(edge_vals[i]), 0, 0)
		push!(zs, edge_vals[i], 0, 0)
	end
	plot(xs, ys, arrows=true, line_z = zs, legend=false, colorbar=true, c=:viridis)
	
end

# ╔═╡ 81b0950d-d5e4-4207-baff-ae8d63dda202
# it works!
edgeplot(
	[1, 2, 4, 5],
	[1, 2, -1, -2],
	[(1,2), (3, 4), (2,4)],
	[1, 3, -4]
)

# ╔═╡ 55a00e1d-0694-4010-b900-1bbd80e1f75d
# define wts in the default orientation; increasing index -> positive edge
begin 
	vars = matread("../data/RGC60_100.mat")
	xs = vars["xyz100"][:,1]
	ys = vars["xyz100"][:,2]
	W = SparseMatrixCSC(UpperTriangular(Diagonal(vars["L100"]) - vars["L100"]))
	heads, tails, wts = findnz(W)
end

# ╔═╡ 40e4b1f0-e3f7-4d21-98bb-7a519e93ce50
# modify the orientation of the base graph: -1 is a flipped direction.

ori = ones(nnz(W))
# ori = rand([1, -1], nnz(W))

# ╔═╡ 881d252d-10ca-4e5d-a6c2-338db7d20ae4
begin
	dirs = wts .* ori
	B = sparse([heads; tails], [1:nnz(W); 1:nnz(W)], [dirs; -dirs])
	L1 = B'*B
	Λ, Φ = eigen(Matrix(L1))
end

# ╔═╡ 8381de83-c3b9-442b-9d4b-404a4f5794ff
begin
	h_ori = [
		ori[i] > 0 ? node : tails[i]
		for (i, node) in enumerate(heads)
	]
	t_ori = [
		ori[i] > 0 ? node : heads[i]
		for (i, node) in enumerate(tails)
	]
end

# ╔═╡ 3ac9ac59-d16c-4ce8-8044-e8315b729d64
edgeplot(xs, ys, zip(heads, tails), ones(nnz(W)))

# ╔═╡ dee70b76-cbea-4c2d-b83c-2f3953c536f5
edgeplot(xs, ys, zip(heads, tails), ori)

# ╔═╡ 67db72fa-1062-4285-a3ff-9e4ed1402f0c
edgeplot(xs, ys, zip(h_ori, t_ori), Φ[:,1])

# ╔═╡ 146d390e-7dcf-4af3-a74f-b341b6c7519f
edgeplot(xs, ys, zip(heads, tails), Φ[:,2] .* ori)

# ╔═╡ f83dfbf7-c366-449f-8b60-c3fa4be76bdf
edgeplot(xs, ys, zip(heads, tails), Φ[:,3] .* ori)

# ╔═╡ 89d473d6-aab0-4025-a5ff-9d2d24a12330
edgeplot(xs, ys, zip(heads, tails), Φ[:,4] .* ori)

# ╔═╡ 4cd3d18a-d1f0-4106-bc00-926532e99eeb
edgeplot(xs, ys, zip(heads, tails), Φ[:,5] .* ori)

# ╔═╡ 542a3e37-e700-4090-97dc-59119cd45eb0
edgeplot(xs, ys, zip(heads, tails), Φ[:,6] .* ori)

# ╔═╡ 0318ca3c-1371-4188-8da0-cd8cffe3c40b
edgeplot(xs, ys, zip(heads, tails), Φ[:,7] .* ori)

# ╔═╡ 38c8c7c6-eafb-4dd7-9f7f-fb828d407423
edgeplot(xs, ys, zip(heads, tails), Φ[:,8] .* ori)

# ╔═╡ 66b69308-1834-41ac-b3b0-891e602c0179
edgeplot(xs, ys, zip(heads, tails), Φ[:,9] .* ori)

# ╔═╡ e0955b56-600e-47c9-9ed9-32ff4f8daa10
edgeplot(xs, ys, zip(heads, tails), Φ[:,10] .* ori)

# ╔═╡ bb051740-abcb-490a-a00d-a0f5f363773b
edgeplot(xs, ys, zip(heads, tails), Φ[:,100] .* ori)

# ╔═╡ de81267d-37cb-4081-8161-020b2a11b892
plot(Λ, legend=:none); plot!([1141, 1141], [0, 5])

# ╔═╡ d10ed3de-528b-4250-9a31-218e0d8485eb
edgeplot(xs, ys, zip(heads, tails), Φ[:,end-1] .* ori)

# ╔═╡ d34bf923-7c34-46dc-ab38-91c8c7c97fc1
edgeplot(xs, ys, zip(heads, tails), Φ[:,end-2] .* ori)

# ╔═╡ af99f0fd-73dc-4d80-8081-357f3088c245
edgeplot(xs, ys, zip(heads, tails), Φ[:,end-3] .* ori)

# ╔═╡ Cell order:
# ╠═3ab55f74-30af-11ec-2bfb-9bbf2581fb13
# ╠═b1e1e11f-ad3f-4995-b9c0-da5c892feff1
# ╠═f90e4ed3-7760-46ce-bdb4-5ebfb6fc1463
# ╠═81b0950d-d5e4-4207-baff-ae8d63dda202
# ╠═0c11b8b1-1734-4918-b5af-d9e7dedf58e5
# ╠═55a00e1d-0694-4010-b900-1bbd80e1f75d
# ╠═40e4b1f0-e3f7-4d21-98bb-7a519e93ce50
# ╠═881d252d-10ca-4e5d-a6c2-338db7d20ae4
# ╠═8381de83-c3b9-442b-9d4b-404a4f5794ff
# ╠═3ac9ac59-d16c-4ce8-8044-e8315b729d64
# ╠═dee70b76-cbea-4c2d-b83c-2f3953c536f5
# ╠═67db72fa-1062-4285-a3ff-9e4ed1402f0c
# ╠═146d390e-7dcf-4af3-a74f-b341b6c7519f
# ╠═f83dfbf7-c366-449f-8b60-c3fa4be76bdf
# ╠═89d473d6-aab0-4025-a5ff-9d2d24a12330
# ╠═4cd3d18a-d1f0-4106-bc00-926532e99eeb
# ╠═542a3e37-e700-4090-97dc-59119cd45eb0
# ╠═0318ca3c-1371-4188-8da0-cd8cffe3c40b
# ╠═38c8c7c6-eafb-4dd7-9f7f-fb828d407423
# ╠═66b69308-1834-41ac-b3b0-891e602c0179
# ╠═e0955b56-600e-47c9-9ed9-32ff4f8daa10
# ╠═bb051740-abcb-490a-a00d-a0f5f363773b
# ╠═de81267d-37cb-4081-8161-020b2a11b892
# ╠═d10ed3de-528b-4250-9a31-218e0d8485eb
# ╠═d34bf923-7c34-46dc-ab38-91c8c7c97fc1
# ╠═af99f0fd-73dc-4d80-8081-357f3088c245
