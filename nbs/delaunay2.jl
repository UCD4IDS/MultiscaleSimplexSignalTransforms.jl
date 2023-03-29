### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ e9b65c32-c95c-11ed-04b4-f11a36d9f557
using Pkg; Pkg.activate(".")

# ╔═╡ cb105fb4-c6bc-4b2c-980b-c04a3e36ac46
using Revise

# ╔═╡ 5f0b0cc3-77b7-4d2f-a8ca-e9a320a379a7
using MultiscaleSimplexSignalTransforms

# ╔═╡ d1032005-6fd2-4b98-bbdc-996d8d1d4a6a
using VoronoiDelaunay, Images

# ╔═╡ 55ec2d35-0e77-444f-9d25-fe82b6811fa7
using SparseArrays, LinearAlgebra, StatsBase, Arpack, Graphs, Plots, Chain, DataStructures, MAT, GraphRecipes, Unzip, Match

# ╔═╡ e9754135-cbb1-4515-9665-9787dee1a74e
using ProfileView

# ╔═╡ caa8621f-21d7-492d-b0ee-0581ff0c81db
import VoronoiDelaunay: from_image, voronoiedges, getplotxy

# ╔═╡ ed016edc-a35b-4b1f-b8b5-f240c6f2430b
import Images: load

# ╔═╡ db0a7566-c177-43fa-a170-692a49576ee2
whichcol = colorimg -> @chain begin
	mean(channelview(colorimg); dims=1)
	dropdims(; dims=1)
	# @. Gray
end

# ╔═╡ b24ad7be-8ba1-4b05-86bc-c4c9fd25cb14
colorimg = load("/Users/eug/Desktop/woman.tiff")

# ╔═╡ 59b83306-64c4-4184-8bad-fa710508ca89
img = whichcol(colorimg)

# ╔═╡ 9469523b-5b50-49d9-b669-03b877cb1c3a
@chain begin
	colorview(Gray, whichcol(colorimg))
	# @aside save("/Users/eug/Desktop/woman-gray.png", _)
end

# ╔═╡ e5d523a1-d65c-4f0f-83e4-594e2bb45c49
n = size(img, 1)

# ╔═╡ 241c62ec-8229-4052-a518-23d52e6ce7ab
npt = 2000

# ╔═╡ bd1603de-2e82-4227-a920-dacddc7725e6
ps = rand(2,npt)

# ╔═╡ f004c8fa-3be3-42f8-952f-79189dcfd59b
pos = Pos(ps[1,:], ps[2,:])

# ╔═╡ 3273f5fd-b4a9-4c41-8738-278aeed37770
function maketess(ps)
	tess = DelaunayTessellation(size(ps, 1))
	push!(tess, Point2D.(collect.(eachrow(ps))...))
	tess
end

# ╔═╡ a902e76f-e071-4ed4-bca1-faf30977e11c
tess = maketess(ps.+1)

# ╔═╡ 6185f91e-4505-435e-8e74-0696cf22073f
vertex_vals = [
	img[n+1-is[2], is[1]]
	for is ∈ eachcol(ceil.(Int, ps .* n))
]

# ╔═╡ 9a643ee0-0d03-463c-8444-3dc116d2125b
tess

# ╔═╡ 79b61df7-9058-4987-8de4-786422de3d5e
@chain begin
scatter(ps[1,:], ps[2,:], mz=vertex_vals, color=:greys, aspect_ratio=:equal, grid=:none, ms=3, msw=0, legend=false, showaxis=false)
# @aside savefig("/Users/eug/Desktop/woman-signal0.png")
end

# ╔═╡ a14e5dbd-7f5a-48ee-9826-39219a8dff35
plot(
	broadcast.(-, getplotxy(delaunayedges(tess)), 1),
	legend=false, aspect_ratio=:equal
)

# ╔═╡ 05f948db-41a3-47e0-9437-34d6b913b7fc
function MultiscaleSimplexSignalTransforms.SimplexTree(tess, ps)
	tree = SimplexTree(size(ps, 2))
	for tri in tess
		insert!(tree, Tuple(
			findfirst(i -> x._x == 1+ps[1,i] && x._y == 1+ps[2,i], 1:size(ps, 2))
			for x ∈ (geta(tri), getb(tri), getc(tri))			
		))
	end
	tree
end

# ╔═╡ b48d78bb-aaa4-45fe-85ee-0c6e7749935b
a = SimplexTree(tess, ps)

# ╔═╡ 0237c4a5-cfcd-4b96-aaa9-a34b1bf2b996
length(simplices(a, 2))

# ╔═╡ 7ede916f-2f62-4d83-a9f4-b28a1c4ef5b2
edge_vals = [mean(vertex_vals[v] for v ∈ e) for e ∈ simplices(a, 1) ]

# ╔═╡ 2bfc5b15-e8d6-4228-8321-d0f24cd1ac66
@chain begin
	splot(KRegion(a, 1), edge_vals, pos; k=1, palette=:greys, size=4)
	# @aside savefig("/Users/eug/Desktop/woman-signal1.png")
end

# ╔═╡ 659dced8-7aea-49da-981f-504aa9dcb643
tri_vals = [mean(vertex_vals[v] for v ∈ t) for t ∈ simplices(a, 2) ]

# ╔═╡ e042f078-234d-4597-9f55-9c5df90ec9e8
@chain begin
splot(
	KRegion(a, 2; weakadjweight=1.0), 
	tri_vals, pos; k=2, palette=:greys, notribds=true)
# @aside savefig("/Users/eug/Desktop/woman-signal2.png")
end

# ╔═╡ 8f55a67e-dbd5-4f45-9858-dc95721a28df
keepfracs = [.01, .05, .1, .25, .5, .75, .9]

# ╔═╡ 1170ebfc-f6de-4f22-9f57-ba82ec2c69a3
function trybases(k)
	reg = KRegion(a, k; weakadjweight=(k==2 ? 1.0 : 0.0))
	ghwt = kGHWT(reg; normalization=Normalization.Combinatorial)
	hglet = kHGLET(reg; normalization=Normalization.Combinatorial)
	(reg, ghwt, hglet)

end

# ╔═╡ 141b5a2a-ed77-429f-8160-0ade48ede451
k=1

# ╔═╡ 60bdcc64-e778-4fb9-ba02-6fb97a294e2c
R, G, H = trybases(k)

# ╔═╡ 51633743-5795-493e-b09a-a3ef4c5da698
gcoef = analyze(G, k==1 ? edge_vals : tri_vals)

# ╔═╡ 75a88580-5800-4450-a90a-f5535fc1f05c
hcoef = analyze(H, k==1 ? edge_vals : tri_vals)

# ╔═╡ f1c6dbf7-ce93-4719-b699-c6380c8cad56
n_simp(R)

# ╔═╡ 4a93684f-fc3a-4b0e-bb19-e0857dc5abe4
n_bd(R)

# ╔═╡ 8798bf3a-12a3-4526-9659-506d68ef55a8
[G[0,0,0] reduce(hcat, [
	G[j-1, k-1, 1] 
	for (j,o) ∈ enumerate(G.inds.offsets) 
	for k ∈ 1:length(o)
	if o[k] < n_simp(R)-1 && (k==length(o) || o[k+1] > o[k]+1)
])]

# ╔═╡ 3711f89d-7249-48db-a70e-d60955e1a10c
length(tri_vals)

# ╔═╡ bfd9054a-7f46-4e10-a7de-cf926d5f9cb2
# gbb = reduce(hcat, [
# 	G[j,k,:]
# 	for (j,k) ∈ best_basis(G, k==1 ? edge_vals : tri_vals).locs
# ])

# ╔═╡ e19dd2b3-189f-473d-839c-5b2891da41ee
function loc_repr(part, bb)
	basis = fill(0.0, partlen(root(part)), 0)
	coefs = Float64[]
	dfs(root(part)) do tsn
		for (tag, fval) ∈ get(bb, tsn.node.sinds, Dict())
			basis = hcat(basis, val(tsn.node.fvals[tag]))
			push!(coefs, val(fval))
		end
	end
	(basis, coefs)
end

# ╔═╡ 2fd2e5bf-bd4c-4680-ad21-7bf45f804539
gbb, gbbcoef = loc_repr(G, alt_best_basis(G, k==1 ? edge_vals : tri_vals; p=0.5))

# ╔═╡ 2768dc6d-2c48-44cb-969a-fe2b16e68640
hbb, hbbcoef = loc_repr(H, alt_best_basis(H, k==1 ? edge_vals : tri_vals; p=0.01))

# ╔═╡ 4c8b8a92-c1f8-4eee-bb03-12b89c5e94ae
getbasis(name) = @match name begin
	:fourier => (H[0,0,:], hcoef[0,0,:])
	:delta => (
		G[lastindex(G.inds.data, 3)-1, :],
		gcoef[lastindex(G.inds.data, 3)-1, :]
	)
	:walsh => (G[0,0,:], gcoef[0,0,:])
	:haar => ([G[0,0,0] reduce(hcat, [
		G[j-1, k-1, 1] 
		for (j,o) ∈ enumerate(G.inds.offsets) 
		for k ∈ 1:length(o)
		if o[k] < n_simp(R)-1 && (k==length(o) || o[k+1] > o[k]+1)
	])], [gcoef[0,0,0]; [
		gcoef[j-1, k-1, 1] 
		for (j,o) ∈ enumerate(G.inds.offsets) 
		for k ∈ 1:length(o)
		if o[k] < n_simp(R)-1 && (k==length(o) || o[k+1] > o[k]+1)
	]])
	:bbghwt => (gbb, gbbcoef)
	:bbhglet => (hbb, hbbcoef)
end

# ╔═╡ 2ea95123-51b5-4a7d-a6aa-36ac2a72d2f2
getbasis(:fourier)

# ╔═╡ a7e59af0-1351-40d6-ba8d-5dbf01f74ad6
function approxerror()
	names = [:delta, :walsh, :haar, :fourier, :bbghwt, :bbhglet]
	[
		@chain (B, cs) begin
			@aside vp  = sortperm(_[2], by=abs, rev=true)
			cumsum(_[1]*Diagonal(_[2])[:, vp]; dims=2) .- (k==1 ? edge_vals : tri_vals)
			mapslices(norm, _; dims=1)[:]
		end
		for (B, cs) ∈ getbasis.(names)
	]
	# vp = sortperm(gcoef[0,0,:], by=abs, rev=true)
	# ip = invperm(vp)
	# [
	# 	G[0,0,:]*[gcoef[0,0,:][vp][1:i]; zeros(n_simp(R)-i)][ip] - tri_vals |> norm
	# 	for i ∈ 1:n_simp(R)
	# ]
end

# ╔═╡ c85fb558-6242-4660-978e-37983b634500
M = approxerror()

# ╔═╡ 21b94e98-2c56-4c9d-89c4-45883de8cafe
MM = reduce(hcat, M) |> m -> m ./ m[1,:]'

# ╔═╡ 9a84f2ad-9dbc-4f79-baf9-724ca253533f
@chain begin
	plot([
		plot(
			MM;
			labels=["delta" "walsh" "haar" "fourier" "κGHWT best basis" "κHGLET best basis"], 
			lw=3, title="κ=$k Approximation Error", xlabel="number of terms", ylabel="L₂ Approximation Error"
		),
		plot(
			log.(MM[1:(size(MM,1)÷2),:]);
			labels=["delta" "walsh" "haar" "fourier" "κGHWT best basis" "κHGLET best basis"], 
			lw=3, title="κ=$k Log Approximation Error", xlabel="number of terms", ylabel="Log L₂ Approximation Error"
		)
	]..., layout=(1,2), size=(1600,600), leftmargin=10Plots.mm)
	# @aside savefig("/Users/eug/Desktop/woman-approx-k$k.png")
end

# ╔═╡ 9b51e143-b432-4783-95fe-b3110ab0d4da
(function()
	x = Vector(sprand(n_simp(R), 0.8))
	bbx = alt_best_basis(G, x)
	@info length(bbx)
	bx, cx = loc_repr(G, bbx)
	norm(bx*cx-x)
end)()

# ╔═╡ 3f011a6d-131e-4c70-b3ff-1a2c2df6c857
function approxsigs()
	names = [:delta, :walsh, :haar, :bbghwt, :fourier, :bbhglet]
	[
		@chain (B, cs) begin
			@aside vp  = sortperm(_[2], by=abs, rev=true)
			@aside ip = invperm(vp)
			[
				_[1]*[_[2][vp][1:i]; zeros(n_simp(R)-i)][ip]
				for i ∈ round.(Int, keepfracs*n_simp(R))
			]
		end
		for (B, cs) ∈ getbasis.(names)
	]
end

# ╔═╡ c7c40326-5723-49e4-bf94-99c7d4ec5654
SS = approxsigs()

# ╔═╡ 58574dec-845a-4a4f-9ec7-9003e2fcf089
@chain begin
	plot([
		splot(R, v, pos; k, palette=:greys, notribds=true) #, margin=-4Plots.mm)
		for S in SS
		for v in S
	# ]..., layout=(4,7), margin=-4Plots.mm, size=(7*600,4*400))
	]..., layout=(6,7), size=(7*400,6*400))
	@aside savefig("/Users/eug/Desktop/woman-tris.png")
end

# ╔═╡ Cell order:
# ╠═e9b65c32-c95c-11ed-04b4-f11a36d9f557
# ╠═cb105fb4-c6bc-4b2c-980b-c04a3e36ac46
# ╠═5f0b0cc3-77b7-4d2f-a8ca-e9a320a379a7
# ╠═d1032005-6fd2-4b98-bbdc-996d8d1d4a6a
# ╠═caa8621f-21d7-492d-b0ee-0581ff0c81db
# ╠═ed016edc-a35b-4b1f-b8b5-f240c6f2430b
# ╠═55ec2d35-0e77-444f-9d25-fe82b6811fa7
# ╠═e9754135-cbb1-4515-9665-9787dee1a74e
# ╠═db0a7566-c177-43fa-a170-692a49576ee2
# ╠═59b83306-64c4-4184-8bad-fa710508ca89
# ╠═b24ad7be-8ba1-4b05-86bc-c4c9fd25cb14
# ╠═9469523b-5b50-49d9-b669-03b877cb1c3a
# ╠═e042f078-234d-4597-9f55-9c5df90ec9e8
# ╠═2bfc5b15-e8d6-4228-8321-d0f24cd1ac66
# ╠═0237c4a5-cfcd-4b96-aaa9-a34b1bf2b996
# ╠═e5d523a1-d65c-4f0f-83e4-594e2bb45c49
# ╠═241c62ec-8229-4052-a518-23d52e6ce7ab
# ╠═bd1603de-2e82-4227-a920-dacddc7725e6
# ╠═f004c8fa-3be3-42f8-952f-79189dcfd59b
# ╠═3273f5fd-b4a9-4c41-8738-278aeed37770
# ╠═a902e76f-e071-4ed4-bca1-faf30977e11c
# ╠═6185f91e-4505-435e-8e74-0696cf22073f
# ╠═9a643ee0-0d03-463c-8444-3dc116d2125b
# ╠═79b61df7-9058-4987-8de4-786422de3d5e
# ╠═a14e5dbd-7f5a-48ee-9826-39219a8dff35
# ╠═05f948db-41a3-47e0-9437-34d6b913b7fc
# ╠═b48d78bb-aaa4-45fe-85ee-0c6e7749935b
# ╠═7ede916f-2f62-4d83-a9f4-b28a1c4ef5b2
# ╠═659dced8-7aea-49da-981f-504aa9dcb643
# ╠═8f55a67e-dbd5-4f45-9858-dc95721a28df
# ╠═1170ebfc-f6de-4f22-9f57-ba82ec2c69a3
# ╠═141b5a2a-ed77-429f-8160-0ade48ede451
# ╠═60bdcc64-e778-4fb9-ba02-6fb97a294e2c
# ╠═51633743-5795-493e-b09a-a3ef4c5da698
# ╠═75a88580-5800-4450-a90a-f5535fc1f05c
# ╠═f1c6dbf7-ce93-4719-b699-c6380c8cad56
# ╠═4a93684f-fc3a-4b0e-bb19-e0857dc5abe4
# ╠═8798bf3a-12a3-4526-9659-506d68ef55a8
# ╠═2ea95123-51b5-4a7d-a6aa-36ac2a72d2f2
# ╠═4c8b8a92-c1f8-4eee-bb03-12b89c5e94ae
# ╠═a7e59af0-1351-40d6-ba8d-5dbf01f74ad6
# ╠═c85fb558-6242-4660-978e-37983b634500
# ╠═21b94e98-2c56-4c9d-89c4-45883de8cafe
# ╠═9a84f2ad-9dbc-4f79-baf9-724ca253533f
# ╠═3711f89d-7249-48db-a70e-d60955e1a10c
# ╠═bfd9054a-7f46-4e10-a7de-cf926d5f9cb2
# ╠═2fd2e5bf-bd4c-4680-ad21-7bf45f804539
# ╠═2768dc6d-2c48-44cb-969a-fe2b16e68640
# ╠═9b51e143-b432-4783-95fe-b3110ab0d4da
# ╠═e19dd2b3-189f-473d-839c-5b2891da41ee
# ╠═3f011a6d-131e-4c70-b3ff-1a2c2df6c857
# ╠═c7c40326-5723-49e4-bf94-99c7d4ec5654
# ╠═58574dec-845a-4a4f-9ec7-9003e2fcf089
