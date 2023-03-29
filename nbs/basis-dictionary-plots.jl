### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ f29d92dc-cc76-11ed-355c-ed15ad184daa
using Pkg; Pkg.activate(".")

# ╔═╡ 4e426395-3235-4af7-8ed1-2c9b5f5a7d15
using Revise

# ╔═╡ 4aa3c240-9496-4ff2-9e14-5e59593189a2
using MultiscaleSimplexSignalTransforms

# ╔═╡ 977797f4-0fbf-4c62-b245-b5f3f7b540d9
using SparseArrays, LinearAlgebra, StatsBase, Arpack, Graphs, Plots, Chain, DataStructures, MAT, GraphRecipes

# ╔═╡ 9bca432a-7e65-40ad-9bd5-36683003a8d8
(dg, pos) = tricycle(5)

# ╔═╡ fa14f601-dd2e-41c6-acde-b574cd5fc7f3
region = KRegion(Graph(dg), 2; weakadjweight=1.0)

# ╔═╡ ca728f61-6f8b-49e4-a324-dd63ddb2c82e
part = kGHWT(region, SubRepresentation.Full; eigenmethod=EigenMethod.Positive)

# ╔═╡ 22fd0a60-46f0-4735-9dce-fcb1ee5e2aed
splot(region, part[0,0,1], pos; palette=:viridis)

# ╔═╡ 0060789d-3e75-4724-96b6-7b4ae7cdeefe
@chain begin
plot([
	splot(
		region, treecolors(root(part); maxlevel=l), pos;
		palette=:viridis, clims = l==0 ? (-1,1) : (0,1)
	)
	for l=0:5
]..., layout=(1,6), size=(6*600,400))
	# @aside savefig("/Users/eug/Desktop/hierarchical-partition-line.png")
end

# ╔═╡ 585cff0e-4b1d-43bf-8794-97a6dc02cdbd
haar = [part[0,0,0] reduce(hcat, [
	part[j-1, k-1, 1] 
	for (j,o) ∈ enumerate(part.inds.offsets) 
	for k ∈ 1:length(o)
	if o[k] < n_simp(region)-1 && (k==length(o) || o[k+1] > o[k]+1)
])]

# ╔═╡ f6234070-24fd-4150-9832-7fa42f2aa0aa
@chain begin
plot([
	splot(
		region, Vector(h), pos; palette= i==1 ? fill(cgrad(:viridis)[1.0], 2) : :viridis
	)
	for (i,h) ∈ enumerate(eachcol(haar))
]..., layout=(@layout [Plots.grid(3,7) ]), size=(7*400,3*400))
# @aside savefig("/Users/eug/Desktop/2HaarBasis.png")
end

# ╔═╡ 776b4de3-f510-4b83-9c28-910d63912cc2
@chain begin
	plot([
		plot([
			splot(region, v, pos; palette= (i,j)==(1,1) ? fill(cgrad(:viridis)[1.0], 2) : :viridis)
			for (j,v) ∈ enumerate(collect.(eachcol(M)))
		]..., layout=(1, 21))
		for (i,M) ∈ enumerate(collect.(eachslice(part[:,:], dims=3)))
	]..., layout=(6,1), size=(400*21, 400*6))
	# @aside savefig("/Users/eug/Desktop/C2F.png")
end

# ╔═╡ 497c7a42-7c8c-465f-a023-ab8abf43bdd6
f2c_ord = [
	sortperm(reduce(vcat,[
		(k, i)
		for (k,lim) ∈ enumerate(diff(v) |> d -> [d; 21-sum(d)])
		for i ∈ 1:lim
	]), by=reverse)
	for v ∈ part.inds.offsets
]

# ╔═╡ 7ddc23f6-3f2b-4a47-b97e-7724ba048e5d
@chain begin
	plot([
		plot([
			splot(region, v, pos; palette= (i,j)==(1,1) ? fill(cgrad(:viridis)[1.0], 2) : :viridis)
			for (j,v) ∈ enumerate(collect.(eachcol(M[:,f2c_ord[i]])))
		]..., layout=(1, 21))
		for (i,M) ∈ reverse(collect(enumerate(collect.(eachslice(part[:,:], dims=3)))))
	]..., layout=(6,1), size=(400*21, 400*6))
	# @aside savefig("/Users/eug/Desktop/F2C.png")
end

# ╔═╡ 52378760-57c7-4bc2-933f-814a83aa3bf7
hgpart = kHGLET(region)

# ╔═╡ fa8b048d-f975-44f6-8f71-40ce171c48b5
@chain begin
	plot([
		plot([
			splot(region, v, pos; palette=:viridis)
			for (j,v) ∈ enumerate(collect.(eachcol(M)))
		]..., layout=(1, 21))
		for (i,M) ∈ enumerate(collect.(eachslice(hgpart[:,:], dims=3)))
	]..., layout=(6,1), size=(400*21, 400*6))
	# @aside savefig("/Users/eug/Desktop/HGLET_combo.png")
end

# ╔═╡ Cell order:
# ╠═f29d92dc-cc76-11ed-355c-ed15ad184daa
# ╠═4e426395-3235-4af7-8ed1-2c9b5f5a7d15
# ╠═4aa3c240-9496-4ff2-9e14-5e59593189a2
# ╠═977797f4-0fbf-4c62-b245-b5f3f7b540d9
# ╠═9bca432a-7e65-40ad-9bd5-36683003a8d8
# ╠═fa14f601-dd2e-41c6-acde-b574cd5fc7f3
# ╠═ca728f61-6f8b-49e4-a324-dd63ddb2c82e
# ╠═22fd0a60-46f0-4735-9dce-fcb1ee5e2aed
# ╠═0060789d-3e75-4724-96b6-7b4ae7cdeefe
# ╠═585cff0e-4b1d-43bf-8794-97a6dc02cdbd
# ╠═f6234070-24fd-4150-9832-7fa42f2aa0aa
# ╠═776b4de3-f510-4b83-9c28-910d63912cc2
# ╠═497c7a42-7c8c-465f-a023-ab8abf43bdd6
# ╠═7ddc23f6-3f2b-4a47-b97e-7724ba048e5d
# ╠═52378760-57c7-4bc2-933f-814a83aa3bf7
# ╠═fa8b048d-f975-44f6-8f71-40ce171c48b5
