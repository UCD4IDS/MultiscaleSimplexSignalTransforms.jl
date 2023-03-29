### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 13c6d514-b80d-11ed-3122-8b24954da7b6
using Pkg;  Pkg.activate(".")

# ╔═╡ 947688eb-2f74-47d3-9af4-bb869d24c0db
using Revise

# ╔═╡ b1a3e2d7-5d18-4b04-9a8c-6c1308012a55
using Downloads, RData

# ╔═╡ b2036b7b-cc77-46a4-88d5-e48346d7ddb7
using Roots, Unzip, Distances, StatsBase, LinearAlgebra, Graphs, SimpleWeightedGraphs, Arpack

# ╔═╡ d17e38e0-c1d9-4bfb-99e7-442eb62dd19c
using MultiscaleSimplexSignalTransforms

# ╔═╡ 38c3ac07-5215-457f-8f38-d5dd3dcc97dd
using Chain

# ╔═╡ 189538b6-6881-418a-9ca3-e4f7bfd12dbe
using Plots

# ╔═╡ 1fe54161-c96c-4dff-8985-89c28e04dfae
using DataStructures

# ╔═╡ 56bbf15b-b27a-4f0f-b904-5d420b67f4c7
using SparseArrays

# ╔═╡ fa7007f1-6898-41b1-9202-e76645e72353
# data = Downloads.download(
# 	"http://scalefreegan.github.io/Teaching/DataIntegration/data/kernel.rda"
# ) |> load

# ╔═╡ d864a54d-2874-441e-af6d-054d8256a7b1
# ns, maxes = dmaxcurve([50, 100, 200, 500, 1000])

# ╔═╡ e922620e-4839-4f46-b187-6fde38b6ebd7
# plot(dmaxcurve([50, 100, 200, 500, 1000])...)

# ╔═╡ 56b1873a-d52a-4c8a-9ed4-eabd32c488c8
struct Poses{R}
	x::Vector{R}
	y::Vector{R}
	c::Vector{Int}
end

# ╔═╡ 670a0b93-fb7d-4257-b6ac-aa3798e9f675
function kernel(pos::Poses; δ = 1, ϵ=1e-6)
	D = pairwise(SqEuclidean(), [pos.x pos.y]; dims=1)
	map(x -> x > ϵ ? x : 0.0, exp.(-D ./ δ) - I)
	# (
	# 	K, # affinity version
	# 	map(x -> x≈0 ? 0.0 : x^-1, K) # distance version
	# )
	# K = map(x -> x > ϵ ? x : 0.0, exp.(-D ./ δ) - I)
end

# ╔═╡ 068464d0-4758-4ee5-8fa6-717c7a59728c
n = 750

# ╔═╡ f1649846-4203-4910-9b58-8f2f81d0dddf
dmax = 80

# ╔═╡ f9f2804f-5167-4767-b916-d5ecf5c13b97
spiralwidth, kernelwidth = 0.6, 0.009

# ╔═╡ e33f1504-b368-4559-a5b6-0ae118b474bb
pos = Poses(spirals(n, 2; δ=spiralwidth)...);

# ╔═╡ f326645b-9b28-43aa-8efe-f8e8b8c1d4e2
W = kernel(pos, δ=kernelwidth);

# ╔═╡ c1eb0535-051e-45e4-b1fb-1e9495b6682f
histogram(sum(W .≉ 0, dims=2)[:], legend=:false)

# ╔═╡ 2236a8df-b152-4b6d-b556-41c91f02fbfe
sum(W.>0) / (size(W,1)*(size(W, 1)-1))

# ╔═╡ 79e99dbf-4688-4dfa-86c2-c18d5526d9ea
g = SimpleWeightedGraph(W);

# ╔═╡ 9c6bfa55-c1fd-41e2-9d10-f50c70bd2558
reg = ZeroRegion(g);

# ╔═╡ 4c4ed5fe-b10d-49c9-8f36-014d1fd90ea9
part = SubmatrixPartition(reg); partition!(part);

# ╔═╡ c0478874-8d4b-45d7-bcb6-d1ca9c1f8c64
begin
	function partscatter(v::Vector{Int}, pos::Poses; pltkwargs...)
		cs = distinguishable_colors(
			length(v)+1, parse.(Colorant, [:white, :blue, :red])
		)[2:end]
		push!(cs, colorant"grey")
		
		scatter(
			pos.x, pos.y; 
			c = getindex.(Ref(cs), v), 
			aspect_ratio=:equal, legend=:none, showaxis=false, grid=false,
			ms=4, msw=0,
			pltkwargs...
		)
	end

	partscatter(v::Vector{<:Real}, pos::Poses; kwargs...) = partscatter(
		map(
			x -> abs(x) < 1e-15 ? length(v)+1 : (x > 0 ? 1 : 2),
			v .* sign(v[1])
		), pos; kwargs...
	)
	partscatter(pos::Poses; kwargs...) = partscatter(pos.c, pos; kwargs...)
end

# ╔═╡ 996b1670-0505-4b28-bc62-92c32406a045
plot(sort(sum(W, dims=2)[:].^0.5), label="sqrt degrees")

# ╔═╡ e94f9356-c8cf-4c12-9ff0-a01175054d02
length(connected_components(g))

# ╔═╡ ca4c4481-61dc-4995-ab23-5773f73fb404
SparseArrays.nnz(M::Symmetric{T, S}) where {T, S<:AbstractSparseArray} = nnz(parent(M))

# ╔═╡ e720366e-8288-421c-96fa-dea4ce0255b2
SparseArrays.nonzeros(M::Symmetric{T, S}) where {T, S<:AbstractSparseArray} = nonzeros(parent(M))

# ╔═╡ 2debe552-860f-4b86-8ce0-7c2c7bcfcc80
function dmaxcurve(ns, dmax=1e12)
	maxes = []
	for n ∈ ns
		@chain n begin
			spirals
			Poses(_...)
			kernel(δ=.02)
			ZeroRegion
			distance_kernel(; dmax)
			2_.M
			nonzeros
			minimum
			dmax-_
			push!(maxes, _)
		end
		
		# pos = Poses(spirals(n)...)
		# reg = ZeroRegion(kernel(pos, δ=.02))
		# push!(maxes, maximum(nonzeros(distance_kernel(reg; dmax=1e12).M)))
	end
	ns, maxes
end

# ╔═╡ 9100491e-5c60-44b0-a6ce-1a6e192cadec
@chain begin
plot(
	# partscatter(pos),
	partscatter(indicator(part), pos),
	# partscatter(dpart, pos)
	title="Laplacian-based partition"
)
	# @aside savefig("/Users/eug/Desktop/lap-spiral.png")
end

# ╔═╡ f6532dad-3d7c-4f77-8e87-890e9f3630af
L = k_laplacian(reg);

# ╔═╡ 702b3b0f-bdd7-452c-aa06-2ca5b235370b
D = Diagonal(sum(W, dims=2)[:]);

# ╔═╡ 3d3df494-8f9f-40f1-a50b-6887ce0bbe42
# plot(
# 	partplot(P₁(length(dpartnop))*dpartnop, title="proj without P"),
# 	partplot(dpart, title="just apply P")
# )

# ╔═╡ 9cc57927-c531-43ee-8f46-b94a678295e1
begin
	P₁(n) = I - ones(n,n)./n
	P₁(M::AbstractMatrix) = P₁(size(M, 1))
end

# ╔═╡ 1a9d58e8-7a27-4307-91c3-4bfd2c4561f1
function component_inds(ccs)
	inds = SortedDict{Int, Int}()
	for (i, comp) ∈ enumerate(ccs), compi ∈ comp
		inds[compi] = i
	end
	collect(values(inds))
end

# ╔═╡ 3cca593f-a581-4b3e-b072-2b73670baba1
ds = [550, 400, 350, 300, 280, 250, 210, 170, 130, 90, 60, 30]
# ds = [300, 30, 20]

# ╔═╡ 9f495949-4c4a-48b9-a96d-fa35a97101eb
Bmult, Kmult = unzip([
	distance_kernel(reg; dmax=d, withdistmat=true)
	for d ∈ ds
]);

# ╔═╡ b6929c09-1743-44e5-9183-229543be3d65
alldists = ds[1] .- filter(<(ds[1]), nonzeros(2Bmult[1]));

# ╔═╡ 065d7b72-7ba0-4ca2-a674-8366756a685b
function interpolate(x, xs, ys)
	i = findfirst(>(x), xs)
	# if the first element is bigger than you, y[1]
	# if no element is bigger than you, y[end]
	# otherwise, you are between two elements x[i-1], x[i]. interpolate.
	if i==1
		ys[1]
	elseif isnothing(i)
		ys[end]
	else
		(xs[i]-x)/(xs[i]-xs[i-1]) * (ys[i-1]-ys[i]) + ys[i]
	end
end

# ╔═╡ 1a5f6166-fe6b-4e40-ba77-d7c10a8f4a36
function markhist(allds, ds)
	hist = fit(Histogram, alldists; nbins=300)
	xmax = maximum(hist.weights)
	dds = filter(≤(xmax), ds)
	ylims = interpolate.(
		dds, Ref(collect(hist.edges[1])[2:end]), Ref(hist.weights)
	)
	plt = plot(
		hist.edges[1] |> x -> [x[1:end-1]';x[2:end]'][:],
		hist.weights |> x -> [x';x'][:];
		legend=false,
		title="Histogram over all distances"
	)
	plot!(
		plt, [dds dds]', [zeros(length(dds)) ylims]', c=:grey, style=:dash
	)

	cdf = cumsum(hist.weights) |> x -> x / x[end]
	plt2 = plot(
		hist.edges[1][2:end], cdf, legend=false,
		title="Sparsity level at given distance cutoff"
	)
	cdflims = interpolate.(
		dds, Ref(collect(hist.edges[1])[2:end]), Ref(cdf)
	)
	plot!(
		plt2, [dds dds]', [zeros(length(dds)) cdflims]', c=:grey, style=:dash
	)
	# plot(plt)
	plot(plt, plt2, layout=(2,1), size=(700, 800))
end

# ╔═╡ 575d4263-3660-4b6f-ac98-f9a5879084d4
@chain begin 
	markhist(alldists, ds)
	# savefig("~/Desktop/spiral-histograms.png")
end

# ╔═╡ 083a53a8-b788-4aff-8073-14657195fa9a
Keigs = [eigs(KK, nev=1, which=:LR)[2][:,1] for KK ∈ Kmult];

# ╔═╡ 7d316a6c-3f9d-46a3-ab80-21419863cda5
@chain (3,4) begin
	plot([
		partscatter(
			KK, pos;
			title="censored at $(ds[i])", ms=8, titlefontsize=40,
		)
		for (i, KK) ∈ enumerate(Keigs)
	]...; layout=_, size=600 .*reverse(_))
	# savefig("~/Desktop/spiral-partitions.png")
end

# ╔═╡ e5c7cf36-fbbe-4b5b-9a67-2f01e00f6728
# plot([
# 	plot(1:10, 1:10)
# 	for _=1:12
# ]...; layout=(6,2), size=(1200, 3600))

# ╔═╡ Cell order:
# ╠═13c6d514-b80d-11ed-3122-8b24954da7b6
# ╠═b1a3e2d7-5d18-4b04-9a8c-6c1308012a55
# ╠═b2036b7b-cc77-46a4-88d5-e48346d7ddb7
# ╠═947688eb-2f74-47d3-9af4-bb869d24c0db
# ╠═d17e38e0-c1d9-4bfb-99e7-442eb62dd19c
# ╠═fa7007f1-6898-41b1-9202-e76645e72353
# ╠═38c3ac07-5215-457f-8f38-d5dd3dcc97dd
# ╠═189538b6-6881-418a-9ca3-e4f7bfd12dbe
# ╠═2debe552-860f-4b86-8ce0-7c2c7bcfcc80
# ╠═d864a54d-2874-441e-af6d-054d8256a7b1
# ╠═e922620e-4839-4f46-b187-6fde38b6ebd7
# ╠═670a0b93-fb7d-4257-b6ac-aa3798e9f675
# ╠═56b1873a-d52a-4c8a-9ed4-eabd32c488c8
# ╠═c1eb0535-051e-45e4-b1fb-1e9495b6682f
# ╠═e33f1504-b368-4559-a5b6-0ae118b474bb
# ╠═f326645b-9b28-43aa-8efe-f8e8b8c1d4e2
# ╠═2236a8df-b152-4b6d-b556-41c91f02fbfe
# ╠═79e99dbf-4688-4dfa-86c2-c18d5526d9ea
# ╠═9c6bfa55-c1fd-41e2-9d10-f50c70bd2558
# ╠═4c4ed5fe-b10d-49c9-8f36-014d1fd90ea9
# ╠═068464d0-4758-4ee5-8fa6-717c7a59728c
# ╠═f1649846-4203-4910-9b58-8f2f81d0dddf
# ╠═f9f2804f-5167-4767-b916-d5ecf5c13b97
# ╠═c0478874-8d4b-45d7-bcb6-d1ca9c1f8c64
# ╠═996b1670-0505-4b28-bc62-92c32406a045
# ╠═e94f9356-c8cf-4c12-9ff0-a01175054d02
# ╠═ca4c4481-61dc-4995-ab23-5773f73fb404
# ╠═e720366e-8288-421c-96fa-dea4ce0255b2
# ╠═9100491e-5c60-44b0-a6ce-1a6e192cadec
# ╠═f6532dad-3d7c-4f77-8e87-890e9f3630af
# ╠═702b3b0f-bdd7-452c-aa06-2ca5b235370b
# ╠═3d3df494-8f9f-40f1-a50b-6887ce0bbe42
# ╠═9cc57927-c531-43ee-8f46-b94a678295e1
# ╠═1fe54161-c96c-4dff-8985-89c28e04dfae
# ╠═1a9d58e8-7a27-4307-91c3-4bfd2c4561f1
# ╠═56bbf15b-b27a-4f0f-b904-5d420b67f4c7
# ╠═3cca593f-a581-4b3e-b072-2b73670baba1
# ╠═9f495949-4c4a-48b9-a96d-fa35a97101eb
# ╠═b6929c09-1743-44e5-9183-229543be3d65
# ╠═065d7b72-7ba0-4ca2-a674-8366756a685b
# ╠═1a5f6166-fe6b-4e40-ba77-d7c10a8f4a36
# ╠═575d4263-3660-4b6f-ac98-f9a5879084d4
# ╠═083a53a8-b788-4aff-8073-14657195fa9a
# ╠═7d316a6c-3f9d-46a3-ab80-21419863cda5
# ╠═e5c7cf36-fbbe-4b5b-9a67-2f01e00f6728
