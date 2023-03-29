### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# ╔═╡ 862160fc-e22a-11ec-2811-6bbfc785b5e9
using MultiscaleSimplexSignalTransforms

# ╔═╡ 48fadd84-dfc1-4b0e-8208-8bba5bf2cec8
using Graphs, Plots, GraphRecipes, LaTeXStrings

# ╔═╡ 5f7f6135-81eb-447f-b64c-b60a18f58a80
using LinearAlgebra, Arpack, Printf

# ╔═╡ 7378739b-4ac4-4e57-b018-45ba63d1f7e3
using MultiscaleGraphSignalTransforms

# ╔═╡ b0f07ece-9a02-43be-aea6-f2960053e841
using StatsBase

# ╔═╡ 0d628f32-46a9-4182-9d61-7ebcde3ac230
using Clustering

# ╔═╡ b0d32412-15a9-488c-affa-2d2e51270f73
using NMF

# ╔═╡ 4fde8488-db77-4cf9-8d52-1508025bbd7d
# automatically use Mac M1 hardware acceleration for linear algebra
# TODO: get SetBlasInt in the registry
# begin
# 	using SetBlasInt
# 	BLAS.lbt_forward("/System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate")
#     setblasint(Int32, :all)
# end

# ╔═╡ d8b0334b-d0fc-4835-a938-df6260a83a2d
function thick_dipath(n)
	g = cartesian_product(path_digraph(n), path_digraph(2))
	for i = 1:n-1
		# add_edge!(g, 2i+2, 2i-1)
		add_edge!(g,2i-1, 2i+2)
	end
	g
end

# ╔═╡ df97f7ef-45cc-49b6-b45c-01be50e0e0f4
function simplex_path(n, k=1)
	@assert k<=1 "orders higher than k=1 not implemented yet"
	k == 0 && return (path_digraph(n), √3*(0:n-1), 1:n)
	g = DiGraph(n+k)
	for i=1:n-1
		add_edge!(g, i+1, i)
		add_edge!(g, i, i+2)
	end
	add_edge!(g, n+1, n)
	(
		g, 
		√3*[0; repeat(1:(n+1)/2, inner=2)][1:end-(n % 2 == 0 ? 0 : 1)],
		collect(Iterators.flatten(zip(1:n/2+1, 0:n/2)))[1:end-(n % 2 == 1 ? 0 : 1)]
	)
end

# ╔═╡ 36f36dcb-b149-48ca-978a-e908a82c505e
n = 9

# ╔═╡ a60ca9f6-166f-493e-ad9a-fc37d9f43a79
gs = simplex_path(n); g = gs[1]

# ╔═╡ a208a712-461d-4223-9ff4-40be5d852c97
_plotg_defaults = Dict{Symbol, Any}(
	:mc => colorant"skyblue",
	:lc => colorant"black",
	:ms => 10,
	:edgemargin => 209,
	:elo => 0.05,
	:ew => 2,
	:ts => 8,
	:lz => nothing
);

# ╔═╡ ed193858-5c71-4476-b564-cf4add2b3159
function _plotg_set_defaults(kwpairs)
	kwargs = Dict{Symbol, Any}(kwpairs)
	for (k, v) in pairs(_plotg_defaults)
		get!(kwargs, k, v)
	end
	kwargs
end

# ╔═╡ b3b02029-cb76-4a01-a29a-93673b404e15
function plotg(
	g::AbstractGraph, x::Array, y::Array; kwpairs...
)
	kwargs = _plotg_set_defaults(kwpairs)
	xplt = []
	yplt = []
	label_pos = Tuple[]
	for (i, e) in enumerate(edges(g))
		s, d = Tuple(e)
		θ = atan(y[d]-y[s], x[d]-x[s])
		
		xₑ = kwargs[:ms]/kwargs[:edgemargin]*cos(θ)
		yₑ = kwargs[:ms]/kwargs[:edgemargin]*sin(θ)

		xₐ = (x[s] + x[d])/2 + kwargs[:elo] * sin(θ)
		yₐ = (y[s] + y[d])/2 - kwargs[:elo] * cos(θ)
		
		push!(xplt, x[s]+xₑ, x[d]-xₑ, NaN)
		push!(yplt, y[s]+yₑ, y[d]-yₑ, NaN)
		push!(label_pos, (xₐ, yₐ, text(L"\mathbf{%$i}", kwargs[:ts])))
	end
	plt = plot(
		xplt, yplt;
		legend=false,
		colorbar=true,
		framestyle=:none,
		arrow=true,
		lc=kwargs[:lc],
		lw=kwargs[:ew],
		lz=isnothing(kwargs[:lz]) ? nothing : repeat(kwargs[:lz], inner=3)
	)
	scatter!(
		x[:], y[:];
		mc = kwargs[:mc],
		ms = kwargs[:ms],
		msw = kwargs[:ew],
		text = [text(L"\mathbf{%$a}", kwargs[:ts]) for a in ('a':'z')[1:ne(g)]]
	)
	annotate!(label_pos)
	plt
end

# ╔═╡ 97c35617-b1fd-4355-8c54-687facbb84c3
L = hodgelaplacian(g; which=:upper)

# ╔═╡ 10ececbc-d408-4dd2-9c35-f0e950dfc231
Ld = L + Diagonal(L)

# ╔═╡ 08e38d75-fbe3-4a55-8846-1a8600edc4c6
eigd = eigs(Ld; which=:SM, nev=size(Ld, 1))

# ╔═╡ eb69fbe8-2b31-427b-bc31-139ddc47cf34
plotg(
	gs...;
	elo= ne(g) > 9 ? 0.12 : 0.05,
	edgemargin = ne(g) > 9 ? 60 : 160,
	lz = eigd[2][:,1],
	lc = :redblue
	# lz = 1:ne(g)
)

# ╔═╡ 62ac414e-05d9-4368-a519-ffd098fee827


# ╔═╡ d7acf7bd-fb26-40df-ae47-6241823c82f6
Ds = Diagonal(Ld)^-.5

# ╔═╡ 031f5d4a-6ed0-4fec-81a6-4d544125916a
Ls = Ds * Ld * Ds

# ╔═╡ 8a4f8001-88e0-4f1d-97b1-fb3442b3924b
Λ, X = eigs(
	# L + Diagonal(L) |> M -> Diagonal(sum(M, dims=2)[:]) \ M; 
	Ls; 
	nev=size(Ls, 1), 
	which=:SM
)

# ╔═╡ 4cc21461-5e9b-475c-9939-b6c6d7719530
plot(Λ)

# ╔═╡ 5b812fc9-d80b-4736-955f-3fa388f178af
use_rw = true

# ╔═╡ 8722a839-492a-40e5-ae60-7525349a2537
Xnorm = use_rw ? Ds*X : X

# ╔═╡ 744c4cd2-08f2-4246-9563-065c88680426
nsame = findlast(≈(Λ[1]), Λ[2:end]) |> i -> isnothing(i) ? 1 : i+1
# nsame = 4

# ╔═╡ d55f96c4-e1b5-4e9b-ad58-59404aa87130
function seminmf(X, k; imax = 100, tol = 1e-8, initconst = 0.2)
	(n, p) = size(X)
	assn = kmeans(X', k).assignments
	G = fill(initconst, n, k)
	for (i, a) ∈ enumerate(assn)
		G[i, a] += 1
	end
	iter = 1
	while true
		F = X'*G/(G'*G)
		XF = X*F
		FF = F'*F
		GFF₊ = G*max.(FF, 0)
		GFF₋ = G*max.(-FF, 0)
		G .*= sqrt.((max.(XF, 0)+GFF₋)./(max.(-XF, 0)+GFF₊))
		res = norm(X - G*F')
		if res < tol || iter >= imax
			println("After $iter iterations, reached residual of $res")
			return (F, G)
		end
		iter += 1
	end
end

# ╔═╡ 6fc84842-d9b2-4c3c-90cd-3f84677c444d
F, G = seminmf(Xnorm[:, 1:nsame], nsame; imax=1000)

# ╔═╡ 404ebacb-0940-4f5f-85c1-6c8594568134
GG = G./sqrt.(sum(G.^2, dims=1))

# ╔═╡ 9523b94f-ccbc-48bc-9423-6c8c64812ecd
plot(
	# plot([GG Xnorm[:,nsame+1]]; title="nonnegative representation"), 
	# plot(Xnorm[:,1:nsame+1]; title="original eigenvectors"),
	plot([GG Xnorm[:,nsame+1]]), 
	plot(Xnorm[:,1:nsame+1]),
	legend=nothing
)

# ╔═╡ d4534817-dd1a-4ac9-98c8-42274531b78e
els = [
	src(e) % 2 == 0 ?
		dst(e) % 2 == 0 ?
			"top" :
			"rung" :
		dst(e) % 2 == 0 ?
			"rung" :
			"bot"
	for e in edges(g)
]

# ╔═╡ bb6f6758-b35d-4a00-aa78-4c3a7dc1e73c
scatter(
	GG[:, 1], X[:, nsame+1];
	series_annotations=els,
	c = [
		a == 1 ? colorant"red" : colorant"blue"
		for a in kmeans([GG[:,1] X[:, nsame+1]]', 2).assignments
	],
	legend=false
)

# ╔═╡ 6963074c-e7a8-4788-b400-ca96d79a1bb7
# XX = StatsBase.transform(fit(UnitRangeTransform, X; dims=1), X)

# ╔═╡ faa2fafe-49f6-4032-af60-c48addf0452e
# cgrad(:viridis, 1:10)

# ╔═╡ 7840d3ba-a256-4a29-a076-3fe842208099
# ╠═╡ disabled = true
#=╠═╡
graphplot(
	g;
	names=[ @sprintf "%2i" i for i in 1:2n],
	nodesize=0.3,
	nodeshape=:circle,
	x = .25*repeat(1:n; inner=2),
	y = .06*[1:2:2n-1 4:2:2n+2]'[:],
	curves=false,
	markerstrokewidth=0,
	# edgelabel=Dict(Tuple(e)=>"$i    \n" for (i, e) in enumerate(edges(g))),
	# edgecolor=Dict(
	# 	Tuple(e) => cgrad(:viridis, XX[i,1])
	# 	for (i, e) in enumerate(edges(g))
	# ),
	# edge_label_box=false
)
  ╠═╡ =#

# ╔═╡ d12dd173-3d1f-4ab2-9181-eb269deee901
assn = kmeans(X[:, nsame .+ (0:1)]', 2).assignments

# ╔═╡ 33e4b782-7967-41e3-855d-cf401b93da4e
# ╠═╡ disabled = true
#=╠═╡
nnmf(X[:,1:3], 3)
  ╠═╡ =#

# ╔═╡ 772fc321-0a2d-45c6-8e00-e943f84f0d53


# ╔═╡ 43c99eca-97cc-4e82-b52d-f934d55c9ebb
scatter(
	X[:, nsame], X[:, nsame+1];
	# X[:, 1], X[:, 4]; 
	series_annotations=Pair.(edges(g)),
	c = [a == 1 ? colorant"red" : colorant"blue" for a in assn],
	legend=false
)

# ╔═╡ 880fd10e-6bda-4c2e-aa30-018ce48681d6
begin
	g1 = copy(g)
	g2 = copy(g)
	es = g |> edges |> collect
	for (i, a) in assn |> enumerate
		rem_edge!(a == 1 ? g2 : g1, es[i])
	end
end

# ╔═╡ e51921b9-b165-4e22-989b-b259bdd73e61
# ╠═╡ disabled = true
#=╠═╡
plot([
	graphplot(
		Graph(gi);
		names=[ @sprintf "%2i" i for i in 1:2n],
		# nodesize=0.06,
		nodeshape=:circle,
		x = .2*repeat(1:n; inner=2),
		y = .04*[1:2:2n-1 4:2:2n+2]'[:],
		curves=false,
		# edgelabel=Dict(Tuple(e)=>"$i    \n" for (i, e) in enumerate(edges(g)))
	)
	for gi in (g1, g2)
]...)
  ╠═╡ =#

# ╔═╡ 117b4b34-c5b3-4c0e-9292-183147b4becd


# ╔═╡ Cell order:
# ╠═862160fc-e22a-11ec-2811-6bbfc785b5e9
# ╠═48fadd84-dfc1-4b0e-8208-8bba5bf2cec8
# ╠═4fde8488-db77-4cf9-8d52-1508025bbd7d
# ╠═d8b0334b-d0fc-4835-a938-df6260a83a2d
# ╠═df97f7ef-45cc-49b6-b45c-01be50e0e0f4
# ╠═36f36dcb-b149-48ca-978a-e908a82c505e
# ╠═a60ca9f6-166f-493e-ad9a-fc37d9f43a79
# ╠═a208a712-461d-4223-9ff4-40be5d852c97
# ╠═ed193858-5c71-4476-b564-cf4add2b3159
# ╠═b3b02029-cb76-4a01-a29a-93673b404e15
# ╠═eb69fbe8-2b31-427b-bc31-139ddc47cf34
# ╠═5f7f6135-81eb-447f-b64c-b60a18f58a80
# ╠═97c35617-b1fd-4355-8c54-687facbb84c3
# ╠═10ececbc-d408-4dd2-9c35-f0e950dfc231
# ╠═08e38d75-fbe3-4a55-8846-1a8600edc4c6
# ╠═62ac414e-05d9-4368-a519-ffd098fee827
# ╠═4cc21461-5e9b-475c-9939-b6c6d7719530
# ╠═d7acf7bd-fb26-40df-ae47-6241823c82f6
# ╠═031f5d4a-6ed0-4fec-81a6-4d544125916a
# ╠═8a4f8001-88e0-4f1d-97b1-fb3442b3924b
# ╠═5b812fc9-d80b-4736-955f-3fa388f178af
# ╠═8722a839-492a-40e5-ae60-7525349a2537
# ╠═7378739b-4ac4-4e57-b018-45ba63d1f7e3
# ╠═744c4cd2-08f2-4246-9563-065c88680426
# ╠═d55f96c4-e1b5-4e9b-ad58-59404aa87130
# ╠═6fc84842-d9b2-4c3c-90cd-3f84677c444d
# ╠═404ebacb-0940-4f5f-85c1-6c8594568134
# ╠═9523b94f-ccbc-48bc-9423-6c8c64812ecd
# ╠═d4534817-dd1a-4ac9-98c8-42274531b78e
# ╠═bb6f6758-b35d-4a00-aa78-4c3a7dc1e73c
# ╠═b0f07ece-9a02-43be-aea6-f2960053e841
# ╠═6963074c-e7a8-4788-b400-ca96d79a1bb7
# ╠═faa2fafe-49f6-4032-af60-c48addf0452e
# ╠═7840d3ba-a256-4a29-a076-3fe842208099
# ╠═0d628f32-46a9-4182-9d61-7ebcde3ac230
# ╠═d12dd173-3d1f-4ab2-9181-eb269deee901
# ╠═b0d32412-15a9-488c-affa-2d2e51270f73
# ╠═33e4b782-7967-41e3-855d-cf401b93da4e
# ╠═772fc321-0a2d-45c6-8e00-e943f84f0d53
# ╠═43c99eca-97cc-4e82-b52d-f934d55c9ebb
# ╠═e51921b9-b165-4e22-989b-b259bdd73e61
# ╠═880fd10e-6bda-4c2e-aa30-018ce48681d6
# ╠═117b4b34-c5b3-4c0e-9292-183147b4becd
