### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ 06a6f062-9936-11ed-0bdf-d7de883b03a4
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 79e982da-6297-4b2e-90d7-4dddcf8d2638
using Revise

# ╔═╡ 6be1844b-aae3-4c83-b3d5-f50d57700346
using MultiscaleSimplexSignalTransforms

# ╔═╡ 7d3af02a-c608-4335-b361-8c1c47ee8a57
using SparseArrays, LinearAlgebra, StatsBase, Arpack, Graphs, Plots, Chain, DataStructures, MAT, GraphRecipes, InvertedIndices, SimpleWeightedGraphs

# ╔═╡ d9d1345f-c9d9-4caa-9c5b-839fcff75ce0
using GraphPlot

# ╔═╡ 78581d3c-48f1-4c59-8024-5e794c39f2be
A = matread("../data/RGC60_100.mat")

# ╔═╡ 596d9956-3cc4-4e0d-a71e-3d848e0dba15
function remove_diag!(M)
	for ij ∈ CartesianIndices(M)
		if ij.I[1] == ij.I[2]
			M[ij] = 0
		end
	end
	M
end

# ╔═╡ 49f92ad5-9f91-418a-8d6b-153cb14db7d3
g = Graph(-remove_diag!(A["L100"] |> copy))

# ╔═╡ 812a9e14-8a7e-47a5-b6cb-6cfc21f36d45
xs, ys = A["xyz100"] |> M -> (M[:,1], M[:,2])

# ╔═╡ 68b4d504-c8d1-4e85-89a2-fe394866a99e
sc = cliquecomplex(g);

# ╔═╡ 5fe9fab9-293b-40d5-95dd-4f6a1acce945
zr = ZeroRegion(sc)

# ╔═╡ 73dcf165-0d63-48de-bb66-48e41c673b52
zpart = kHGLET(zr)

# ╔═╡ f8a04b4a-4ccc-4cd1-8d58-fdccad1f21a8
@chain splot(zr, zpart[0,0,1], pos) _ #@aside savefig("/Users/eug/Desktop/fiedler-neuron.png")

# ╔═╡ 3ca7d5cb-19a1-4b29-b113-dfc3ef9b52d7


# ╔═╡ b3851a78-fd4d-40d4-96d5-6f6ae801b0d2
k = 1

# ╔═╡ e58b6dc0-2a39-4970-acf8-8317187ae822
zreg = KRegion(sc, 0);

# ╔═╡ 0f040b63-9c2d-4806-b1bf-b43dcb30894d
reg = KRegion(sc, k; weakadjweight = (k==0 ? 0. : 1.));

# ╔═╡ ee08ed68-8468-4d1b-b535-a8ee819b0b85
# # part = kHGLET(reg);
# part = kGHWT(
# 	reg;
# 	# normalization=Normalization.Combinatorial,
# 	input=PartitionInput.Identity,
# 	method=PartitionMethod.Kmeans
# );
partvec = eigs(
	k_laplacian(
		reg;
	)+I;
	nev=1, which=:SM, maxiter=10000
)[2][:,1]
# part = SubmatrixPartition(
# 	reg;
# 	normalization=Normalization.Combinatorial,
# 	input=PartitionInput.Identity,
# 	method=PartitionMethod.FiedlerSign
# ); MultiscaleSimplexSignalTransforms.partition!(part); partvec = indicator(part)

# ╔═╡ 2f636633-1973-46ce-9c42-cee6394e55e9
l, fourier = eigen(Matrix(k_laplacian(reg; normalization=Normalization.Combinatorial)))

# ╔═╡ 00a15a36-c113-40a5-9103-3e6e9f9c38c3
plot(l)

# ╔═╡ e5ada405-8d0d-4bb1-96bf-a7bcf538afe4
pos = Pos(xs, ys)

# ╔═╡ f483ea63-5858-4cb8-af27-38e019bc1994
splot(reg, sign.(fourier[:,1]).*fourier[:,4], pos; k=1)

# ╔═╡ 9016916b-1574-4a60-af2d-c8d18e79da48
# part2 = kHGLET(reg);

# ╔═╡ d64c419d-f8af-4c14-afc2-60ab357b635d
@chain begin
	plot([
		splot(reg, sign.(fourier[:,i]); k, xs, ys, simpsize=4, palette=:viridis)
		for i=2:5
	]..., layout=(1,4), size=(1600,400))
	# @aside savefig("/Users/eug/Desktop/dendrite-higher-eigs.png")
end

# ╔═╡ 7d4b1f53-ac94-4a40-9dde-27c5824336b0
hubs = @chain g begin
	edges
	collect
	[
		findall(e-> dst(e) ≠ src(e)+1, _)[2:end];
		# findall(e -> Graphs.degree(g, dst(e))==3, _)
	]
	unique
	sort
end
# 	filter(≤(nv(g)), sort(unique([
# 	findall(>(2), Graphs.degree(g));
# 	findall(==(1), Graphs.degree(g)) .+ 1
# ])))

# ╔═╡ 70f6e2c2-d526-440a-b8cc-471d6d0af707
findall(e-> dst(e) ≠ src(e)+1, collect(edges(g)))

# ╔═╡ 380acbc4-ce91-4ea1-9bd7-491b42dfc6fe
asunit(v) = extrema(v) |> mM -> 2(v .- mM[1]) ./ (mM[2]-mM[1]) .- 1

# ╔═╡ 27f71d57-d13a-4c15-a29b-fb2722ec4e57
branchcols = reduce(vcat, [
	asunit(
		(partvec[v] > 0 ? identity : reverse)(
		# (
			v:(i ≤ length(hubs) ? hubs[i] : ne(g)+1)-1
		)
	)
	for (i, v) ∈ enumerate([1; hubs])
])

# ╔═╡ 4b8a5946-c6e4-4ba6-826e-083e0a4fcba3
[
	asunit(v:(i ≤ length(hubs) ? hubs[i] : ne(g)+1)-1)
	for (i, v) ∈ enumerate([1; hubs])
]

# ╔═╡ 865b74e2-8ec9-4ba3-ad79-80acefac0019


# ╔═╡ e7aaef06-a81b-44c7-b8ec-ef05ecfd2001
# hub_id = sparsevec(findall(>(2), Graphs.degree(g)), 1.0, n)

# ╔═╡ a35b6304-6342-4cb7-b709-06981c5174f8
p = Pos(xs, ys)

# ╔═╡ fdf317eb-17eb-40b1-b6d0-f342f806b4c7
@chain begin
	splot(reg, branchcols, p; k, simpsize=4, palette=:viridis)
	# @aside savefig("/Users/eug/Desktop/dendrite-reorient.png")
end

# ╔═╡ e564a47e-77e4-4599-9ac6-8aef5b78fd94


# ╔═╡ b5cbba82-a40e-42c0-9cba-9444596c04cc
# splot(reg, part2[1, 0, 0]; k, xs, ys, size=4, palette=:redblue)


# ╔═╡ af54a055-54ef-4196-8cf4-9fae59996dbe
splot(
	zreg,
	Vector(sparsevec(27*40 .+(1:40), 1:40, n_simp(zreg)));
	k=0, xs, ys, simpsize=3,
	palette=:redblue,
	clims=(1, 40)
)

# ╔═╡ d1b59bb4-47b9-4b4f-9674-fd4ed9759996
# part3 = kGHWT(zreg; representation=Representation.Distance)

# ╔═╡ efb2b701-e055-4320-8ad1-853cb36a1cdf
t = randomtree(10)

# ╔═╡ ce488d89-d599-4d51-a51c-6b50a763b9e4
treg = ZeroRegion(t)

# ╔═╡ 6e351d96-5259-4608-b3b3-a6dae0d0c241
# tpart = kGHWT(
# 	ZeroRegion(t), SubRepresentation.Full;
# 	representation=Representation.Distance,
# 	input=PartitionInput.Identity
# )

# ╔═╡ b7ed2cc8-3fde-4489-b8a8-01395e6ed6d7
# tpart = SubmatrixPartition(
# 	treg;
# 	representation=Representation.Distance,
# 	input=PartitionInput.Identity
# ); partition!(tpart);

# ╔═╡ 570fe0e7-9e1f-4dd6-a669-b7fe235f12d4
x, y = GraphRecipes.tree_graph(adjacency_matrix(t))

# ╔═╡ b869ba37-694b-46ae-936b-af2c882c8fc7
# indicator(tpart)

# ╔═╡ 93dcee4a-0fa8-4f0d-adcb-6990d312f9f0
#Λ, X = distance_kernel(treg) |> eigen

# ╔═╡ de4dc2a4-c9c5-4182-a906-be668b6f3f48
# splot(treg, indicator(tpart); k=0, xs=x, ys=y)

# ╔═╡ 21e22972-59b2-4f49-8328-aaa64ac45c87
# SortedDict(
# 	l => sum(vs)
# 	for (l, vs) ∈ bin(
# 		x -> x[1],
# 		bfs(tsn -> (tsn.level, length(tsn.node.fvals)), part.part.root);
# 		valf = x -> x[2]
# 	)
# )

# ╔═╡ 4449be15-32ba-4ae7-8bdf-72e6ec4567f7
Base.return_types(sum, (eltype([1,2,3]),))[1]

# ╔═╡ 306b57a3-760d-4b81-b3dc-4aa3e025d2f2
function bin(f, xs; valf=identity)
	keytype = Base.return_types(f, (eltype(xs),))[1]
	valtype = Vector{Base.return_types(valf, (eltype(xs),))[1]}
	res = DefaultDict{keytype, valtype}(() -> [])
	for x ∈ xs
		push!(res[f(x)], valf(x))
	end
	res
end

# ╔═╡ c801d3e0-4d82-4b0a-b368-c59f208e25cf
pg = path_graph(10)

# ╔═╡ 9823485e-fd71-4e59-bc56-29c704a56e65
B = distmat(pg)

# ╔═╡ d26a7cec-f20b-4676-a5d2-0fe568cca3fc
K = -B/2

# ╔═╡ ee9418c6-ac2e-42a7-986b-a2a04cb384d5
L = laplacian_matrix(pg)

# ╔═╡ 636b9f68-c25c-434a-916b-e632d6931dd8
Lpinv = pinv(Matrix(L))

# ╔═╡ 0cf0ac1f-32b3-4ff8-bed2-d1c3304fbde6
Pone = I - ones(10, 10)/10

# ╔═╡ ebf69a0a-7b27-400a-8e99-8af4ca202570
Pone*K*Pone - Lpinv |> norm

# ╔═╡ 1e51e750-5271-46c8-afec-a98e440b3d16


# ╔═╡ aa5df510-3aa4-43ae-8f43-033f75e812b9
splot(zreg, 1:n_simp(zreg), pos; c=:viridis)

# ╔═╡ bf726737-561f-4ca3-9986-03b4066e7004
M = adjacency_matrix(g) |> UpperTriangular

# ╔═╡ 4c3d2136-7491-4615-85ec-b6d9b9dda45d
Ldir = Diagonal(sum(adjacency_matrix(g), dims=2)[:]) - M

# ╔═╡ 33d9063a-8afd-4c07-8d26-2873bf6a9d48
U, S, V = svd(Ldir)

# ╔═╡ d73b3587-53e2-4ba6-875c-ecee5210c7dd
plot(S); plot!([1141,1141], [1.0,3.5])

# ╔═╡ 2f96f9fd-dc69-4ad8-b5d6-7b093d6cf8ad
S[1140:end]

# ╔═╡ c3d5ca3c-7254-444f-9a36-d60062a890f5
splot(zreg, U[:,end-18], pos; c=:viridis)

# ╔═╡ 94d98795-a309-4d45-b06a-ee269a3df37b
begin
	dmax = 10000
	Δ = fill(dmax, size(M))
	for (i,j,v) in zip(findnz(distmat(DiGraph(M')))...)
		Δ[i,j] = v
	end
	Δ -= dmax*I
end

# ╔═╡ 60eb389d-9427-4178-9611-59dec0ec63e6
function Kharmonic()
	D = distmat(DiGraph(M'))
	II, JJ, V = findnz(D)
	K = sparse(vcat(II,  [1154]), vcat(JJ , [1154]), vcat(1.0 ./V, 0))
	# K -= I
end

# ╔═╡ 283456dc-2415-45c8-a811-cfd93ab23beb
Kh = Kharmonic()

# ╔═╡ 8e1ca445-a738-40be-ab4c-f948c18e2f41
KhU, KhS, KhV = svd(Matrix(Kh))

# ╔═╡ cdc36c27-df16-4448-9c7e-1ad22cfa4fd2
plot(KhS)

# ╔═╡ b3b2dc90-edc9-4200-8459-734a2c1dfcda
14 |> n-> splot_many(zreg, pos, eachcol(KhV[:,end-n:end])...; size=(600, 400*n), c=:seismic, simpsize=4)

# ╔═╡ ce6a1868-f8dc-49a9-99cf-9add95abd51e
SortedDict(Base.Order.Reverse, countmap(Δ))

# ╔═╡ 6081e9ba-72df-4ade-801d-d2434de14e57
DU, DS, DV = svd(Δ)

# ╔═╡ f99dcaff-8367-4a64-be52-7979ab003c20
plot(DS[1140:end], yscale=:log10)

# ╔═╡ 6b7d9922-0b34-4dba-9148-be13a92726d7
8 |> n-> splot_many(zreg, pos, eachcol(DV[:,1:n])...; size=(600, 400*n))

# ╔═╡ d977411f-dccd-4aad-8251-c570a9ffec32
begin
	# splot(zreg, ones(n_simp(zreg)), pos; c=:grays, simpsize=1)
	gplot(zreg.weights, pos.x, -pos.y)
	# scatter(pos.x, pos.y; marker_z=DV[:,5], c=:viridis)
end

# ╔═╡ 287a403e-ce53-472f-9995-db091d2d21a8
pos

# ╔═╡ 98b7e234-194e-44c1-b78e-3507e688eaea
# splot_many(zreg, pos, eachcol(DU[:,end-15:end])...; size=(600, 400*8))

# ╔═╡ 50096ed8-e785-4e90-8fa8-1637cec1be43
Plots.default(:size)

# ╔═╡ 933f7893-81ea-49e2-8c30-0bc6866ce7f5
Kdir = distance_kernel(DiGraph(M'); dmax=1000, issymmetric=false)

# ╔═╡ 45d26735-360f-4b30-a0eb-437280138d69
Kdense = Kdir*I

# ╔═╡ 0c8141b3-8a18-4189-a756-4e925ca28a05
KU, KS, KV = svd(Kdense)

# ╔═╡ cdfe43e2-b6f9-417e-b0a0-a07756e5e8a9
plot(KS[1150:end], yscale=:log10)

# ╔═╡ 013c164a-f4fc-4452-a2ad-d52f4b3c5cd2
plot(KS, yscale=:log10)

# ╔═╡ 77ea6942-4053-4fa9-b9c5-e60ae853847f
8 |> n -> splot_many(zreg, pos, eachcol(KV[:,1:n])...; size=(600, 400*n))

# ╔═╡ 957aece0-fc23-40f3-a6cf-cfb37ea6051a


# ╔═╡ Cell order:
# ╠═06a6f062-9936-11ed-0bdf-d7de883b03a4
# ╠═79e982da-6297-4b2e-90d7-4dddcf8d2638
# ╠═6be1844b-aae3-4c83-b3d5-f50d57700346
# ╠═7d3af02a-c608-4335-b361-8c1c47ee8a57
# ╠═78581d3c-48f1-4c59-8024-5e794c39f2be
# ╠═596d9956-3cc4-4e0d-a71e-3d848e0dba15
# ╠═49f92ad5-9f91-418a-8d6b-153cb14db7d3
# ╠═812a9e14-8a7e-47a5-b6cb-6cfc21f36d45
# ╠═68b4d504-c8d1-4e85-89a2-fe394866a99e
# ╠═5fe9fab9-293b-40d5-95dd-4f6a1acce945
# ╠═73dcf165-0d63-48de-bb66-48e41c673b52
# ╠═f8a04b4a-4ccc-4cd1-8d58-fdccad1f21a8
# ╠═3ca7d5cb-19a1-4b29-b113-dfc3ef9b52d7
# ╠═b3851a78-fd4d-40d4-96d5-6f6ae801b0d2
# ╠═e58b6dc0-2a39-4970-acf8-8317187ae822
# ╠═0f040b63-9c2d-4806-b1bf-b43dcb30894d
# ╠═ee08ed68-8468-4d1b-b535-a8ee819b0b85
# ╠═2f636633-1973-46ce-9c42-cee6394e55e9
# ╠═00a15a36-c113-40a5-9103-3e6e9f9c38c3
# ╠═e5ada405-8d0d-4bb1-96bf-a7bcf538afe4
# ╠═f483ea63-5858-4cb8-af27-38e019bc1994
# ╠═9016916b-1574-4a60-af2d-c8d18e79da48
# ╠═d64c419d-f8af-4c14-afc2-60ab357b635d
# ╠═7d4b1f53-ac94-4a40-9dde-27c5824336b0
# ╠═70f6e2c2-d526-440a-b8cc-471d6d0af707
# ╠═380acbc4-ce91-4ea1-9bd7-491b42dfc6fe
# ╠═27f71d57-d13a-4c15-a29b-fb2722ec4e57
# ╠═4b8a5946-c6e4-4ba6-826e-083e0a4fcba3
# ╠═865b74e2-8ec9-4ba3-ad79-80acefac0019
# ╠═e7aaef06-a81b-44c7-b8ec-ef05ecfd2001
# ╠═a35b6304-6342-4cb7-b709-06981c5174f8
# ╠═fdf317eb-17eb-40b1-b6d0-f342f806b4c7
# ╠═e564a47e-77e4-4599-9ac6-8aef5b78fd94
# ╠═b5cbba82-a40e-42c0-9cba-9444596c04cc
# ╠═af54a055-54ef-4196-8cf4-9fae59996dbe
# ╠═d1b59bb4-47b9-4b4f-9674-fd4ed9759996
# ╠═efb2b701-e055-4320-8ad1-853cb36a1cdf
# ╠═ce488d89-d599-4d51-a51c-6b50a763b9e4
# ╠═6e351d96-5259-4608-b3b3-a6dae0d0c241
# ╠═b7ed2cc8-3fde-4489-b8a8-01395e6ed6d7
# ╠═570fe0e7-9e1f-4dd6-a669-b7fe235f12d4
# ╠═b869ba37-694b-46ae-936b-af2c882c8fc7
# ╠═93dcee4a-0fa8-4f0d-adcb-6990d312f9f0
# ╠═de4dc2a4-c9c5-4182-a906-be668b6f3f48
# ╠═21e22972-59b2-4f49-8328-aaa64ac45c87
# ╠═4449be15-32ba-4ae7-8bdf-72e6ec4567f7
# ╠═306b57a3-760d-4b81-b3dc-4aa3e025d2f2
# ╠═c801d3e0-4d82-4b0a-b368-c59f208e25cf
# ╠═9823485e-fd71-4e59-bc56-29c704a56e65
# ╠═d26a7cec-f20b-4676-a5d2-0fe568cca3fc
# ╠═ee9418c6-ac2e-42a7-986b-a2a04cb384d5
# ╠═636b9f68-c25c-434a-916b-e632d6931dd8
# ╠═0cf0ac1f-32b3-4ff8-bed2-d1c3304fbde6
# ╠═ebf69a0a-7b27-400a-8e99-8af4ca202570
# ╠═1e51e750-5271-46c8-afec-a98e440b3d16
# ╠═aa5df510-3aa4-43ae-8f43-033f75e812b9
# ╠═bf726737-561f-4ca3-9986-03b4066e7004
# ╠═4c3d2136-7491-4615-85ec-b6d9b9dda45d
# ╠═33d9063a-8afd-4c07-8d26-2873bf6a9d48
# ╠═d73b3587-53e2-4ba6-875c-ecee5210c7dd
# ╠═2f96f9fd-dc69-4ad8-b5d6-7b093d6cf8ad
# ╠═c3d5ca3c-7254-444f-9a36-d60062a890f5
# ╠═94d98795-a309-4d45-b06a-ee269a3df37b
# ╠═60eb389d-9427-4178-9611-59dec0ec63e6
# ╠═283456dc-2415-45c8-a811-cfd93ab23beb
# ╠═8e1ca445-a738-40be-ab4c-f948c18e2f41
# ╠═cdc36c27-df16-4448-9c7e-1ad22cfa4fd2
# ╠═b3b2dc90-edc9-4200-8459-734a2c1dfcda
# ╠═ce6a1868-f8dc-49a9-99cf-9add95abd51e
# ╠═6081e9ba-72df-4ade-801d-d2434de14e57
# ╠═f99dcaff-8367-4a64-be52-7979ab003c20
# ╠═6b7d9922-0b34-4dba-9148-be13a92726d7
# ╠═d9d1345f-c9d9-4caa-9c5b-839fcff75ce0
# ╠═d977411f-dccd-4aad-8251-c570a9ffec32
# ╠═287a403e-ce53-472f-9995-db091d2d21a8
# ╠═98b7e234-194e-44c1-b78e-3507e688eaea
# ╠═50096ed8-e785-4e90-8fa8-1637cec1be43
# ╠═933f7893-81ea-49e2-8c30-0bc6866ce7f5
# ╠═45d26735-360f-4b30-a0eb-437280138d69
# ╠═0c8141b3-8a18-4189-a756-4e925ca28a05
# ╠═cdfe43e2-b6f9-417e-b0a0-a07756e5e8a9
# ╠═013c164a-f4fc-4452-a2ad-d52f4b3c5cd2
# ╠═77ea6942-4053-4fa9-b9c5-e60ae853847f
# ╠═957aece0-fc23-40f3-a6cf-cfb37ea6051a
