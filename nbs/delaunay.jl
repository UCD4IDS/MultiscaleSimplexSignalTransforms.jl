### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 5196acf0-9172-11ed-3894-83abfa3491a1
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ be2fc790-f600-433b-bf3b-a122684a88b9
using Revise

# ╔═╡ f1c6a19e-2879-4b55-8989-9f2ca1f3721a
using MultiscaleSimplexSignalTransforms

# ╔═╡ 9bb0c156-1024-465f-825e-e05dfc7c66cf
using VoronoiDelaunay

# ╔═╡ bf2a8642-4876-4694-9e86-d9f70a4aa563
using HTTP

# ╔═╡ e92c0793-3e2c-40ab-b838-b89317f53f35
using Plots

# ╔═╡ 38568ca2-2cbd-47ff-aa82-66d7d11139de
using Distances

# ╔═╡ 4e5c62bd-bb53-432b-b486-c73b5d999d61
using Graphs

# ╔═╡ 499124ad-8ce0-4eda-b02e-9e3e19979419
using Combinatorics, LinearAlgebra, Arpack

# ╔═╡ f0780b5e-4dec-45b5-9304-e1a28b5bcc78
using GraphRecipes

# ╔═╡ 8f9a970e-fb95-46f4-ad4a-e4e4f821d9cd
using SimpleWeightedGraphs

# ╔═╡ 72ca230f-a40f-4cca-bf9b-329d0ebfacba
using MAT

# ╔═╡ 4c16a1cc-63ea-41a4-ae6d-b1a5c0fcb6d1
using StatsBase

# ╔═╡ 2ff3a4a7-b7e1-487c-816b-31492b6c01f0
begin
	import VoronoiDelaunay: from_image, voronoiedges, getplotxy
	import Images: load
end

# ╔═╡ 50ee8de1-1e27-4559-81a6-74df60da7472
img = HTTP.get("https://raw.githubusercontent.com/JuliaGeometry/VoronoiDelaunay.jl/master/examples/julia.png").body |> IOBuffer |> load

# ╔═╡ 5939630b-e551-453e-b03d-0df76a2b26e6
emojimg = load("/Users/eug/Pictures/emoji/emoji-1000.png")

# ╔═╡ 2342740d-1cf7-4f48-bd2b-6bac2c52dc70
begin
	intensity(c::RGB)  = c.b
	intensity(c::RGBA) = c.b
	intensity(c)       = getfield(c, 1) # Workaround. Gray needs to be imported from images, which would take to long.
end

# ╔═╡ 3e7c27e3-01b2-4145-9b9f-10162522e10d
# Create DelaunayTessellation with npts points from an image
function my_from_image(img, npts)

    # placing points in places that represent the image
    pts = Point2D[]
    for i in 1:npts
        x = rand()
        y = rand()
        if intensity(img[Int64(floor(x * size(img)[1])) + 1, Int64(floor(y * size(img)[2])) + 1]) > 0.1
            # if rand() < 0.100
            #     push!(pts, Point2D(1.0 + rand(), 1.0 + rand()))
            # end
            continue
        end
        push!(pts, Point2D(1.0 + x, 2.0 - y))
    end

    tess = DelaunayTessellation(npts)
    push!(tess, pts)

    tess
end

# ╔═╡ 2d4d1fdc-e741-48a4-abdb-dc271bdf50cf
tess = my_from_image(img, 3000)

# ╔═╡ 241431ef-f8ca-4c69-91ae-a24a3e7cb859
begin
	x, y = getplotxy(delaunayedges(tess))
	xpt = filter(!isnan, unique(x))
	ypt = filter(!isnan, unique(y))
	plot(-y, -x, legend=false, ylims=(-2,-1), xlims=(-2,-1))
	scatter!(-ypt, -xpt, ms=2)
	# p = plot(x=x, y=y, Geom.path, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.25, maxvalue=1.75))
end

# ╔═╡ b554ff69-dc4b-4997-908f-3e78a1d6b6a8
ds = pairwise(Euclidean(), [xpt ypt]; dims=1)

# ╔═╡ 1ac4d417-2f6f-489a-9025-f82350018845
first(tess)

# ╔═╡ a14bfd1d-d04a-4e09-8f92-cbbe8b469e46
Base.length(q::DelaunayTessellation2D) = length(q._trigs)

# ╔═╡ 7dc88a13-7f59-4418-b8e1-de5f958bd17f
begin
	wt(simp) = length(simp) == 1 ?
		1.0 :
		1 / sum(ds[ab...] for ab ∈ combinations(simp, 2))
	wt(e::Edge) = wt(Tuple(e))
end

# ╔═╡ 49daee56-fdb5-491c-961b-9e807108eea9
begin
	tris = NTuple{3, Int}[]
	for tri in tess
		push!(tris, Tuple(
			findfirst(i -> x._x == xpt[i] && x._y == ypt[i], 1:length(xpt))
			for x ∈ (geta(tri), getb(tri), getc(tri))			
		))
	end
	tris
end

# ╔═╡ 8dbb29bf-a69a-4d8d-b96b-149aed483bc5
begin
	a = SimplexTree(length(xpt))
	for tri in tris
		insert!(a, tri)
	end
	# a = cliquecomplex(a, 3);
	nothing
end

# ╔═╡ 94e62d94-4cb2-425a-9b18-41fb06f18b66
k = 2

# ╔═╡ 2eb4bd06-2f02-40f8-b6be-66abe9c287aa
edgereg = KRegion(a, 1);

# ╔═╡ d97f62bc-7f54-4743-936e-812e2ca9d336
DiGraph(triu(adjacency_matrix(a)))

# ╔═╡ 42eedc71-b109-407d-ae22-4ef2e685376f
function oriented(a, oris)
	g = DiGraph(length(MultiscaleSimplexSignalTransforms.simplices(a, 0)))
	for (i, e) ∈ enumerate(simplices(a, 1))
		add_edge!(
			g,
			(oris[i] ≥ 0 ? e : reverse(e))...
		)
	end
	g
end

# ╔═╡ e9fce938-aa53-4d30-8441-731f13506ec7
partvec = eigs(
	k_laplacian(
		edgereg;
	)+I;
	nev=1, which=:SM, maxiter=10000
)[2][:,1]

# ╔═╡ 98edb45b-9e42-4558-88ea-106e718c6131
splot(edgereg, sign.(partvec); k=1, xs=-ypt, ys=-xpt, size=1); scatter!(
	-ypt, -xpt, mc=1:length(xpt), palette=:redblue, ms=3
)

# ╔═╡ a37820a0-f5da-4688-8d4f-9b7da51608e9
[sum(adjacency_matrix(oriented(a, partvec)), dims=i)[:] for i=1:2] |> ((a,b),) -> a ./ (a .+ b) |> plot

# ╔═╡ ccc24fce-9b2d-48d5-a0a6-ca693b07b68e
plot(oriented(a, partvec); x=-ypt, y=-xpt, curves=false)

# ╔═╡ da631ece-e78d-49d4-ae43-f44896fa07d9
juliareg = KRegion(a, k; weakadjweight=(k==1 ? 0.0 : 1.0));

# ╔═╡ 0b6f0b5a-f8d5-41b6-a2c1-1431164304c9
jpart = kGHWT(
	juliareg;
	normalization=Normalization.Symmetric
);

# ╔═╡ 0dc92a9e-50de-40ba-be2a-28d385c6a9c2
fpart = kGHWT(
	juliareg, SubRepresentation.Full;
	normalization=Normalization.Symmetric
);

# ╔═╡ 85ef6cc8-4cfa-4a15-98d3-53717de5a349
# plot([
# 	splot(
# 		juliareg, 
# 		treecolors(jpart.part.root, maxlevel=i);
# 		# ones(length(n_simp(juliareg))) .|> Float64; 
# 		k, xs = -ypt, ys = -xpt,
# 		palette = :viridis,
# 		size=1
# 	)
# 	for i=0:3
# ]...)

# ╔═╡ 2b22f86e-1712-4559-9978-fd8b91ab1180
plot([
	splot(
		juliareg,
		v;
		k, xs = -ypt, ys = -xpt,
		palette = :redblue,
		size=1
	)
	for Is ∈ [(0,0,10),]
	for v ∈ [jpart[Is], fpart[Is]]
]...)

# ╔═╡ 5c8aefe5-834b-4fb6-88c3-17a9ef05f2c2
jpart[0,0,40]

# ╔═╡ 5e8fc1ba-622d-47d0-8ecb-4f68e549cd9a
map(
	x -> x[2],
	filter(xx -> xx[1]%3==0, collect(enumerate(x)))
) .|> isnan |> all

# ╔═╡ ada98a3e-e864-4609-b973-b0b8ff53d892
filter(!isnan, unique(x))

# ╔═╡ 53f2ad37-e676-4feb-b26d-28a5470e6834
pts = [
	x
	for e ∈ delaunayedges(tess)
	for x ∈ (geta(e), getb(e))
] |> unique .|> pt -> (-pt._y, -pt._x)

# ╔═╡ 260ee4df-8162-46db-89d8-630bfaca30b1
# g = SimpleGraphFromIterator(Edge.([(1,2), (2, 3), (3, 4)]))
g = SimpleWeightedGraph([1,2,3], [2, 3, 4], [5, 6, 7])

# ╔═╡ f8980019-124d-4bb3-9c80-05609e4f153e
# ZeroRegion(g; hullweightfn = _ -> rand()).weights.weights
Graphs.weights(ZeroRegion(g).weights)

# ╔═╡ ca998b74-4056-479b-8830-bc099669e3c7
ifelse.(Graphs.weights(ZeroRegion(g).weights).==0.0, 0.0, 1.0)

# ╔═╡ 5d4b20f5-127a-468c-b446-fa96ea69087e
# groundtruth = matread("/Users/eug/Downloads/Indian_pines_gt.mat")["indian_pines_gt"]

# ╔═╡ e4717f61-8a9c-4750-91d6-7d63bf61bcd6
# dt = matread("/Users/eug/Downloads/Indian_pines_corrected.mat")["indian_pines_corrected"]

# ╔═╡ 69ff5869-8677-41da-8321-f9d830d51073
# function tensor2sc(mat)
# 	szlims = (0, cumsum(size(mat))[1:end-1]...)
# 	sc = SimplexTree()
# 	for is ∈ CartesianIndices(mat)
# 		insert!(sc, Tuple(is) .+ szlims)
# 	end
# 	KRegion(sc, length(szlims)-2; hullweightfn = s -> mat[(s .- szlims)...]
# 	)
# 	# KRegion(sc, length(szlims)-1)
# end

# ╔═╡ 6f1d47f1-dcdc-4626-a011-6ed63d900f50
# dt_region = tensor2sc(dt)

# ╔═╡ 0afd5e87-8941-4ead-b605-e3a23ebebb47
# dt_part = SubmatrixPartition(dt_region)

# ╔═╡ 76d43003-0666-4019-a8d5-76c91f6af567
# basis_dictionary!(kGHWT(dt_part))

# ╔═╡ Cell order:
# ╠═5196acf0-9172-11ed-3894-83abfa3491a1
# ╠═be2fc790-f600-433b-bf3b-a122684a88b9
# ╠═f1c6a19e-2879-4b55-8989-9f2ca1f3721a
# ╠═9bb0c156-1024-465f-825e-e05dfc7c66cf
# ╠═2ff3a4a7-b7e1-487c-816b-31492b6c01f0
# ╠═bf2a8642-4876-4694-9e86-d9f70a4aa563
# ╠═50ee8de1-1e27-4559-81a6-74df60da7472
# ╠═5939630b-e551-453e-b03d-0df76a2b26e6
# ╠═2d4d1fdc-e741-48a4-abdb-dc271bdf50cf
# ╠═e92c0793-3e2c-40ab-b838-b89317f53f35
# ╠═2342740d-1cf7-4f48-bd2b-6bac2c52dc70
# ╠═3e7c27e3-01b2-4145-9b9f-10162522e10d
# ╠═241431ef-f8ca-4c69-91ae-a24a3e7cb859
# ╠═38568ca2-2cbd-47ff-aa82-66d7d11139de
# ╠═b554ff69-dc4b-4997-908f-3e78a1d6b6a8
# ╠═4e5c62bd-bb53-432b-b486-c73b5d999d61
# ╠═7dc88a13-7f59-4418-b8e1-de5f958bd17f
# ╠═499124ad-8ce0-4eda-b02e-9e3e19979419
# ╠═1ac4d417-2f6f-489a-9025-f82350018845
# ╠═a14bfd1d-d04a-4e09-8f92-cbbe8b469e46
# ╠═49daee56-fdb5-491c-961b-9e807108eea9
# ╠═8dbb29bf-a69a-4d8d-b96b-149aed483bc5
# ╠═94e62d94-4cb2-425a-9b18-41fb06f18b66
# ╠═2eb4bd06-2f02-40f8-b6be-66abe9c287aa
# ╠═d97f62bc-7f54-4743-936e-812e2ca9d336
# ╠═42eedc71-b109-407d-ae22-4ef2e685376f
# ╠═e9fce938-aa53-4d30-8441-731f13506ec7
# ╠═98edb45b-9e42-4558-88ea-106e718c6131
# ╠═f0780b5e-4dec-45b5-9304-e1a28b5bcc78
# ╠═a37820a0-f5da-4688-8d4f-9b7da51608e9
# ╠═ccc24fce-9b2d-48d5-a0a6-ca693b07b68e
# ╠═da631ece-e78d-49d4-ae43-f44896fa07d9
# ╠═0b6f0b5a-f8d5-41b6-a2c1-1431164304c9
# ╠═0dc92a9e-50de-40ba-be2a-28d385c6a9c2
# ╠═85ef6cc8-4cfa-4a15-98d3-53717de5a349
# ╠═2b22f86e-1712-4559-9978-fd8b91ab1180
# ╠═5c8aefe5-834b-4fb6-88c3-17a9ef05f2c2
# ╠═5e8fc1ba-622d-47d0-8ecb-4f68e549cd9a
# ╠═ada98a3e-e864-4609-b973-b0b8ff53d892
# ╠═53f2ad37-e676-4feb-b26d-28a5470e6834
# ╠═8f9a970e-fb95-46f4-ad4a-e4e4f821d9cd
# ╠═260ee4df-8162-46db-89d8-630bfaca30b1
# ╠═f8980019-124d-4bb3-9c80-05609e4f153e
# ╠═ca998b74-4056-479b-8830-bc099669e3c7
# ╠═72ca230f-a40f-4cca-bf9b-329d0ebfacba
# ╠═5d4b20f5-127a-468c-b446-fa96ea69087e
# ╠═4c16a1cc-63ea-41a4-ae6d-b1a5c0fcb6d1
# ╠═e4717f61-8a9c-4750-91d6-7d63bf61bcd6
# ╠═69ff5869-8677-41da-8321-f9d830d51073
# ╠═6f1d47f1-dcdc-4626-a011-6ed63d900f50
# ╠═0afd5e87-8941-4ead-b605-e3a23ebebb47
# ╠═76d43003-0666-4019-a8d5-76c91f6af567
