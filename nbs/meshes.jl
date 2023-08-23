### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 0e9f569e-be69-11ed-2f95-4ba551eb69fc
begin
	using Pkg; Pkg.activate(".")
	using FileIO, MeshIO
end

# ╔═╡ 4281dfdd-91d6-4bdf-8016-420b072dd21d
using Chain

# ╔═╡ f8cee7ed-cec2-41d2-94c9-9e9ade4e3048
using MultiscaleSimplexSignalTransforms

# ╔═╡ 217833cc-dae0-428c-a5fd-be6231ca4100
using Makie, WGLMakie

# ╔═╡ 6ebaf401-b1ba-45af-bba6-b9cb5026bc79
using Meshes, MeshViz

# ╔═╡ 6a42dbd3-ae89-4d08-a289-c6f324241125
using StatsBase

# ╔═╡ 0b29f261-a28e-43f2-bf4f-6da2f366b1fa
import GeometryBasics: value, Mesh, coordinates, faces

# ╔═╡ 450c5a82-182d-404b-b65e-9cac2f28c5a4
# mymesh = load("/Users/eug/Downloads/meshes/mesh_0165.obj")
mymesh = load("/Users/eug/Downloads/Toysmith_Windem_Up_Flippin_Animals_Dog/meshes/model.obj")

# ╔═╡ 42cff982-f93a-40a5-9208-b54648dffbef
function MultiscaleSimplexSignalTransforms.SimplexTree(mesh::Mesh)
	a = SimplexTree(length(coordinates(mesh)))
	for face ∈ faces(mesh)
		insert!(a, Tuple(value.(face)))
	end
	a
end

# ╔═╡ 2a26fb4e-c1d8-40aa-b008-9642709acf62


# ╔═╡ b12e0d38-f9ce-41d6-b7f8-5fc53a4254fe
tree = SimplexTree(mymesh);

# ╔═╡ 901e0164-47d4-4e44-a385-507f2512dcdc
reg = KRegion(tree, 2; weakadjweight=1.0);

# ╔═╡ 4d80c98c-b81b-4c05-882b-d54f84a1bc60
submpart = kGHWT(reg; normalization=Normalization.Symmetric);

# ╔═╡ af207c7e-8ba6-42ca-8186-9ed1e46b0654
fullpart = kGHWT(reg, SubRepresentation.Full; normalization=Normalization.Symmetric);

# ╔═╡ 878ad117-eac5-4229-bbc6-8ae05cbb0e78
[
	value.(ind)
	for f in faces(mymesh)
	for ind in f
] |> unique |> sort

# ╔═╡ 09cf64b9-6abf-40b1-9041-3c89f94bd7a1
# findfirst(==(mesh[faces(mesh)[1][1]][1]), coordinates(mesh))
@chain mymesh begin
	faces
	_[1]
	@. value
	# coordinates(mesh)[_]
	# @aside @info mesh[1]
end

# ╔═╡ fa319490-e13a-473c-886a-d60df6154802


# ╔═╡ aa64f0b1-2c4d-4a98-b4c1-6f582ff38a49
smesh = Meshes.SimpleMesh(
	Meshes.Point3.(coordinates(mymesh)),
	connect.(Tuple.(broadcast.(value, faces(mymesh))))
)

# ╔═╡ 4ac5aaf6-9d75-4835-b023-75677d26cbdf
Is = (4,7,17)

# ╔═╡ c76431ae-8f1a-4ae1-8853-c8822e0b383b
Makie.loadasset

# ╔═╡ ab7b185f-25e6-4979-a4d9-13b0051e1c78


# ╔═╡ 5e30a6d0-0a09-4cc9-a12c-d3ddece4aca4
subfig = @chain begin
	viz(smesh, color=submpart[Is], axis=(; show_axis=false))
	first
	@aside Label(_[1,1,Top()], "Submatrix Hierarchical Partition (j=$(Is[1]), k=$(Is[2]), l=$(Is[3]))", padding=(0,0,0,0))
	_
end

# ╔═╡ c1288024-92ac-40f9-90f9-58618da52929
dog = load("/Users/eug/Downloads/texture.png");

# ╔═╡ 24aee0bb-aa59-4c7d-98fe-7becbca901ff
fullfig = @chain begin
	Makie.mesh(smesh, color=dog, axis=(; show_axis=false))
	# first
	# @aside Label(_[1,1,Top()], "Full Hierarchical Partition (j=$(Is[1]), k=$(Is[2]), l=$(Is[3]))", padding=(0,0,0,0))
	_
end

# ╔═╡ 70e826d3-72b8-46dc-adb7-253ea22945f6
meshparts(filename) = @chain filename begin
	load
	@aside pltmesh = Meshes.SimpleMesh(
		Meshes.Point3.(coordinates(_)),
		connect.(Tuple.(broadcast.(value, faces(_))))
	)
	SimplexTree
	KRegion(2; weakadjweight=1.0)
	@aside submpart = kGHWT(_; normalization=Normalization.Symmetric)
	kGHWT(_, SubRepresentation.Full; normalization=Normalization.Symmetric)
	(; submpart, fullpart=_, pltmesh)
end

# ╔═╡ d97937cb-7234-4963-9b44-be68dc972306
unitmap(x) = extrema(x) |> mM -> (x.-mM[1]) ./ (mM[2]-mM[1])


# ╔═╡ Cell order:
# ╠═0e9f569e-be69-11ed-2f95-4ba551eb69fc
# ╠═4281dfdd-91d6-4bdf-8016-420b072dd21d
# ╠═f8cee7ed-cec2-41d2-94c9-9e9ade4e3048
# ╠═217833cc-dae0-428c-a5fd-be6231ca4100
# ╠═6ebaf401-b1ba-45af-bba6-b9cb5026bc79
# ╠═0b29f261-a28e-43f2-bf4f-6da2f366b1fa
# ╠═450c5a82-182d-404b-b65e-9cac2f28c5a4
# ╠═42cff982-f93a-40a5-9208-b54648dffbef
# ╠═2a26fb4e-c1d8-40aa-b008-9642709acf62
# ╠═b12e0d38-f9ce-41d6-b7f8-5fc53a4254fe
# ╠═901e0164-47d4-4e44-a385-507f2512dcdc
# ╠═4d80c98c-b81b-4c05-882b-d54f84a1bc60
# ╠═af207c7e-8ba6-42ca-8186-9ed1e46b0654
# ╠═878ad117-eac5-4229-bbc6-8ae05cbb0e78
# ╠═09cf64b9-6abf-40b1-9041-3c89f94bd7a1
# ╠═fa319490-e13a-473c-886a-d60df6154802
# ╠═aa64f0b1-2c4d-4a98-b4c1-6f582ff38a49
# ╠═4ac5aaf6-9d75-4835-b023-75677d26cbdf
# ╠═24aee0bb-aa59-4c7d-98fe-7becbca901ff
# ╠═c76431ae-8f1a-4ae1-8853-c8822e0b383b
# ╠═ab7b185f-25e6-4979-a4d9-13b0051e1c78
# ╠═5e30a6d0-0a09-4cc9-a12c-d3ddece4aca4
# ╠═c1288024-92ac-40f9-90f9-58618da52929
# ╠═70e826d3-72b8-46dc-adb7-253ea22945f6
# ╠═6a42dbd3-ae89-4d08-a289-c6f324241125
# ╠═d97937cb-7234-4963-9b44-be68dc972306
