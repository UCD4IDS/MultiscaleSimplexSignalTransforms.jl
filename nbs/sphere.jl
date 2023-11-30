### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 1bf3bd34-f03d-11ed-31d2-d5935d3e207d
using Pkg; Pkg.activate(".")

# ╔═╡ 27a51fc5-e9cd-4194-83c2-0f59ae68ea07
using Graphs, MultiscaleSimplexSignalTransforms, Makie, Unzip, Meshes, MeshViz, Arpack

# ╔═╡ 4c524f1b-971c-46e0-9713-056459ecdb1f
using WGLMakie

# ╔═╡ 633a515b-c76d-4662-8eb6-93968aa5eab7
using Chain: @chain

# ╔═╡ cbdd573d-8007-4f64-aa5f-ff949117d4e5
begin
using GeometryTypes
using LinearAlgebra

"""
    subdivide(msh::HomogenousMesh,f::Function)

Returns a subdived triangular mesh from passed mesh `msh` and interpolator function `f`. The interpolator `f` is expected to accept two vertex indicies of the edge and to return a tuple of coordinates of the middle point (as one wishes to define it).
"""
function subdivide(msh::HomogenousMesh,f::Function)
    edges = filter(x->x[1]<x[2],GeometryTypes.decompose(GeometryTypes.Face{2,Int},msh))
    epoint(v1,v2) = findfirst(x->x==(v2>v1 ? GeometryTypes.Face(v1,v2) : GeometryTypes.Face(v2,v1)), edges) + length(msh.vertices)

    newfaces = GeometryTypes.Face{3,Int}[]

    newvertices = copy(msh.vertices)
    resize!(newvertices,length(msh.vertices) + length(edges))

    for (v1,v2,v3) in msh.faces

        ev3 = epoint(v1,v2)
        ev1 = epoint(v2,v3)
        ev2 = epoint(v3,v1)

        ### Usually, does assignment twice but is important if the surface is not tightly connected
        newvertices[ev3] = GeometryTypes.Point(f(v1,v2)...)
        newvertices[ev1] = GeometryTypes.Point(f(v2,v3)...)
        newvertices[ev2] = GeometryTypes.Point(f(v3,v1)...)

        push!(newfaces,GeometryTypes.Face(v1,ev3,ev2))
        push!(newfaces,GeometryTypes.Face(v2,ev1,ev3))
        push!(newfaces,GeometryTypes.Face(v3,ev2,ev1))
        push!(newfaces,GeometryTypes.Face(ev1,ev2,ev3))
    end

    return HomogenousMesh(newvertices,newfaces)
end

"""
    unitsphere(n)

Returns a sphere mesh made out of icosahedron subdivided `n` times.
"""
function unitsphere(subdivissions::Int64)
    t = ( 1 + sqrt( 5 ) ) / 2;

    vertices = GeometryTypes.Point{3,Float64}[
        [ -1,  t,  0 ], [  1, t, 0 ], [ -1, -t,  0 ], [  1, -t,  0 ],
        [  0, -1,  t ], [  0, 1, t ], [  0, -1, -t ], [  0,  1, -t ],
        [  t,  0, -1 ], [  t, 0, 1 ], [ -t,  0, -1 ], [ -t,  0,  1 ]
    ]

    faces = GeometryTypes.Face{3,Int64}[
        [1, 12, 6], [1, 6, 2], [1, 2, 8], [1, 8, 11], [1, 11, 12], [2, 6, 10], [6, 12, 5],
        [12, 11, 3], [11, 8, 7], [8, 2, 9], [4, 10, 5], [4, 5, 3], [4, 3, 7], [4, 7, 9],
        [4, 9, 10], [5, 10, 6], [3, 5, 12], [7, 3, 11], [9, 7, 8], [10, 9, 2]
    ]

    msh = HomogenousMesh(vertices ./ sqrt(1+t^2),faces)

    function sf(msh,v1,v2)
        p1 = msh.vertices[v1]
        p2 = msh.vertices[v2]
        return (p1 + p2)/norm(p1+p2)
    end


    for i in 1:subdivissions
        msh = subdivide(msh,(v1,v2)->sf(msh,v1,v2))
    end

    return msh
end
end

# ╔═╡ 531ba901-27cf-4992-8006-e3e92b2c1f3c
# using GLMakie

# ╔═╡ ec72e9fb-a7ae-4fa6-be05-e39dc81e1f96
# ts = [
# 	Shape(collect(collect.(S.vertices[f])))
# 	for f in S.faces
# ]

# ╔═╡ 8f7ee235-5c06-4d3b-b851-d526e05c6c7b
# plotlyjs()
# # gr()

# ╔═╡ 1e5329e8-7a3e-4a29-8966-24c87dbd6ecd
# @chain begin
# 	S.vertices[S.faces[1]]
# 	@. collect
# 	reduce(hcat, _)
# 	[
# 		_*[a,b,1-(a+b)]
# 		for a ∈ 0:.1:1
# 		for b ∈ 0:.1:1-a
# 	]
# 	@. Tuple
# 	unzip
# 	surface(_..., mz=2)
# end

# ╔═╡ 4b2afde4-b19b-49ef-bee0-7bc22c16ec74
@chain begin
	# [[[1,0,0],[0,1,0],[0,0,1]], [[1,0,0],[1,1,0],[1,0,1]]]
	# [[[1,0,0],[0,1,0],[0,0,1]]]
	[[1,0,0],[0,1,0],[0,0,1]]
	reduce(hcat, _)
	[
		Tuple(_*[a,b,1-(a+b)])
		for a ∈ 0:.01:1
		for b ∈ 0:.01:1-a
	]
	unzip
	surface(_..., )
end

# ╔═╡ 05984623-5fb6-4a2c-8462-818c9ccb703a
@chain begin
	[[1,0,0],[0,1,0],[0,0,1], [1, 1, 1]]
	@. Tuple
	unzip
	surface(_...)
end

# ╔═╡ e989aef5-72a8-4758-a720-10eecc7d1b07


# ╔═╡ 4bbd0db0-bf34-49de-9de7-33662aa00bf1


# ╔═╡ 3ed9a4a4-dff6-4d7f-96d9-21f6e49f8fe5


# ╔═╡ 321d56ca-b530-4078-8606-bb1e287925c9
S = unitsphere(2)

# ╔═╡ c50271e1-d29f-498d-9120-e566ae17bdb3
typeof(S)

# ╔═╡ 685a1e66-3dea-49bc-8ec6-5ab0b514ad22
complex = @chain begin
	Set(
		sort(f[i])
		for f ∈ S.faces
		for i ∈ [[1,2], [1,3], [2,3]]
	)
	collect
	@. Tuple
	@. Edge
	Graph
	cliquecomplex(4)
end

# ╔═╡ fa41e924-649c-412a-ab8f-9ed0d2696c20
region = KRegion(complex, 2; weakadjweight=1.0)

# ╔═╡ 669cf3a3-c9bf-4687-96dc-98d50d44bc96
part = kGHWT(region)

# ╔═╡ d759307e-0270-4a12-bb9b-b5f2294ce6af
Λ, X = eigs(k_laplacian(region), nev=100, which=:SM)

# ╔═╡ 886c6342-f114-40bb-9413-e62d81e44c7a
Λm, Xm = eigs(k_laplacian(region), nev=100)

# ╔═╡ 94eb788c-e4f6-42ea-8d09-4a90cfc14e76
# plot(X[:,3] .* sign.(X[:,2]))

# ╔═╡ 616020d7-c99d-41b2-b946-c9b0fa0cda29
# plot((X[:,1]))

# ╔═╡ ee5eda6e-cf98-4d2a-a17f-8137d26be02d
# plot(Xm[:,6])

# ╔═╡ a7e2a5bd-9a87-4cd8-a140-38f9c8e260a5
# plot(Λ)

# ╔═╡ cb48ea58-29e9-4dc4-8053-152e5f3cad36
# scatter3d([getindex.(S.vertices, i) for i in 1:3]..., ms=1, mz=part[1,0,1], legend=false)

# ╔═╡ d9c6dc5e-3c07-4c83-bd67-5e098d46a61f
# plot(L)

# ╔═╡ Cell order:
# ╠═1bf3bd34-f03d-11ed-31d2-d5935d3e207d
# ╠═27a51fc5-e9cd-4194-83c2-0f59ae68ea07
# ╠═4c524f1b-971c-46e0-9713-056459ecdb1f
# ╠═531ba901-27cf-4992-8006-e3e92b2c1f3c
# ╠═633a515b-c76d-4662-8eb6-93968aa5eab7
# ╠═cbdd573d-8007-4f64-aa5f-ff949117d4e5
# ╠═ec72e9fb-a7ae-4fa6-be05-e39dc81e1f96
# ╠═8f7ee235-5c06-4d3b-b851-d526e05c6c7b
# ╠═1e5329e8-7a3e-4a29-8966-24c87dbd6ecd
# ╠═4b2afde4-b19b-49ef-bee0-7bc22c16ec74
# ╠═05984623-5fb6-4a2c-8462-818c9ccb703a
# ╠═e989aef5-72a8-4758-a720-10eecc7d1b07
# ╠═4bbd0db0-bf34-49de-9de7-33662aa00bf1
# ╠═3ed9a4a4-dff6-4d7f-96d9-21f6e49f8fe5
# ╠═321d56ca-b530-4078-8606-bb1e287925c9
# ╠═c50271e1-d29f-498d-9120-e566ae17bdb3
# ╠═685a1e66-3dea-49bc-8ec6-5ab0b514ad22
# ╠═fa41e924-649c-412a-ab8f-9ed0d2696c20
# ╠═669cf3a3-c9bf-4687-96dc-98d50d44bc96
# ╠═d759307e-0270-4a12-bb9b-b5f2294ce6af
# ╠═886c6342-f114-40bb-9413-e62d81e44c7a
# ╠═94eb788c-e4f6-42ea-8d09-4a90cfc14e76
# ╠═616020d7-c99d-41b2-b946-c9b0fa0cda29
# ╠═ee5eda6e-cf98-4d2a-a17f-8137d26be02d
# ╠═a7e2a5bd-9a87-4cd8-a140-38f9c8e260a5
# ╠═cb48ea58-29e9-4dc4-8053-152e5f3cad36
# ╠═d9c6dc5e-3c07-4c83-bd67-5e098d46a61f
