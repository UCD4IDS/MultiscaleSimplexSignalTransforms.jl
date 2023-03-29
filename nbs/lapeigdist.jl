### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ ce59fbbf-c630-4ce9-8a7e-4249e9e3575b
begin
	using Pkg; Pkg.activate(".")
end

# ╔═╡ 4c5b1af6-52f4-4c2c-bc46-9c5bc77c5f88
using Revise

# ╔═╡ 6bf34b72-bdbd-11ed-393c-1194451fa347
using MultiscaleSimplexSignalTransforms

# ╔═╡ d3410e3c-d1ef-407a-afd8-1d17556d5d36
using Graphs, Plots, Chain, Arpack, LinearAlgebra, StatsBase, Distances, Unzip, SimpleWeightedGraphs, Roots, Polynomials

# ╔═╡ 5f9b8c79-6c7d-491f-a6c5-874bcbe68ca7
using LaTeXStrings

# ╔═╡ d94137f1-91c8-4c7f-b1cc-1e4424ee540c
function manifold(n, ϵ=.05, a=0.1, b=4)
	t = rand(n)*2π
	ss, cs = unzip(sincos.(t))
	r = sqrt.(maximum.(roots.([
		Polynomial([cs[i]^4, ss[i]^2 - (b-a)*cs[i]^2, -a*b])
		for i ∈ 1:n
	])))
	x = r .* cs
	y = r .* ss
	# x = sqrt(b) * (2rand(n) .- 1)
	# s = rand((-1, 1), n)
	# y = @. s * sqrt((a+x^2)*(b-x^2))
	noise = rand(n,2) * ϵ
	[x y] + noise
end

# ╔═╡ db2e9f6e-fe6f-4cd1-b3ca-2cf4fa1bdf6c
x, y = collect.(eachcol(manifold(1500)))

# ╔═╡ 9a0cadd7-3c2c-4153-851c-ef3fd500ee09
θ = atan.(y,x)

# ╔═╡ fe618802-df36-45be-a316-80d7d57272d9
embedmanifold = scatter(
	x,y, 
	marker_z=θ,
	color=:rainbow,
	legend=false, title="Original manifold samples"
)

# ╔═╡ 1907060f-8d6c-4732-a87d-4048a84803ea
# savefig(embedmanifold, "/Users/eug/Desktop/embed-manifold.png")

# ╔═╡ 91e97c52-6a05-45ea-9f50-b577245934fe
inds = sortperm(atan.(y, x))

# ╔═╡ 16143238-c8ea-4e69-9d2e-aec5f9f29951
D = pairwise(SqEuclidean(), [x y], dims=1)

# ╔═╡ 0afb2716-60d9-404e-ab38-69ed0f575b4b
W = map(x -> x > 1e-6 ? x : 0.0, exp.(-D ./ .01) - I)

# ╔═╡ 7b467d99-600e-494e-88c2-53b4fc247545
Ω = @chain W begin
	sum(dims=2)
	_[:]
	Diagonal
	_^-0.5
end

# ╔═╡ f2c1eb07-64bd-48cb-8834-36cb06404931
L = Diagonal(sum(W, dims=2)[:])^-0.5 |> D -> I - D*W*D |> Symmetric;

# ╔═╡ 19d4c312-d75d-47dd-a0f2-c5392433be6a
X = eigs(I+L, nev=3, which=:SM, maxiter=10000)[2][:,2:3]

# ╔═╡ 765d241e-f40d-46f0-9ea1-1d329e274a96
Q = eigen(L).vectors

# ╔═╡ dc3915ba-c919-4859-884c-3df4a0fe22b1
K = distance_kernel(ZeroRegion(W); normalization=Normalization.Symmetric, dmax=1e2)

# ╔═╡ be2a4c89-2791-43d7-aaef-70104cd56a78
Λ, Y = eigs(K, nev=3, which=:LR)

# ╔═╡ 3c869f09-cfd5-4447-a08b-83d02977f430
usesym = false

# ╔═╡ 17c9e98d-e096-4c10-befc-e532b138d8a7
titleword, normM = usesym ?
	("symmetric", I) : 
	("random walk", Ω)

# ╔═╡ 60e45070-1aa6-448f-9bd5-a2bceebf4111
lapfig = scatter(eachcol(normM*X[:,1:2])..., marker_z=θ, color=:rainbow, colorbar=false, legend=false, title="L-based embedding, $titleword normalization", ms=4, msw=0, aspect_ratio=:equal)

# ╔═╡ 8da36a35-5a38-49b1-88e7-4d251336b6dd
distfig = scatter(eachcol(normM*Y[:,1:2])..., marker_z=θ, color=:rainbow, colorbar=false, legend=false, title="H-based embedding, $titleword normalization", ms=2.5, msw=0, aspect_ratio=:equal)

# ╔═╡ 5faecb52-013d-42de-aa1e-82aa06d7aa76
scatter3d(
	eachcol(normM*Y[:,1:3])...,
	marker_z=θ, color=:rainbow, colorbar=false, legend=false, title="H-based embedding, $titleword normalization", ms=2, msw=0,
	aspect_ratio=:equal
)

# ╔═╡ e1d0f1b7-198c-4716-98e5-3f103a2636e0
# savefig(distfig, "~/Desktop/embed-rwdist.png")

# ╔═╡ Cell order:
# ╠═ce59fbbf-c630-4ce9-8a7e-4249e9e3575b
# ╠═4c5b1af6-52f4-4c2c-bc46-9c5bc77c5f88
# ╠═6bf34b72-bdbd-11ed-393c-1194451fa347
# ╠═d3410e3c-d1ef-407a-afd8-1d17556d5d36
# ╠═5f9b8c79-6c7d-491f-a6c5-874bcbe68ca7
# ╠═d94137f1-91c8-4c7f-b1cc-1e4424ee540c
# ╠═db2e9f6e-fe6f-4cd1-b3ca-2cf4fa1bdf6c
# ╠═9a0cadd7-3c2c-4153-851c-ef3fd500ee09
# ╠═fe618802-df36-45be-a316-80d7d57272d9
# ╠═1907060f-8d6c-4732-a87d-4048a84803ea
# ╠═91e97c52-6a05-45ea-9f50-b577245934fe
# ╠═16143238-c8ea-4e69-9d2e-aec5f9f29951
# ╠═0afb2716-60d9-404e-ab38-69ed0f575b4b
# ╠═7b467d99-600e-494e-88c2-53b4fc247545
# ╠═f2c1eb07-64bd-48cb-8834-36cb06404931
# ╠═19d4c312-d75d-47dd-a0f2-c5392433be6a
# ╠═765d241e-f40d-46f0-9ea1-1d329e274a96
# ╠═dc3915ba-c919-4859-884c-3df4a0fe22b1
# ╠═be2a4c89-2791-43d7-aaef-70104cd56a78
# ╠═3c869f09-cfd5-4447-a08b-83d02977f430
# ╠═17c9e98d-e096-4c10-befc-e532b138d8a7
# ╠═60e45070-1aa6-448f-9bd5-a2bceebf4111
# ╠═8da36a35-5a38-49b1-88e7-4d251336b6dd
# ╠═5faecb52-013d-42de-aa1e-82aa06d7aa76
# ╠═e1d0f1b7-198c-4716-98e5-3f103a2636e0
