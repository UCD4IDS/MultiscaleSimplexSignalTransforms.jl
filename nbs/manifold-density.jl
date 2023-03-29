### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 9de11aba-ca92-11ed-1474-192dd67c97c1
using Pkg; Pkg.activate(".")

# ╔═╡ 90f1912c-b304-4312-9177-ad4faaa01cc8
using Revise

# ╔═╡ c9360999-a739-4bd6-b34e-0af5ca5ca5a8
using MultiscaleSimplexSignalTransforms

# ╔═╡ 88f5f36d-bb48-47de-a088-1da5cf023959
using Plots, ForwardDiff, QuadGK, Roots

# ╔═╡ c2a89404-212f-40f1-b4b7-edbae06e956a
2sin.(π.*(rand(1000).-0.5)) |> x -> scatter(
	[x;-x],
	[sqrt.(4 .- x.^2);-sqrt.(4 .- x.^2)],
	aspect_ratio=:equal,
	ms=0.5
)

# ╔═╡ 7e4d01fd-e4f5-47af-b808-bc2f8839936b
rand(1000).*4 .-2 |> x -> scatter(
	[x;-x],
	[sqrt.(4 .- x.^2);-sqrt.(4 .- x.^2)],
	aspect_ratio=:equal,
	ms=0.5
)

# ╔═╡ 716a52fd-7ee6-40f1-b5eb-b50ac3c407e6
a, b = sqrt(0.1), 2

# ╔═╡ a859cb3f-781e-4efa-bc61-294802a1009b
f(x;a=sqrt(0.1), b=2) = sqrt.((a^2 .+ x.^2) .* (b^2 .- x.^2))

# ╔═╡ e7d87f7b-f489-4a1c-9e5c-20d16260e6ff
rand(1000) .* 2b .- b |> x -> scatter(x, f, aspect_ratio=:equal)

# ╔═╡ 072345ae-0213-4369-a3a2-9e56211f6d6a
nx = 100000

# ╔═╡ 1efd89f4-ce4d-4da3-8995-05b0f75cd50a
X = (-b:1/nx:b)[2:end-1]

# ╔═╡ 45d3d135-bfc8-43e2-8edd-e7a619341458
length(X)

# ╔═╡ 29dd81e4-4a50-4cbe-a06e-01d3f4d9baf2
floor(Int, rand()*2b*nx)

# ╔═╡ db990f83-2198-4aeb-84e1-46c4f5a491f3
F = ForwardDiff.derivative.(f, X)

# ╔═╡ b31b7dac-9875-41a7-89a4-dcece4a710eb
Ff(x) = floor(Int, (x+b)*2b*nx) |> i -> i == 0 ? F[1] : i > lastindex(F) ? F[end] : F[i]

# ╔═╡ 706ab776-5b31-47f3-b21b-e7a3e7c5ae1d
ds(x) = sqrt.(1 .+ Ff.(x).^2)

# ╔═╡ 3486abe3-a2b5-47da-9b11-c5260f6c9fe6
s(x) = quadgk(ds, -b, x)[1]

# ╔═╡ 7fd5e93f-45a6-48aa-bc74-75933f1b0c52
sm = s(b)

# ╔═╡ 0541734b-2cdd-4328-ac82-2ca39ac71f55
Finv(x) = find_zero(u -> s(u)-x, (-b, b))

# ╔═╡ d3b054ce-965a-4efc-ac05-6e35e8e9bd83
# seems very slow, 14K seconds wasnt enough
# samples = Finv.(rand(1000))

# ╔═╡ Cell order:
# ╠═9de11aba-ca92-11ed-1474-192dd67c97c1
# ╠═90f1912c-b304-4312-9177-ad4faaa01cc8
# ╠═c9360999-a739-4bd6-b34e-0af5ca5ca5a8
# ╠═88f5f36d-bb48-47de-a088-1da5cf023959
# ╠═c2a89404-212f-40f1-b4b7-edbae06e956a
# ╠═7e4d01fd-e4f5-47af-b808-bc2f8839936b
# ╠═716a52fd-7ee6-40f1-b5eb-b50ac3c407e6
# ╠═a859cb3f-781e-4efa-bc61-294802a1009b
# ╠═e7d87f7b-f489-4a1c-9e5c-20d16260e6ff
# ╠═072345ae-0213-4369-a3a2-9e56211f6d6a
# ╠═1efd89f4-ce4d-4da3-8995-05b0f75cd50a
# ╠═45d3d135-bfc8-43e2-8edd-e7a619341458
# ╠═29dd81e4-4a50-4cbe-a06e-01d3f4d9baf2
# ╠═db990f83-2198-4aeb-84e1-46c4f5a491f3
# ╠═b31b7dac-9875-41a7-89a4-dcece4a710eb
# ╠═706ab776-5b31-47f3-b21b-e7a3e7c5ae1d
# ╠═3486abe3-a2b5-47da-9b11-c5260f6c9fe6
# ╠═7fd5e93f-45a6-48aa-bc74-75933f1b0c52
# ╠═0541734b-2cdd-4328-ac82-2ca39ac71f55
# ╠═d3b054ce-965a-4efc-ac05-6e35e8e9bd83
