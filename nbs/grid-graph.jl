### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 954c2b26-c649-11ed-153e-d30afc8bd148
using Pkg; Pkg.activate(".")

# ╔═╡ 4ede4804-50cc-4b7c-b2bc-e508f4ca168c
using Revise

# ╔═╡ c1555624-a4c5-4a62-b446-2f646a82bf38
using MultiscaleSimplexSignalTransforms

# ╔═╡ 085588fd-852d-4559-9375-7c90fc4019c5
using NetworkLayout, Graphs, LinearAlgebra, SparseArrays, Unzip, Plots, Chain, GraphRecipes

# ╔═╡ 4d8b5fcd-c752-427d-91ef-01d370e2b9b5
r, c = 11, 17

# ╔═╡ aaf7ae27-a757-4774-8dba-b98b2d56d89c
n = r*c

# ╔═╡ dddf050a-dcd8-4501-9abb-30deeeaa6cb5
G = Graphs.grid([r,c])

# ╔═╡ c15695d2-6762-479d-a802-bd2def986401
vs = eigen(Matrix(distmat(G))).values

# ╔═╡ a0b9d35f-1638-499d-a496-eecd8ec3fa51
(sum(abs.(vs) .< 1e-10 ), r*c, (r-1)*(c-1)+1)

# ╔═╡ 5fd7d44f-2b58-4bdf-9292-fe15e3ce3eef
plot(eigen(distmat(G)|>Matrix).values, legend=false)

# ╔═╡ e04c921a-a7fe-406a-bfd4-f22f9c6d3b1b
D = distance_kernel(G)

# ╔═╡ 43396ea3-0cd9-4bc0-a3e9-e51b7c2e4134
L = k_laplacian(ZeroRegion(G));

# ╔═╡ 59851beb-9f94-4b52-ad16-60c453d429fa
Leig = eigen(Matrix(L));

# ╔═╡ 36b9f698-977a-43b8-80ac-6ffb266046d6
LX = Leig.vectors .* sign.(Leig.vectors[1,:])';

# ╔═╡ 2e03e373-bbae-4733-8978-ebcd59ad290a
Deig = eigen(Symmetric(D));

# ╔═╡ 14b2e5f1-4cad-48e3-a349-4febc6ef2b69
DX = Deig.vectors .* sign.(Deig.vectors[1,:])';

# ╔═╡ 23a8ba98-e0e1-401c-8384-0ac6fb2bbd92
plot([Leig.values (Deig.values[[1; n:-1:2]] .|> v -> abs(v)<1e-1 ? v : 1/v)])

# ╔═╡ 4ad20b3c-b8cc-4aa8-b7c3-1fea90975a4e
Deig.values |> reverse

# ╔═╡ 2cf1c772-f542-4565-b87b-4a355ed29a3b
eigind = 1

# ╔═╡ 32a249cc-9a27-4831-8e24-3475ba049d85
plot([LX[:,2+eigind] DX[:,end-eigind]])
# plot([Deig.vectors[:,end]])

# ╔═╡ 0d59bad4-6726-46f1-82d3-3cb778467481
asunit(v) = extrema(v) |> mM -> (v .- mM[1]) ./ (mM[2]-mM[1])

# ╔═╡ 2154c61a-e8cd-47f3-a510-217314a5c612
plot(G; curves=false, nodecolor=cgrad(:redblue)[asunit(DX[:,end-eigind])], nodesize=0.4, msw=0, shape=:circle)

# ╔═╡ 3679705a-0983-4cb9-b4d6-6a7ca0dbf081
asunit(LX[:, 2+eigind])

# ╔═╡ Cell order:
# ╠═954c2b26-c649-11ed-153e-d30afc8bd148
# ╠═4ede4804-50cc-4b7c-b2bc-e508f4ca168c
# ╠═c1555624-a4c5-4a62-b446-2f646a82bf38
# ╠═085588fd-852d-4559-9375-7c90fc4019c5
# ╠═aaf7ae27-a757-4774-8dba-b98b2d56d89c
# ╠═4d8b5fcd-c752-427d-91ef-01d370e2b9b5
# ╠═dddf050a-dcd8-4501-9abb-30deeeaa6cb5
# ╠═c15695d2-6762-479d-a802-bd2def986401
# ╠═a0b9d35f-1638-499d-a496-eecd8ec3fa51
# ╠═5fd7d44f-2b58-4bdf-9292-fe15e3ce3eef
# ╠═e04c921a-a7fe-406a-bfd4-f22f9c6d3b1b
# ╠═43396ea3-0cd9-4bc0-a3e9-e51b7c2e4134
# ╠═59851beb-9f94-4b52-ad16-60c453d429fa
# ╠═36b9f698-977a-43b8-80ac-6ffb266046d6
# ╠═2e03e373-bbae-4733-8978-ebcd59ad290a
# ╠═14b2e5f1-4cad-48e3-a349-4febc6ef2b69
# ╠═23a8ba98-e0e1-401c-8384-0ac6fb2bbd92
# ╠═4ad20b3c-b8cc-4aa8-b7c3-1fea90975a4e
# ╠═2cf1c772-f542-4565-b87b-4a355ed29a3b
# ╠═32a249cc-9a27-4831-8e24-3475ba049d85
# ╠═2154c61a-e8cd-47f3-a510-217314a5c612
# ╠═0d59bad4-6726-46f1-82d3-3cb778467481
# ╠═3679705a-0983-4cb9-b4d6-6a7ca0dbf081
