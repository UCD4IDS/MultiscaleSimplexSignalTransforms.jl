### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ d8cba18c-2f63-11ee-062b-337a52a86175
using Pkg;
Pkg.activate(".");

# ╔═╡ 8b08d2ea-9c62-40a4-b1d7-f2fd86a6118e
using Revise

# ╔═╡ 2d1ddd7f-f3dc-4310-bfb1-84a197e32dc9
using LinearAlgebra

# ╔═╡ 04f67a57-3a5d-479b-85d2-f2f3deac35e3
using MultiscaleSimplexSignalTransforms, Graphs, Arpack

# ╔═╡ bead31a0-1589-4517-bea4-943128676de7
n = 100

# ╔═╡ 0a4c2fb0-dbf8-440a-9d9e-a0dd666377f2


# ╔═╡ ef5f18e2-8448-4edf-b0aa-2d9679bf3109
g = path_graph(n)

# ╔═╡ c8a2b78e-5025-4f4e-b7ca-1fa3e0606957


# ╔═╡ 4d1715c9-e7d1-4b8f-8f0c-c43a949eb80d
region = ZeroRegion(g)

# ╔═╡ d26bf1d9-15a9-4486-959d-ca0cc65b9c84
L = k_laplacian(region)

# ╔═╡ 7b8f5f94-3c8b-4b35-babe-d51f871fcefc
pl = f -> splot(
    region, f;
    xs=collect(Float64.(1:n)),
    ys=zeros(n),
    palette=:redblue,
    size=10
)

# ╔═╡ b8791d74-0008-455b-aaef-99118f0571be
pls = M -> splot_many(
    region, collect(eachcol(M))...;
    xs=collect(Float64.(1:n)),
    ys=zeros(n),
    palette=:redblue,
    size=10
)

# ╔═╡ a9102099-b4eb-4842-8c18-bb1761243c7f
Λ, X = eigs(L, which=:SM)

# ╔═╡ 8e134977-e7b8-4adb-8a0b-780b08d6a687
pls(X)

# ╔═╡ 6a38b860-c16e-4329-bb74-9c8ea319ce14
ntree = 4

# ╔═╡ 4495899b-3bb6-4691-91c8-9fbc94df8910
function tree_coor(n)
    xs = Iterators.flatten([[
            range(start=0.5 + 2^(n - k - 1), step=2^(n - k), length=2^(k - 1))
            for k in 1:(n-1)
        ]..., 1:2^(n-1)])
    ys = Iterators.flatten([
        fill(Float64(n - k + 1), 2^(k - 1))
        for k in 1:n
    ])
    collect.([xs, ys])
end

# ╔═╡ 01dabf78-fcaa-4cb8-940d-fd787f551d66
bt = binary_tree(ntree)

# ╔═╡ 33ec4aa6-d075-422d-8ebf-df81ee6a68c7
btr = KRegion(bt, 1, weakadjweight=1.0)

# ╔═╡ e538c15a-5b04-4089-8a64-fa9d347831e4
bt_msst = kGHWT(btr)

# ╔═╡ a96f5ac8-b02c-49e6-97fa-e5b6287113b4
splot_many(btr, Pos(tree_coor(ntree)...), [bt_msst[0, 0, i] for i in 0:5]..., palette=:redblue)

# ╔═╡ 42cb7f68-3483-4041-b9b8-d3c6f292f30c
eigen(Matrix(k_laplacian(btr)))

# ╔═╡ c38ab2cf-b7d1-4f3e-8aa7-360cd1aae8fb


# ╔═╡ 960e3093-bd76-4063-88f5-b53619b886f4


# ╔═╡ Cell order:
# ╠═d8cba18c-2f63-11ee-062b-337a52a86175
# ╠═2d1ddd7f-f3dc-4310-bfb1-84a197e32dc9
# ╠═bead31a0-1589-4517-bea4-943128676de7
# ╠═8b08d2ea-9c62-40a4-b1d7-f2fd86a6118e
# ╠═04f67a57-3a5d-479b-85d2-f2f3deac35e3
# ╠═0a4c2fb0-dbf8-440a-9d9e-a0dd666377f2
# ╠═ef5f18e2-8448-4edf-b0aa-2d9679bf3109
# ╠═c8a2b78e-5025-4f4e-b7ca-1fa3e0606957
# ╠═4d1715c9-e7d1-4b8f-8f0c-c43a949eb80d
# ╠═d26bf1d9-15a9-4486-959d-ca0cc65b9c84
# ╠═7b8f5f94-3c8b-4b35-babe-d51f871fcefc
# ╠═b8791d74-0008-455b-aaef-99118f0571be
# ╠═a9102099-b4eb-4842-8c18-bb1761243c7f
# ╠═8e134977-e7b8-4adb-8a0b-780b08d6a687
# ╠═6a38b860-c16e-4329-bb74-9c8ea319ce14
# ╠═4495899b-3bb6-4691-91c8-9fbc94df8910
# ╠═01dabf78-fcaa-4cb8-940d-fd787f551d66
# ╠═33ec4aa6-d075-422d-8ebf-df81ee6a68c7
# ╠═e538c15a-5b04-4089-8a64-fa9d347831e4
# ╠═a96f5ac8-b02c-49e6-97fa-e5b6287113b4
# ╠═42cb7f68-3483-4041-b9b8-d3c6f292f30c
# ╠═c38ab2cf-b7d1-4f3e-8aa7-360cd1aae8fb
# ╠═960e3093-bd76-4063-88f5-b53619b886f4
