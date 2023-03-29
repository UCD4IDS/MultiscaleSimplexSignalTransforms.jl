### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 36182664-b839-11ed-36a5-7b82aaee1af4
using Pkg; Pkg.activate(".")

# ╔═╡ df3cd885-b12f-4453-8c8e-6ffced6cb3da
using Revise

# ╔═╡ 429269de-ba3e-463d-92cf-303a648af161
using MultiscaleSimplexSignalTransforms

# ╔═╡ 17ae611b-733a-4c21-960d-f41f624c8d12
using Graphs, SimpleWeightedGraphs, Plots, LinearAlgebra, Arpack, SparseArrays, Clustering, InvertedIndices, DataStructures, MAT, Chain, StatsBase, GraphRecipes

# ╔═╡ f73b1925-8fef-468e-9782-0b32e4849887
begin
	minnmat = matread("../data/minnesota/minnesota.mat")
	g = minnmat["Problem"]["A"] |> Graph
	x, y = minnmat["Problem"]["aux"]["coord"] |> eachcol .|> Array
	
	# the data includes a disconnected component!
	for i ∈ [349, 348] # descending to avoid index counting issues
		rem_vertex!(g, i) # remove the node and all associated info in the graph
		# imitate this process for the coordinates -- the graph method is to switch
		# the position with the final vertex and the deleted one.
		x[i] = x[end]
		pop!(x)
		y[i] = y[end]
		pop!(y)
	end

	@assert length(connected_components(g)) == 1 "graph is disconnected"
end

# ╔═╡ 2e291a45-4d11-4049-9523-f0f022bd9825
k=0

# ╔═╡ 31c90731-e9b6-44f5-ae1a-c8037f7f97ce
reg = KRegion(cliquecomplex(g), k);

# ╔═╡ 8cb26c61-a5cd-483d-b33d-c84f3cb373c5
part = kGHWT(reg; representation=Representation.Distance)

# ╔═╡ ec6c7bee-45f0-4feb-8d5d-105438d5fab8


# ╔═╡ 4e9f94ef-296a-4969-8791-acfa79396362
n = 30

# ╔═╡ 0cd5d453-2a56-489a-b2d7-e6d4dbe267be
t = randomtree(n)

# ╔═╡ 92fce89b-1af6-4359-865b-081f4fbd8224
treg = ZeroRegion(t);

# ╔═╡ 477e06fe-adf9-4e02-ad04-fb83e7462445
L = k_laplacian(treg)

# ╔═╡ 1f99d5bf-f20c-401e-b37a-619141cb79c0
P = I - ones(n,n)./n;

# ╔═╡ 2f9563cc-bd1d-4344-ac7c-d4c3b32bac7f
K = -0.5*floyd_warshall_shortest_paths(t).dists

# ╔═╡ 77884cb2-1c7e-4e83-8cde-700417d2ca47
Lsym = k_laplacian(treg; normalization=Normalization.Symmetric)

# ╔═╡ 5ce5fa81-49fc-4b05-a255-5608bcb50234
d = diag(L).^0.5

# ╔═╡ 957ffebc-5ef9-46d9-9db6-5e892e47556e
Ksym = Diagonal(d) |> D -> D*K*D

# ╔═╡ ee29547b-ac13-4503-bafa-29758ec3a66e
Pd = I - d*d' / (d'*d)

# ╔═╡ 221f2957-91e6-4c2f-885e-c3ee7ac3de51
P*L - L |> norm

# ╔═╡ 03200d61-7c1b-483a-8489-d66750f642c0
Pd * Lsym - Lsym |> norm

# ╔═╡ 786d873b-f458-43a3-988c-dce0c0ec5254
P*K*L - L*K*P |> norm

# ╔═╡ 106318d8-7740-4d01-8a07-ff1534ce33f9
Pd * Ksym * Lsym - Lsym * Ksym * Pd |> norm

# ╔═╡ 6a1f6034-c5dd-4105-adc0-4d4f45af33d1
distance_kernel(treg) - P*K*P |> norm

# ╔═╡ 1424a262-2339-4c25-a02d-f4452a816d70
distance_kernel(treg; normalization=Normalization.Symmetric) - Pd*Ksym*Pd |> norm

# ╔═╡ Cell order:
# ╠═36182664-b839-11ed-36a5-7b82aaee1af4
# ╠═df3cd885-b12f-4453-8c8e-6ffced6cb3da
# ╠═429269de-ba3e-463d-92cf-303a648af161
# ╠═17ae611b-733a-4c21-960d-f41f624c8d12
# ╠═f73b1925-8fef-468e-9782-0b32e4849887
# ╠═2e291a45-4d11-4049-9523-f0f022bd9825
# ╠═31c90731-e9b6-44f5-ae1a-c8037f7f97ce
# ╠═8cb26c61-a5cd-483d-b33d-c84f3cb373c5
# ╠═ec6c7bee-45f0-4feb-8d5d-105438d5fab8
# ╠═4e9f94ef-296a-4969-8791-acfa79396362
# ╠═0cd5d453-2a56-489a-b2d7-e6d4dbe267be
# ╠═92fce89b-1af6-4359-865b-081f4fbd8224
# ╠═477e06fe-adf9-4e02-ad04-fb83e7462445
# ╠═1f99d5bf-f20c-401e-b37a-619141cb79c0
# ╠═2f9563cc-bd1d-4344-ac7c-d4c3b32bac7f
# ╠═77884cb2-1c7e-4e83-8cde-700417d2ca47
# ╠═5ce5fa81-49fc-4b05-a255-5608bcb50234
# ╠═957ffebc-5ef9-46d9-9db6-5e892e47556e
# ╠═ee29547b-ac13-4503-bafa-29758ec3a66e
# ╠═221f2957-91e6-4c2f-885e-c3ee7ac3de51
# ╠═03200d61-7c1b-483a-8489-d66750f642c0
# ╠═786d873b-f458-43a3-988c-dce0c0ec5254
# ╠═106318d8-7740-4d01-8a07-ff1534ce33f9
# ╠═6a1f6034-c5dd-4105-adc0-4d4f45af33d1
# ╠═1424a262-2339-4c25-a02d-f4452a816d70
