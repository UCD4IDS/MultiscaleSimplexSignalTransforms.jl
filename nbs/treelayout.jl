### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ e45f7a26-7341-4178-9a19-8ca25b4b39bd
using Pkg; Pkg.activate(".")

# ╔═╡ cae1ba45-0991-49dd-a947-2e1c2cbaa495
using Revise

# ╔═╡ 243921ee-0d8c-40a0-9763-040b579f3161
using MultiscaleSimplexSignalTransforms

# ╔═╡ 812e4a50-c5dd-11ed-0220-152dfc8cbf8b
using NetworkLayout, Graphs, LinearAlgebra, SparseArrays, Unzip, Plots, Chain

# ╔═╡ 225831c6-e59a-40e2-8cb7-ac8b2a9645b1
G = randomtree(50)

# ╔═╡ a3e0f127-0704-4ac3-9033-8c5b85762e30
usedfs = false

# ╔═╡ 39cd2b9e-3233-4e88-922f-eb961d80f908
treg, pos = tree_layout(G, usedfs)

# ╔═╡ fcfcec59-7d29-416b-b981-052d586cc072
begin
	a = SimplexTree()
	insert!(a, [1,2,3,5])
	insert!(a, [4,5,6])
	insert!(a, [5,6])
end

# ╔═╡ 513093f7-05e2-4526-8c83-59d1d3c32430
areg, apos = tree_layout(a)

# ╔═╡ 078e1bea-8751-498f-b850-9650e966583e
@chain begin
plot(
	splot(
		areg, ones(n_simp(areg)), apos;
		k=0, size=30, palette=palette([:skyblue, :skyblue], 2)
	);
	size=(1000, 1000),
	annotations = [
		(p.+[-0,0] ..., (
			simplices(a)[i] |> s -> isempty(s) ? "∅" : s[end], 
		20))
		for (i,p) in enumerate(zip(apos.x, apos.y))
	]
)
# @aside savefig("/Users/eug/Desktop/simplextree.png")
end

# ╔═╡ 9cb30bd3-c09f-4127-9d16-37666fd8c6b6
tpart = kGHWT(treg)

# ╔═╡ 6adf754f-abb2-4040-8520-26ba22ec6d4f
plot(
	splot(treg, tpart[0, 0, 0], pos; k=0, size=30, palette=:redblue);
	size=(1000, 1800),
	annotations = [
		(p.+[-0,0] ..., (string(i), 20))
		for (i,p) in enumerate(zip(pos.x, pos.y))
	]
)

# ╔═╡ 81e2647e-58eb-4ff9-be04-0b2e3c2744a9
plot(
	[
		splot(
			treg,
			tpart[1,0,i],
			pos; 
			k=0, size=20, palette=:redblue
		)
		for i=0:11
		# for v ∈ eachcol(tpart[tpart.inds.sindmap[[6,11,15,16,17,21,22,27]]...,:])
	]...;
	size=(1000,2000)
)

# ╔═╡ 2c98692d-5eee-4442-90e4-01281020d563
# plot([
# 	splot(treg, v, pos; k=0)
# 	for v ∈ eachcol(analyze(tpart, ones(50))[:,:][:,1:10])
# ]...)

# ╔═╡ e14b6658-441f-43e7-b51d-0f20b4368b90
B = best_basis(tpart, rand(50))

# ╔═╡ bc051c48-a2ef-4d36-8101-b69ef75c681e
# tpart[B.locs[1]...,:]'*ones(50)

# ╔═╡ 0fec6b62-aa60-4486-a42a-b7f719ebdc6c


# ╔═╡ Cell order:
# ╠═e45f7a26-7341-4178-9a19-8ca25b4b39bd
# ╠═cae1ba45-0991-49dd-a947-2e1c2cbaa495
# ╠═243921ee-0d8c-40a0-9763-040b579f3161
# ╠═812e4a50-c5dd-11ed-0220-152dfc8cbf8b
# ╠═225831c6-e59a-40e2-8cb7-ac8b2a9645b1
# ╠═a3e0f127-0704-4ac3-9033-8c5b85762e30
# ╠═39cd2b9e-3233-4e88-922f-eb961d80f908
# ╠═6adf754f-abb2-4040-8520-26ba22ec6d4f
# ╠═fcfcec59-7d29-416b-b981-052d586cc072
# ╠═513093f7-05e2-4526-8c83-59d1d3c32430
# ╠═078e1bea-8751-498f-b850-9650e966583e
# ╠═9cb30bd3-c09f-4127-9d16-37666fd8c6b6
# ╠═81e2647e-58eb-4ff9-be04-0b2e3c2744a9
# ╠═2c98692d-5eee-4442-90e4-01281020d563
# ╠═e14b6658-441f-43e7-b51d-0f20b4368b90
# ╠═bc051c48-a2ef-4d36-8101-b69ef75c681e
# ╠═0fec6b62-aa60-4486-a42a-b7f719ebdc6c
