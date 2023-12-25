### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ e4618ea1-7d1c-43aa-a7d7-f848c4418e50
using Pkg; Pkg.activate(".")

# ╔═╡ 8fcdef9d-2116-40b2-8bcc-a34c6c0ddbe7
using Revise

# ╔═╡ 582cfd91-faa0-48c8-b0d7-307f5718d53b
using LinearAlgebra, Chain, Graphs, SimpleWeightedGraphs, MultiscaleSimplexSignalTransforms, SparseArrays

# ╔═╡ b5012d98-8e25-4602-b6c9-25708f7b7817
using Plots, GraphRecipes

# ╔═╡ 85b06160-9ded-11ee-3e2e-a5ef3b35e93f
a1, a2 = rand(2)

# ╔═╡ 112e3f92-ffd5-4b98-8c4e-f6c529e2803b
b1, b2 = rand(2)

# ╔═╡ 23ac1f26-4e96-44b0-9a98-ba8bd937454b
B1 = [
	0 a1
	b1 0
]

# ╔═╡ ca24f9e8-3eb3-43e7-b9b5-b1244e616929
B2 = [
	0 a1 a1+a2
	b1 0 a2
	b1+b2 b2 0
]

# ╔═╡ bbf1eb37-ac22-475c-9c7f-3bd9a7a4b92e
B1i = inv(B1)

# ╔═╡ d85bf01f-537d-4aec-9294-be458a2e649a
B2i = inv(B2)

# ╔═╡ 1940b00a-bd5f-4759-9141-4bcb2c67ab83
u1 = B1i*ones(2)

# ╔═╡ 93bb3277-af6d-4987-8392-d45c075efa88
u2 = B2i*ones(3)

# ╔═╡ d8f5aed6-ef14-419c-82ac-906d10dbb7f4
v1 = B1i'*ones(2)

# ╔═╡ 1cf0e998-13d3-4bb9-aaa1-e28c0c1da1a7
v2 = B2i'*ones(3)

# ╔═╡ 55dfa8cf-0ab4-494a-ab80-de133fa00213
r1 = sum(u1)

# ╔═╡ c6063a4e-3a23-459c-8485-c930d7277748
norm(sum(v1).-r1)

# ╔═╡ 11070e57-5f1b-4c68-8c7e-c26f3584fde6
r2 = sum(u2)

# ╔═╡ 01e6adb2-1f15-43c2-9212-49200a25de74
norm(sum(v2).-r2)

# ╔═╡ 11dae83f-b391-4f05-9370-c7b01b48f680
s2 = a2+b2+a2*b2*r1

# ╔═╡ 25af3359-4ebb-429a-93f0-289c138d8ffc
B1iun = B1i[:,end]

# ╔═╡ cf61357a-ea16-4125-b001-127dbd42ee7a
B1ivn = B1i'[:,end]

# ╔═╡ 3fdf5ba3-68e4-43b6-a9fb-82746c1c9197
B1inn = B1iun[end]

# ╔═╡ abf09987-0430-4cda-9f73-8dfb60a2dc19
norm(B1inn - B1ivn[end])

# ╔═╡ 69d38e94-810a-4983-85bd-ac53b706799d
B2iun = B2i[:, end]

# ╔═╡ c8a86c1f-3cf0-410a-9e7e-b56d2d93926c
B2ivn = B2i'[:, end]

# ╔═╡ be62fe4a-d9be-45d6-bab1-64dd4df55412
B2inn = B2iun[end]

# ╔═╡ 1449c3ed-619e-4437-93af-8dfbc91e0afc
norm(B2inn - B2ivn[end])

# ╔═╡ d0c585ac-812c-4984-9105-3fd1fd9bcc07
R(n, pos=true) = [
	I(n-2) 0 0
	0 1 (pos ? 1 : -1)
	0 0 1
]

# ╔═╡ d665aa82-e862-4f49-b482-7d132de4ef17
function pad(M)
	n = size(M, 1)
	[
		M zeros(n)
		zeros(n)' 0
	]
end

# ╔═╡ 7d4759a5-58d2-46ba-9ea9-149a440b1fe6
function Rsandwich(M)
	n = size(M, 1)
	R(n+1)' * pad(M) * R(n+1)
end

# ╔═╡ 26856de0-3ca1-4e83-94e0-6d79bc38ee59
function Risandwich(M)
	n = size(M, 1)
	R(n+1, false) * pad(M) * R(n+1, false)'
end

# ╔═╡ ec7bc395-5e79-4998-ada6-5666b75190c8
function padded!(v, amt=1)
	v[end] += amt
	v
end

# ╔═╡ cd16cd34-cef7-44ea-9234-24777b07261f
norm(B2i - (
	pad(B1i) - 1/s2 .* [padded!(a2*u1); -1] * [padded!(b2*v1); -1]'
))

# ╔═╡ ae74dd3b-d2fc-4e9a-b75a-4e5f80c66c22
(
	norm(u2 - (
		[u1; 0] - b2*r1/s2 .* [padded!(a2.*u1);-1]
	)),
	norm(u2 - (
		[padded!((a2+b2).*u1, -b2*r1); b2*r1]./s2
	))
)

# ╔═╡ b3144c6e-e1aa-4106-9037-8f5faaf221fd
(
	norm(v2 - (
		[v1; 0] - a2*r1/s2 .* [padded!(b2.*v1); -1]
	)),
	norm(v2 - (
		[padded!((a2+b2).*v1, -a2*r1); a2*r1]./s2
	))
)

# ╔═╡ 186a42f7-d5c3-4908-844d-d4b8d1110122
norm(1/r2 - (
	1/r1 + 1/(1/a2+1/b2)
))

# ╔═╡ 62d73438-8ad9-47e7-8ab5-d60b3342c54d
md"scratch below here"

# ╔═╡ 7a07e29b-27ab-4c75-bba6-525179f6396b
@chain begin 
	[
		B1 a2*ones(2)
		b2*ones(2)' -(a2+b2)
	]
	R(3)'*_*R(3)
	_ - B2
end

# ╔═╡ bb054392-49d1-42c0-93f5-07c8cbd000aa
@chain begin
	[
		B1 a2*ones(2)
		b2*ones(2)' -(a2+b2)
	]
	inv
	B2i - R(3)*_*R(3)'
	_ * s2
end

# ╔═╡ 66a5fa4b-60e5-4754-9df9-f8f43054210b
W = diagm(1 => rand(10), -1 => rand(10))

# ╔═╡ a3448bae-501d-4d4e-b1a8-5874f1616a95
B = distmat(SimpleWeightedDiGraph(W))

# ╔═╡ 8ec4b42e-a539-4469-946d-824632d88adf
1 ./ sum(1 ./ (1 ./diag(W, 1) .+ 1 ./ diag(W, -1))) - sum(inv(Matrix(B)))

# ╔═╡ ec855c85-1310-4eae-88a5-ddf05b8ae380
L1 = diagm(
	2,2,
	1  => 1 ./ [a1],
	-1 => 1 ./ [b1],
) |> (W -> Diagonal(sum(W, dims=2)[:]) - W)

# ╔═╡ 03f9e4a0-c4af-4d58-877c-ca669feae000
L2 = diagm(
	3,3,
	1  => 1 ./ [a1, a2],
	-1 => 1 ./ [b1, b2],
) |> (W -> Diagonal(sum(W, dims=2)[:]) - W)

# ╔═╡ 74fc19c1-3d50-4868-a067-251b36298bb7
B1i

# ╔═╡ cfa74425-a07e-4904-b27a-85b550873679
B1i + (a1*b1)/(a1^2+b1^2) * L1 |> rank

# ╔═╡ 89b38b7d-0189-48ff-8063-107bbcf6c93a
B2i + (a2*b2)/(a2^2+b2^2) * L2

# ╔═╡ 754123bb-9979-4183-8a5a-828e9ca25349
(B2i - u2 * v2' ./ r2)

# ╔═╡ 37c1542c-645c-483b-89a1-e47783c37f8f
(1 ./[a1, a2], 1 ./[b1, b2])

# ╔═╡ 4b2406b3-5ac9-4fbb-962c-8c8788161f6b
g = @chain begin
	randomtree(20)
	adjacency_matrix(_) .* rand(nv(_), nv(_))
	SimpleWeightedDiGraph
end

# ╔═╡ 4e1083a9-7bd2-4a8c-a82c-bdc73d32b8f7
is_connected(g)

# ╔═╡ a0b9b881-ae12-4447-ae3b-7218aac0c4a4
begin
	gLs = weights(g)+weights(g)'
	map!(x -> 1 ./ x, gLs.nzval, gLs.nzval)
	gLs = Diagonal(sum(gLs, dims=2)[:])-gLs
end

# ╔═╡ c9162c85-3e17-4cfe-a1b5-94db4283da3c
gwi = inv(Matrix(distmat(g)))

# ╔═╡ 35a686de-bfe1-4532-938a-d608ab96398d
gwil = gwi - sum(gwi, dims=2) * sum(gwi, dims=1) / sum(gwi)

# ╔═╡ b892c89f-9a38-4f3c-8ac1-3b1c2b25e1b4
sum(gwi)

# ╔═╡ 74afc94b-4d6e-454e-82ed-488a5282f4db
norm(gwil+gLs)

# ╔═╡ 3a2e08cc-a559-4cc1-9d53-4fb7d7373439
P(M) = I-ones(size(M))/size(M,1) |> PP -> PP*M*PP

# ╔═╡ 06e68dfc-4cfd-49c5-b4bb-dca564f85c0b
P(distmat(g)) |> M -> norm(M-M')

# ╔═╡ a1e5ffd7-b58e-4bb3-86cb-d488f213c24e
pinv(Matrix(gLs)) + P(distmat(g)) |> norm

# ╔═╡ ab7cc904-cd1c-48bc-81fe-174ff4ec178b
spgwil = (@. ifelse(abs(gwil) > 1e-12, gwil, 0)) |> sparse

# ╔═╡ 858e5a0a-008c-45fb-ad87-562689eac34b
function makegwr()
	jj,ii,vv = findnz(weights(g)')
	sparse(ii, jj, 1 ./vv, nv(g), nv(g))
end

# ╔═╡ 2f72b739-2c62-4f19-a689-35e40882922e
gwr = makegwr()

# ╔═╡ 2e935757-b7a7-4be7-9c31-f1d5bc6f4b75
gL = gwr |> W -> Diagonal(sum(W, dims=2)[:]) - W

# ╔═╡ a4b368ff-3587-4d6c-9ace-9e5a839ca103


# ╔═╡ 205b4773-607d-47f0-9219-dd77a8ea70c3
@chain begin
	findnz(spgwil)
	sparse(_[1], _[2], _[3] ./ -nonzeros(gL))
	_[3,:]
	nonzeros
end

# ╔═╡ 2bbf8ae0-5fa7-4d21-b9b5-604618cd3f83
BB = [
	# 1  -1  0  0
	# -1  1  1 -1
	# 0   0 -1  1
	1  0
	-1 1
	0 -1
]

# ╔═╡ fe0c3aff-664e-4338-9b0a-da57c2348dc6
d = [
	(a1-b1)/(a1+b1)
	(b1-a1)/(a1+b1)
	# (a2-b2)/(a2+b2)
	# (b2-a2)/(a2+b2)
]

# ╔═╡ 2ce741cd-8a90-4e60-9a78-4115e88078ff
ones(3) - BB*d

# ╔═╡ 28a8bc01-05ef-4beb-874d-2487bae5d4f4
[a1/(a1+b1),a2/(a2+b2)-a1/(a2+b2),1-a2/(a2+b2)]

# ╔═╡ 569cec6f-99d1-4179-b138-947c3f01a26d
eg = erdos_renyi(100,0.5)

# ╔═╡ 53509fc5-f8eb-4599-ab6d-728106c1da00
is_connected(eg)

# ╔═╡ d5e8fc9a-fb4a-43f6-8004-053ab661b424
egw = SimpleWeightedDiGraph(adjacency_matrix(eg) .* rand(100,100))

# ╔═╡ e5b4e901-e7dd-4987-b67f-c166e539bd06
P(distmat(egw)) |> M -> norm(M-M')

# ╔═╡ 84271ee4-5a1f-4812-870d-a204592a38f9
function spmap(f, M)
	i, j, v = findnz(M)
	sparse(i, j, map(f, v), size(M)...)
end

# ╔═╡ c46a86f7-34d8-4960-9318-cf10fb4a7665
function check_grid()
	# g = Graphs.grid((20,25))
	g = randomtree(20)
	W = UpperTriangular(adjacency_matrix(g))
	B = distmat(DiGraph(W))
	L = spmap(x->1/x, W+W') |> WW -> Diagonal(sum(WW, dims=2)[:]) - WW
	P(B)*L-L*P(B) |> norm
end

# ╔═╡ dd451116-f477-4ce8-966a-cdb366eed45e
check_grid()

# ╔═╡ Cell order:
# ╠═e4618ea1-7d1c-43aa-a7d7-f848c4418e50
# ╠═8fcdef9d-2116-40b2-8bcc-a34c6c0ddbe7
# ╠═582cfd91-faa0-48c8-b0d7-307f5718d53b
# ╠═85b06160-9ded-11ee-3e2e-a5ef3b35e93f
# ╠═112e3f92-ffd5-4b98-8c4e-f6c529e2803b
# ╠═23ac1f26-4e96-44b0-9a98-ba8bd937454b
# ╠═ca24f9e8-3eb3-43e7-b9b5-b1244e616929
# ╠═bbf1eb37-ac22-475c-9c7f-3bd9a7a4b92e
# ╠═d85bf01f-537d-4aec-9294-be458a2e649a
# ╠═1940b00a-bd5f-4759-9141-4bcb2c67ab83
# ╠═93bb3277-af6d-4987-8392-d45c075efa88
# ╠═d8f5aed6-ef14-419c-82ac-906d10dbb7f4
# ╠═1cf0e998-13d3-4bb9-aaa1-e28c0c1da1a7
# ╠═55dfa8cf-0ab4-494a-ab80-de133fa00213
# ╠═c6063a4e-3a23-459c-8485-c930d7277748
# ╠═11070e57-5f1b-4c68-8c7e-c26f3584fde6
# ╠═01e6adb2-1f15-43c2-9212-49200a25de74
# ╠═11dae83f-b391-4f05-9370-c7b01b48f680
# ╠═25af3359-4ebb-429a-93f0-289c138d8ffc
# ╠═cf61357a-ea16-4125-b001-127dbd42ee7a
# ╠═3fdf5ba3-68e4-43b6-a9fb-82746c1c9197
# ╠═abf09987-0430-4cda-9f73-8dfb60a2dc19
# ╠═69d38e94-810a-4983-85bd-ac53b706799d
# ╠═c8a86c1f-3cf0-410a-9e7e-b56d2d93926c
# ╠═be62fe4a-d9be-45d6-bab1-64dd4df55412
# ╠═1449c3ed-619e-4437-93af-8dfbc91e0afc
# ╠═d0c585ac-812c-4984-9105-3fd1fd9bcc07
# ╠═d665aa82-e862-4f49-b482-7d132de4ef17
# ╠═7d4759a5-58d2-46ba-9ea9-149a440b1fe6
# ╠═26856de0-3ca1-4e83-94e0-6d79bc38ee59
# ╠═ec7bc395-5e79-4998-ada6-5666b75190c8
# ╠═cd16cd34-cef7-44ea-9234-24777b07261f
# ╠═ae74dd3b-d2fc-4e9a-b75a-4e5f80c66c22
# ╠═b3144c6e-e1aa-4106-9037-8f5faaf221fd
# ╠═186a42f7-d5c3-4908-844d-d4b8d1110122
# ╟─62d73438-8ad9-47e7-8ab5-d60b3342c54d
# ╠═7a07e29b-27ab-4c75-bba6-525179f6396b
# ╠═bb054392-49d1-42c0-93f5-07c8cbd000aa
# ╠═66a5fa4b-60e5-4754-9df9-f8f43054210b
# ╠═a3448bae-501d-4d4e-b1a8-5874f1616a95
# ╠═8ec4b42e-a539-4469-946d-824632d88adf
# ╠═ec855c85-1310-4eae-88a5-ddf05b8ae380
# ╠═03f9e4a0-c4af-4d58-877c-ca669feae000
# ╠═74fc19c1-3d50-4868-a067-251b36298bb7
# ╠═cfa74425-a07e-4904-b27a-85b550873679
# ╠═89b38b7d-0189-48ff-8063-107bbcf6c93a
# ╠═754123bb-9979-4183-8a5a-828e9ca25349
# ╠═37c1542c-645c-483b-89a1-e47783c37f8f
# ╠═4b2406b3-5ac9-4fbb-962c-8c8788161f6b
# ╠═4e1083a9-7bd2-4a8c-a82c-bdc73d32b8f7
# ╠═a0b9b881-ae12-4447-ae3b-7218aac0c4a4
# ╠═c9162c85-3e17-4cfe-a1b5-94db4283da3c
# ╠═35a686de-bfe1-4532-938a-d608ab96398d
# ╠═b892c89f-9a38-4f3c-8ac1-3b1c2b25e1b4
# ╠═74afc94b-4d6e-454e-82ed-488a5282f4db
# ╠═3a2e08cc-a559-4cc1-9d53-4fb7d7373439
# ╠═06e68dfc-4cfd-49c5-b4bb-dca564f85c0b
# ╠═a1e5ffd7-b58e-4bb3-86cb-d488f213c24e
# ╠═ab7cc904-cd1c-48bc-81fe-174ff4ec178b
# ╠═858e5a0a-008c-45fb-ad87-562689eac34b
# ╠═2f72b739-2c62-4f19-a689-35e40882922e
# ╠═2e935757-b7a7-4be7-9c31-f1d5bc6f4b75
# ╠═a4b368ff-3587-4d6c-9ace-9e5a839ca103
# ╠═205b4773-607d-47f0-9219-dd77a8ea70c3
# ╠═2bbf8ae0-5fa7-4d21-b9b5-604618cd3f83
# ╠═fe0c3aff-664e-4338-9b0a-da57c2348dc6
# ╠═2ce741cd-8a90-4e60-9a78-4115e88078ff
# ╠═28a8bc01-05ef-4beb-874d-2487bae5d4f4
# ╠═569cec6f-99d1-4179-b138-947c3f01a26d
# ╠═53509fc5-f8eb-4599-ab6d-728106c1da00
# ╠═d5e8fc9a-fb4a-43f6-8004-053ab661b424
# ╠═e5b4e901-e7dd-4987-b67f-c166e539bd06
# ╠═b5012d98-8e25-4602-b6c9-25708f7b7817
# ╠═84271ee4-5a1f-4812-870d-a204592a38f9
# ╠═c46a86f7-34d8-4960-9318-cf10fb4a7665
# ╠═dd451116-f477-4ce8-966a-cdb366eed45e
