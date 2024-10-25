### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 6b6b348c-b9b5-11ed-2dd1-811b6386fcab
using Pkg; Pkg.activate(".");

# ╔═╡ b2d5ecf8-66c7-476c-9c19-aa250f1e2247
using Revise

# ╔═╡ 169da0e5-dd5d-4ccf-b45d-1ff130d4ee59
using MultiscaleSimplexSignalTransforms

# ╔═╡ 921e885d-c55f-4608-a1b9-c88d7199a881
using Graphs, SimpleWeightedGraphs, Plots, LinearAlgebra, Arpack, SparseArrays, Clustering, InvertedIndices, DataStructures, MAT, Chain, StatsBase, GraphRecipes

# ╔═╡ ca68daff-c171-4d07-ac51-b54b731e9d62
using EnumX, Match

# ╔═╡ 64e8bf92-36d1-4a0a-8e62-b2b49723d336
using Unzip

# ╔═╡ bf93b262-c100-43a2-b56e-59a673bc0c05
using Roots

# ╔═╡ f2669aa6-c442-415e-bfd7-22b18d698e0f
function import_minnesota(weighted=false)
    minnmat = matread("../data/minnesota/minnesota.mat")

    x, y = minnmat["Problem"]["aux"]["coord"] |> eachcol .|> Array
    g = if !weighted
        minnmat["Problem"]["A"] |> Graph
    else
        adj = minnmat["Problem"]["A"]
        wts = findnz(adj) |> ((inz, jnz, vnz),) -> sparse(inz, jnz, [
            1 ./ norm([x[i], y[i]] - [x[j], y[j]])
            for (i, j) ∈ zip(inz, jnz)
        ])
        SimpleWeightedGraph(wts)
    end

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
    (g, Pos(x, y))
end

# ╔═╡ a406bd25-ca0a-4c16-82bd-4a54d643be11
g, p = import_minnesota(false)

# ╔═╡ dfa29fbc-3857-4764-8b5e-763f64e99563
minreg = KRegion(g, 1; weakadjweight=1.0)

# ╔═╡ f73267d6-6b2c-4aca-bf38-7cd298cb5c1d
minlap = k_laplacian(
    minreg;
    normalization=Normalization.Combinatorial
)

# ╔═╡ 82c35c67-9c23-42d4-a95b-da275c20d44a
mineigs = eigen(Matrix(minlap))

# ╔═╡ 65f41964-2476-428d-b814-d9ae8d94cab4
minpos = kGHWT(
    minreg,
    SubRepresentation.Full;
    eigenmethod=EigenMethod.Positive
)

# ╔═╡ c1ffc86e-149f-404f-8787-01d25a8c8708
@chain (3, 4) begin
    plot(
        [
            splot(
                minreg, minpos[0, 0, i], p; k=1, size=8,
                palette=:redblue
            )
            for i ∈ 0:prod(_)-1
        ]...,
        layout=reverse(_), size=(800 * _[1], 600 * _[2])
    )
end

# ╔═╡ 19cb36ba-e7a8-4a85-9af0-ad4eae492b72


# ╔═╡ c2b9103b-f47c-4963-9b82-6adec1f5d0f3
findfirst(>(1e-10), mineigs.values)

# ╔═╡ 3a581b39-c28b-434f-8366-9cc23f0a4127
plot(mineigs.values)

# ╔═╡ 0c644bb4-c6c3-4983-842d-8b752b81f484
# partition!(minip, 664, 666)

# ╔═╡ f8a784d3-ec7b-4a86-9bbd-6db529f3aae5
# mvpointer!(minip, Direction.Right)

# ╔═╡ fbb92e88-fab5-4730-85ca-e7c57edf17fd
# setbasis!(minip)

# ╔═╡ 7c8e0eb7-35ca-4eac-b0c0-7acbc5c70d79
@enumx Direction Up Left Right

# ╔═╡ 5b85c010-9a3f-4e42-abc2-1c3a614c7c3a
begin
    mutable struct InteractivePartition
        complex::SimplexTree
        k::Int
        part::PartitionTree
        pointer::Vector{Bool}
        function InteractivePartition(complex::SimplexTree, k::Int)
            ip = new(
                complex, k, PartitionTree(length(simplices(complex, k))), []
            )
            setbasis!(ip)
            ip
        end
    end

    function getnode(ip::InteractivePartition)
        ret = ip.part
        for goright ∈ ip.pointer
            ret = goright ? ret.rchild : ret.lchild
        end
        ret
    end

    function setbasis!(ip::InteractivePartition)
        node = getnode(ip)

        eigs = @chain ip.complex begin
            KRegion(
                ip.k;
                weakadjweight=1.0, subregion_inds=node.sinds
            )
            k_laplacian(; normalization=Normalization.Combinatorial)
            Matrix
            eigen
        end

        node.fvals = TagDict(
            i - 1 => Vector(sparsevec(node.sinds, c, partlen(ip.part)))
            for (i, c) = enumerate(eachcol(eigs.vectors))
        )
        node.fvals[-1] = eigs.values
    end

    function mvpointer!(ip::InteractivePartition, dir::Direction.T)
        @match dir begin
            $(Direction.Up) => !isempty(ip.pointer) && pop!(ip.pointer)
            $(Direction.Left) => push!(ip.pointer, false)
            $(Direction.Right) => push!(ip.pointer, true)
        end
        ip.pointer
    end

    Base.getindex(ip::InteractivePartition, i::Int) = @chain ip begin
        getnode
        _.fvals[i-1]
        val
    end

    vals(ip::InteractivePartition) = ip[0]
end

# ╔═╡ 898f6858-f699-4ced-ae15-0ef4cee33ec0
minip = InteractivePartition(cliquecomplex(g), 1)

# ╔═╡ a22a85de-7cc8-4a97-9c78-e720b46b7935
@chain (10, 1e-4) begin
    @aside dcind = findfirst(>(_[2]), vals(minip))
    @aside @info dcind
    plot([
            splot(
                minreg,
                sign.(minip[dcind]) .* minip[dcind+i],
                p; k=1, size=8
            )
            for i ∈ 0:_[1]-1
        ]..., layout=(ceil(Int, _[1] / 2), 2), size=(1600, 300 * _[1]))

    # @aside savefig("/Users/eug/Desktop/minnesota-edge-part3.png")
end

# ╔═╡ e21a0dc3-e20c-4a15-9358-64aa415a8960
plot(vals(minip))

# ╔═╡ 92073ca9-82be-4d27-8c1f-4566ac409f90
findfirst(>(1e-4), minip[0])

# ╔═╡ 1cebc527-2c5a-4d92-a5c5-373b444f38ba
function MultiscaleSimplexSignalTransforms.partition!(ip::InteractivePartition, dc_i, fiedler_i)
    node = getnode(ip)
    if length(children(node)) > 0
        @info "Tried to partition an already-partitioned node, aborting"
        return
    end

    if MultiscaleSimplexSignalTransforms.isleaf(node)
        return
    end

    if partlen(node) == 2
        MultiscaleSimplexSignalTransforms.lchild!(node, PartitionTree(; sinds=[node.sinds[1]]))
        MultiscaleSimplexSignalTransforms.rchild!(node, PartitionTree(; sinds=[node.sinds[2]]))
        return
    end

    (lpart, rpart) = part_by_sign(sign.(ip[dc_i]) .* ip[fiedler_i])
    any(isempty.((lpart, rpart))) && error("partition was trivial")

    MultiscaleSimplexSignalTransforms.lchild!(node, PartitionTree(; sinds=node.sinds[lpart]))
    MultiscaleSimplexSignalTransforms.rchild!(node, PartitionTree(; sinds=node.sinds[rpart]))

    nothing
end

# ╔═╡ 35075de0-f6f2-4e97-8511-a7f5a68c3938
k = 0

# ╔═╡ e5e2f528-c8ef-452d-8b34-f676be04ad65
reg = ZeroRegion(g);

# ╔═╡ 1313b07a-436d-4565-84d2-897d97cdcf1f
part = kGHWT(
    reg, SubRepresentation.Full;
    normalization=Normalization.Combinatorial
)

# ╔═╡ 50d68e4a-deae-42e6-9039-3d75dbd8d3d4
plot([
    splot(reg, ppart, p; k=0)
    for ppart ∈ eachcol(part[1, 0, :][:, 1:10])
]...)

# ╔═╡ 8f9eda5c-1f14-4447-a47c-ccc9f5456072
M = k_laplacian(reg);

# ╔═╡ aa8bd062-076d-4a3d-8a36-18ae229f856a
@chain begin
    [treecolors(part.part.root; maxlevel=i) for i ∈ 0:5]
    [
        splot(reg, _[i], p; k, palette=:viridis, size=5, msw=0, clims=(minimum(_[i]), 1))
        for i ∈ eachindex(_)
    ]
    plot(_..., layout=(2, 3), size=(1800, 1300))
    # @aside savefig("/Users/eug/Desktop/hier-part-minnesota.png")
end

# ╔═╡ b4321dc0-5b9a-444f-8eb9-094927cb6f49
lpart = kGHWT(reg; normalization=Normalization.Combinatorial);

# ╔═╡ 89bf067c-ac02-4bf3-891e-456e02e9e435
compare(Is; kwargs...) = plot(
    splot(reg, sign.(part[Is]), p; k, title="Full Partition $(Is[1:2])", kwargs...),
    splot(reg, sign.(lpart[Is]), p; k, title="Submatrix Partition $(Is[1:2])", kwargs...),
    size=(800, 400)
)

# ╔═╡ 052d17dd-310e-4cbd-93e1-7161e3d6d16b
@chain begin
    Is = (1, 1, 1)
    compare(; palette=:redblue, size=5, msw=0)
    # @aside savefig("/Users/eug/Desktop/subfull-min0-$(Is[1:2])")
end

# ╔═╡ 1c900f73-07ea-4cad-bffe-167e168c3120
function wiggles(n=50, k=1, lower=true, fiedlerflip=true)
    k = lower ? k + 1 : k
    (pathg, pathpos) = lower ?
                       simplex_path_lower(n, k) :
                       simplex_path_lower(n, k + 1)
    pathregion = KRegion(Graph(pathg), k; weakadjweight=lower ? 1.0 : 0.0)
    pathpart = kHGLET(
        pathregion;
        normalization=Normalization.Combinatorial,
        input=fiedlerflip ? PartitionInput.DCOrientation : PartitionInput.Identity
    )
    (pathregion, pathpart, pathpos, k)
end

# ╔═╡ 620dd37f-adac-45ba-b1e1-ae810b5d8e5a
nlin = 50

# ╔═╡ ff2bfce6-481f-4b1d-82d0-f1ebb9df615b
preg, ppart, ppos, pk = wiggles(nlin, 1, true, true)

# ╔═╡ 21c2e289-41d9-4391-b38e-04c41a24849e
@chain begin
    # maximum(abs.(Xpath[:,1:15]))
    plot(
        # splot_many(ZeroRegion(preg.tree), ppos, [
        splot_many(preg, ppos, [
                # Xpath[:,i+1]
                ppart[0, 0, i]
                for i ∈ 0:14
            ] .* sign.([ppart[0, 0, i][1] for i ∈ 0:14])...;
            # k=0,
            k=pk,
            size=3, palette=:viridis);
        size=(600, 1000)
        # clims=0.5.*(-_,_)
    )
    # @aside savefig("/Users/eug/Desktop/path-two.png")
end

# ╔═╡ d5cb6bb8-80be-4e92-9471-4dfed8ec0978
Lpath, Xpath = @chain begin
    eigen(Matrix(k_laplacian(ZeroRegion(preg.tree))))
    (_.values, _.vectors .* _.vectors[1, :]')
end

# ╔═╡ 2839684e-3bec-4201-a225-5fdb267de92d
@chain begin
    plot([ppart[0, 0, i] for i = 1:3], labels=["ϕ₁" "ϕ₂" "ϕ₃"])
    # , title="Eigenvectors of P_{50,1}")
    # @aside savefig("/Users/eug/Desktop/path-upper-plot.png")
end

# ╔═╡ 41bf07bd-cda5-46f8-8093-854850611b16
Ib, Js = unzip([
    (i, j)
    for (i, b) ∈ enumerate(boundaries(preg))
    for j ∈ values(preg.hingemap[b].leaves)
])

# ╔═╡ c2904a25-4169-4e54-91de-6a98c09aa9ef
B = sparse(Ib, Js, 1.0, n_bd(preg), n_simp(preg))

# ╔═╡ d5cd7efe-3d94-4cd1-af20-e8247ee01aa4
bddegs = @chain begin
    preg.simplexweights.weights
    values
    collect
    _ .^ -1
    sparse(B') \ _
    _ .^ -1
end

# ╔═╡ f77a6c31-4167-40b4-a013-9e9f3ff4cfb3
B' * (bddegs .^ -1) - collect(values(preg.simplexweights.weights)) .^ -1

# ╔═╡ 48e75a9b-cda0-4c61-8bca-261526a2a538
origbd = @chain begin
    preg.bdweights.weights
    values
    collect
end

# ╔═╡ 7913dc93-3a13-4a71-989c-da6e977b9542
plot([bddegs origbd])

# ╔═╡ cce684af-020d-415e-93bb-1aa9d9f3effc
[
    s .+ (i - 1) * (length(simplices(preg.tree, 0)))
    for i ∈ 1:2
    for s ∈ first.(simplices(preg.tree, 0))
]

# ╔═╡ f07d27da-756d-49bf-b6a6-bcca0b31f5ad
repeat(rand(5), 3)

# ╔═╡ 1523c30b-d46d-4f24-aa8d-c26383746115
qlower = false

# ╔═╡ 24217822-d536-437a-b838-7aed80cccd56
qfiedler = false

# ╔═╡ 5e75c2ef-369e-4131-a7bf-56a46291c46e
qreg, qpart, qpos, qk = wiggles(nlin, 8, qlower, qfiedler)

# ╔═╡ 52453948-00df-4435-9287-266fcda5d2d4
nn = qlower ? nlin - 1 : nlin * (qk + 1)

# ╔═╡ f10372c1-35b8-49d0-a525-c4c6ab28c2a2
# plot(qpart[0,0,qk%2==0 ? nlin-2 : 1])
plot(qpart[0, 0, nn])

# ╔═╡ 2cf56e36-c87f-4861-9a5a-360adc114dbd
qpart[0, 0, 1]

# ╔═╡ 82195c7f-74a8-467c-b001-dcec35c54574
l0 = find_zero(x -> coth(sqrt(-x) / 2) - sqrt(-x) / 2, -5)

# ╔═╡ cb4788e3-9d28-48b4-b44f-0638557d98ae
plot(0:0.01:1, x -> cosh(sqrt(-l0) * (x - 0.5)) * (0.5 + sinh(sqrt(-l0)) / (2sqrt(-l0)))^-0.5, ylims=[0, 1.5])

# ╔═╡ Cell order:
# ╠═6b6b348c-b9b5-11ed-2dd1-811b6386fcab
# ╠═b2d5ecf8-66c7-476c-9c19-aa250f1e2247
# ╠═169da0e5-dd5d-4ccf-b45d-1ff130d4ee59
# ╠═921e885d-c55f-4608-a1b9-c88d7199a881
# ╠═f2669aa6-c442-415e-bfd7-22b18d698e0f
# ╠═a406bd25-ca0a-4c16-82bd-4a54d643be11
# ╠═dfa29fbc-3857-4764-8b5e-763f64e99563
# ╠═f73267d6-6b2c-4aca-bf38-7cd298cb5c1d
# ╠═82c35c67-9c23-42d4-a95b-da275c20d44a
# ╠═65f41964-2476-428d-b814-d9ae8d94cab4
# ╠═c1ffc86e-149f-404f-8787-01d25a8c8708
# ╠═19cb36ba-e7a8-4a85-9af0-ad4eae492b72
# ╠═c2b9103b-f47c-4963-9b82-6adec1f5d0f3
# ╠═3a581b39-c28b-434f-8366-9cc23f0a4127
# ╠═a22a85de-7cc8-4a97-9c78-e720b46b7935
# ╠═0c644bb4-c6c3-4983-842d-8b752b81f484
# ╠═f8a784d3-ec7b-4a86-9bbd-6db529f3aae5
# ╠═fbb92e88-fab5-4730-85ca-e7c57edf17fd
# ╠═e21a0dc3-e20c-4a15-9358-64aa415a8960
# ╠═92073ca9-82be-4d27-8c1f-4566ac409f90
# ╠═898f6858-f699-4ced-ae15-0ef4cee33ec0
# ╠═ca68daff-c171-4d07-ac51-b54b731e9d62
# ╠═7c8e0eb7-35ca-4eac-b0c0-7acbc5c70d79
# ╠═5b85c010-9a3f-4e42-abc2-1c3a614c7c3a
# ╠═1cebc527-2c5a-4d92-a5c5-373b444f38ba
# ╠═35075de0-f6f2-4e97-8511-a7f5a68c3938
# ╠═e5e2f528-c8ef-452d-8b34-f676be04ad65
# ╠═1313b07a-436d-4565-84d2-897d97cdcf1f
# ╠═50d68e4a-deae-42e6-9039-3d75dbd8d3d4
# ╠═8f9eda5c-1f14-4447-a47c-ccc9f5456072
# ╠═aa8bd062-076d-4a3d-8a36-18ae229f856a
# ╠═b4321dc0-5b9a-444f-8eb9-094927cb6f49
# ╠═89bf067c-ac02-4bf3-891e-456e02e9e435
# ╠═052d17dd-310e-4cbd-93e1-7161e3d6d16b
# ╠═1c900f73-07ea-4cad-bffe-167e168c3120
# ╠═620dd37f-adac-45ba-b1e1-ae810b5d8e5a
# ╠═ff2bfce6-481f-4b1d-82d0-f1ebb9df615b
# ╠═21c2e289-41d9-4391-b38e-04c41a24849e
# ╠═d5cb6bb8-80be-4e92-9471-4dfed8ec0978
# ╠═2839684e-3bec-4201-a225-5fdb267de92d
# ╠═64e8bf92-36d1-4a0a-8e62-b2b49723d336
# ╠═41bf07bd-cda5-46f8-8093-854850611b16
# ╠═c2904a25-4169-4e54-91de-6a98c09aa9ef
# ╠═d5cd7efe-3d94-4cd1-af20-e8247ee01aa4
# ╠═f77a6c31-4167-40b4-a013-9e9f3ff4cfb3
# ╠═48e75a9b-cda0-4c61-8bca-261526a2a538
# ╠═7913dc93-3a13-4a71-989c-da6e977b9542
# ╠═cce684af-020d-415e-93bb-1aa9d9f3effc
# ╠═f07d27da-756d-49bf-b6a6-bcca0b31f5ad
# ╠═5e75c2ef-369e-4131-a7bf-56a46291c46e
# ╠═1523c30b-d46d-4f24-aa8d-c26383746115
# ╠═24217822-d536-437a-b838-7aed80cccd56
# ╠═52453948-00df-4435-9287-266fcda5d2d4
# ╠═f10372c1-35b8-49d0-a525-c4c6ab28c2a2
# ╠═2cf56e36-c87f-4861-9a5a-360adc114dbd
# ╠═bf93b262-c100-43a2-b56e-59a673bc0c05
# ╠═82195c7f-74a8-467c-b001-dcec35c54574
# ╠═cb4788e3-9d28-48b4-b44f-0638557d98ae
