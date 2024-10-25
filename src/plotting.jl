export splot, splot_many,
    tree_layout

function splot(
    c_vertices::Vector, c_edges::Vector, c_triangles::Vector, f::AbstractVector{R};
    k::Int, xs::Vector{S}, ys::Vector{S},
    clims=nothing,
    simpsize=nothing,
    palette=:seismic,
    tol=1e-12,
    notribds=true,
    scatterkwargs...
) where {R<:Real,S<:Real}
    v_xs = xs[c_vertices]
    v_ys = ys[c_vertices]

    e_coor = [[(v_xs[i], v_ys[i]) for i ∈ e] for e ∈ c_edges]
    t_coor = [[(v_xs[i], v_ys[i]) for i ∈ t] for t ∈ c_triangles]

    simpsize = ifelse(
        isnothing(simpsize),
        (lw=6, ms=4),
        (lw=simpsize, ms=simpsize)
    )
    isnode, isedge, istri = k .== 0:2
    isconst = diff(collect(extrema(f)))[1] < tol

    balancedextrema(f) = extrema(f) .|> abs |> maximum |> m -> (-m, m)
    asunit(f, mM=balancedextrema(f)) = mM |> ((m, M),) -> (f .- m) ./ (M - m)

    c = isconst ?
        cgrad(palette)[0.5] :
        (cgrad(palette)[isnothing(clims) ? asunit(f) : asunit(f, clims)] |> permutedims)

    plot(
        Shape.(t_coor);
        legend=:none,
        framestyle=:none,
        aspect_ratio=:equal,
        fillalpha=istri ? 1 : 0,
        linealpha=0,
        c=istri ? c : nothing
    )
    plot!(
        Shape.(e_coor);
        lw=isedge ?
           permutedims(map(v -> v == 0 ? simpsize.lw / 2 : simpsize.lw, f)) :
           (notribds ? 0 : 1),
        lc=isedge ? c : :black
    )
    scatter!(
        v_xs,
        v_ys;
        ms=isnode ? map(v -> v == 0 ? simpsize.ms / 2 : simpsize.ms, f) : simpsize.ms,
        alpha=isnode ? 1 : 0,
        c=cgrad(palette),
        msw=0,
        mz=isnode ? (isconst ? fill(0.0, length(f)) : f) : nothing,
        clims=isnothing(clims) ? balancedextrema(f) : clims,
        scatterkwargs...
        # mz = isnode ? asunit(vals) : nothing
    )
end

function splot(
    c_edges::Vector, f::AbstractVector{R};
    k::Int=0, xs::Vector{S}, ys::Vector{S},
    clims=nothing,
    simpsize=nothing,
    palette=:seismic,
    tol=1e-12,
    scatterkwargs...
) where {R<:Real,S<:Real}
    @assert k == 0 "can only plot nodes for a 0-complex"

    e_coor = [[(xs[i], ys[i]) for i ∈ e] for e ∈ c_edges]

    simpsize = ifelse(
        isnothing(simpsize),
        (lw=1, ms=4),
        (lw=1, ms=simpsize)
    )
    isconst = diff(collect(extrema(f)))[1] < tol

    balancedextrema(f) = extrema(f) .|> abs |> maximum |> m -> (-m, m)

    plot(
        Shape.(e_coor);
        legend=:none,
        framestyle=:none,
        aspect_ratio=:equal,
        simpsize.lw,
        lc=:black
    )
    scatter!(
        xs,
        ys;
        ms=map(v -> v == 0 ? simpsize.ms / 2 : simpsize.ms, f),
        alpha=1,
        c=cgrad(palette),
        msw=0,
        mz=(isconst ? fill(0.0, length(f)) : f),
        clims=isnothing(clims) ? balancedextrema(f) : clims,
        scatterkwargs...
        # mz = isnode ? asunit(vals) : nothing
    )
end

splot(
    region::R, f::AbstractVector{S}, pos::Pos; kwargs...
) where {R<:Region,S<:Real} = splot(region, f; xs=pos.x, ys=pos.y, kwargs...)

function splot(
    region::KRegion{KL,K}, f::AbstractVector{R}; k=KL, kwargs...
) where {KL,K,R<:Real}
    @assert K ∈ (2, 3) "need both edges and triangles in the complex"
    c_edges, c_triangles = K == 2 ?
                           (simplices(region), hulls(region)) :
                           (boundaries(region), simplices(region))
    c_vertices = first.(simplices(region.tree, 0))
    splot(c_vertices, c_edges, c_triangles, f; k, kwargs...)
end

splot(
    region::ZeroRegion, f::AbstractVector{R}; k=0, kwargs...
) where {R<:Real} = splot(hulls(region), f; k, kwargs...)

splot(
    region::BipartiteRegion, f::AbstractVector{R}; k=0, kwargs...
) where {R<:Real} = splot(hulls(region), f; k, kwargs...)

function splot_many(
    region::KRegion{KL,K}, fs::AbstractVector{R}...; xs, ys, k=KL, kwargs...
) where {KL,K,R<:Real}
    @assert K ∈ (2, 3) "need both edges and triangles in the complex"
    c_edges, c_triangles = K == 2 ?
                           (simplices(region), hulls(region)) :
                           (boundaries(region), simplices(region))
    c_vertices = first.(simplices(region.tree, 0))
    nf = length(fs)
    nvert = length(simplices(region.tree, 0))
    repeatsimp(cs, repval=nvert) = [c .+ (i - 1) * repval for i ∈ 1:nf for c ∈ cs]
    rv = diff(collect(extrema(ys)))[1]
    rv = rv == 0 ?
        :simpsize ∈ keys(kwargs) ?
        	kwargs[:simpsize] :
        	1 :
        rv

    splot(
        repeatsimp.([c_vertices, c_edges, c_triangles])...,
        [fi for f ∈ fs for fi ∈ f];
        xs=repeat(xs, nf),
        ys=repeatsimp(ys, -1.1 * rv),
        k, kwargs...
    )
end

function splot_many(
    region::ZeroRegion, fs::AbstractVector{R}...; xs, ys, k=0, kwargs...
) where {R<:Real}
    c_edges = map(e -> Tuple(e)[1:2], edges(region.weights))
    nf = length(fs)
    repeatsimp(cs, repval=nv(region.weights)) = [c .+ (i - 1) * repval for i ∈ 1:nf for c ∈ cs]

    rv = diff(collect(extrema(ys)))[1]
    rv = rv == 0 ?
    	:simpsize ∈ keys(kwargs) ?
        	kwargs[:simpsize] :
        	1 :
        rv

    splot(
        repeatsimp(c_edges),
        [fi for f ∈ fs for fi ∈ f];
        xs=repeat(xs, nf),
        ys=repeatsimp(ys, -1.1 * rv),
        k, kwargs...
    )
end

splot_many(
    region::R, pos::Pos, fs::AbstractVector{S}...; kwargs...
) where {R<:Region,S<:Real} = splot_many(region, fs...; xs=pos.x, ys=pos.y, kwargs...)

tree_layout(g::AbstractGraph, usedfs=true) = rooted(g, usedfs) |> W -> (
    ZeroRegion(W),
    Pos(unzip(Tuple.(buchheim(triu(W))))...)
)

function tree_layout(t::SimplexTree)
    ss = simplices(t)
    vs = first.(filter(s -> length(s) == 1, ss))
    d = Dict(s => i for (i, s) ∈ enumerate(ss))
    adjlist = [
        filter(!isnothing, [
            get(d, vcat(s, v), nothing)
            for v ∈ vs if (isempty(s) || v > s[end])
        ])
        for s ∈ ss
    ]
    pos = Pos(unzip(Tuple.(buchheim(adjlist)))...)
    reg = @chain adjlist begin
        sparse(
            reduce(vcat, [fill(i, length(v)) for (i, v) ∈ enumerate(_)]),
            reduce(vcat, _),
            1.0, length(_), length(_)
        )
        Symmetric
        Graph
        ZeroRegion
    end
    (reg, pos)
end
