alphabet_stream = Iterators.map(
  cs -> reverse(String([c for c in cs])),
  Iterators.flatten(
    Iterators.map(
      n -> Iterators.product(fill('a':'z', n)...),
      Iterators.countfrom()
    )
  )
)

number_stream = Iterators.countfrom()

_plotg_defaults = Dict{Symbol, Any}(
  :mc => colorant"skyblue",   # node color
  :lc => colorant"black",     # line color
  :ms => 10,                  # node size
  :edgemargin => 209,         # margin betwen edges and nodes
  :elo => 0.05,               # edge label offset
  :ew => 2,                   # edge width
  :ts => 8,                   # text size
  :lz => nothing,             # z-level for edge color
  :nls => alphabet_stream,    # labels for nodes (as iterator)
  :els => number_stream       # labels for edges (as iterator)
);

function _plotg_set_defaults(kwpairs)
  kwargs = Dict{Symbol, Any}(kwpairs)
  for (k, v) in pairs(_plotg_defaults)
    get!(kwargs, k, v)
  end
  kwargs
end

# plot a graph. Work in progress until I replace this directly with a Plots recipe
function plotg(
  g::AbstractGraph, x::Array, y::Array; kwpairs...
)
  kwargs = _plotg_set_defaults(kwpairs)
  xplt = []
  yplt = []
  els = isnothing(kwargs[:els]) ? fill("", ne(g)) : first(number_stream, ne(g))
  label_pos = Tuple[]
  for (i, e) in enumerate(edges(g))
    s, d = Tuple(e)
    any(isnan.([x[s], y[s], x[d], y[d]])) && continue
    sθ, cθ = sincos(atan(y[d]-y[s], x[d]-x[s]))
    
    xₑ = kwargs[:ms]/kwargs[:edgemargin] * cθ
    yₑ = kwargs[:ms]/kwargs[:edgemargin] * sθ

    xₐ = (x[s] + x[d])/2 + kwargs[:elo] * sθ
    yₐ = (y[s] + y[d])/2 - kwargs[:elo] * cθ
    
    push!(xplt, x[s]+xₑ, x[d]-xₑ, NaN)
    push!(yplt, y[s]+yₑ, y[d]-yₑ, NaN)
    ei = els[i]; !isnothing(kwargs[:els]) && 
      push!(label_pos, (xₐ, yₐ, text(L"\mathbf{%$ei}", kwargs[:ts])))
  end
  plt = plot(
    xplt, yplt;
    legend=false,
    colorbar=true,
    framestyle=:none,
    arrow=true,
    lc=kwargs[:lc],
    lw=kwargs[:ew],
    lz=isnothing(kwargs[:lz]) ? nothing : repeat(kwargs[:lz], inner=3),
    clims = isnothing(kwargs[:lz]) ? nothing : (
      maximum(abs.(extrema(kwargs[:lz]))) |> x -> (-x, x)
    )
  )
  scatter!(
    x[.!isnan.(x)], y[.!isnan.(y)];
    mc = kwargs[:mc],
    ms = kwargs[:ms],
    msw = kwargs[:ew],
    text = isnothing(kwargs[:nls]) ?
      nothing :
      [text(L"\mathbf{%$a}", kwargs[:ts]) for a in first(kwargs[:nls], nv(g))]
  )
  !isnothing(kwargs[:els]) && annotate!(label_pos)
  plt
end
