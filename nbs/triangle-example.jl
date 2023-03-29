### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ db92dbbb-b8fb-4811-add1-700118bdcdd0
using Pkg; Pkg.activate(".")

# ╔═╡ 59dc1312-a723-48d8-af81-00ce04031243
using Revise

# ╔═╡ 7a3a232e-c941-11ed-29ca-d9152e333d5d
using Plots, GraphRecipes, Graphs, LaTeXStrings, Unzip, Chain

# ╔═╡ b61067d9-e6f2-488e-9a94-472cc219b6cb
using MultiscaleSimplexSignalTransforms

# ╔═╡ 70902d55-29cb-4dae-b541-345149c0587c
margin = 0.02

# ╔═╡ a5496645-4f45-4b6c-ba70-38e69d952b3b
function triplot(margin, textmargin=margin; plotkwargs...)
	m = margin
	ym, xm = m.*sincos(π/3)
	tx, ty = textmargin.*sincos(π/3)
	lw = 3
	plot(
		[
			m; 1-m; NaN;
			0.5+xm; 1-xm; NaN;
			xm; 0.5-xm; NaN
		],
		[
			0; 0; NaN;
			1/sqrt(2)-ym; ym; NaN; 
			ym; 1/sqrt(2)-ym; NaN
		];
		arrow=Plots.arrow(:simple, :head),
		grid=:hide,
		showaxis=:no,
		legend=false,
		color=:black,
		linewidth=lw,
		annotations=[
			(0.25-tx, 1/(2sqrt(2))+ty, (L"e", 30)),
			(0.75+tx, 1/(2sqrt(2))+ty, (L"f", 30)),
			(0.5, -sqrt(2)ty, (L"g", 30)),
			(0.5, 0.25, (L"s", 30))
		],
		ylims=(-.1,.7),
		xlims=(-.05, .755).* (3/2)
	)
	cy, cx = 0.1 .* collect.(unzip(sincos.(3π/4:.01:9π/4)))
	plot!(
		cx .+ 0.5, cy .+ 0.25;
		color=:black,
		linewidth=lw,
		arrow=(:simple, :tail)
	)
end

# ╔═╡ aca59676-10fe-4e32-880a-bd3cf90e678b
@chain begin
	triplot(margin, 5margin)
	# @aside savefig("/Users/eug/Desktop/tri-consistency.png")
end

# ╔═╡ 2c033c1b-2820-434b-9c69-fc99638dd3ea
g = randomtree(20)

# ╔═╡ a2476fb6-46c0-4ab2-88dc-b04789d1b867
reg = KRegion(g, 1; weakadjweight=1.0)

# ╔═╡ cdb7459a-7fdb-422e-9b03-46cfb8a993ea
part = kGHWT(reg)

# ╔═╡ 586ea482-306e-4245-9332-fd4a02f42907
part[0,0,20]

# ╔═╡ Cell order:
# ╠═db92dbbb-b8fb-4811-add1-700118bdcdd0
# ╠═7a3a232e-c941-11ed-29ca-d9152e333d5d
# ╠═70902d55-29cb-4dae-b541-345149c0587c
# ╠═a5496645-4f45-4b6c-ba70-38e69d952b3b
# ╠═aca59676-10fe-4e32-880a-bd3cf90e678b
# ╠═59dc1312-a723-48d8-af81-00ce04031243
# ╠═b61067d9-e6f2-488e-9a94-472cc219b6cb
# ╠═2c033c1b-2820-434b-9c69-fc99638dd3ea
# ╠═a2476fb6-46c0-4ab2-88dc-b04789d1b867
# ╠═cdb7459a-7fdb-422e-9b03-46cfb8a993ea
# ╠═586ea482-306e-4245-9332-fd4a02f42907
