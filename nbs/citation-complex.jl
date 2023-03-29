### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 67fb4bd4-cd0b-11ed-331e-09af9880a257
using Pkg; Pkg.activate(".")

# ╔═╡ 1d53191c-4880-4b00-bdd7-52eab9b4eb93
using Revise

# ╔═╡ a7cebdbb-6f29-4259-b655-d36f38fdaaa0
using MultiscaleSimplexSignalTransforms

# ╔═╡ de547e0b-ec1f-445e-8b77-2c2ed8ea2679
using SparseArrays, LinearAlgebra, StatsBase, Arpack, Graphs, Plots, Chain, DataStructures, MAT, GraphRecipes, Unzip, Match, JSON3

# ╔═╡ 5eae16db-b13f-463c-b753-0e4f1c6e6c98
using Combinatorics

# ╔═╡ ceaa14a7-8397-4127-bcc5-1e5450ca2da3
# base = "/Users/eug/Downloads/cit-HepTh-abstracts/"

# ╔═╡ 1c4719da-f3a7-4338-baee-1f71b1daa7aa
# authors = [
# 	@chain begin
# 		readlines(joinpath(base, "1992", f))
# 		filter(s -> startswith(s, "Author"), _)
# 		# isempty(_) ? ": NA" : first(_)
# 		first
# 		split(':'; limit=2)
# 		getindex(2)
# 		strip
# 		split("and")
# 		@. strip
# 		split(f, '.')[1] => _
# 	end
# 	for f ∈ readdir(joinpath(base, "1992"))
# ]

# ╔═╡ 24d8d502-7326-443d-a297-4187fa1bd46d
# readdir(joinpath(base, "1992"))

# ╔═╡ 54afad08-cf51-40ab-9245-8a7e829c79fa
# md"------------------------------------------------"

# ╔═╡ c37d9c30-2b30-414c-81b5-b29ca0d18330
papers = "/Users/eug/dev/python/citation/simplicial_neural_networks/data/s2_1_raw"

# ╔═╡ d3840f2a-bc19-49fd-b46c-dae52d1cb782
struct Paper
	authors::Vector{Int}
	citations::Int
	Paper(obj::JSON3.Object) = new(
		[parse(Int, i) for a ∈ obj.authors for i ∈ a.ids],
		length(obj.inCitations)
	)
end

# ╔═╡ effd9508-0c9c-42d7-a9df-7680b6087a88
function readpapers()
	authors = 0
	authormap = DefaultDict{Int, Int}(() -> (authors += 1))
	signal = DefaultDict{Vector{Int}, Int}(()->0)
	for i ∈ 0:0
		for line in eachline(joinpath(papers, "s2-corpus-$(lpad(i,2,"0")).json"))
			p = Paper(JSON3.read(line))
			simplex::Vector{Int} = getindex.(Ref(authormap), p.authors)
			for v ∈ powerset(simplex)
				signal[v] += p.citations
			end
		end
	end
	signal, authormap
end

# ╔═╡ 19c7f9d9-c3c3-4a33-9992-edfccd625539
# q = [
# 	@chain line begin
# 		JSON3.read
# 		Paper
# 	end
# 	for line in eachline(joinpath(papers, "s2-corpus-00.json"))
# ]

# ╔═╡ c4fa7a23-fdaa-4563-bf6d-5719c2810bd3
readpapers()

# ╔═╡ Cell order:
# ╠═67fb4bd4-cd0b-11ed-331e-09af9880a257
# ╠═1d53191c-4880-4b00-bdd7-52eab9b4eb93
# ╠═a7cebdbb-6f29-4259-b655-d36f38fdaaa0
# ╠═de547e0b-ec1f-445e-8b77-2c2ed8ea2679
# ╠═5eae16db-b13f-463c-b753-0e4f1c6e6c98
# ╠═ceaa14a7-8397-4127-bcc5-1e5450ca2da3
# ╠═1c4719da-f3a7-4338-baee-1f71b1daa7aa
# ╠═24d8d502-7326-443d-a297-4187fa1bd46d
# ╠═54afad08-cf51-40ab-9245-8a7e829c79fa
# ╠═c37d9c30-2b30-414c-81b5-b29ca0d18330
# ╠═d3840f2a-bc19-49fd-b46c-dae52d1cb782
# ╠═effd9508-0c9c-42d7-a9df-7680b6087a88
# ╠═19c7f9d9-c3c3-4a33-9992-edfccd625539
# ╠═c4fa7a23-fdaa-4563-bf6d-5719c2810bd3
