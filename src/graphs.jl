export
    randomtree,
	spirals

function randomtree(n)
	g = Graph(n)
	seq = rand(1:n-1, n-2)
	degrees = @chain seq begin
		countmap
		[get(_, i, 0) for i=1:n]
		_ .+ 1
	end
	
	for seqi ∈ seq
		j = findfirst(==(1), degrees)
		add_edge!(g, seqi, j)
		degrees[seqi] -= 1
		degrees[j] -= 1
	end
	add_edge!(g, findall(==(1), degrees)...)
	g
end

function spirals(n, k=2; a=1/(2π), b=1, l=2π, δ=0.5)
	θ₀ = find_zero((x -> a-(a*x + b)*tan(x)), .75π)
	arclen(t) = a/2 * (t*sqrt(1+t^2) + asinh(t))

	θ = [find_zero(t -> arclen(t) - ri, 0.5) for ri ∈ rand(k*n)] .* l/2 .+ θ₀
	r = a .* θ .+ b .+ (rand(k*n) .* δ .* 2π/k .* a .- 1)
	y, x = unzip(broadcast.(*, r, sincos.(θ)))

	for i = 2:k
		(i-1)*n+1:i*n |>
			I -> (
				(x[I], y[I]) = LinearAlgebra.rotate!(x[I], y[I], reverse(sincos(2π*(i-1)/k))...)
			)
	end

	x, y, reduce(vcat, [i*ones(Int, n) for i ∈ 1:k])
end
