export 	
  simplex_path,
  simplex_path_lower,
  triangle_stack,
  tricycle

function simplex_path(n, k=1; flip=false)
  n = n+1
  @assert k<=1 "orders higher than k=1 not implemented yet"
  k == 0 && return (path_digraph(n), √3*(0:n-1), 1:n)
  g = DiGraph(n+k)
  for i=1:n-1
    i % 2 == 1 && flip ?
      add_edge!(g, i, i+1) :
      add_edge!(g, i+1, i)
    i % 2 == 0 && flip ?
      add_edge!(g, i+2, i) :
      add_edge!(g, i, i+2)
  end
  n % 2 == 1 && flip ?
    add_edge!(g, n, n+1) :
    add_edge!(g, n+1, n)
  (
    g, 
    # √3*[0; repeat(1:(n+1)/2, inner=2)][1:end-(n % 2 == 0 ? 0 : 1)],
    # collect(Iterators.flatten(zip(1:n/2+1, 0:n/2)))[1:end-(n % 2 == 1 ? 0 : 1)]
    Pos(
		collect(Iterators.flatten(zip(
			# fill(1:ceil((n+k)/2), 2)...
			(1:ceil((n+k)/2) , (1:ceil((n+k)/2)) .+ 0.5 )...
			)))[1:end-(n+k)%2],
    	repeat([1;1+√3/2], ceil(Int, (n+k)/2))[1:end-(n+k)%2]
	)
  )
end

# an oriented 1-complex, of a stack of triangles.
function triangle_stack(n)
  g = DiGraph((n+1)*(n+2)÷2)
  for l=1:n
    for i=1:l
      add_edge!(g, l*(l-1)÷2+i, l*(l+1)÷2+i)
      add_edge!(g, l*(l+1)÷2+i, l*(l+1)÷2+i+1)
      add_edge!(g, l*(l+1)÷2+i+1, l*(l-1)÷2+i)
    end
  end

  (
    g,
    [(n+l)/2+i-l+1 for l=1:n+1 for i=0:l-1],
    [Float64(n-l) for l=1:n+1 for i=1:l]
  )
end

function simplex_path_lower(n, k=1)
  n == 0 && return (DiGraph(), [], []) # a trivial k-complex with 0 simplices
  k == 0 && return (DiGraph(n), 1:n |> collect, ones(n)) # n disconnected points
  g = DiGraph(n+k)
  for i=1:k, j=i+1:k
    add_edge!(g, i, j)
  end
  for i=1:n, j=1:k
    add_edge!(g, i+k-j, i+k)
  end
  k == 1 && return (g, collect(1:n+k), ones(n+k))

  (
    g,
    Pos(
		collect(Iterators.flatten(zip(
			(1:ceil((n+k)/2) , (1:ceil((n+k)/2)) .+ 0.5 )...
		)))[1:end-(n+k)%2],
    	repeat([1;1+√3/2], ceil(Int, (n+k)/2))[1:end-(n+k)%2]
	)
  )
end

function tricycle(n)
	@assert n ≥ 4 "need larger n"
	g = DiGraph(6n-9)
	
	xind = 0
	yind = n
	x = Float64[]
	y = Float64[]

	# top except last corner
	for j=1:2n
		xind += 1
		yind += j%2==0 ? -1 : 1
		push!(x, xind)
		push!(y, yind)
		add_edge!(g, j, j+1)
	end

	# last corner up top (2n+1)
	xind += 1
	yind += 1
	push!(x, xind)
	push!(y, yind)

	# first outer right side (2n+2)
	xind -= 2
	yind -= 2
	push!(x, xind)
	push!(y, yind)
	
	# connect top to right side
	add_edge!(g, 2n-2, 2n+2)
	add_edge!(g, 2n-2, 2n+3)
	add_edge!(g, 2n, 2n+2)
	add_edge!(g, 2n+2, 2n+3)

	# rest of right side except bottom corner
	for j=2n+3:4n-3
		xind += j%2==0 ? 1 : -2
		yind += j%2==0 ? -1 : 0
		push!(x, xind)
		push!(y, yind)
		add_edge!(g, j, j+1)
	end

	# bottom corner (4n-2)
	xind += 1
	yind += -1
	push!(x, xind)
	push!(y, yind)

	yind += 2
	# connect right side to left side
	add_edge!(g, 4n-5, 4n-1)
	add_edge!(g, 4n-5, 4n)
	add_edge!(g, 4n-3, 4n-1)

	# left side
	for j=4n-1:6n-9
		xind += j%2==0 ? 1 : -2
		yind += j%2==0 ? 1 : 0
		push!(x, xind)
		push!(y, yind)
		add_edge!(g, j, j+1)
	end

	# connect top to left side
	add_edge!(g, 2, 6n-9)
	add_edge!(g, 4, n==4 ? 11 : 6n-10) #edge-case since n=4 has minimum hole
	add_edge!(g, 4, 6n-9)

	# inner and outer paths
	[ add_edge!(g, j, j+2) for j = 1:2n-1 ]
	[ add_edge!(g, j, j+2) for j = 2n+2:4n-4 ]
	[ add_edge!(g, j, j+2) for j = 4n-1:6n-11 ]
	
	(g, Pos(x, √3y))
end
