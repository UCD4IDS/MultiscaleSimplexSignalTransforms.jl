export
    Pos

"""
shorthand for Type Unions
"""
Base.:(|)(A::Type, B::Type) = Union{A,B}
Base.:(|)(A::Type, B::TypeVar) = Union{A,B}
Base.:(|)(B::TypeVar, A::Type) = A|B

spdelta(k, n, v=1.0) = sparsevec([k], [v], n)
delta(k, n, v=1.0) = Vector(spdelta(k, n, v))

struct Pos{R<:Real}
    x::Vector{R}
    y::Vector{R}
    c::Vector{Int}
end

Pos(x, y) = Pos(x, y, Int[])
