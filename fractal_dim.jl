using StaticArrays

import Base: in, length, split, intersect

const Point = SVector{3}
const TOL = 1E-13

struct Segment
    A::Point
    B::Point
end

struct Box
    X::Point
    Y::Point

    function Box(X, Y)
        @assert X[1] < Y[1] && X[2] < Y[2] && X[3] < Y[3]
        new(X, Y)
    end
end

function split(box::Box)::NTuple{8, Box}
    x, y = box.X, box.Y
    m = 0.5*(x+y)

    return (Box(Point(x[1], x[2], x[3]), Point(m[1], m[2], m[3])),
            Box(Point(x[1], m[2], x[3]), Point(m[1], y[2], m[3])),
            Box(Point(m[1], x[2], x[3]), Point(y[1], m[2], m[3])),
            Box(Point(m[1], m[2], x[3]), Point(y[1], y[2], m[3])),
            Box(Point(x[1], x[2], m[3]), Point(m[1], m[2], y[3])),
            Box(Point(x[1], m[2], m[3]), Point(m[1], y[2], y[3])),
            Box(Point(m[1], x[2], m[3]), Point(y[1], m[2], y[3])),
            Box(Point(m[1], m[2], m[3]), Point(y[1], y[2], y[3])))
end

# Bounding boxes of tdim > 0 entities
bounding_box(thing::Box) = thing

function bounding_box(line::Segment)
    X = Point((min(line.A[i], line.B[i]) for i in 1:3)...)
    Y = Point((max(line.A[i], line.B[i]) for i in 1:3)...)

    Box(X, Y)
end

# collision of bounding boxes / things
function collides(b::Box, B::Box)
    !(b.Y[1]+TOL < B.X[1] || B.Y[1] < b.X[1]-TOL ||                                                                         
      b.Y[2]+TOL < B.X[2] || B.Y[2] < b.X[2]-TOL ||                                                                         
      b.Y[3]+TOL < B.X[3] || B.Y[3] < b.X[3]-TOL)   
end

collides(a::Segment, b::Box) = collides(bounding_box(a), b)
collides(a::Box, b::Segment) = collides(a, bounding_box(b))

# point in objects
function in(p::Point, line::Segment)
    t = dot(p-line.A, line.B-line.A)/length(line)
    -TOL < t < 1. + TOL
end

function in(p::Point, box::Box)
    ((box.X[1] - TOL < p[1] < box.Y[1] + TOL) &&
     (box.X[2] - TOL < p[2] < box.Y[2] + TOL) &&
     (box.X[3] - TOL < p[3] < box.Y[3] + TOL))
end

function intersect(line::Segment, box::Box)
    # No intersection if both outside
    (line.A ∉ box && line.B ∉ box) && return false
    # Otherwise ok
    true
end

# using PyCall

# unshift!(PyVector(pyimport("sys")["path"]), "/usr/lib/python2.7/dist-packages/")
# @pyimport fenics

# mesh = fenics.Mesh("vasc_mesh.xml.gz")
# @assert py"($mesh).topology().dim()" == 1
# @assert py"($mesh).geometry().dim()" == 3

# X = py"$mesh.coordinates()"

# segments = Vector{Segment}()
# # This is a very clumsy way of reading in the mesh
# cells = fenics.cells(mesh)
# for c in cells
#     v0, v1 = py"$c.entities(0)"
#     v0 += 1
#     v1 += 1
#     segment = Segment(Point(X[v0, 1], X[v0, 2], X[v0, 3]),
#                       Point(X[v1, 1], X[v1, 2], X[v1, 3]))
#     push!(segments, segment)
# end

# xmin = minimum(X, 1)
# xmax = maximum(X, 1)

xmin = [0, 0, 0.]
xmax = [1, 1, 1.]

segments = [Segment(Point(xmin...), Point(xmax...))]

root = Box(Point(xmin...), Point(xmax...))
# Root pixel
size = reduce(*, xmax-xmin)^(1./3)

oldgen = [root]
level = 0
while level < 10
    tic()
    count = 0
    newgen = Vector{Box}()
    # Loop over current generation
    for box in oldgen
        println(box)
        # Box is intesected by any segment
        for line in segments
            if intersect(line, box)
                count += 1
                println("isect")
                # Childen will be visitid in the next round
                for child in split(box)
                    push!(newgen, child)
                end
                break
            end
        end
    end
    dt = toq()
    # Repoints
    @show (level, count, size, dt, 8^level, count/8^level)

    oldgen = newgen
    @show length(oldgen)
    size /= 2
    level += 1
end
