using FractionalDim

# Should be 1
l = Segment(Point(0, 0), Point(1, 1))
l = Segment(Point(0., 0.3), Point(1., 0.3))

counter = BoxCounter([l])

sizes, Ns = Vector{Float64}(), Vector{Int}()

#state = start(counter)

# Note this is necessary because straigh line does not have well defined
# bbox
state = (1., 1, [Box(Point(0, 0), Point(1, 1))])
while !done(counter, state)
    state = next(counter, state)

    size, N = state[1], state[2]
    push!(sizes, size)
    push!(Ns, N)

    if length(sizes) > 1
        estimate = -first([log.(sizes) ones(sizes)]\log.(Ns))
    else
        estimate = NaN
    end
    println((size, N, estimate))
end

# FIXME: Koch - does it work
#        XML reader
#        3d algorithms
#        the real deal
