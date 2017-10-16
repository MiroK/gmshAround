using FractionalDim

function KochFractal(niters::Int)
    @assert niters > 0
    
    #segments = Vector{Segment{2}}([Segment(Point(0, 0), Point(0.5, sqrt(3)/2)),
    #                               Segment(Point(0.5, sqrt(3)/2), Point(1, 0)),
    #                               Segment(Point(1, 0), Point(0, 0))])

    segments = Vector{Segment{2}}([Segment(Point(0, 0), Point(1, 0))])
    level = 0
    while level < niters
        new_segments = Vector{Segment{2}}()
        for segment in segments
            append!(new_segments, KochFractal(segment))
        end
        segments = new_segments
        level += 1
    end
    segments
end

function KochFractal(segment::Segment{2})
    # Split the thirds
    p0, p1, p2, p3 = segment.O, segment(1./3), segment(2./3), segment.B
    L = segment.length/3

    midpoint = 0.5*(p1 + p2)
    lift = Point(-segment.t[2], segment.t[1])
    tip = midpoint + (sqrt(3)/2)*L*lift
    
    (Segment(p0, p1), Segment(p1, tip), Segment(tip, p2), Segment(p2, p3))
end



# Should be 1
#l = [Segment(Point(0, 0), Point(1, 1))]
#l = [Segment(Point(0., 0.3), Point(1., 0.3))]

koch = KochFractal(7)

using PyCall
@pyimport matplotlib.pyplot as plt


counter = BoxCounter(koch)

sizes, Ns = Vector{Float64}(), Vector{Int}()

state = start(counter)

# Note this is necessary because straigh line does not have well defined
# bbox
# state = (1., 1, [Box(Point(0, 0), Point(1, 1))])

while !done(counter, state)
    tic()
    state = next(counter, state)
    dt = toq()
    
    size, N = state[1], state[2]
    push!(sizes, size)
    push!(Ns, N)

    if length(sizes) > 1
        estimate = -first([log.(sizes) ones(sizes)]\log.(Ns))
        estimate2 = -log(Ns[end]/Ns[end-1])/log(sizes[end]/sizes[end-1])
    else
        estimate = NaN
        estimate2 = NaN
    end
    println((size, N, estimate, estimate2, dt, length(last(state))))

    # plt.figure()
    # for line in koch
    #     plt.plot([line.O[1], line.B[1]], [line.O[2], line.B[2]], "k")
    # end

    # for box in last(state)
    #     for line in surface(box)
    #         plt.plot([line.O[1], line.B[1]], [line.O[2], line.B[2]], "b")
    #     end
    # end
    # plt.show()
end

# FIXME: Koch - does it work
#        XML reader
#        3d algorithms
#        the real deal
