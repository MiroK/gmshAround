using StaticArray

import Base: length, split, in

const Point = SVector

struct Segment{D}
    A::Point{D}
    B::Point{D}

    function S
end

length{D}(line::Segment{D}) = norm(line.B - line.A)
tangent{D}(line::Segment{D}) = line.B - line.A

struct Rectangle{D}
    X::Point{D}
    Y::Point{D}

    function Rectangle(X, Y)
        # Don't allow for degenerate
        @assert X[1] < Y[1] && X[2] < Y[2]
        
    end
end
