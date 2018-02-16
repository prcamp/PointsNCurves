using PointsNCurves
fn = "test.svg"
C = ClosedCurves()

cents = linspace(0, 50, 17)
for c in cents
   push!(C, gensquare(19,Point(c,c)))
end

m,o = masked_occluded(C[1],C[2])


DrawSVGCurves(fn, Curves(vcat(o,[C[2]])))