module PointsNCurves

export Point,
       Polar,
       Segment,
       Curve,
       Curves,
       ClosedCurve,
       ClosedCurves,
       curve,
       closedcurve,
       polar,
       closed2curve,
       cleancurve,
       curves2array,
       d,
       x,
       y,
       point2vec,
       norm,
       minx,
       miny,
       maxx,
       maxy,
       genrandpoint,
       genrandcrv,
       genrandcirc,
       gencirc,
       genstar,
       gentestsquare,
       gensquare,
       boxit,
       orientation,
       mask_sequence,
       occlude_sequence,
       masked_occluded,
       inshape,
       strictlyinshape,
       mightintersect,
       segmentintersect,
       onsegment,
       rgbconv,
       MakeSVGShape,
       MakeSVGCirc,
       MakeSVGGradCircs,
       MakeSVGCurve,
       DrawSVGCurves,
       MakeSVGFileCentered,
       simplerecursor,
       trnsfrm1,
       combinator1,
       simplerecursor2,
       trnsfrm2,
       combinator2,
       whichInterval,
       GroupCurvesBy,
       SortCurvesForDrawing,
       DrawPath,
       parsecnc




include("points_and_curves_interfaces.jl")
include("points_and_curves_utils.jl")
include("p&c_recursors.jl")
include("MakeSVGDrawingCmds.jl")
include("XYPlotterFcnsfor3.0.jl")
end