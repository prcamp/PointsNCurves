function simplerecursor(trnsfrm,shp::ClosedCurve,combinator,n::Int,H::Curves=Curves())
  (n==0) && (return H)
  shp2 = trnsfrm(shp)
  H0 = combinator(shp,shp2)
  H = vcat(H,H0)
  simplerecursor(trnsfrm,shp2,combinator,n-1,H)
end

function trnsfrm1(shp::ClosedCurve)
    shp2 = .7*(shp+(-1)*(mean(shp)))+mean(shp)
    shp2
end
combinator1(shp1,shp2) = [shp2]
