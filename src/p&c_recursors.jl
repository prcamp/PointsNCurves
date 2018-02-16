function simplerecursor(trnsfrm::Function,shp::ClosedCurve,combinator,n::Int,H::Curves=Curves())
  (n==0) && (return H)
  shp2 = trnsfrm(shp)
  H0 = combinator(shp,shp2)
  H = vcat(H,H0)
  simplerecursor(trnsfrm,shp2,combinator,n-1,H)
end

trnsfrm1(shp::ClosedCurve) = .7*(shp+(-1)*(mean(shp)))+mean(shp)
combinator1(shp1,shp2) = [shp2]
