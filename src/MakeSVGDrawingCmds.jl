#using Color
# using PointsNCurves
# using Pkg
dir = pwd()
pkgdir = "/Users/prc/github" # Pkg.dir()
#cd(WD)
#println("Made it this far!")

function rgbconv(rgb)
  # rgb should be a vector of values between 0 and 1
  if rgb[1] <= 1
    rgb=map(t -> round(Int,255*t),rgb)
  else
    rgb=round(Int,rgb)
  end
  h=""
  for k in 1:3
    #println(k)
    hk = string(rgb[k], base = 16)
    if length(hk)==1
      hk="0"*hk
    end
    h=h*hk
  end
  h
end
#println("rgbconv done")

function MakeSVGShape(c::Curve,clr,sclr)
  sw = .3
  txt = "<polyline points=\""
  for p in c
    txt = string(txt,"$(x(p)),$(-y(p)) ")
  end
    txt = string(txt,"\" style=\"fill:#$(rgbconv(clr));stroke:#$(rgbconv(sclr));stroke-width:$sw\" />")
end
MakeSVGShape(c::Curve,clr) = MakeSVGShape(c,clr,clr)

# function MakeSVGShape(Coords::Curve,Clr::Array)
#   Txt=['<polyline points="',sprintf('%0.2f,%0.2f  ',Coords'),'" ',...
#       'style="stroke-linecap:round;stroke-linejoin: round;fill:#',Clr,';stroke:#',StrokeClr,';stroke-width:',num2str(width),';fill-opacity:',num2str(Opacity),';stroke-opacity:',num2str(StrokeOpacity),';stroke-position:inside" />'];
#   Txt = "<polyline points=\""
#   for c in Coords
#     Txt = string(Txt,"$(c.x), $(c.y)  ")
#   end
#     Txt
#   ',sprintf('%0.2f,%0.2f  ',Coords'),'" ',...
#           'style="stroke-linecap:round;stroke-linejoin: round;fill:#',Clr,';stroke:#',StrokeClr,';stroke-width:',num2str(width),';fill-opacity:',num2str(Opacity),';stroke-opacity:',num2str(StrokeOpacity),';stroke-position:inside" />'];
#   return Txt
# end

function MakeSVGCirc(Center::Point = Point(0.,0.),
                     Rad::Number=25,
                     Clr=rand(1,3);
                     nofill::Bool=true,
                     strkWdth::Number=1/100*Rad)
  Clr=rgbconv(Clr)
  if nofill
    Txt=string("<circle cx=\"",Center.x,"\" cy=\"",Center.y,"\" r=\"",Rad,"\" style=\"stroke:#",Clr,";stroke-width:",string(strkWdth),";fill:none \" />\n")
  else
    Txt=string("<circle cx=\"",Center.x,"\" cy=\"",Center.y,"\" r=\"",Rad,"\" style=\"stroke-width:",string(strkWdth),";fill:#",Clr," \" />\n")
  end
  return Txt
end
#println("MakeSVGCirc done")

function MakeSVGGradCircs(Cents::Array=[-50 -50;50  50],
  Rads::Array = [80 10],
  Clr::Array = rand(2,3),
  n=50,
  nf::Bool=false)
  # If Cents has more than two entries, then size(Cents)[2] will be taken as n
  if size(Cents)[1]>2
    n=size(Cents)[1]
  else
    Cents=[linspace(Cents[1,1],Cents[2,1],n) linspace(Cents[1,2],Cents[2,2],n)]
    Clr=[linspace(Clr[1,1],Clr[2,1],n) linspace(Clr[1,2],Clr[2,2],n) linspace(Clr[1,3],Clr[2,3],n)]
    Rads=linspace(Rads[1],Rads[2],n)
  end
  Txt=[]
  for i in 1:n
    crc=MakeSVGCirc(Cents[i,:],Rads[i],Clr[i,:];nofill=nf)
    Txt=[Txt;crc]
  end
  return Txt
end

function MakeSVGCurve(c::Union{Curve,ClosedCurve},sw::Number=2,clr=[0,0,0])
  if typeof(c)==ClosedCurve
    c = closed2curve(c)
  end
  txt = "<polyline points=\""
  for p in c
    txt = string(txt,"$(x(p)),$(-y(p)) ")
  end
  txt = string(txt,"\" style=\"fill:none;stroke:#$(rgbconv(clr));stroke-width:$sw\" />")
end

function MakeSVGCurve(c::ClosedCurve)
  txt = MakeSVGCurve(closed2curve(c))
end

function DrawSVGCurves(NewFile,c::Curve)
  txt = MakeSVGCurve(c)
  MakeSVGFileCentered(NewFile,txt)
  return txt
end

function DrawSVGCurves(NewFile,c::ClosedCurve)
  txt = MakeSVGCurve(c)
  MakeSVGFileCentered(NewFile,txt)
  return txt
end


function DrawSVGCurves(NewFile,c::ClosedCurves)
  # C=Curves()
  # for cc in c
  #   push!(C,cc.vertices)
  # end
  txt = DrawSVGCurves(NewFile,map(cc->Curve(cc),c))
  return txt
end

function DrawSVGCurves(NewFile,crvs::Curves)
  println("Drawing $(length(crvs)) curves!")
  txt = ""
  for c in crvs
    if length(c)>0
      txtc = MakeSVGCurve(c)
      txt = string(txt,"\n",txtc)
    end
  end
  MakeSVGFileCentered(NewFile,txt)
  return txt
end

function MakeSVGFileCentered(NewFile,TXT)
  f=open(pkgdir*"/PointsNCurves/src/TestingCentered.svg")
  lines=readlines(f)
  close(f)
  hdr = ""
  for k in 1:6
    hdr = hdr*lines[k]*"\n"
  end
  ftr = string(lines[end-1],"\n",lines[end])
  txt = string(hdr,"\n",TXT,"\n",ftr)
  filestring=dir*"/"*NewFile
  println("drawing the svg to $filestring")
  clipboard(filestring)
  NF=open(filestring,"w")
  write(NF,txt)
  close(NF)
end

ispoly(l::String) = l[2:9] == "polyline"

trnsfrm_line = r".*\<g transform=\"matrix\(([^,]*),([^,]*),([^,]*),([^,]*),([^,]*),([^,]*)\)"
istransform(l::Union{SubString,String}) = occursin(trnsfrm_line,l)
parsetransform(l::Union{SubString,String}) = map(x -> tryparse(Float64,x),match(trnsfrm_line,l).captures)
affine(vec::Vector{Float64}) = [vec[1] vec[3] vec[5]; vec[2] vec[4] vec[6]; 0 0 1]

apply_affine(a::Array{Float64,2},p::Point) =

function parseSVG(fn::String)::Curves
  f = open(fn)
  lines = redlines(f)
  close(f)
  c = Curves()
  transform = LinearAlgebra.I
  for l in lines
    if ispoly(l)

    else if istransform(l)

    else if ispath(l)

    end
end
