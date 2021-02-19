#using Color
# using PointsNCurves
# using Pkg
using Random
dir = pwd()
pkgdir = "/Users/prc/github" # Pkg.dir()
using Colors
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
    MakeSVGCurve(closed2curve(c))
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
    Mx = maxx(crvs)
    My = maxy(crvs)
    mx = minx(crvs)
    my = miny(crvs)
    cntr = Point(.5*(mx+Mx),.5*(my+My))
    for c in crvs
        if length(c)>0
            txtc = MakeSVGCurve(c + (-1*cntr))
            txt = string(txt,"\n",txtc)
        end
    end
    MakeSVGFileCentered(NewFile,txt)
    return txt
end

function MakeSVGFileCentered(NewFile,TXT)
    f=try
        open(pkgdir*"/PointsNCurves/src/TestingCentered.svg")
    catch
        open("src/TestingCentered.svg")
    end
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

trnsfrm_regex = r".*\<g transform=\"matrix\(([^,]*),([^,]*),([^,]*),([^,]*),([^,]*),([^,]*)\)"
issvgtransform(l::Union{SubString,String}) = occursin(trnsfrm_regex,l)
parsetransform(l::Union{SubString,String}) = affine(map(x -> tryparse(Float64,x),match(trnsfrm_regex,l).captures))
affine(vec::Vector{Float64}) = [vec[1] vec[3] vec[5]; vec[2] vec[4] vec[6]; 0 0 1]
apply_affine(a::Array{Float64,2},p::Point) = Point(a[1,1]*p.x + a[1,2]*p.y + a[1,3], a[2,1]*p.x + a[2,2]*p.y + a[2,3])

path_regex = r"\<path d=\"M([^\"]*)"
issvgpath(l::Union{String,SubString}) = occursin(path_regex,l)
flipy(p::Point) = Point(p.x,-p.y)
function parse_path(path::Union{String,SubString},transform::Array{Float64,2}=diagm([1.,1.,1.]))
    try
        pthstr = match(path_regex,path).captures[1]
        isclosed = occursin("Z",pthstr)
        pthstr = replace(pthstr,Pair("Z",""))
        parts = split(pthstr,"L")
        c = Curve(map(v -> flipy(apply_affine(transform,Point(tryparse(Float64,v[1]),tryparse(Float64,v[2])))), map(xy -> split(xy,","),parts)))
        isclosed && (push!(c,c[1]))
        c
    catch e
        println(path)
        error(e)
    end
end

function parse_paths(text::Union{String,SubString},transform::Array{Float64,2}=diagm([1.,1.,1.]))
    lines = split(text,"\n")
    c = Curves()
    for l in lines
        if issvgpath(l)
            push!(c,parse_path(l,transform))
        end
    end
    c
end

polyline_regex = r"\<polyline points=\"([^\"]*)\""
issvgpolyline(l::Union{String,SubString}) = occursin(polyline_regex,l)
function parse_polyline(polyline::Union{String,SubString},transform::Array{Float64,2}=diagm([1.,1.,1.]))
    plystr = match(polyline_regex,polyline).captures[1]
    parts = filter(s -> !isempty(s),split(plystr," "))
    c = Curve(map(v -> apply_affine(transform,Point(tryparse(Float64,v[1]),tryparse(Float64,v[2]))), map(xy -> split(xy,","),parts)))
    c
end

path_regex = r"\<path "
issvgpath(l::Union{String,SubString}) = occursin(path_regex,l)
# function parse_path(polyline::Union{String,SubString},transform::Array{Float64,2}=diagm([1.,1.,1.]))
#   plystr = match(polyline_regex,polyline).captures[1]
#   parts = filter(s -> !isempty(s),split(plystr," "))
#   c = Curve(map(v -> apply_affine(transform,Point(tryparse(Float64,v[1]),tryparse(Float64,v[2]))), map(xy -> split(xy,","),parts)))
#   c
# end


closeblock_regex = r"\<\/g\>"
iscloseblock(l::Union{String,SubString}) = occursin(closeblock_regex,l)

function parseSVG(fn::String)::Curves
    f = open(fn)
    lines = readlines(f)
    close(f)
    cs = Curves()
    transform = diagm([1.,1.,1.])
    # assuming depth 3 or less transform
    invtrnsform = zeros(Float64,(3,3,3))
    invtrnsform[:,:,1] = transform
    depth = 1
    for l in lines
        if issvgpolyline(l)
            push!(cs,parse_polyline(l,transform))
        elseif issvgpath(l)
            println("parsing svg paths is not generally supported yet")
        elseif issvgtransform(l)
            ntransform = parsetransform(l)
            invtrnsform[:,:,depth] = inv(ntransform)
            transform = transform * ntransform
            depth += 1
        elseif issvgpath(l)
            push!(cs,parse_path(l,transform))
        elseif iscloseblock(l)
            transform = transform * invtrnsform[:,:,depth]
            depth -= 1
        end
    end
    cs
end

function scribit_draw(Crvs::Union{Curves,ClosedCurves}, fn::String, drawingdir::String="")
    AllCrvs = Vector{typeof(Crvs)}([Crvs])
    scribit_draw(AllCrvs,fn, drawingdir)
end

function draw_many_scribits(AllCrvs::Union{Vector{Curves},Vector{ClosedCurves}}, fn::String, drawingdir::String)
    if check_scribit_bounds(AllCrvs)
        error("curves must have coordinates between 0 and 1024")
    end
    try
        out = run(`ls $drawingdir`);
        @show drawingdir
    catch
        out = run(`mkdir $drawingdir`);
    end
    for (i,C) in enumerate(AllCrvs)
        scribit_draw(C,drawingdir*"/"*fn*"_color_$i.svg","")
    end
    txt = drawmultiscolorscribit(AllCrvs)
    scribit_create_drawing(txt,fn*".svg",drawingdir)
end

function drawmultiscolorscribit(AllCrvs::Union{Vector{Curves},Vector{ClosedCurves}})
        numcolors = length(AllCrvs)
        stylesection = makestyles(numcolors)

        hdr, ftr = scribit_template_header_footer(stylesection)

        body_txt = ""
        for (i,C) in enumerate(AllCrvs)
            subbody = """<g id="M$i">"""
            class_str = """class="st$i" """
            crvs_txt = curves_to_svg_paths(C,class_str)
            subbody = subbody * "\n" * crvs_txt * "\n" * "</g>" * "\n"
            body_txt = body_txt * "\n" * subbody
        end
        txt = string(hdr,"\n",body_txt,"\n",ftr)
        txt
end
prisma_colors = Dict{Int,String}(
    4 => "E20206",
    175 => "FFDB0D",
    55 => "F505FC",
    40 => "4DB6FF",
    88 => "6F753C",
    6 => "FF578C",
    39 => "69D4FF",
    177 => "FF8AD1",
    151 => "AE0765",
    179 => "149ACA",
    104 => "A6ADB5",
    214 => "595126",
    211 => "131817",
    174 => "FCE925",
    98 => "060808",
    181 => "A8C860",
    37 => "43C8B5"
)
function makestyles(numcolors)
    # if numcolors > length(prisma_colors)
    #     error("add more colors to prisma_colors to do more than $(length(prisma_colors)) colors")
    # end
    i = 1
    style = ""

    clrs = map(hex,Colors.sequential_palette(200,numcolors+2;c=1,s=1,w=1,b=1))[2:end-1]
    idx = randperm(numcolors)
    clrs = clrs[idx]
    # for i in 1:numcolors
    #     rvs = rand(3)
    #     rvs = rvs .- .8*minimum(rvs)
    #     rvsn = scls[i] * rvs./maximum(rvs)
    #     push!(clrs,hex(RGB(rvsn...)))
    # end

    for v in clrs
       style *= ".st$i{fill:none;stroke:#$v;stroke-miterlimit:16;}"*"\n\t"
       i += 1
    end
    style
end

function scribit_draw(AllCrvs0::Union{Vector{Curves},Vector{ClosedCurves}}, fn::String="testing.svg", drawingdir::String="")
    if length(AllCrvs0) > 1
        draw_many_scribits(AllCrvs0,fn,drawingdir)
        return
    else
        AllCrvs = Vector{typeof(AllCrvs0[1])}([[],[],[],[]])
        AllCrvs[1] = AllCrvs0[1]

        if check_scribit_bounds(AllCrvs)
            error("curves must have coordinates between 0 and 1024")
        end
        hdr, ftr = scribit_template_header_footer()

        states = [0,1,4,5]
        ms = [1,2,3,4]
        body_txt = ""
        for (i,C) in enumerate(AllCrvs)
            subbody = """<g id="M$(ms[i])">"""
            class_str = """class="st$(states[i])" """
            crvs_txt = curves_to_svg_paths(C,class_str)
            subbody = subbody * "\n" * crvs_txt * "\n" * "</g>" * "\n"
            body_txt = body_txt * "\n" * subbody
        end
        txt = string(hdr,"\n",body_txt,"\n",ftr)
        @show drawingdir
        scribit_create_drawing(txt,fn,drawingdir)
    end
end

function scribit_create_drawing(txt::String,fn::String,drawingdir::String)
    filestring = if drawingdir != ""
        try
            out = run(`ls $drawingdir`)
        catch
            out = run(`mkdir $drawingdir`)
        end
        dir*"/"*drawingdir*"/"*fn
    else
        dir*"/"*fn
    end
    println("drawing the svg to $filestring")
    clipboard(filestring)
    NF=open(filestring,"w")
    write(NF,txt)
    close(NF)
end

function scribit_template_header_footer(stylesection::String="")
    f=try
        open(pkgdir*"/PointsNCurves/SCR_template/scribit_template.svg")
    catch
        println(pwd())
        open("SCR_template/scribit_template.svg")
    end
    lines=readlines(f)
    close(f)
    hdr = ""
    for k in 1:5
        hdr = hdr*lines[k]*"\n"
    end

    styles = stylesection == "" ? join(lines[6:15],"\n") : stylesection

    hdr *= styles*"\n"
    hdr *= lines[16]*"\n"

    ftr = string(lines[end-1],"\n",lines[end])
    hdr, ftr
end

function check_scribit_bounds(AllCrvs::Union{Vector{Curves},Vector{ClosedCurves}})
    bxmax = bymax = 1024
    bxmin = bymin = 0

    allcrvsxmin = bxmax
    allcrvsxmax = bxmin
    allcrvsymin = bymax
    allcrvsymax = bymin
    for C in AllCrvs
        subemp = true
        for c in C
            subemp = subemp && isempty(c)
        end
        if !isempty(C) && subemp
            println(isempty(C))
            println(subemp)

            allcrvsxmin = min(allcrvsxmin, minx(C))
            allcrvsxmax = max(allcrvsxmax, maxx(C))
            allcrvsymin = min(allcrvsymin, miny(C))
            allcrvsymax = max(allcrvsymax, maxy(C))

        end
    end
    allcrvsxmin < bxmin || allcrvsxmax > bxmax || allcrvsymin < bymin || allcrvsymax > bymax
end


function curves_to_svg_paths(C::Union{Curves,ClosedCurves}, class_str::String)
    txt = ""
    for c in C
        txt *= "\n" * curve_to_svg_path(c, class_str)
    end
    txt
end


function curve_to_svg_path(c::Union{Curve,ClosedCurve}, class_str::String = "")
    #<path class="st0" d="M189.5,378.2c-73.6-65.4,13.8-175.1-11.2-199.9c25,25,134.5-62.4,200,11.3"/>
    crv = "<path $class_str d=\"M"
    isfirst = true
    for pnt in c
        if isfirst
            isfirst = false
            crv *= "$(pnt.x),$(pnt.y)"
        else
            crv *= "L$(pnt.x),$(pnt.y)"
        end

    end
    crv *= typeof(crv) == ClosedCurve ? " Z" : ""
    crv *= "\"/>"
    crv
end



function hexfor(code::Int)::String
    clr = prisma_colors[code]
    # clipboard(clr)
    clr
end

function scribitprismarecolor(fn::String,clrs::Vector{Int})
# assuming svg generated using the scribit template
# numbers correspond to the prisma marker numbers
    f = open(fn)
    lns = readlines(f)
    close(f)
    lns[6] = "\t.st0{fill:none;stroke:#"*hexfor(clr1)*";stroke-miterlimit:10;}"
    lns[7] = "\t.st1{fill:none;stroke:#"*hexfor(clr2)*";stroke-miterlimit:10;}"
    lns[10] = "\t.st4{fill:none;stroke:#"*hexfor(clr3)*";stroke-miterlimit:10;}"
    lns[11] = "\t.st5{fill:none;stroke:#"*hexfor(clr4)*";stroke-miterlimit:10;}"
    f = open(fn,"w")
    write(f,join(lns,"\n"))
    close(f)
end
