global xypathspeed = 15000

MkCrv(N)=begin
  c = Any[]
  map(i->push!( c, rand(100,2) ),1:N)
  return c
end

function whichInterval(a,b,numSubIntervals,val)
  dt = (b-a)/(numSubIntervals)
  idx = 1
  if val <= a
    return idx
  elseif val >= b
    return numSubIntervals
  else
    idx = min(floor(Int, (val - a)/dt ) + 1, numSubIntervals)
  end
  if idx < 1
    println("boooooo")
    @show a
    @show b
    @show dt
    @show numSubIntervals
    @show val
    @show idx
  end
  return idx
end

function GroupCurvesBy(meas::Function,crvs::Curves,numSlices::Int=18) #::Vector{Set{Int}}()
  sortedcurvesindices = Vector{Set{Int}}(numSlices)
  if !isempty(crvs)
    vals = map(meas,crvs)
    M = maximum(vals)
    m = minimum(vals)
    for (i,v) in enumerate(vals)
      idx = whichInterval(m,M,numSlices,v)
      try
      typeof(sortedcurvesindices[idx]) == Set{Int}
      catch e
        sortedcurvesindices[idx] = Set{Int}()
      end
      push!(sortedcurvesindices[idx],i)
    end
  end
  return sortedcurvesindices
end

function SortCurvesForDrawing(crvs::Curves)
  sortedcurvestodraw = Curves()
  drawboxes = 18
  if !isempty(crvs)
    xind_slices = GroupCurvesBy(c -> x(c[1]),crvs,drawboxes)
    yind_slices = GroupCurvesBy(c -> y(c[1]),crvs,drawboxes)
    #
    # @show (typeof(xind_slices),size(xind_slices))
    # @show xind_slices
    # @show yind_slices

    for j in 1:drawboxes
      for i in 1:drawboxes

        idx = Set{Int}()
        try
          idx = intersect(xind_slices[i],yind_slices[j])
        catch
        end
          for ix in idx
            push!(sortedcurvestodraw,crvs[ix])
          end
      end
    end
  end
  sortedcurvestodraw
end
#
# function SortCurvesForDrawing(crvs::Curves)
#   sortedcurvestodraw = Curves()
#   if !isempty(crvs)
#     xpos = map(c -> x(c[1]),crvs)
#     ypos = map(c -> y(c[1]),crvs)
#
#     drawboxes = 18
#
#     Mx = maximum(xpos)
#     mx = minimum(xpos)
#     My = maximum(ypos)
#     my = minimum(ypos)
#
#     boxbounds_x = linspace(mx,Mx,drawboxes+1)
#     boxbounds_y = linspace(my,My,drawboxes+1)
#
#     numcurvestodraw = length(crvs)
#
#     xord = Vector{Int}(numcurvestodraw)
#     yord = Vector{Int}(numcurvestodraw)
#
#
#
#     for j in 1:drawboxes
#       for i in 1:drawboxes
#         idxinxbox = (i < drawboxes) ? find(x -> (boxbounds_x[i] <= x < boxbounds_x[i+1]), xpos) : find(x -> (boxbounds_x[i] <= x <= boxbounds_x[i+1]), xpos)
#
#         idxinybox = (j < drawboxes) ?  find(y -> (boxbounds_y[j] <= y < boxbounds_y[j+1]), ypos) : find(y -> (boxbounds_y[j] <= y <= boxbounds_y[j+1]), ypos)
#         idx = intersect(idxinxbox,idxinybox)
#         sortedcurvestodraw = vcat(sortedcurvestodraw,crvs[idx])
#       end
#     end
#   end
#   sortedcurvestodraw
# end


function DrawPath(c::Curves, NewFileName::AbstractString; sortcurves = true)
  crvs = c
  if sortcurves
    crvs = SortCurvesForDrawing(c)
  end
  arr = curves2array(crvs)
  DrawPath(arr, NewFileName)

  # Go ahead and make an svg as well:
  NewFileNameSVG = NewFileName[1:end-4]*".svg"
  DrawSVGCurves(NewFileNameSVG, crvs)
end

function DrawPath{T<:Number}( c::Array{Array{T,2},1}, NewFileName::AbstractString )
  # c is a vector of nx2 paths, nx2xNumPaths
  # For example: see MkCrv.
  NPaths = length( c )
  Dirs = map(i -> c[i][2:end,:]-c[i][1:end-1,:],1:NPaths)
  Starts = map(i -> c[i][1,:], 1:NPaths)
  Ends = map(i -> c[i][end,:], 1:NPaths)
  Transitions=map(i -> Starts[i+1]-Ends[i],1:NPaths-1)
  @show size(Transitions[1])
  if length(c)>1
    push!(Transitions,Starts[1]-Ends[NPaths-1])
  end
  Txt=Any[]
  push!(Txt,"G21") # sets units to mm.
  push!(Txt,"G91") # sets to relative positioning
  push!(Txt,string("G1 X",c[1][1,1]," Y",c[1][1,2]," F$(xypathspeed) "))
  for i in 1:NPaths
    push!( Txt, "G1 Z40.0 F2500.0 ")
    for j in 1:size(Dirs[i])[1]
      push!( Txt, string("G1 X",Dirs[i][j,1]," Y",Dirs[i][j,2]," F$(xypathspeed)"))
    end
    push!( Txt, "G1 Z-40.0 F2500.0 ")
    if i < NPaths
        push!(Txt,string("G1 X",Transitions[i][1]," Y",Transitions[i][2]," F$(xypathspeed)"))
    end
  end
  #push!( Txt, "G1 Z-40.0 F2500.0 ")
  fn = open(NewFileName, "w")
  for i in 1:length(Txt)
    write(fn,Txt[i],"\n")
  end
  close(fn)
end

function makecncfile( c::Curves, NewFileName )
  # c is a vector of nx2 paths, nx2xNumPaths
  # For example: see MkCrv.
  NPaths = length(c)
  crvs = Any[]
  for i in 1:NPaths
    crv = hcat(x(c[i]),y(c[i]))
    push!(crvs,crv)
  end

  Dirs = map(i -> crvs[i][2:end,:]-crvs[i][1:end-1,:],1:NPaths)
  Starts = map(i -> crvs[i][1,:], 1:NPaths)
  Ends = map(i -> crvs[i][end,:], 1:NPaths)
  Transitions=map(i -> Starts[i+1]-Ends[i],1:NPaths-1)
  if length(crvs)>1
    push!(Transitions,Starts[1]-Ends[NPaths-1])
  end
  Txt=Any[]
  push!(Txt,"G21") # sets units to mm.
  push!(Txt,"G91") # sets to relative positioning
  push!(Txt,string("G1 X",crvs[1][1,1]," Y",crvs[1][1,2]," F$(xypathspeed) "))
  for i in 1:NPaths
    push!( Txt, "G1 Z40.0 F2500.0 ")
    for j in 1:size(Dirs[i])[1]
      push!( Txt, string("G1 X",Dirs[i][j,1]," Y",Dirs[i][j,2]," F$(xypathspeed)"))
    end
    push!( Txt, "G1 Z-40.0 F2500.0 ")
    if i < NPaths
      push!(Txt,string("G1 X",Transitions[i][1,1]," Y",Transitions[i][1,2]," F$(xypathspeed)"))
    end
  end
  #push!( Txt, "G1 Z-40.0 F2500.0 ")
  fn = open(NewFileName, "w")
  for i in 1:length(Txt)
    write(fn,Txt[i],"\n")
  end
  close(fn)
end
DrawPath( c )=DrawPath( c, "DrawPathTry.cnc" )
