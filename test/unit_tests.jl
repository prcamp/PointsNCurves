test("basic point tests") do 
	p = Point(0,0)
	q = Point(1.,1.)
	@test norm(p) == 0
	@test isequal(q + (-1)*q,p)
end

test("closed curve tests") do 
	sq = gensquare()
	@test isequal(sq[0],Point(-.5,-.5))
	@test isequal(sq[1],Point(-.5,0.5))
	@test isequal(sq[0],sq[4])
end

test("basic masking/occluding") do
	sq0 = gensquare(1,Point(0.,0.))
	sq1 = gensquare(1,Point(0.5,0.5))
	@time (m,o) = masked_occluded(sq1,sq0)
	@test isequal(typeof(m) , Curves)
	@test isequal(typeof(o) , Curves)
	@test isequal(length(m) , 2)
	@test isequal(m[1] , [Point(0.0, 0.0),
 				   Point(0.0, 0.5)])
    @test isequal(m[2] , [Point(0.5, 0.0),
 				   Point(0, 0)])
    @test isequal(o[1] , [Point(0.0, 0.5),
					Point(0.0, 1.0),
					Point(1.0, 1.0),
					Point(1.0, 0.0),
					Point(0.5, 0.0)])	
end
