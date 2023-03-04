#=
TestingTestTesting.jl
Jacob Child
Mar 3, 2023
Pseudocode: Use the Test package to run tests on a file and learn how 
testing works generally.

=#

#Packages to use
using Test

#Super simple example the format is generally @test statement == expected outcome
@test 1 == 1 #if running just this line it will show Test Passed, running more than 1 line it just continues
#@test 1 == 2 #Test Failed-> the code stops here and doesn't run the rest of the file
@test 2+2 == 4
#println("did we get here?") #? Ans: No, it looks like it stops once a test fails 
#any print statement will cover up whatever ran most recently, but if it gets to that statement, it means it passed the test

#Further formating @test f(args...) key=val... or @test f(args..., key=val...)
#you can see if things are approximatley close with @test a ≈ b (≈ is \approx + tab) or say isapprox
#you can set an absolute tolerance with atol= and relative tolerance with rtol=
@test pi ≈ 3.14 atol=.1
@test isapprox(pi, 3.14, atol=.1)
@test isapprox(pi, 3.14, rtol=.1) #? what is relative tolerance vs absolute?
#test relative tolerance vs absolute to see if relative is a percent 
@test isapprox(2, 2.02, rtol=.01) #true because it is within 1% of 2
#@test isapprox(2, 2.03, rtol=.01) #false because it is not within 1% of 2
@test isapprox(2, 2.01, atol=.01) #false because 2.02 is not within .02 of 2
#? Ans: relative tolerance is a percent absolute tolerance can be considered less than or equal to the difference or a bounds

#?what are the default okay bounds?
@test 1 ≈ .99999999 #8 decimal places and up of similiarity passes
#@test 1 ≈ .9999999 #below 8 decimal places fails without setting a tolerance
#test if it is within a number of places
#@test 1000000000 ≈ 1000000099 # this fails even though 8/10 digits match 
@test 1000000000 ≈ 1000000000.000001 # this passes even though it is different at 6 decimal places
#?see if the default is a relative difference?
percentpass = (1 - .99999999) #according to above
#println("Percent to Pass is approximate by default: ", percentpass*100, "%") #0.000001% difference
percentfail = (1 - .9999999) #according to above
@test 1000000000 ≈ 1000000000*(1+percentpass) #passes
#@test 1000000000 ≈ 1000000000*(1+percentfail) #fails, so the default is a relative difference
#? Ans: Yes, the default is a relative difference 

#Other conditions: broken and skip 
#broken -> the format is @test a = b broken=true 
@test 1 == 2 broken=true #this should fail, but because of the broken=true condition it is flagged as broken and runs the next line
#@test 1 == 1 broken=true #this throws an error and does not run the next line because it is an "Unexpected Pass"
@test 1 == 1
#skip -> the format is @test a = b skip=true
@test 1 == 2 skip=true #this should fail, but because of the skip=true condition it is flagged as broken and runs the next line
@test 1 == 1 

#Another Macro: @test_throws looks for specific error messages  
#the format is @test_throws ExceptionType statement
@test_throws DomainError sqrt(-1) #this passes as it expects a DomainError and gets one
#note we don't see the error message as Test Passed hides it
#you can still test full expressions with @test_throws and see if the expected error happens in the function 
foo(x) = length(x)^2 #inline function foo that squares the length of the input x
#@test_throws DomainError foo("hello") #this fails as it expected a domain error and there was no error at all 
@test_throws MethodError foo(:cat) #this passes as length() doesn't work with symbols and throws a method error 

#Working with @testset -> all the tests in the set are run and then a summary of results is printed after 
#format for @testset is @testset "name" begin ... end;
@testset "Testing sqrt" begin
    @test sqrt(4) == 2
    @test sqrt(9) == 3
    @test sqrt(16) == 4
end #this all passes and a summary is given 
#you can nest testsets
@testset "Testing" begin
    @testset "ints" begin
        #@test 1 == 2 #fail
        @test 3 == 3 #pass
        @test 1 == 7 broken=true #Broken
    end
    @testset "floats" begin
        #@test 1.0 == 1.1 #fail
        @test 2.2222 == 2.2222 #pass
        @test_throws DomainError sqrt(-1) #pass
        #@test sqrt(-1) #Error
        @test sqrt(-1) skip=true #Broken
    end  
end #everything runs and a summary is given, if anything in the summary fails or errors it does not run past this set
#testset can run for loops 
@testset "Is blah even? -> see it in two ways!" begin
    for i in 1:5
        #@test iseven(i)
        @test (mod(1,2) == 0) #!there is an error in here and I can't find out why
    end #this all runs and a summary is given
end