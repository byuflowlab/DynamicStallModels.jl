testupdown = [1 1; 2 2; 3 3; 4 4; 5 5; 4 3; 3 2; 2 1; 1 0; 0 -1;]
testup = zeros(length(testupdown[:,1]), 2)
testdown = zeros(length(testupdown[:,1]), 2)
#testup = [Float64[] Float64[];]
#testdown = [Float64[] Float64[];]
#push!(testup[1,:], testupdown[1,:])
for i in 2:length(testupdown[:,1])
    if testupdown[i,1] > testupdown[i-1,1]
        testup[i,:] = testupdown[i,:]
        println("top if")
    elseif testupdown[i,1] < testupdown[i-1,1]
        testdown[i,:] = testupdown[i,:]
        println("elseif")
    else 
        println("I am here")
    end
end