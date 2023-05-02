#=
Functions to make testing more convinient. 
=#

err(x, xt) = x.-xt
relerr(x, xt) = (x.-xt)./xt
function RMS(x, xt)
    diff = @. (x -xt)^2
    return sqrt(sum(diff)/length(x))
end