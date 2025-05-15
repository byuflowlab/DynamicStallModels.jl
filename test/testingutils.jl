#=
Functions to make testing more convinient. 
=#

calculate_error(x, xt) = x.-xt
calculate_relerr(x, xt) = (x.-xt)./xt
function RMS(x, xt)
    diff = @. (x -xt)^2
    return sqrt(sum(diff)/length(x))
end