using DelimitedFiles, FLOWMath

lift = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/static.csv", ',')
drag = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/staticdrag.csv", ',')
moment = readdlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/staticmoment.csv", ',')

liftfit = Akima(lift[:,1], lift[:,2])
dragfit = Akima(drag[:,1], drag[:,2])
momentfit = Akima(moment[:,1], moment[:,2])

alphavec = collect(-2.7:0.1:30)
liftvec = liftfit.(alphavec)
dragvec = dragfit.(alphavec)
momentvec = momentfit.(alphavec)

polar = hcat(alphavec, liftvec, dragvec, momentvec)

# writedlm("/Users/adamcardoza/Box/research/FLOW/bladeopt/experimentaldata/Hansen2004/figure9a_nonlinearsolver/polar.csv", polar, ',')