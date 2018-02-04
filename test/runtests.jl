using SANA
using Base.Test

function test1()
    n = 10; G = sprand(n,n,0.3); G = G - Diagonal(G); G = G' + G; G.nzval[:] = 1.0; 
    sanas3(G,G,maxtime=10)
    true
end

@test test1()

