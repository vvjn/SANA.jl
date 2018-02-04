# SANA

[![Build Status](https://travis-ci.org/vvjn/SANA.jl.svg?branch=master)](https://travis-ci.org/vvjn/SANA.jl) [![codecov.io](http://codecov.io/github/vvjn/SANA.jl/coverage.svg?branch=master)](http://codecov.io/github/vvjn/SANA.jl?branch=master) 

SANA performs network alignment using simulated annealing.
SANA can optimize the S3[1] and WEC network alignment measures, the
DS3[2] and DWEC[2] dynamic network alignment measures, and node conservation
(i.e. the assignment problem), as well as convex combinations of any
of these measures.

SANA.jl is a partial implementation of this paper [3]. SANA.jl
is similar to the paper in that it
uses the same swapping mechanism to find neighbors of an alignment,
and it optimizes some of the same network alignment measures. It is
different from the paper in that it calculates initial temperature
(`k`) and temperature decay rate (`lambda`) using different
heuristics, since the heuristics used in the paper weren't published.
Simulated annealing is highly dependent on the `k` and
`lambda` parameters; so this will not have similar results to that of
the paper.

[1] MAGNA++: Maximizing Accuracy in Global Network Alignment via both
node and edge conservation, V. Vijayan, V. Saraph, and T. Milenkovic,
Bioinformatics, Volume 31, Issue 14, 15 July 2015, Pages 2409-2411.
[2] Alignment of dynamic networks, V. Vijayan, D. Critchlow, and
T. Milenkovic, Bioinformatics, Volume 33, Issue 14, 15 July 2017,
Pages i180-i189.
[3] SANA: simulated annealing far outperforms many other search
algorithms for biological network alignment" Nil Mamano Wayne B. Hayes
Bioinformatics, Volume 33, Issue 14, 15 July 2017, Pages 2156-2164

# Example usage

```julia
using SANA, LightGraphs

g1 = erdos_renyi(200,0.1)
g2 = g1
G1 = adjacency_matrix(g1)
G2 = adjacency_matrix(g2)

res = sanas3(G1,G2,maxtime=600,tempschedule=TempAdaptive())
f = res.f
mean(1:length(f) .== f[1:length(f)])

res = sanas3(G1,G2,maxtime=600,tempschedule=TempDecay())
f = res.f
mean(1:length(f) .== f[1:length(f)])
```

SANA can be installed as follows.

```julia
Pkg.clone("https://github.com/vvjn/SANA.jl")
```

# Documentation

Available [here](https://vvjn.github.io/SANA.jl/latest).
