"""
This uses alignment measures from the NetalignMeasures package

 sagenaux generates an auxiliary structure
 for a given NetalignMeasure, which is used as
 temporary scratch space when running the
 measureswap function

  measureswap calculates alignment score of nf,
 where nf is f with two nodes in the alignment swapped,
 i.e. nf[l],nf[k] = f[k],f[l], nf[i] = f[i] otherwise
"""

immutable PermState
    f :: Vector{Int} # current alignment
    inv :: Vector{Int} # inverse map of f
    nf :: Vector{Int} # new alignment, f plus a transposition
    ninv :: Vector{Int} # inverse map of nf
end

sagenaux(::NetalignMeasure) = nothing

sagenaux(meas::ConvexCombMeasure) = (sagenaux(meas.S),sagenaux(meas.T))

function measureswap(meas::ConvexCombMeasure,
                     es::ConvexCombScore,
                     st::PermState,
                     k::Int,l::Int,aux) # f[k],f[l] = f[l],f[k]
    a = measureswap(meas.S, es.s, st, k, l, aux[1])
    b = measureswap(meas.T, es.t, st, k, l, aux[2])
    s = meas.alpha * a.score + (1-meas.alpha) * b.score
    ConvexCombScore(a,b,s)
end

function measureswap(meas::NodeSimMeasure,
                     es::NodeSimScore,
                     st::PermState,
                     k::Int,l::Int, aux) # f[k],f[l] = f[l],f[k]
    S = meas.S
    m = dim(meas,1)
    simorig = (k>m ? 0.0 : S[k,st.f[k]]) + (l>m ? 0.0 : S[l,st.f[l]])
    simnew = (k>m ? 0.0 : S[k,st.f[l]]) + (l>m ? 0.0 : S[l,st.f[k]])
    score = es.score + (simnew - simorig)/size(S,1)
    NodeSimScore(score)
end

sagenaux(meas::S3Measure) =
    Vector{Int}(maximum(sum(meas.G1,1)) + maximum(sum(meas.G2,1)))

function measureswap(meas::S3Measure,
                     es::S3Score,
                     st::PermState,
                     k::Int,l::Int, aux) # f[k],f[l] = f[l],f[k]
    G1 = meas.G1
    G2 = meas.G2
    Nc,Nn,newNc,newNn = 0,0,0,0
    if k <= size(G1,1)
        nc,nn = nodecounthoodce(G1,G2,st.f,st.inv,k, aux)
        Nc += nc
        Nn += nn
        nc,nn = nodecounthoodce(G1,G2,st.nf,st.ninv,k, aux)
        newNc += nc
        newNn += nn
    end
    if l <= size(G1,1)
        nc,nn = nodecounthoodce(G1,G2,st.f,st.inv,l, aux)
        Nc += nc
        Nn += nn
        nc,nn = nodecounthoodce(G1,G2,st.nf,st.ninv,l, aux)
        newNc += nc
        newNn += nn
    end
    Nc = es.Nc - Nc + newNc
    Nn = es.Nn - Nn + newNn
    score = Nc/(Nc + Nn)
    S3Score(Nc,Nn,score)
end

function nodecounthoodce(G1::SparseMatrixCSC,
                         G2::SparseMatrixCSC,
                         f::Vector{Int},
                         finv::Vector{Int},
                         k::Int, aux)
    m = size(G1,1)
    adj1ix = nzrange(G1,k)
    nb1 = length(adj1ix)
    for i = 1:nb1
        aux[i] = G1.rowval[adj1ix[i]]
    end
    adj2ix = nzrange(G2,f[k])
    nb2 = length(adj2ix)
    nbc = nb1
    for i = 1:nb2
        x = finv[G2.rowval[adj2ix[i]]]
        if x <= m
            nbc += 1
            aux[nbc] = x
        end
    end
    #println("$(nbc) : $(aux[1:nbc])")
    sort!(@view aux[1:nbc])
    nc = 0
    nn = 0
    i = 1
    while i < nbc
        if aux[i+1] == aux[i]
            nc += 1
            i += 2
        else
            nn += 1
            i += 1
        end
    end
    if i == nbc
        nn += 1
    end
    nc,nn
end

function s3test(G1::SparseMatrixCSC,
                G2::SparseMatrixCSC,
                f::Vector{Int})
    aux = Vector{Int}(maximum(sum(G1,1)) + maximum(sum(G2,1)))
    finv = invperm(f)
    Nc = 0; Nn = 0
    for k = 1:size(G1,1)
        nc,nn = nodecounthoodce(G1,G2,f,finv,k,aux)
        #println((k,nc,nn))
        Nc += nc; Nn += nn
    end
    S3Score(Nc,Nn, Nc/(Nc+Nn))
end

# Counts conserved and non-conserved edges
# adjacent to node k
function nodecounthoodceslow(G1::SparseMatrixCSC,
                             G2::SparseMatrixCSC,
                             f::Vector{Int},
                             finv::Vector{Int},
                             k::Int, aux)
    m = size(G1,1)
    a1 = adjnodes(G1,k)
    a2 = finv[adjnodes(G2,f[k])]
    a2 = a2[a2 .<= m]
    Ak = sparsevec(vcat(a1,a2), 1, m)
    nc = 0
    nn = 0
    #println("$(length(Ak)) : $(vcat(a1,a2))")
    for x in nonzeros(Ak)
        nc += Int(x==2)
        nn += Int(x==1)
    end
    nc,nn
end

function s3testslow(G1::SparseMatrixCSC,
                    G2::SparseMatrixCSC,
                    f::Vector{Int})
    aux = nothing
    finv = invperm(f)
    Nc = 0; Nn = 0
    for k = 1:size(G1,1)
        nc,nn = nodecounthoodceslow(G1,G2,f,finv,k,nothing)
        #println((nc,nn))
        Nc += nc; Nn += nn
    end
    S3Score(Nc,Nn, Nc/(Nc+Nn))
end

function measureswap(meas::DS3Measure,
                     es::DS3Score,
                     st::PermState,
                     k::Int,l::Int, aux) # f[k],f[l] = f[l],f[k]
    G1 = meas.G1
    G2 = meas.G2
    Tc,Tn,newTc,newTn = 0.0,0.0,0.0,0.0
    if k <= size(G1,1)
        tc,tn = nodecounthoodcet(G1,G2,st.f,st.inv,k, aux)
        Tc += tc
        Tn += tn
        tc,tn = nodecounthoodcet(G1,G2,st.nf,st.ninv,k, aux)
        newTc += tc
        newTn += tn
    end
    if l <= size(G1,1)
        tc,tn = nodecounthoodcet(G1,G2,st.f,st.inv,l, aux)
        Tc += tc
        Tn += tn
        tc,tn = nodecounthoodcet(G1,G2,st.nf,st.ninv,l, aux)
        newTc += tc
        newTn += tn
    end
    Tc = es.Tc - Tc + newTc
    Tn = es.Tn - Tn + newTn
    score = Tc/(Tc + Tn)
    DS3Score(Tc,Tn,score)
end

immutable DS3Aux
    ix :: Vector{Int}
    rowval :: Vector{Int}
    nzval :: Vector{Events}
end

function sagenaux(meas::DS3Measure)
    F1 = flatten(meas.G1)
    F2 = flatten(meas.G2)
    n = maximum(sum(F1,1)) + maximum(sum(F2,1))
    DS3Aux(Vector{Int}(n), Vector{Int}(n), Vector{Events}(n))
end

function nodecounthoodcet(G1::SparseMatrixCSC,
                          G2::SparseMatrixCSC,
                          f::Vector{Int},
                          finv::Vector{Int},
                          k::Int, aux)
    m = size(G1,1)
    adj1ix = nzrange(G1,k)
    nb1 = length(adj1ix)
    for i = 1:nb1
        aix = adj1ix[i]
        aux.rowval[i] = G1.rowval[aix]
        aux.nzval[i] = G1.nzval[aix]
    end
    adj2ix = nzrange(G2,f[k])
    nb2 = length(adj2ix)
    nbc = nb1
    for i = 1:nb2
        aix = adj2ix[i]
        x = finv[G2.rowval[aix]]
        if x <= m
            nbc += 1
            aux.rowval[nbc] = x
            aux.nzval[nbc] = G2.nzval[aix]
        end
    end
    #println("$(nbc) : $(aux[1:nbc])")
    sortperm!(view(aux.ix,1:nbc), view(aux.rowval,1:nbc))
    tc = 0.0
    tn = 0.0
    i = 1
    while i < nbc
        ixc = aux.ix[i]
        ixn = aux.ix[i+1]
        aric = aux.rowval[ixc]
        arin = aux.rowval[ixn]
        if aric == arin
            tcp,tnp = cet_ncet(aux.nzval[ixc],aux.nzval[ixn])
            tc += tcp
            tn += tnp
            i += 2
        else
            tcp,tnp = cet_ncet(aux.nzval[ixc])
            tn += tnp
            i += 1
        end
    end
    if i == nbc
        ixc = aux.ix[i]
        tcp,tnp = cet_ncet(aux.nzval[ixc])
        tn += tnp
    end
    tc,tn
end

function ds3test(G1::SparseMatrixCSC,
                 G2::SparseMatrixCSC,
                 f::Vector{Int})
    aux = sagenaux(DS3Measure(G1,G2))
    finv = invperm(f)
    Nc = 0; Nn = 0
    for k = 1:size(G1,1)
        nc,nn = nodecounthoodcet(G1,G2,f,finv,k,aux)
        #println((k,nc,nn))
        Nc += nc; Nn += nn
    end
    Nc,Nn, Nc/(Nc+Nn)
end

# Counts conserved and non-conserved event time
# adjacent to node k
function nodecounthoodcetslow(G1::SparseMatrixCSC,
                              G2::SparseMatrixCSC,
                              f::Vector{Int},
                              finv::Vector{Int},
                              k::Int, aux)
    m = size(G1,1)
    a1 = adjnodes(G1,k)
    a2 = finv[adjnodes(G2,f[k])]
    a2ix = a2 .<= m
    a2 = a2[a2ix]
    v1 = adjnzval(G1,k)
    v2 = adjnzval(G2,f[k])[a2ix]
    Ak = sparsevec(vcat(a1,a2), vcat(v1,v2), m, mergesorted)
    tc = 0.0
    tn = 0.0
    for x in nonzeros(Ak)
        tcp,tnp = cet_ncet(x)
        tc += tcp
        tn += tnp
    end
    tc,tn
end

function ds3testslow(G1::SparseMatrixCSC,
                     G2::SparseMatrixCSC,
                     f::Vector{Int})
    aux = nothing
    finv = invperm(f)
    Nc = 0; Nn = 0
    for k = 1:size(G1,1)
        nc,nn = nodecounthoodcetslow(G1,G2,f,finv,k,aux)
        #println((k,nc,nn))
        Nc += nc; Nn += nn
    end
    Nc,Nn, Nc/(Nc+Nn)
end

sagenaux(meas::WECMeasure) =
    Vector{Int}(maximum(sum(meas.G1,1)) + maximum(sum(meas.G2,1)))

function measureswap(meas::WECMeasure,
                     es::WECScore,
                     st::PermState,
                     k::Int,l::Int, aux) # f[k],f[l] = f[l],f[k]
    G1 = meas.G1
    G2 = meas.G2
    S = meas.S
    Nc, newNc = 0.0, 0.0
    if k <= size(G1,1)
        nc = nodecounthoodwec(G1,G2,S,st.f,st.inv,k, aux)
        Nc += nc
        nc = nodecounthoodwec(G1,G2,S,st.nf,st.ninv,k, aux)
        newNc += nc
    end
    if l <= size(G1,1)
        nc = nodecounthoodwec(G1,G2,S,st.f,st.inv,l, aux)
        Nc += nc
        nc = nodecounthoodwec(G1,G2,S,st.nf,st.ninv,l, aux)
        newNc += nc
    end
    WECScore(es.score-Nc+newNc)
end

function nodecounthoodwec(G1::SparseMatrixCSC,
                          G2::SparseMatrixCSC,
                          S::AbstractMatrix,
                          f::Vector{Int},
                          finv::Vector{Int},
                          k::Int, aux)
    m = size(G1,1)
    adj1ix = nzrange(G1,k)
    nb1 = length(adj1ix)
    for i = 1:nb1
        aux[i] = G1.rowval[adj1ix[i]]
    end
    adj2ix = nzrange(G2,f[k])
    nb2 = length(adj2ix)
    nbc = nb1
    for i = 1:nb2
        x = finv[G2.rowval[adj2ix[i]]]
        if x <= m
            nbc += 1
            aux[nbc] = x
        end
    end
    sort!(@view aux[1:nbc])
    sw = 0.0
    i = 1
    while i < nbc
        ai = aux[i]
        if aux[i+1] == ai
            sw += S[ai,f[ai]]
            i += 2
        else
            i += 1
        end
    end
    sw / min(nnz(G1),nnz(G2))
end

function wectest(G1::SparseMatrixCSC,
                 G2::SparseMatrixCSC,
                 S::AbstractMatrix,
                 f::Vector{Int})
    aux = Vector{Int}(maximum(sum(G1,1)) + maximum(sum(G2,1)))
    finv = invperm(f)
    Sc = 0.0
    for k = 1:size(G1,1)
        sc = nodecounthoodwec(G1,G2,S,f,finv,k,aux)
        Sc += sc
    end
    Sc
end

# Counts weighted edge conservation
# adjacent to node k
function nodecounthoodwecslow(G1::SparseMatrixCSC,
                              G2::SparseMatrixCSC,
                              S::AbstractMatrix,
                              f::Vector{Int},
                              finv::Vector{Int},
                              k::Int, aux)
    m = size(G1,1)
    a1 = adjnodes(G1,k)
    a2 = finv[adjnodes(G2,f[k])]
    a2 = a2[a2 .<= m]
    Ak = sparsevec(vcat(a1,a2), 1, m)
    #sk = S[k,f[k]]
    sw = 0.0
    for i in 1:length(Ak.nzind)
        aix = Ak.nzind[i]
        sw += ifelse(Ak.nzval[i] == 2, S[aix,f[aix]], 0.0)
        #sw += ifelse(Ak.nzval[i] == 2, 0.5 * (sk + S[aix,f[aix]]), 0.0)
    end
    sw / min(nnz(G1),nnz(G2))
end

function wectestslow(G1::SparseMatrixCSC,
                     G2::SparseMatrixCSC,
                     S::AbstractMatrix,
                     f::Vector{Int})
    aux = Vector{Int}(maximum(sum(G1,1)) + maximum(sum(G2,1)))
    finv = invperm(f)
    Sc = 0.0
    for k = 1:size(G1,1)
        Sc += nodecounthoodwecslow(G1,G2,S,f,finv,k,aux)
    end
    Sc
end

immutable DWECAux
    ix :: Vector{Int}
    rowval :: Vector{Int}
    nzval :: Vector{Events}
end

function measureswap(meas::DWECMeasure,
                     es::DWECScore,
                     st::PermState,
                     k::Int,l::Int, aux) # f[k],f[l] = f[l],f[k]
    G1 = meas.G1
    G2 = meas.G2
    S = meas.S
    Nc, newNc = 0.0, 0.0
    if k <= size(G1,1)
        nc = nodecounthoodwec(G1,G2,S,st.f,st.inv,k,
                              meas.activitysum1,meas.activitysum2,aux)
        Nc += nc
        nc = nodecounthoodwec(G1,G2,S,st.nf,st.ninv,k,
                              meas.activitysum1,meas.activitysum2,aux)
        newNc += nc
    end
    if l <= size(G1,1)
        nc = nodecounthoodwec(G1,G2,S,st.f,st.inv,l,
                              meas.activitysum1,meas.activitysum2,aux)
        Nc += nc
        nc = nodecounthoodwec(G1,G2,S,st.nf,st.ninv,l,
                              meas.activitysum1,meas.activitysum2,aux)
        newNc += nc
    end
    DWECScore(es.score-Nc+newNc)
end

function sagenaux(meas::DWECMeasure)
    F1 = flatten(meas.G1)
    F2 = flatten(meas.G2)
    n = maximum(sum(F1,1)) + maximum(sum(F2,1))
    DWECAux(Vector{Int}(n), Vector{Int}(n), Vector{Events}(n))
end

function nodecounthooddwec(G1::SparseMatrixCSC,
                           G2::SparseMatrixCSC,
                           S::AbstractMatrix,
                           f::Vector{Int},
                           finv::Vector{Int},
                           k::Int, activitysum1, activitysum2, aux)
    m = size(G1,1)
    adj1ix = nzrange(G1,k)
    nb1 = length(adj1ix)
    for i = 1:nb1
        aix = adj1ix[i]
        aux.rowval[i] = G1.rowval[aix]
        aux.nzval[i] = G1.nzval[aix]
    end
    adj2ix = nzrange(G2,f[k])
    nb2 = length(adj2ix)
    nbc = nb1
    for i = 1:nb2
        aix = adj2ix[i]
        x = finv[G2.rowval[aix]]
        if x <= m
            nbc += 1
            aux.rowval[nbc] = x
            aux.nzval[nbc] = G2.nzval[aix]
        end
    end
    #println("$(nbc) : $(aux[1:nbc])")
    sortperm!(view(aux.ix,1:nbc), view(aux.rowval,1:nbc))
    sw = 0.0
    i = 1
    while i < nbc
        ixc = aux.ix[i]
        ixn = aux.ix[i+1]
        aric = aux.rowval[ixc]
        arin = aux.rowval[ixn]
        if aric == arin
            tcp,tnp = cet_ncet(aux.nzval[ixc],aux.nzval[ixn])
            sw += tcp * S[aric,f[aric]]
            i += 2
        else
            i += 1
        end
    end
    sw / min(activitysum1,activitysum2)
end

function dwectest(G1::SparseMatrixCSC,
                  G2::SparseMatrixCSC,
                  S::AbstractMatrix,
                  f::Vector{Int})
    meas = DWECMeasure(G1,G2,S)
    aux = sagenaux(meas)
    finv = invperm(f)
    Sc = 0.0
    for k = 1:size(G1,1)
        Sc += nodecounthooddwec(G1,G2,S,f,finv,k,
                                meas.activitysum1,meas.activitysum2,aux)
    end
    Sc
end

# Counts conserved and non-conserved event time
# adjacent to node k
function nodecounthooddwecslow(G1::SparseMatrixCSC,
                               G2::SparseMatrixCSC,
                               S::AbstractMatrix,
                               f::Vector{Int},
                               finv::Vector{Int},
                               k::Int, activitysum1, activitysum2, aux)
    m = size(G1,1)
    a1 = adjnodes(G1,k)
    a2 = finv[adjnodes(G2,f[k])]
    a2ix = a2 .<= m
    a2 = a2[a2ix]
    v1 = adjnzval(G1,k)
    v2 = adjnzval(G2,f[k])[a2ix]
    Ak = sparsevec(vcat(a1,a2), vcat(v1,v2), m, mergesorted)
    sw = 0.0
    for i = 1:length(Ak.nzind)
        aix = Ak.nzind[i]
        tcp,tnp = cet_ncet(Ak.nzval[i])
        sw += tcp * S[aix,f[aix]]
    end
    sw / min(activitysum1,activitysum2)
end

function dwectestslow(G1::SparseMatrixCSC,
                      G2::SparseMatrixCSC,
                      S::AbstractMatrix,
                      f::Vector{Int})
    meas = DWECMeasure(G1,G2,S)
    aux = nothing
    finv = invperm(f)
    Sc = 0.0
    for k = 1:size(G1,1)
        Sc += nodecounthooddwecslow(G1,G2,S,f,finv,k,
                                    meas.activitysum1,meas.activitysum2,aux)
    end
    Sc
end

