export SanaResult, SanaParams, sanads3, sanads3nodesim,
    sanas3, sanas3nodesim,sananodesim, calculate_k_lambda, sana, sanawec

immutable SanaParams
    k :: Float64
    lambda :: Float64
    maxtime :: Float64 # stop if over this many seconds
    tempthresh :: Float64 # " below this temperature
    tempschedule :: TempSchedule
    maxiter :: Float64 # " " many iterations
    rng :: MersenneTwister
    stepiter :: Int
    function SanaParams(k::Number,lambda::Number,maxtime::Number,tempthresh::Number,tempschedule::TempSchedule,
                        maxiter::Number,rng=Base.GLOBAL_RNG,stepiter::Integer=100_000)
        new(k,lambda,maxtime,tempthresh,tempschedule,maxiter,rng,stepiter)
    end
end

immutable SanaResult{M}
    f::Vector{Int} # alignment
    score::Float64 # alignment score (bigger is better)
    es::M # alignment score details
    minscorediff::Float64 # min/max difference of scores between iterations
    maxscorediff::Float64
    time::Float64 # total time to run
    temp::Float64 # final temperature
    iter::Int # total iterations
    numswaps::Int # number of iterations where f was modified
    numupswaps::Int # number of iterations where score increased
    params::SanaParams
end

"""
Main simulated annealing algorithm
Input:
meas -- alignment measure to optimize
f -- initial alignment
params -- includes initial temperature, decay rate, max. iterations, max. time etc.
verbose -- if true, prints every params.stepiter iterations
"""
function sana(meas::NetalignMeasure,f::Vector{Int},params::SanaParams,verbose=true)
    st = PermState(copy(f),invperm(f),copy(f),invperm(f))
    m = dim(meas,1)
    n = length(f)
    n != dim(meas) && error("Bad args")

    es = measure(meas,st.f)
    verbose && println("Initial score: $(score(es))")
    minscorediff = Inf
    maxscorediff = -Inf
    tstart = usercputime_us() # measures time in cpu time
    tcur = tstart
    tdiff = 0.0 # current time - start time, in seconds
    tempcur = params.k
    numswaps = 0
    numupswaps = 0
    numconsnegdiffe = 0 # number of consecutive iterations w/ no improvement

    aux = sagenaux(meas) # auxiliary space for calculating delta scores

    tempcur = 0.0
    probcur = 0.0

    iter = 0
    while true
        # During each iteration, we apply a transposition to f, creating nf.
        # That is, we swap two nodes in f, creating nf.
        # If nf is better than f, then swap. Otherwise swap with some probability.

        # We randomly select a transposition by selecting k,l randomly.
        # If k>m && l>m then it doesn't change the alignment since f[m+1:n] are dummy nodes.
        k = 0; l = 0
        while k==l || (k>m && l>m)
            k = rand(params.rng, 1:n); l = rand(params.rng, 1:n)
        end
        # Apply transposition.
        st.nf[k] = st.f[l]; st.nf[l] = st.f[k]
        st.ninv[st.nf[k]] = k; st.ninv[st.nf[l]] = l
        # Calculate alignment score of nf.
        # We use measureswap to calculate score of nf efficiently:
        # measureswap returns the score of nf given the score of f,
        # the swapped nodes (k and l), and f, inv(f), nf, inv(nf).
        # Ideally, measureswap should have O(d), d = avg. degree of node, or
        # better complexity in order to be efficient.
        nes = measureswap(meas,es,st,k,l, aux)
        diffe = score(nes) - score(es)

        if diffe < 0.0 numconsnegdiffe += 1 end
        if diffe > 0.0 numconsnegdiffe = 0 end
        if diffe >= 0.0
            doswap = true
        elseif (tempcur = sanatemp(params.tempschedule, iter,params.k,params.lambda,numconsnegdiffe);
                probcur = exp(diffe/tempcur);
                rand(params.rng) < probcur)
            doswap = true
            numupswaps += 1
        else
            doswap = false
        end

        if doswap # Change f to nf
            numswaps += 1
            es = nes
            st.f[k] = st.nf[k]; st.f[l] = st.nf[l];
            st.inv[st.f[k]] = k; st.inv[st.f[l]] = l
        else # Revert nf to f
            st.nf[k] = st.f[k]; st.nf[l] = st.f[l]
            st.ninv[st.nf[k]] = k; st.ninv[st.nf[l]] = l
        end
        absdiffe = abs(diffe)
        if absdiffe > 0.0 minscorediff = min(absdiffe, minscorediff) end
        maxscorediff = max(absdiffe, maxscorediff)
        iter += 1

        # Logging
        if mod(iter,params.stepiter)==0
            tcur = usercputime_us()
            tdiff = ((tcur-tstart)/1e6)
            if verbose
                @printf "\rTime %.1fs temp %.5e prob %.5e score %.5f iter %d nswaps %d nupswaps %d  numconsnegdiffe %d  " tdiff tempcur probcur score(es) iter numswaps numupswaps numconsnegdiffe
                flush(STDOUT)
            end
            if (iter >= params.maxiter) || (tempcur <= params.tempthresh) ||
                (tdiff >= params.maxtime)
                break
            end
        end
    end
    if verbose println("\nFinal score: $(es.score)") end
    SanaResult(st.f,score(es),es,minscorediff,maxscorediff,tdiff,tempcur,iter,numswaps,numupswaps,params)
end

function sana(meas::NetalignMeasure;maxtime::Number=30.0,
              k::Number=NaN,lambda::Number=NaN,
              maxiter::Number=Inf,
              tempthresh::Number=0.0,
              tempschedule::TempSchedule=TempDecay(),
              rng=Base.GLOBAL_RNG,
              stepiter::Number=100_000,
              pstart=1e-3,pend=1e-40,f=randperm(dim(meas)))
    if isnan(k) || isnan(lambda)
        k,lambda = calculate_k_lambda_2(meas,f,maxtime,k=k,lambda=lambda,tempschedule=tempschedule,
                                        rng=rng,pstart=pstart,pend=pend)
    end
    sana(meas,f,SanaParams(k,lambda,maxtime,tempthresh,tempschedule,maxiter,rng,stepiter))
end

"Optimize the assignment problem (node similarity)"
sananodesim(b::AbstractMatrix; args...) =
    sana(NodeSimMeasure(b); stepiter=100_000, args...)
"Optimize DS3"
sanads3(G1::SparseMatrixCSC,G2::SparseMatrixCSC; args...) =
    sana(DS3Measure(G1,G2); args...)
"Optimize convex combination of DS3 and node similarity"
sanads3nodesim(G1::SparseMatrixCSC,G2::SparseMatrixCSC,S::AbstractMatrix,alpha::Float64; args...) =
    sana(ConvexCombMeasure(DS3Measure(G1,G2),NodeSimMeasure(S),alpha); args...)
"Optimize S3"
sanas3(G1::SparseMatrixCSC,G2::SparseMatrixCSC; args...) =
    sana(S3Measure(G1,G2); args...)
"Optimize convex combination of S3 and node similarity"
sanas3nodesim(G1::SparseMatrixCSC,G2::SparseMatrixCSC,S::AbstractMatrix,alpha::Float64; args...) =
    sana(ConvexCombMeasure(S3Measure(G1,G2),NodeSimMeasure(S),alpha); args...)
"Optimize WEC"
sanawec(G1::SparseMatrixCSC,G2::SparseMatrixCSC, S::AbstractMatrix; args...) =
    sana(WECMeasure(G1,G2,S); args...)
