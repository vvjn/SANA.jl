    export TempDecay, TempAdaptive, TempSchedule

abstract type TempSchedule end
immutable TempDecay <: TempSchedule end
immutable TempAdaptive <: TempSchedule end

# the lambda log(1+r_i) part is from the paper
# "Hybrid Simulated Annealing Algorithm Based on Adaptive
# Cooling Schedule for TSP" by Liu et al (2009)
# Just the adaptive cooling schedule part, not the hybrid or other parts
# the adaptive cooling shedule increases the temperature
# based on how long there has not been any improvement in the objective function

# k is initial temperature
# lambda is the cooling rate
# with time, the temperature decreases exponentially
# numconsnegdiffe is the # of consecutive iterations with no improvement in objective function
sanatemp(::TempDecay, i::Integer,k::Real,lambda::Real,numconsnegdiffe::Integer=0) =
    k * exp(-lambda*i)

sanatemp(::TempAdaptive, i::Integer,k::Real,lambda::Real,numconsnegdiffe::Integer=0) =
    k * exp(-lambda*i) + lambda * log(1 + numconsnegdiffe)

sanaprob(de::Real,i::Integer,k::Real,lambda::Real,numconsnegdiffe::Integer=0) =
    exp(de/sanatemp(i,k,lambda,numconsnegdiffe))
sanaprob(de::Real,temp::Real) = exp(de/temp)

function calculate_k_lambda_2(meas::NetalignMeasure,f::Vector{Int},
                              tmax::Number; # tmax in seconds
                              k::Number=NaN,lambda::Number=NaN,
                              tempschedule::TempSchedule=TempDecay(),
                              rng=Base.GLOBAL_RNG,
                              pstart=1e-3,pend=1e-40)
    println("********\nEstimate k and lambda parameters pstart $pstart pend $pend")
    updatek = isnan(k)
    updatelam = isnan(lambda)
    function calc_step(k::Number,lambda::Number,tmax_step::Number)
        if updatek || updatelam
            res = sana(meas,f,SanaParams(k,lambda,tmax_step,0.0,tempschedule,Inf,rng))
            #minscore,maxscore,titer = sana_stats(m)
            iter = Int(ceil(tmax/(res.time/res.iter)))
        end
        kfinal = 0.0
        if updatek
            k = -res.maxscorediff/log(pstart)
        end
        if updatelam
            kfinal = -res.minscorediff/log(pend)
            lambda = -log(kfinal/k)/iter
        end
        println("Using params k $k kfinal $kfinal lambda $lambda (iter $iter tmax $tmax mindiffe $(res.minscorediff) maxdiffe $(res.maxscorediff) ")
        k,lambda
    end
    k,lambda = calc_step(10.0,0.01,5) # before was 5,10,250
    k,lambda = calc_step(k,lambda,5)
    k,lambda = calc_step(k,lambda,10)
    println("End estimation k and lambda parameters with k = $k and lambda = $lambda\n*******")
    k,lambda
end

function calculate_k_lambda_1(meas::NetalignMeasure,f::Vector{Int},
                              tmax::Number; # tmax in seconds
                              k::Number=NaN,lambda::Number=NaN,
                              tempschedule::TempSchedule=TempDecay(),
                              rng=Base.GLOBAL_RNG,
                              pstart=1e-3,pend=1e-40)
    res = sana(meas,f,SanaParams(1.0,0.01,tmax,0.0,tempschedule,dim(meas)*Int(ceil(sqrt(dim(meas))))))
    iter = Int(ceil(tmax/(res.time/res.iter)))
    if isnan(lambda)
        lamvec = logspace(-20,-3,10)
    else
        lamvec = [lambda]
    end
    if isnan(k)
        kvec = flipdim(logspace(-5,1,10),1)
    else
        kvec = [k]
    end
    # find biggest k and smallest lamda s.t. pstart <= p <= pend or
    # if that doesn't work s.t. p <= pend
    kf = 0.0
    lamf = 0.0
    p = 1.0
    for i = 1:length(kvec), j = 1:length(lamvec)
        kf = kvec[i]
        lamf = lamvec[j]
        p = sanaprob(-res.minscorediff, iter, kf, lamf)
        if pstart <= p <= pend
            break
        end
    end
    for i = 1:length(kvec), j = 1:length(lamvec)
        kf = kvec[i]
        lamf = lamvec[j]
        p = sanaprob(-res.minscorediff, iter, kf, lamf)
        if p <= pend
            break
        end
    end
    println("Using params k $kf lambda $lamf for p $p")
    kf,lamf
end
