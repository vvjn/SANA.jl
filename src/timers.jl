
function usercputime_us()
    rusage = Libc.malloc(4*sizeof(Clong) + 14*sizeof(UInt64))  # sizeof(uv_rusage_t); this is different from sizeof(rusage)
    ccall(:uv_getrusage, Cint, (Ptr{Void},), rusage)
    utime = UInt64(1000000)*unsafe_load(convert(Ptr{Clong}, rusage + 0*sizeof(Clong))) +    # user CPU time
    unsafe_load(convert(Ptr{Clong}, rusage + 1*sizeof(Clong)))
    #stime = UInt64(1000000)*unsafe_load(convert(Ptr{Clong}, rusage + 2*sizeof(Clong))) +    # system CPU time
    #unsafe_load(convert(Ptr{Clong}, rusage + 3*sizeof(Clong)))
    #ttime = utime + stime  # total CPU time
    Libc.free(rusage)
    return utime
end
function usercputime_s()
    rusage = Libc.malloc(4*sizeof(Clong) + 14*sizeof(UInt64))  # sizeof(uv_rusage_t); this is different from sizeof(rusage)
    ccall(:uv_getrusage, Cint, (Ptr{Void},), rusage)
    utime = UInt64(unsafe_load(convert(Ptr{Clong}, rusage + 0*sizeof(Clong))))     # user CPU time
    Libc.free(rusage)
    return utime
end
