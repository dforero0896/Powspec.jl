# Boxes, periodic
function CATA(data::Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, data_w::AbstractArray{T},
    ) where T <: AbstractFloat

    cat = CATA(1)
    data_ptr, sumw2_data, wdata = build_catalog(data, data_w, true)
    unsafe_store!(cat.data, data_ptr, 1)
    unsafe_store!(cat.ndata,length(data_w),1)
    unsafe_store!(cat.wdata,wdata[1],1)
    cat
end #func

# Boxes, cross periodic
function CATA(data::Tuple{Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}}, 
            data_w::Tuple{AbstractArray{T}, AbstractArray{T}},
    ) where T <: AbstractFloat
    cat = CATA(length(data))
    for i in 1:2
        data_ptr, sumw2_data, wdata = build_catalog(data[i], data_w[i], true)
        unsafe_store!(cat.data, data_ptr, i)
        unsafe_store!(cat.ndata,length(data_w[i]),i)
        unsafe_store!(cat.wdata,wdata[1],i)
    end #for
    cat
end #func


function power_spectrum(data::Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, 
                        data_w::AbstractArray{T},
                        powspec_conf_file::AbstractString; output_file = nothing, precision = :double) where T <: AbstractFloat
    save_out = output_file != nothing
    
    test_output = save_out ? "--auto=$output_file" : "--auto=test.dat"
    test_output_ptr = Base.unsafe_convert(Cstring, test_output)


    conf_file = "--conf=$powspec_conf_file"
    conf_file_ptr = Base.unsafe_convert(Cstring, conf_file)

    argc::Cint = 3
    argv_vec = [Base.unsafe_convert(Cstring, "POWSPEC"), conf_file_ptr, test_output_ptr]
    argv = Base.unsafe_convert(Ptr{Cstring}, argv_vec)
    @show argv
    int_cache = Ptr{Cint}(Base.Libc.calloc(2, sizeof(Cint)))
    cat = CATA(data, data_w)
    if precision == :single
        pk = ccall((:compute_pk, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec_f.so"), Ptr{PK}, (Ref{CATA}, Cint, Cint, Ptr{Cint}, Cint, Ptr{Cstring}), cat, save_out, false, int_cache, argc, argv)
        
    elseif precision == :double
        pk = ccall((:compute_pk, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec.so"), Ptr{PK}, (Ref{CATA}, Cint, Cint, Ptr{Cint}, Cint, Ptr{Cstring}), cat, save_out, false, int_cache, argc, argv)
        
    end #if
    nbin = unsafe_load(int_cache, 1)
    nl = unsafe_load(int_cache, 2)

    k = zeros(T, nbin)
    kavg = zeros(T, nbin)
    kmin = zeros(T, nbin)
    kmax = zeros(T, nbin)
    nmodes = zeros(T, nbin)
    multipoles = zeros(T, nl, nbin)

    #copy_pk_results(PK *pk, double *k, double, *kavg, double *kedge, int *nmodes, double *auto_multipole1, double *auto_multipole2, double *cross_multipole)
    ccall((:copy_pk_results, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec.so"), 
          Ptr{Cvoid}, 
          (Ptr{PK}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          pk, vec(k), vec(kavg), vec(kmin), vec(kmax), vec(nmodes), vec(multipoles), C_NULL, C_NULL)
     
    out = Dict()
    out[:k] = k
    out[:kavg] = kavg
    out[:kmin] = kmin
    out[:kmax] = kmax
    out[:nmodes] = nmodes
    out[:multipoles] = multipoles
    out[:n_data_objects] = unsafe_load(cat.ndata, 1)
    out[:w_data_objects] = unsafe_load(cat.wdata, 1)
    out[:shot_noise] = unsafe_load(cat.shot, 1)
    out[:normalisation] = unsafe_load(cat.norm, 1)

    ccall((:powspec_destroy, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec.so"), 
            Ptr{Cvoid},
            (Ptr{PK},),
            pk)

    out
end #func


function power_spectrum(data::Tuple{Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}}, 
                        data_w::Tuple{AbstractArray{T}, AbstractArray{T}},
                        powspec_conf_file::AbstractString; output_auto = nothing, output_cross = nothing, precision = :double) where T <: AbstractFloat
    save_auto = output_auto != nothing
    save_cross = output_auto != nothing
    
    auto_output = save_auto ? "--auto=[$(output_auto[1]), $(output_auto[2])]" : "--auto=[test1.dat, test1.dat]"
    @show auto_output
    #auto_output_ptr = Base.unsafe_convert(Cstring, auto_output)
    auto_output_ptr = Cstring(pointer(deepcopy(auto_output)))
    cross_output = save_cross ? "--cross=$output_cross" : "--cross=cross.dat"
    cross_output_ptr = Cstring(pointer(deepcopy(cross_output)))


    conf_file = "--conf=$powspec_conf_file"
    conf_file_ptr = Cstring(pointer(deepcopy(conf_file)))

    argc::Cint = 4
    argv_vec = [Base.unsafe_convert(Cstring, "POWSPEC"), conf_file_ptr, auto_output_ptr, cross_output_ptr]
    argv = Base.unsafe_convert(Ptr{Cstring}, argv_vec)
    
    int_cache = Ptr{Cint}(Base.Libc.calloc(2, sizeof(Cint)))
    cat = CATA(data, data_w)
    
    if precision == :single
        pk = ccall((:compute_pk, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec_f.so"), Ptr{PK}, (Ref{CATA}, Cint, Cint, Ptr{Cint}, Cint, Ptr{Cstring}), cat, save_auto & save_cross, false, int_cache, argc, argv)
        
    elseif precision == :double
        pk = ccall((:compute_pk, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec.so"), Ptr{PK}, (Ref{CATA}, Cint, Cint, Ptr{Cint}, Cint, Ptr{Cstring}), cat, save_cross & save_cross, false, int_cache, argc, argv)
        
    end #if
    if pk == C_NULL
        error("C library returned NULL.")
    end #if
    nbin = unsafe_load(int_cache, 1)
    nl = unsafe_load(int_cache, 2)

    k = zeros(T, nbin)
    kavg = zeros(T, nbin)
    kmin = zeros(T, nbin)
    kmax = zeros(T, nbin)
    nmodes = zeros(T, nbin)
    auto_multipoles1 = zeros(T, nl, nbin)
    auto_multipoles2 = zeros(T, nl, nbin)
    cross_multipoles = zeros(T, nl, nbin)
    

    #copy_pk_results(PK *pk, double *k, double, *kavg, double *kedge, int *nmodes, double *auto_multipole1, double *auto_multipole2, double *cross_multipole)
    ccall((:copy_pk_results, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec.so"), 
          Ptr{Cvoid}, 
          (Ptr{PK}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          pk, vec(k), vec(kavg), vec(kmin), vec(kmax), vec(nmodes), vec(auto_multipoles1), vec(auto_multipoles2), vec(cross_multipoles))
     
    out = Dict()
    out[:k] = k
    out[:kavg] = kavg
    out[:kmin] = kmin
    out[:kmax] = kmax
    out[:nmodes] = nmodes
    out[:auto_multipoles] = (auto_multipoles1, auto_multipoles2)
    out[:cross_multipoles] = cross_multipoles
    out[:n_data_objects] = unsafe_wrap(Vector{Csize_t}, cat.ndata, 2; own = true)
    out[:w_data_objects] = unsafe_wrap(Vector{Cdouble}, cat.wdata, 2; own = true)
    out[:shot_noise] = unsafe_wrap(Vector{Cdouble}, cat.shot, 2, own = true)
    out[:normalisation] = unsafe_wrap(Vector{Cdouble}, cat.norm, 2, own = true)


    ccall((:powspec_destroy, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec.so"), 
            Ptr{Cvoid},
            (Ptr{PK},),
            pk)


    out
end #func
