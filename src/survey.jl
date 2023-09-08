# Survey, randoms
function CATA(data::Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, data_w::AbstractArray{T}, data_fkp::AbstractArray{T}, data_nz::AbstractArray{T},
    rand::Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, rand_w::AbstractArray{T}, rand_fkp::AbstractArray{T}, rand_nz::AbstractArray{T},
    ) where T <: AbstractFloat
    cat = CATA(1)
    data_ptr, sumw2_data, wdata = build_catalog(data, data_w, false; data_fkp = data_fkp, data_nz = data_nz)
    rand_ptr, sumw2_rand, wrand = build_catalog(rand, rand_w, false; data_fkp = rand_fkp, data_nz = rand_nz)
    unsafe_store!(cat.data, data_ptr, 1)
    unsafe_store!(cat.ndata,length(data_w),1)
    unsafe_store!(cat.wdata,wdata[1],1)

    unsafe_store!(cat.rand, rand_ptr, 1)
    unsafe_store!(cat.nrand,length(rand_w),1)
    unsafe_store!(cat.wrand,wrand[1],1)

    unsafe_store!(cat.alpha, wdata[1] / wrand[1], 1)
    unsafe_store!(cat.shot, sumw2_data[1] + unsafe_load(cat.alpha, 1)^2  * sumw2_rand[1], 1)
    if sumw2_data[1] == 0
        unsafe_store!(cat.norm, unsafe_load(cat.alpha, 1) * sumw2_rand[2], 1)
    elseif sumw2_rand[1] == 0
        unsafe_store!(cat.norm, sumw2_data[2], 1)
    else
        unsafe_store!(cat.norm, unsafe_load(cat.alpha, 1) * sumw2_rand[2], 1)
    end #if
    
    cat
end #func


# Survey cross, randoms
function CATA(data::Tuple{Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}}, 
            data_w::Tuple{AbstractArray{T}, AbstractArray{T}},
            data_fkp::Tuple{AbstractArray{T}, AbstractArray{T}},
            data_nz::Tuple{AbstractArray{T}, AbstractArray{T}},
            rand::Tuple{Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}}, 
            rand_w::Tuple{AbstractArray{T}, AbstractArray{T}},
            rand_fkp::Tuple{AbstractArray{T}, AbstractArray{T}},
            rand_nz::Tuple{AbstractArray{T}, AbstractArray{T}},
    ) where T <: AbstractFloat
    cat = CATA(length(data))
    for i in 1:2
        data_ptr, sumw2_data, wdata = build_catalog(data[i], data_w[i], false, data_fkp = data_fkp[i], data_nz = data_nz[i])
        rand_ptr, sumw2_rand, wrand = build_catalog(rand[i], rand_w[i], false, data_fkp = rand_fkp[i], data_nz = rand_nz[i])
        unsafe_store!(cat.data, data_ptr, i)
        unsafe_store!(cat.ndata,length(data_w[i]),i)
        unsafe_store!(cat.wdata,wdata[1],i)


        unsafe_store!(cat.rand, rand_ptr, i)
        unsafe_store!(cat.nrand,length(rand_w[i]),i)
        unsafe_store!(cat.wrand,wrand[1],i)

        unsafe_store!(cat.alpha, wdata[1] / wrand[1], i)
        unsafe_store!(cat.shot, sumw2_data[1] + unsafe_load(cat.alpha, i)^2  * sumw2_rand[1], i)
        if sumw2_data[1] == 0
            unsafe_store!(cat.norm, unsafe_load(cat.alpha, i) * sumw2_rand[2], i)
        elseif sumw2_rand[1] == 0
            unsafe_store!(cat.norm, sumw2_data[2], i)
        else
            unsafe_store!(cat.norm, unsafe_load(cat.alpha, i) * sumw2_rand[2], i)
        end #if
    end #for
    cat
end #func


function power_spectrum(data::Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, 
                        data_w::AbstractArray{T}, data_fkp::AbstractArray{T}, data_nz::AbstractArray{T},
                        rand::Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, 
                        rand_w::AbstractArray{T}, rand_fkp::AbstractArray{T}, rand_nz::AbstractArray{T},
                        powspec_conf_file::AbstractString; output_file = nothing, precision = :double) where T <: AbstractFloat
    save_out = output_file != nothing
    
    test_output = save_out ? "--auto=$output_file" : "--auto=test.dat"
    test_output_ptr = Base.unsafe_convert(Cstring, deepcopy(test_output))


    conf_file = "--conf=$powspec_conf_file"
    conf_file_ptr = Base.unsafe_convert(Cstring, deepcopy(conf_file))

    argc::Cint = 3
    argv_vec = [Base.unsafe_convert(Cstring, deepcopy("POWSPEC")), conf_file_ptr, test_output_ptr]
    argv = Base.unsafe_convert(Ptr{Cstring}, argv_vec)
    
    int_cache = Ptr{Cint}(Base.Libc.calloc(2, sizeof(Cint)))
    cat = CATA(data, data_w, data_fkp, data_nz, rand, rand_w, rand_fkp, rand_nz)
    if precision == :single
        pk = ccall((:compute_pk, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec_f.so"), 
                    Ptr{PK}, 
                    (Ref{CATA}, Cint, Cint, Ptr{Cint}, Cint, Ptr{Cstring}), 
                    cat, save_out, true, int_cache, argc, argv)
        
    elseif precision == :double
        pk = ccall((:compute_pk, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec.so"), 
                    Ptr{PK}, 
                    (Ref{CATA}, Cint, Cint, Ptr{Cint}, Cint, Ptr{Cstring}), 
                    cat, save_out, true, int_cache, argc, argv)
        
    end #if
    if pk == C_NULL
        error("Could not compute power spectra")
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
    out[:n_rand_objects] = unsafe_load(cat.nrand, 1)
    out[:w_rand_objects] = unsafe_load(cat.wrand, 1)
    out[:shot_noise] = unsafe_load(cat.shot, 1)
    out[:normalisation] = unsafe_load(cat.norm, 1)
    out
end #func



function power_spectrum(data::Tuple{Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}}, 
                        data_w::Tuple{AbstractArray{T}, AbstractArray{T}},
                        data_fkp::Tuple{AbstractArray{T}, AbstractArray{T}},
                        data_nz::Tuple{AbstractArray{T}, AbstractArray{T}},
                        rand::Tuple{Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}}, 
                        rand_w::Tuple{AbstractArray{T}, AbstractArray{T}},
                        rand_fkp::Tuple{AbstractArray{T}, AbstractArray{T}},
                        rand_nz::Tuple{AbstractArray{T}, AbstractArray{T}},
                        powspec_conf_file::AbstractString; output_auto = nothing, output_cross = nothing, precision = :double) where T <: AbstractFloat
    save_auto = output_auto != nothing
    save_cross = output_auto != nothing
    
    auto_output = save_auto ? "--auto=[$(output_auto[1]), $(output_auto[2])]" : "--auto=[test1.dat, test1.dat]"
    auto_output_ptr = Base.unsafe_convert(Cstring, deepcopy(auto_output))
    cross_output = save_cross ? "--cross=$output_cross" : "--cross=cross.dat"
    cross_output_ptr = Base.unsafe_convert(Cstring, deepcopy(cross_output))


    conf_file = "--conf=$powspec_conf_file"
    conf_file_ptr = Base.unsafe_convert(Cstring, deepcopy(conf_file))

    argc::Cint = 4
    argv_vec = [Base.unsafe_convert(Cstring, deepcopy("POWSPEC")), conf_file_ptr, auto_output_ptr, cross_output_ptr]
    argv = Base.unsafe_convert(Ptr{Cstring}, argv_vec)
    
    int_cache = Ptr{Cint}(Base.Libc.calloc(2, sizeof(Cint)))
    cat = CATA(data, data_w, data_fkp, data_nz, rand, rand_w, rand_fkp, rand_nz)
    if precision == :single
        pk = ccall((:compute_pk, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec_f.so"), Ptr{PK}, (Ref{CATA}, Cint, Cint, Ptr{Cint}, Cint, Ptr{Cstring}), cat, save_auto & save_cross, true, int_cache, argc, argv)
        
    elseif precision == :double
        pk = ccall((:compute_pk, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec.so"), Ptr{PK}, (Ref{CATA}, Cint, Cint, Ptr{Cint}, Cint, Ptr{Cstring}), cat, save_cross & save_cross, true, int_cache, argc, argv)
        
    end #if
    if pk == C_NULL
        error("Could not compute power spectra")
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
    out[:n_rand_objects] = unsafe_wrap(Vector{Csize_t}, cat.nrand, 2; own = true)
    out[:w_rand_objects] = unsafe_wrap(Vector{Cdouble}, cat.wrand, 2; own = true)
    out[:shot_noise] = unsafe_wrap(Vector{Cdouble}, cat.shot, 2, own = true)
    out[:normalisation] = unsafe_wrap(Vector{Cdouble}, cat.norm, 2, own = true)
    out
end #func

