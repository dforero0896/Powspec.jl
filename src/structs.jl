@with_kw mutable struct DATA
    x::Ptr{Float64} = C_NULL
    w::Float64 = zero(Float64)
end #DATA

@with_kw mutable struct CATA{T}
    num::Cint = 0 #number of catalogs
    data::Ptr{Ptr{DATA}} = C_NULL
    rand::Ptr{Ptr{DATA}} = C_NULL
    ndata::Ptr{Csize_t} = C_NULL
    nrand::Ptr{Csize_t} = C_NULL
    wdata::Ptr{T} = C_NULL
    wrand::Ptr{T} = C_NULL
    alpha::Ptr{T} = C_NULL
    shot::Ptr{T} = C_NULL
    norm::Ptr{T} = C_NULL
    
end #CATA

@with_kw mutable struct PK{T}
    issim::Bool = false
    log::Bool = false         
    isauto::Ptr{Bool} = C_NULL 
    iscross::Bool = false
    nl::Cint = 0
    nbin::Cint = 0      
    nmu::Cint = 0    
    poles::Ptr{Cint} = C_NULL
    los::Ptr{T} = C_NULL
    dk::T = zero(T)
    kedge::Ptr{T} = C_NULL
    k::Ptr{T} = C_NULL
    km::Ptr{T} = C_NULL
    cnt::Ptr{Csize_t} = C_NULL
    lcnt::Ptr{T} = C_NULL
    pl::Ptr{Ptr{Ptr{T}}} = C_NULL
    xpl::Ptr{Ptr{T}} = C_NULL
    nomp::Cint = 0
    pcnt::Ptr{T} = C_NULL
    plcnt::Ptr{T} = C_NULL
end #PK


function build_catalog(data::Tuple{AbstractVector{T}, AbstractVector{T},AbstractVector{T}}, 
                        data_w::AbstractArray{T}, is_sim::Bool; data_fkp = nothing, data_nz = nothing) where T <: Float64
    ndata = length(data_w)
    sumw2_vec = zeros(T, 2)
    wdata = zeros(T,1)
    if is_sim
        out_ptr = ccall((:build_catalog_sim, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec.so"), 
                Ptr{DATA},
                (Csize_t, Ptr{T}, Ptr{T}, Ptr{T}, Ptr{T}, Ptr{T}, Ptr{T}),
                ndata, sumw2_vec, wdata, data[1], data[2], data[3], data_w)
    else
        out_ptr = ccall((:build_catalog_lc, "$(ENV["LIBPOWSPEC_PATH"])/libpowspec.so"), 
                Ptr{DATA},
                (Csize_t, Ptr{T}, Ptr{T}, Ptr{T}, Ptr{T}, Ptr{T}, Ptr{T}, Ptr{T}, Ptr{T}),
                ndata, sumw2_vec, wdata, data[1], data[2], data[3], data_w, data_fkp, data_nz)
    end #if
    out_ptr, sumw2_vec, wdata
end #func




function CATA(num::Int)
    cat = CATA{Float64}()
    cat.num = num

    cat.data = Ptr{Ptr{DATA}}(Base.Libc.malloc(cat.num * sizeof(DATA)))
    cat.rand = Ptr{Ptr{DATA}}(Base.Libc.malloc(cat.num * sizeof(DATA)))
    
    cat.ndata = Ptr{Csize_t}(Base.Libc.calloc(cat.num, sizeof(Csize_t)))
    cat.wdata = Ptr{Float64}(Base.Libc.calloc(cat.num, sizeof(Float64)))  

    cat.nrand = Ptr{Csize_t}(Base.Libc.calloc(cat.num, sizeof(Csize_t)))
    cat.wrand = Ptr{Float64}(Base.Libc.calloc(cat.num, sizeof(Float64)))

    cat.alpha = Ptr{Float64}(Base.Libc.calloc(cat.num, sizeof(Float64)))
    cat.norm = Ptr{Float64}(Base.Libc.calloc(cat.num, sizeof(Float64)))
    cat.shot = Ptr{Float64}(Base.Libc.calloc(cat.num, sizeof(Float64)))

    for i in 1:num
        unsafe_store!(cat.rand,C_NULL,i)
        unsafe_store!(cat.data,C_NULL,i)
    end #for    
    cat
end #func