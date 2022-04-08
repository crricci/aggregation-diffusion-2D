
@with_kw struct incomeReallocation{T}
    # aggregation diffusion
    γa::T = - 0.1       # aggregation (< 0)
    γd::T = 0.025       # diffusion
    h::T = 0.3          # bandwith  
    γ₁::T = 0.01        # growth linear
    γ₂::T = -0.02       # growth quadratic
    γᵥ::T = 1.0         # potential strength
    
    # reallocation gain
    V_GRD::T = 0.0      # reallocation gain physical capital
    V_GRA::T = 0.0      # reallocation gain labor 

    # domain
    L::T = 2.0
    T_end::T = 1.0
    borderLength::T = 0.5
    @assert h < borderLength

    # numerical
    Δx::T = 1e-2
    Nx::Int = Int(L/Δx)
    Ny::Int = Int(L/Δx)
    x::LinRange{T, Int64} = LinRange(0,L,Nx)
    y::LinRange{T, Int64} = LinRange(0,L,Nx)
    Wd::Matrix{T} = make_WDiscrete(Δx,h)

    # confining potential around border
    V::Matrix{T} = make_Vflat(Nx,Ny,Δx,h/2,Wd)

    # initial condition
    u₀::Matrix{T} = make_u₀(Nx,Ny,Δx)
    # NBumps::Int = 1
    # bumpCenter::Vector{T} = [(0,0)]
    # bumpRadius::Vector{T} = [L/4]
    # y₀::Vector{T} = initialCondition(NBumps,bumpCenter,bumpRadius,Δx,Nx,x)
    # centerOfMass::T = sum(y₀ .* x) * Δx
    # @assert length(bumpCenter) == NBumps
    # @assert length(bumpRadius) == NBumps

    # human capital (optional)
    β::Float64 = 1.0
    γ::Float64 = 1.0
    δ::Float64 = 1.0
    ϕ::Float64 = 1.0
    ψ::Float64 = 1.0

end

function make_u₀(Nx,Ny,Δx)
    u₀ = zeros(Nx,Ny)

    u₀[Int(Nx/2):Int(Nx/2)+10,Int(Nx/2):Int(Nx/2)+10] .= 1

    return u₀
end

function make_Vflat(Nx,Ny,Δx,borderLength,Wd)
    V = -ones(Nx,Ny)
    borderNpts = ceil(Int,borderLength/Δx)
    V[1:borderNpts,:] .= 0.0
    V[:,1:borderNpts] .= 0.0
    V[end-borderNpts:end,:] .= 0.0
    V[:,end-borderNpts:end] .= 0.0

    V = imfilter(V,Wd,Fill(0,V))
    return V
end

function W(x)
    """ smoothing kernel """
    return norm(x) <= 1 ? 1-norm(x) : 0.0
end

function W(x,h)
    """ rescaled kernel """
    return 1/h*W(x/h)
end

function make_WDiscrete(Δx,bandwith)

    Npt = ceil(Int,bandwith/Δx)
    Wd = zeros(2Npt+1,2Npt+1)
    x = -Npt*Δx:Δx:Npt*Δx
    y = -Npt*Δx:Δx:Npt*Δx
    for i in 1:length(x), j in 1:length(y)
        Wd[i,j] = W([x[i],y[j]],bandwith)
    end
    return normalize(Wd,1)
end

