
function df!(du,u,p,t)
    @unpack_incomeReallocation p 
    Wu = imfilter(u,Wd,Fill(0,u))    
    
    du .=   γ₁ * u + γ₂ * u.^2 + 
    γd * Δ(u.^δ,p) + γa * (∂x(u .* ∂x(Wu,p),p) +  ∂y(u .* ∂y(Wu,p),p)) +
    V_GRD * ( (∂x(u,p)).^2 + (∂y(u,p)).^2 ) + 
    V_GRA * ( (∂x(Wu,p)).^2 + (∂y(Wu,p)).^2 ) + 
    γᵥ * (∂x(u .* ∂xV,p) +  ∂y(u .* ∂yV,p))

    # @show t
    # du .= γd * Δ(u.^δ,p) + p.γᵥ * (∂x(u .* ∂xV,p) +  ∂y(u .* ∂yV,p))
end

function Δ(u,p) return Δ(u,p.Nx,p.Δx) end

function Δ(u,Nx,Δx)

    Δu = similar(u)
    for i in 2:Nx-1, j in 2:Nx-1
        Δu[i,j] = 1/(Δx)^2 * (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] - 4*u[i,j])
    end

    # dirichlet 0
    for i in 2:Nx-1   
        Δu[i,1] = 1/(Δx)^2 * (u[i-1,1] + u[i+1,1] + u[i,2] - 4*u[i,1])
        Δu[i,Nx] = 1/(Δx)^2 * (u[i-1,Nx] + u[i+1,Nx] + u[i,Nx-1] - 4*u[i,Nx])
    end
    for j = 2:Nx-1
        Δu[1,j] = 1/(Δx)^2 * (u[2,j] + u[1,j-1] + u[1,j+1] - 4*u[1,j])
        Δu[Nx,j] = 1/(Δx)^2 * (u[Nx-1,j] + u[Nx,j-1] + u[Nx,j+1] - 4*u[Nx,j])
    end
    Δu[1,1] = 1/(Δx)^2 * (u[2,1] + u[1,2] - 4*u[1,1])
    Δu[1,Nx] = 1/(Δx)^2 * (u[2,Nx] + u[1,Nx-1] - 4*u[1,Nx])
    Δu[Nx,1] = 1/(Δx)^2 * (u[Nx-1,1] + u[Nx,2] - 4*u[Nx,1])
    Δu[Nx,Nx] = 1/(Δx)^2 * (u[Nx-1,Nx] + u[Nx,Nx-1] - 4*u[Nx,Nx])

    return Δu
end

function ∂x(u,p) return ∂x(u,p.Nx,p.Δx) end

function ∂x(u,Nx,Δx)

    ∂x = similar(u)
    for i in 1:Nx, j in 2:Nx-1
        ∂x[i,j] = 1/(2*Δx) * (u[i,j+1] - u[i,j-1])
    end

    # dirichlet 0
    for i in 1:Nx
        ∂x[i,1] =  1/(2*Δx) * u[i,2] 
        ∂x[i,Nx] = - 1/(2*Δx) * u[i,Nx-1]
    end
    return ∂x
end

function ∂y(u,p) return ∂y(u,p.Nx,p.Δx) end

function ∂y(u,Nx,Δx)

    ∂y = similar(u)
    for i in 2:Nx-1, j in 1:Nx
        ∂y[i,j] = 1/(2*Δx) * (u[i-1,j] - u[i+1,j])
    end

    # dirichlet 0
    for j in 1:Nx
        ∂y[1,j] =  - 1/(2*Δx) * u[2,j] 
        ∂y[Nx,j] = 1/(2*Δx) * u[Nx-1,j]
    end
    return ∂y
end



# function conv(u,p,bandwith)
#     """ 2-D convolution between wₕ and u
#         if u (Nx x Ny) return same size matrix 
#     """
#     @unpack_incomeReallocation p;

#     Wu = zeros(size(u))
#     @showprogress for k in 1:Nx, k′ in 1:Nx
#         for i in 1:Nx, j in 1:Nx
#             Wu[k,k′] += W([x[k]-x[i],y[k′]-y[j]],bandwith) * u[i,j] * (Δx)^2
#         end
#     end
#     return Wu            
# end


