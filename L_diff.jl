
function df!(du,u,p,t)
    (∂x,∂y,Δ,∂xV,∂yV),p = q
    @unpack_incomeReallocation p 
    
    # Wu = imfilter(u,Wd,Fill(0,u))    
    

    
    # du .=   γ₁ * u + γ₂ * u.^2 + 
    # γd * (Δ * u.^δ) +  p.γᵥ * (∂x * (u .* ∂xV) +  ∂y * (u .* ∂yV)) + p.γa * ∇⋅(u.^p.γ.*∇wu,p) + 
    # + p.V_GRD * abs.(∇u).^2 + p.V_GRA * abs.(∇(Wu,p)).^2 + 
    # + γᵥ*∇⋅(u.*∇V,p) 

    du .= p.γd * (Δ * u) 
    # du .= p.γd * (Δ * u) + p.γᵥ * (∂x * (u .* ∂xV) +  ∂y * (u .* ∂yV))
    # @show t
end

function Δ_hand(u,p)
    @unpack_incomeReallocation p 

    Δu = similar(u)
    for i in 2:Nx-1, j in 2:Nx-1
        Δu[i,j] = 1/(Δx)^2 * (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] - 4*u[i,j])
    end

    for i in 2:Nx-1   # dirichlet 0
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
    Δu[Nx,Nx] = 1/(Δx)^2 * (u[2,1] + u[1,2] - 4*u[1,1])

    return Δu
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


