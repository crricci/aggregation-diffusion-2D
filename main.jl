# Started on April 6th 2022
# Cristiano Ricci - cristiano.ricci6@gmail.com

include("L_LoadAll.jl")


function main(p)
    
    T = eltype(p.γa)
    # 1st order
    ∂x = CenteredDifference{1}(1, p.approx_order, p.Δx, p.Nx)
    ∂y = CenteredDifference{2}(1, p.approx_order, p.Δx, p.Ny)

    # 2nd order
    ∂xx = CenteredDifference{1}(2, p.approx_order, p.Δx, p.Nx)
    ∂yy = CenteredDifference{2}(2, p.approx_order, p.Δx, p.Ny)
    Δ = ∂xx + ∂yy
    # boundary conditions
    # Qx, Qy = Neumann0BC(T, (p.Δx, p.Δx), p.approx_order, (p.Nx, p.Nx))
    Qx, Qy = DirichletBC(zero(T), zero(T),(p.Nx,p.Nx))
    Q = compose(Qx, Qy)
    
    ∂x = ∂x * Q
    ∂y = ∂y * Q
    Δ = Δ * Q
    ∂xV = ∂x * p.V
    ∂yV = ∂y * p.V

    T_span = (0.0,p.T_end)
    q = ((∂x,∂y,Δ,∂xV,∂yV),p);
    prob = ODEProblem(df!, p.u₀, T_span, q)
    sol = solve(prob,Tsit5())


end
