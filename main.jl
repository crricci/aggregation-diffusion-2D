# Started on April 6th 2022
# Cristiano Ricci - cristiano.ricci6@gmail.com

include("L_LoadAll.jl")


function main(p)
    
    T_span = (0.0,p.T_end)
    prob = ODEProblem(df!, p.u₀, T_span, p)
    sol = solve(prob,ROCK2(),progress=true)
    return sol,p
end
