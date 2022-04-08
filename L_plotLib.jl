

function plotSol(sol,p; t_step = 0.1)
    
    fig = figure()
    t_plot = sol.t[1]:t_step:sol.t[end]

    for t in t_plot
        fig.clear()
        imshow(log.(sol(t).+1),origin="lower")
        title("T = $t")
        colorbar()
        pause(0.1)
    end
end

function plotSurfSol(sol,p; t_step = 0.1)
    fig = figure()
    t_plot = sol.t[1]:t_step:sol.t[end]

    X = repeat(p.x,1,p.Ny)
    Y = repeat(p.y',p.Nx,1)
    for t in t_plot
        fig.clear()
        surf(X,Y,log.(sol(t).+1),facecolors=get_cmap("jet")(log.(sol(t).+1)))
        title("T = $t")
        pause(0.1)
    end
end