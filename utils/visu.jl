using Graphs, MetaGraphsNext
using GraphRecipes
using Plots
using Colors

function visu_graph(g)
    w = []
    for i_node in 1:nv(g)
        if i_node < 10
            push!(w, 10)
        else
            push!(w, 1)
        end
    end
    graphplot(g, 
        node_weights = w,
        names=string.(1:nv(g)),
        curvature_scalar=0.01, 
        node_size = 0.2)
end



function visu_partitioning(g, partitionning)
    w = []
    for i_node in 1:nv(g)
        if i_node < 10
            push!(w, 10)
        else
            push!(w, 1)
        end
    end

    colors = distinguishable_colors(length(partitionning), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)

    marker_cols = []
    for i_node in 1:nv(g)
        push!(marker_cols, colors[partitionning[i_node]])
    end

    p = graphplot(g, 
        node_weights=w,
        names=string.(1:nv(g)),
        markercolor=marker_cols,
        curvature_scalar=0.01, 
        node_size=0.2)
    display(p) 

end


function write_single_partitionning(g, nodes)
    w = []
    for i_node in 1:nv(g)
        if i_node < 10
            push!(w, 10)
        else
            push!(w, 1)
        end
    end

    colors = distinguishable_colors(2, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)

    marker_cols = []
    for i_node in 1:nv(g)
        if i_node ∈ nodes
            push!(marker_cols, colors[2])
        else
            push!(marker_cols, colors[1])
        end
    end

    p = graphplot(g, 
        node_weights=w,
        names=string.(1:nv(g)),
        markercolor=marker_cols,
        curvature_scalar=0.01, 
        node_size=0.2)
    savefig(p, "") 

end



function write_added_nodes(g, nodes, added, name)
    w = []
    for i_node in 1:nv(g)
        if i_node < 10
            push!(w, 10)
        else
            push!(w, 1)
        end
    end

    colors = distinguishable_colors(3, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)

    marker_cols = []
    for i_node in 1:nv(g)
        if i_node ∈ added
            push!(marker_cols, colors[2])
        elseif i_node ∈ nodes
            push!(marker_cols, colors[3])
        else
            push!(marker_cols, colors[1])
        end
    end

    p = graphplot(g, 
        node_weights=w,
        names=string.(1:nv(g)),
        markercolor=marker_cols,
        curvature_scalar=0.01, 
        node_size=0.2)
    savefig(p, "$name") 

end


# Visu de la solution !



# visu perfo profile

function performance_profile(times, names_algo)
    time_max = maximum([sum(times[i]) for i in 1:length(times)])
    plot(fmt = :png)

    for (i_algo, times) in enumerate(times)

        x = [0.]
        y = [0.]

        nb_instance = 0
        nb_total_instances = length(times)
        time_total = 0.

        for time in times
            push!(y, nb_instance/nb_total_instances)
            nb_instance +=1
            push!(y, nb_instance/nb_total_instances)

            time_total += time
            push!(x, time_total)
            push!(x, time_total)
        end
        push!(x, time_max)
        push!(y, 1.)

        plot!(x, y, label = names_algo[i_algo], legend = :bottomright, xaxis = "Time (s)", yaxis = "% solved instances",linewidth=3)
    end

    savefig("benchmark_results.pdf" )# Display the plot
end