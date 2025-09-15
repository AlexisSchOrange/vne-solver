using Graphs, MetaGraphsNext
using GraphRecipes
using Plots
using Colors
using  NetworkLayout


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
        node_size=0.1)
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




function visu_graph_fixed(g)
    pos = spring(g; iterations=2000)
    x = [pos[i][1] for i in 1:nv(g)]
    y = [pos[i][2] for i in 1:nv(g)]

    graphplot(g, 
        x=x,
        y=y,
        node_shape=:circle, 
        node_size = 2.)
end


function visu_partitionning_fixed(g, partitionning)
    pos = spring(g; iterations=2000)
    x = [pos[i][1] for i in 1:nv(g)]
    y = [pos[i][2] for i in 1:nv(g)]

    colors = distinguishable_colors(length(partitionning), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    marker_cols = []
    for i_node in 1:nv(g)
        push!(marker_cols, colors[partitionning[i_node]])
    end


    graphplot(g, 
        x=x,
        y=y,
        node_shape=:circle, 
        markercolor=marker_cols,
        node_size = 2.)
end


function visu_pos(g, pos)

    x = [pos[i][1] for i in 1:nv(g)]
    y = [pos[i][2] for i in 1:nv(g)]

    graphplot(g, 
        x=x,
        y=y,
        node_shape=:circle, 
        node_size = 1)
end
