using Graphs, MetaGraphsNext
using GraphRecipes
using Plots

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

function print_graph(mg)
    println("Metagraph ")
    print("Nodes:")
    for node in vertices(mg)
        print("\nNode " * string(node) * " with ")
        for k in keys(mg[node])
            print(string(k) * ": " * string(mg[node][k]))
        end
    end
    print("\nEdges:")
    for edge in edges(mg)
        println("\nEdge " * string(edge))
        for k in keys(mg[src(edge), dst(edge)])
            print(" with " * string(k) * ": " * string(mg[src(edge), dst(edge)][k]))
        end
    end
end

function print_substrate(mg)
    println("Metagraph ")
    for node in vertices(mg)
        println("Node " * string(node) * " with cap " * string(mg[node][:cap]) * " and cost " * string(mg[node][:cost]))
    end
    for edge in edges(mg)
        println("Edge " * string(edge) * " with cap " * string(mg[src(edge), dst(edge)][:cap]) * " and cost " * string(mg[src(edge), dst(edge)][:cost]) )
    end
end

function print_virtual(mg)
    println("Metagraph ")
    for node in vertices(mg)
        println("Node " * string(node) * " with dem " * string(mg[node][:dem]))
    end
    for edge in edges(mg)
        println("Edge " * string(edge) * " with dem " * string(mg[src(edge), dst(edge)][:dem]) )
    end
end



# Visu de la solution !