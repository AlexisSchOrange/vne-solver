using Graphs, MetaGraphsNext


### === Structs
struct NetworkDecompositionOverlapping
    subgraphs
    v_nodes_assignment
    v_edges_master
    overlapping_nodes
end

struct Subgraph
    graph
    nodes_of_main_graph
end


function print_stuff_subgraphs(original_graph, subgraphs)

    println("For $original_graph, there is $(length(subgraphs)) subgraphs:")
    for subgraph in subgraphs
        println("       $(subgraph.graph[][:name]) with $(nv(subgraph.graph)) nodes and $(ne(subgraph.graph)) edges")
    end
    
end


function set_up_decompo_overlapping(instance, node_partitionning)

    vn = instance.v_network

    # Here, node assigment is a dict, I don't know why. And node assigment ofa  node is alsoa dic,
    # because it can be in several sub-networks.
    node_assignment = Dict()
    for v_node in vertices(vn)
        node_assignment[v_node] = Dict()
    end

    # getting the subgraphs and the node assignment
    subgraphs = []
    for (i_subgraph, v_nodes) in enumerate(node_partitionning)
        subgraph = Subgraph(my_induced_subgraph(vn, v_nodes, "subgraph_$i_subgraph"), v_nodes)
        for (i_node, v_node) in enumerate(v_nodes)
            node_assignment[v_node][subgraph] = i_node
        end
        push!(subgraphs, subgraph)
    end


    # finding out the master virtual edges
    v_edge_master = [] 
    for v_edge in edges(vn)
        in_master = true
        for subgraph_src in keys(node_assignment[src(v_edge)])
            for subgraph_dst in keys(node_assignment[dst(v_edge)])
                if subgraph_src == subgraph_dst
                    in_master = false
                end
            end
        end
        if in_master
            push!(v_edge_master, v_edge)
        end
    end


    # finding and removing overlapping edges!
    for v_edge in edges(vn)
        common_subgraph = collect(keys(node_assignment[src(v_edge)]) ∩ keys(node_assignment[dst(v_edge)]))
        if length(common_subgraph) ≥ 2
            println("A virtual edge is in two subgraph! $common_subgraph")
            for subgraph in common_subgraph[2:end]
                src_in_subgraph = node_assignment[src(v_edge)][subgraph]
                dst_in_subgraph = node_assignment[dst(v_edge)][subgraph]
                rem_edge!(subgraph.graph, src_in_subgraph, dst_in_subgraph)
                if !is_connected(subgraph.graph)
                    println("Well, the graph is not connected anymore... That sucks... Fix this or include overlapping edges...")
                end
            end
            println("I removed the edge: $common_subgraph")
        end
    end


    # finding overlapping nodes
    v_node_overlapping = Dict()
    for v_node in vertices(vn)
        if length(keys(node_assignment[v_node])) > 1
            v_node_overlapping[v_node] = keys(node_assignment[v_node])
            println("$v_node is overlapping !")
        end
    end

    vn_decompo = NetworkDecompositionOverlapping(subgraphs, node_assignment, v_edge_master, v_node_overlapping)

    return vn_decompo
end

