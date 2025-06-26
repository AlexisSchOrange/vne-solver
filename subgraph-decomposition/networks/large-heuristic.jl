
using Revise
using Statistics
using Graphs, MetaGraphsNext

includet("base-ghostz/utils/utils-subgraphdecompo.jl")
includet("../utils/metis_wrapper.jl")
includet("../heuristics/uepso.jl")




function large_heuristic(instance) 

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    # Decompose the virtual network into subgraph of 10 nodes (approx)
    nb_part = floor(Int, nv(v_network.graph)/10)+1
    v_node_partitionning = partition_graph_metis(v_network.graph, nb_part)
    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs
    println("Decomposition set: ")
        println("For $v_network, there is $(length(vn_decompo.subgraphs)) subgraphs:")

    for subgraph in vn_decompo.subgraphs
        println("       $(subgraph.graph[][:name]) with $(nv(subgraph.graph)) nodes")
    end
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")


    # Decompose the substrate network into the same number of subgraph
    sn_subgraphs = []
    s_node_partitionning = partition_graph_metis(instance.s_network.graph, nb_part)
    for (i_cluster, cluster) in enumerate(s_node_partitionning)
        sub_s_network = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
        push!(sn_subgraphs, Subgraph(sub_s_network, cluster))
    end


    # Map the v_subgraph into s_subgraph
    # Not taking the previous mapping for now, even if that would be probably a bit better...
    overall_cost = 0
    node_placement = zeros(Integer, nv(instance.v_network))
    edge_routing = Dict()

    for i_subgraph in 1:nb_part

        v_subgraph = vn_subgraphs[i_subgraph]
        s_subgraph = sn_subgraphs[i_subgraph]
        sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)
            
        sub_mapping, cost = solve_UEPSO(sub_instance; nb_particle=25, nb_iter=50, time_max=0.1, print_things=false)
        overall_cost += cost



        for v_node in vertices(v_subgraph.graph)

            original_v_node = v_subgraph.nodes_of_main_graph[v_node]
            original_s_node =  s_subgraph.nodes_of_main_graph[sub_mapping.node_placement[v_node]]
            node_placement[original_v_node] = original_s_node

        end


        for v_edge in edges(v_subgraph.graph)
            original_v_edge = get_edge(v_network, v_subgraph.nodes_of_main_graph[src(v_edge)], v_subgraph.nodes_of_main_graph[dst(v_edge)])
            used_edges = []
            for s_edge in sub_mapping.edge_routing[v_edge].edges
                real_s_edge = get_edge(s_network_dir, s_subgraph.nodes_of_main_graph[src(s_edge)], s_subgraph.nodes_of_main_graph[dst(s_edge)])
                push!(used_edges, real_s_edge)
            end
            edge_routing[original_v_edge] = order_path(s_network_dir, used_edges, node_placement[src(original_v_edge)], node_placement[dst(original_v_edge)]) 
        end
    end        

    edge_routing, additional_routing_cost = route_cut_edges(instance, vn_decompo, node_placement, edge_routing)
    overall_cost += additional_routing_cost

    println("We obtained a mapping of cost $overall_cost !")
    return Mapping(v_network, s_network, node_placement, edge_routing)
end







function partition_graph_metis(graph, nb_clusters)
    
    # Partitionning. Since connectivity is enforced, sometime, it will not the best
    println("$nb_clusters clusters to do... Partitionning done with METIS!")
    best_clusters = nothing
    best_imb = 10000
    imb = [1.01, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
    for imbalance in imb
        partition = partition_metis(graph, nb_clusters, imbalance)

        clusters = [Vector{Int64}() for i in 1:nb_clusters]
        for s_node in vertices(graph)
            push!(clusters[partition[s_node]], s_node)
        end

        moyenne = mean([length(cluster) for cluster in clusters])
        current_imb = maximum([length(cluster) / moyenne for cluster in clusters])
        if current_imb < 1.10
            best_clusters = clusters
            best_imb = current_imb
            break
        end
        if current_imb < best_imb
            best_imb = current_imb
            best_clusters = clusters
        end
    end
    
    println("Best partition found has imbalance of $best_imb.")

    return best_clusters
end





function route_cut_edges(instance, vn_decompo, v_node_placement, edge_routing)

    s_network_dir = instance.s_network_dir
    s_network_dir_copy = deepcopy(instance.s_network_dir)
    additional_cost = 0

    # tackle the v_edges already done
    for v_edge in keys(edge_routing)

        path = edge_routing[v_edge]
        
        for s_edge in path.edges
            set_attribute_edge(s_network_dir_copy, s_edge, :cap,  get_attribute_edge(s_network_dir_copy, s_edge, :cap)-1)
            set_attribute_edge(s_network_dir_copy, get_reverse_edge(s_network_dir_copy, s_edge), :cap,  get_attribute_edge(s_network_dir_copy, get_reverse_edge(s_network_dir_copy, s_edge), :cap)-1)
            if get_attribute_edge(s_network_dir_copy, s_edge, :cap) <= 0
                #println("Well it is time to stop using $s_edge kids")
                distmx[src(s_edge), dst(s_edge)] = 0
                distmx[dst(s_edge), src(s_edge)] = 0
                rem_edge!(s_network_dir_copy, src(s_edge), dst(s_edge))
                rem_edge!(s_network_dir_copy, dst(s_edge), src(s_edge))
            end
        end

    end


    # dstmax matrix
    distmx = zeros(Int, nv(s_network_dir_copy), nv(s_network_dir_copy))
    for s_edge in edges(s_network_dir_copy)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir_copy, s_edge, :cost)
    end
    
    for v_edge in vn_decompo.v_edges_master

        s_src = v_node_placement[src(v_edge)]
        s_dst = v_node_placement[dst(v_edge)]
        shortest_path = a_star(s_network_dir_copy, s_src, s_dst, distmx)


        if shortest_path == []
            #println("No shortest path found: the graph is full!")
            #println("I had the following routing: $edge_routing")
            return Dict(), 99999999
        end

        edge_routing[v_edge] = order_path(s_network_dir, shortest_path, s_src, s_dst) 

        for s_edge in shortest_path
            set_attribute_edge(s_network_dir_copy, s_edge, :cap,  get_attribute_edge(s_network_dir_copy, s_edge, :cap)-1)
            set_attribute_edge(s_network_dir_copy, get_reverse_edge(s_network_dir_copy, s_edge), :cap,  get_attribute_edge(s_network_dir_copy, get_reverse_edge(s_network_dir_copy, s_edge), :cap)-1)
            if get_attribute_edge(s_network_dir_copy, s_edge, :cap) <= 0
                distmx[src(s_edge), dst(s_edge)] = 0
                distmx[dst(s_edge), src(s_edge)] = 0
                rem_edge!(s_network_dir_copy, src(s_edge), dst(s_edge))
                rem_edge!(s_network_dir_copy, dst(s_edge), src(s_edge))
            end
        end

        additional_cost += edge_routing[v_edge].cost
    end


    return edge_routing, additional_cost

end



