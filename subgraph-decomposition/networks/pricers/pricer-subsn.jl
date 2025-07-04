includet("../../../utils/import_utils.jl")
includet("../../../utils/metis_wrapper.jl")

using CPLEX, JuMP



struct PricerSubSubstrate
    vn_subgraph
    sn_subgraph
    sub_instance
    original_instance
    model
end

function Base.show(io::IO, pricer::PricerSubSubstrate)
    println(io, "Some pricer for subvn $(pricer.vn_subgraph.graph[][:name]) on subsn $(pricer.sub_instance.s_network[][:name])")
end



function set_up_pricer_sn_decompo(instance, vn_subgraph, clusters, type_pricer ="")

    s_network = instance.s_network
    pricers = []
    for (i_cluster, cluster) in enumerate(clusters)
        sub_s_network = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")

        sn_subgraph = Subgraph(sub_s_network, cluster)

        sub_instance = Instance(vn_subgraph.graph, sub_s_network)


        model = Model(CPLEX.Optimizer)

        if type_pricer == "ghost"
            set_up_model_pricer_subsn_ghost(model, sub_instance, instance, vn_subgraph, sn_subgraph)
        elseif type_pricer == "constraint"
            set_up_model_pricer_subsn_cons(model, sub_instance, instance, vn_subgraph, sn_subgraph)
        else
            set_up_model_pricer_subsn_basic(model, sub_instance, instance, vn_subgraph, sn_subgraph)
        end



        push!(pricers,  PricerSubSubstrate(vn_subgraph, sn_subgraph, sub_instance, instance, model))
    end

    return pricers
end



# computing the sn subgraphs, using a partition of the sn network again..
function get_sn_decompo(s_network, nb_clusters, nb_nodes_per_clusters)


    graph = instance.v_network.graph
    
    #1) Partitionning. Since connectivity is enforced, sometime, it will not the best and quite unbalanced,
    # thus is do some different imbalance. It's very fast, and sometime, it's more balanced... Go understand why!
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




    # 2) adding nodes. For now, any adjacent nodes will do.
    clusters = best_clusters
    i_cluster = 1
    subgraphs = []

    for cluster in clusters
        all_neighbors = Dict()
        for s_node in cluster
            for neigh in neighbors(s_network, s_node)
                if neigh ∉ cluster
                    if neigh ∉ keys(all_neighbors)
                        all_neighbors[neigh] = 1
                    else
                        all_neighbors[neigh] += 1
                    end
                end
            end
        end

        added = []

        nb_nodes_with_capacity=0
        for s_node in cluster
            if s_network[s_node][:cap] >= 1
                nb_nodes_with_capacity+=1
            end
        end

        while nb_nodes_with_capacity < nb_nodes_per_clusters

            # ranking the neighbors
            ranking = sort(collect(keys(all_neighbors)), by = x->all_neighbors[x], rev=true)
            # add the most connected neighbor
            new_s_node = ranking[1]
            if s_network[new_s_node][:cap] >= 1
                nb_nodes_with_capacity+=1
            end
            push!(cluster, new_s_node)
            push!(added, new_s_node)
            delete!(all_neighbors, new_s_node)

            for neigh in neighbors(s_network, new_s_node)
                if neigh ∉ cluster
                    if neigh ∉ keys(all_neighbors)
                        all_neighbors[neigh] = 1
                    else
                        all_neighbors[neigh] += 1
                    end
                end
            end
        end

        #println("Well I have $(length(cluster)) whereas I need $(nb_nodes_per_clusters)")
        sub_s_network = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
        i_cluster += 1
        push!(subgraphs, sub_s_network)

        #write_added_nodes(s_network.graph, cluster, added, sub_s_network[][:name])
    end



    return clusters

end


# computing the sn decompo. It depends on the subvirtual network...
function get_sn_decompo_kahip(s_network, nb_clusters, nb_nodes_per_clusters)


    # 1 : Partitionner
    inbalance = 0.10
    partition = partition_kahip(s_network.graph, nb_clusters, inbalance)
    clusters = [Vector{Int64}() for i in 1:nb_clusters]
    for s_node in vertices(s_network)
        push!(clusters[partition[s_node]], s_node)
    end

    # 2 : Corriger

    # 2a) correction by removing unconnected sets
    for cluster in clusters
        simple_subgraph, vmap = induced_subgraph(s_network.graph, cluster)
        if !is_connected(simple_subgraph)
            println("One of the subgraph is not connected ! :(")
            #println("At the beginning the cluster is : $cluster")
            components = connected_components(simple_subgraph)
            component_sorted = sort(components, by=x->length(x), rev=true)
            new_cluster = [vmap[i] for i in component_sorted[1]]
            for subcluster in component_sorted[2:length(component_sorted)]
                #Let's add all those nodes to a (most) connected subgraph
                nodes_original = [vmap[node] for node in subcluster]
                subgraph_neighbors = zeros(Int, nb_clusters)
                for node in nodes_original
                    for neighbor in neighbors(s_network, node)
                        if neighbor ∉ cluster
                            subgraph_neighbors[partition[neighbor]] += 1
                        end
                    end
                end
                most_connected_subgraph = sortperm(subgraph_neighbors, rev=true)
                append!(clusters[most_connected_subgraph[1]], nodes_original)
                #println("Well let's add $nodes_original to cluster $(clusters[most_connected_subgraph[1]])")
                filter!(e->e∉nodes_original, cluster)
            end
        end
    end

    # 2b) adding nodes. For now, any adjacent nodes will do.

    i_cluster = 1
    subgraphs = []


    for cluster in clusters

        all_neighbors = Dict()
        for s_node in cluster
            for neigh in neighbors(s_network, s_node)
                if neigh ∉ cluster
                    if neigh ∉ keys(all_neighbors)
                        all_neighbors[neigh] = 1
                    else
                        all_neighbors[neigh] += 1
                    end
                end
            end
        end

        added = []

        nb_nodes_with_capacity=0
        for s_node in cluster
            if s_network[s_node][:cap] >= 1
                nb_nodes_with_capacity+=1
            end
        end

        while nb_nodes_with_capacity < nb_nodes_per_clusters

            # ranking the neighbors
            ranking = sort(collect(keys(all_neighbors)), by = x->all_neighbors[x], rev=true)
            # add the most connected neighbor
            new_s_node = ranking[1]
            if s_network[new_s_node][:cap] >= 1
                nb_nodes_with_capacity+=1
            end
            push!(cluster, new_s_node)
            push!(added, new_s_node)
            delete!(all_neighbors, new_s_node)

            for neigh in neighbors(s_network, new_s_node)
                if neigh ∉ cluster
                    if neigh ∉ keys(all_neighbors)
                        all_neighbors[neigh] = 1
                    else
                        all_neighbors[neigh] += 1
                    end
                end
            end
        end

        #println("Well I have $(length(cluster)) whereas I need $(nb_nodes_per_clusters)")
        #sub_s_network = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
        i_cluster += 1
        #push!(subgraphs, sub_s_network)

        #write_added_nodes(s_network.graph, cluster, added, sub_s_network[][:name])
    end

    println("Length of each substrate subgraph:")
    for cluster in clusters
        print(" $(length(cluster))")
    end
    return clusters

end


function set_up_model_pricer_subsn_basic(model, sub_instance, original_instance, vn_subgraph, sn_subgraph)


    s_network = sub_instance.s_network
    s_network_dir = sub_instance.s_network_dir
    v_network = sub_instance.v_network


    original_v_network = original_instance.v_network
    original_s_network = original_instance.s_network


    set_silent(model)
    ### Variables
    @variable(model, x[vertices(v_network), vertices(s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);
 

    ### Objective
    placement_cost = @expression(model, sum( s_network[s_node][:cost] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);

    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
    end

    # if one to one : one virtual node per substrate node
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end



    # node capacity
    for s_node in vertices(s_network)
        @constraint(model, 
            sum( x[v_node, s_node] 
                for v_node in vertices(v_network) ) 
            <= 
            s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity (undirected version)
    for s_edge in edges(s_network)
        @constraint(model, 
        sum( (y[v_edge, s_edge] + y[v_edge, get_reverse_edge(s_network_dir, s_edge)]  )
                for v_edge in edges(vn_subgraph.graph)) 
            <= 
            s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(s_network)
        for v_edge in edges(v_network)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end


    ## Departure cst : Node + Edge
    for s_node in vertices(s_network)
        for v_node in vertices(v_network)
            for v_edge in get_out_edges(v_network, v_node)
                @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) >= x[v_node, s_node])
            end
        end
    end
    
    # Removing some x, due to sub_vn and sub_sn connectivity & capacities. THIS IS LOCAL TO THE PRICER!
    for v_node in vertices(v_network)
        necessary_bw = degree(v_network, v_node)
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end

    # Removing some x, due to the overall v_network and s_network
    for v_node in vertices(v_network)
        necessary_bw = degree(original_v_network, vn_subgraph.nodes_of_main_graph[v_node])
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(original_s_network, sn_subgraph.nodes_of_main_graph[s_node], neighbor) for neighbor in neighbors(original_s_network, sn_subgraph.nodes_of_main_graph[s_node])]
            available_bw = sum(original_s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end
    

end



function set_up_model_pricer_subsn_cons(model, sub_instance, original_instance, vn_subgraph, sn_subgraph)


    s_network = sub_instance.s_network
    s_network_dir = sub_instance.s_network_dir
    v_network = sub_instance.v_network


    original_v_network = original_instance.v_network
    original_s_network = original_instance.s_network


    set_silent(model)
    ### Variables
    @variable(model, x[vertices(v_network), vertices(s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);
 

    ### Objective
    placement_cost = @expression(model, sum( s_network[s_node][:cost] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);

    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
    end

    # if one to one : one virtual node per substrate node
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end



    # node capacity
    for s_node in vertices(s_network)
        @constraint(model, 
            sum( x[v_node, s_node] 
                for v_node in vertices(v_network) ) 
            <= 
            s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity (undirected version)
    for s_edge in edges(s_network)
        @constraint(model, 
        sum( (y[v_edge, s_edge] + y[v_edge, get_reverse_edge(s_network_dir, s_edge)]  )
                for v_edge in edges(vn_subgraph.graph)) 
            <= 
            s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(s_network)
        for v_edge in edges(v_network)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end


    ## Departure cst : Node + Edge
    for s_node in vertices(s_network)
        for v_node in vertices(v_network)
            for v_edge in get_out_edges(v_network, v_node)
                @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) >= x[v_node, s_node])
            end
        end
    end
    
    # Removing some x, due to sub_vn and sub_sn connectivity & capacities. THIS IS LOCAL TO THE PRICER!
    for v_node in vertices(v_network)
        necessary_bw = degree(v_network, v_node)
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end

    # Removing some x, due to the overall v_network and s_network
    for v_node in vertices(v_network)
        necessary_bw = degree(original_v_network, vn_subgraph.nodes_of_main_graph[v_node])
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(original_s_network, sn_subgraph.nodes_of_main_graph[s_node], neighbor) for neighbor in neighbors(original_s_network, sn_subgraph.nodes_of_main_graph[s_node])]
            available_bw = sum(original_s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end


    # Additional constraint ! BUT must be related to original substrate network !
    for s_node in vertices(s_network)
        s_edges_incident_original_sn = [get_edge(original_s_network, sn_subgraph.nodes_of_main_graph[s_node], neighbor) for neighbor in neighbors(original_s_network, sn_subgraph.nodes_of_main_graph[s_node])]
        available_bw = sum(original_s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident_original_sn; init=0.0)

        @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node) for v_edge in edges(v_network))
                                + sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node) for v_edge in edges(v_network)) 
                                + sum(x[v_node, s_node] * (degree(original_v_network, vn_subgraph.nodes_of_main_graph[v_node]) - degree(v_network, v_node)) for v_node in vertices(v_network)) 
                                    <= available_bw)
    end

end



function set_up_model_pricer_subsn_ghost(model, sub_instance, original_instance, vn_subgraph, sn_subgraph)

    s_network = sub_instance.s_network
    s_network_dir = sub_instance.s_network_dir
    v_network = sub_instance.v_network


    original_v_network = original_instance.v_network
    original_s_network = original_instance.s_network

    set_silent(model)

    ghost_nodes = []
    ghost_nodes_appearances = Dict() # For each ghost node, we record how many neighbor it has in the subgraph.
    # ^ useful for additional constraints ?
    ghost_edges = []
    
    for v_node in vertices(v_network)
        for v_neighbor in neighbors(original_v_network, vn_subgraph.nodes_of_main_graph[v_node])
            if (v_neighbor ∉ ghost_nodes) && (v_neighbor ∉ vn_subgraph.nodes_of_main_graph)
               push!(ghost_nodes, v_neighbor) 
               ghost_nodes_appearances[v_neighbor] = 0
            end

            if v_neighbor ∉ vn_subgraph.nodes_of_main_graph
                edge_neigh = Dict()
                edge_neigh[:src] = v_node
                edge_neigh[:dst] = v_neighbor
                push!(ghost_edges, edge_neigh)
                ghost_nodes_appearances[v_neighbor] = ghost_nodes_appearances[v_neighbor] + 1
            end
        end
    end

    #println("Well, subgraph vn : $(vn_subgraph.nodes_of_main_graph), and the neighbors are: $(ghost_nodes), and ghost edges : $ghost_edges")
    
    ### Variables
    @variable(model, x[vertices(v_network), vertices(s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);
    @variable(model, x_ghost[ghost_nodes, vertices(s_network)], binary=true);
    @variable(model, y_ghost[ghost_edges, edges(s_network_dir)], binary=true);


    ### Objective
    placement_cost = @expression(model, sum( s_network[s_node][:cost] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);


    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
    end

    # for ghost nodes!
    for v_neighbor in ghost_nodes
        @constraint(model, sum(x_ghost[v_neighbor, s_node] for s_node in vertices(s_network)) == 1)
    end
    
    

    # if one to one : one virtual node per substrate node
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end


    # node capacity
    for s_node in vertices(s_network)
        @constraint(model, 
            sum( x[v_node, s_node] 
                for v_node in vertices(v_network))
            + sum( x_ghost[v_neighbor, s_node]
                for v_neighbor in ghost_nodes)
            <= 
            s_network[s_node][:cap] )
    end



    ## Edges 
    
    # edge capacity (undirected version)
    for s_edge in edges(s_network)
        @constraint(model, 
            sum( (y[v_edge, s_edge] + y[v_edge, get_reverse_edge(s_network_dir, s_edge)] ) 
                for v_edge in edges(v_network)) 
            + sum((y_ghost[v_adjacent_edge, s_edge] + y_ghost[v_adjacent_edge, get_reverse_edge(s_network_dir, s_edge)]  )
                for v_adjacent_edge in ghost_edges) 
            <= 
            s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    

    # Flow conservation
    for s_node in vertices(s_network)
        for v_edge in edges(v_network)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end

    # for ghosts edges!
    for s_node in vertices(s_network)
        for v_adjacent_edge in ghost_edges
            @constraint(model, 
                x[v_adjacent_edge[:src], s_node] - x_ghost[v_adjacent_edge[:dst], s_node] 
                ==
                sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end


    ## Departure cst : Node + Edge
    for s_node in vertices(s_network)
        for v_node in vertices(v_network)
            for v_edge in get_out_edges(v_network, v_node)
                @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) >= x[v_node, s_node])
            end
        end
    end

    
    # for ghost edges!
    for s_node in vertices(s_network)
        for v_adjacent_edge in ghost_edges
            @constraint(model, sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[v_adjacent_edge[:src], s_node])
        end
    end
    


    # Removing some x, due to sub_vn and sub_sn connectivity & capacities. THIS IS LOCAL TO THE PRICER!
    for v_node in vertices(v_network)
        necessary_bw = degree(v_network, v_node)
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end

    # Removing some x, due to original vn AND sn connectivity & capacities. THIS IS OVERALL FOR THE PRICER!
    for v_node in vertices(v_network)
        necessary_bw = degree(original_v_network, vn_subgraph.nodes_of_main_graph[v_node])
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(original_s_network, sn_subgraph.nodes_of_main_graph[s_node], neighbor) for neighbor in neighbors(original_s_network, sn_subgraph.nodes_of_main_graph[s_node])]
            available_bw = sum(original_s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end
        
    # ALSO ON GHOST NODES: Removing some x, due to sub_vn and sub_sn connectivity & capacities. THIS IS LOCAL TO THE PRICER!
    for v_neighbor in ghost_nodes
        necessary_bw = ghost_nodes_appearances[v_neighbor]
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, x_ghost[v_neighbor, s_node] == 0)
            end 
        end
    end

    # ALSO ON GHOST NODES: Removing some x, due to original vn AND sn connectivity & capacities. THIS IS OVERALL FOR THE PRICER!
    for v_neighbor in ghost_nodes
        necessary_bw = degree(original_v_network, v_neighbor)
        for s_node in vertices(s_network)
            s_edges_incident = [get_edge(original_s_network, sn_subgraph.nodes_of_main_graph[s_node], neighbor) for neighbor in neighbors(original_s_network, sn_subgraph.nodes_of_main_graph[s_node])]
            available_bw = sum(original_s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.0)
            if necessary_bw > available_bw
                @constraint(model, x_ghost[v_neighbor, s_node] == 0)
            end 
        end
    end

    #= Hum, at most 3 ?
    print("Hum ")
    for v_adjacent_edge in ghost_edges
        @constraint(model, sum(y_ghost[v_adjacent_edge, s_edge] for s_edge in edges(s_network_dir)) <= 3)
    end
    =#

end







function update_pricer_sn_decompo(vn_decompo, pricer, dual_costs)

    model = pricer.model

    sn_subgraph = pricer.sn_subgraph
    sub_s_network_dir = pricer.sub_instance.s_network_dir
    original_s_network = pricer.original_instance.s_network
    sub_s_network = pricer.sub_instance.s_network
    v_subgraph = pricer.vn_subgraph

    ### Objective
    placement_cost = @expression(model, 
        sum( ( sub_s_network_dir[s_node][:cost] - dual_costs.capacity_s_node[sn_subgraph.nodes_of_main_graph[s_node]] )  * model[:x][v_node, s_node] 
            for v_node in vertices(v_subgraph.graph) for s_node in vertices(sub_s_network_dir) ))

    routing_cost = @expression(model, sum( 
        ( sub_s_network[src(s_edge), dst(s_edge)][:cost] - dual_costs.capacity_s_edge[get_edge(original_s_network, sn_subgraph.nodes_of_main_graph[src(s_edge)], sn_subgraph.nodes_of_main_graph[dst(s_edge)])] ) 
        * v_subgraph.graph[src(v_edge), dst(v_edge)][:dem] * (model[:y][v_edge, get_edge(sub_s_network_dir, src(s_edge), dst(s_edge))] + model[:y][v_edge, get_edge(sub_s_network_dir, dst(s_edge), src(s_edge))])
                for v_edge in edges(v_subgraph.graph) for s_edge in edges(sub_s_network) ))


            
    # flow conservation
    flow_conservation_cost = AffExpr(0.)

    for s_node in vertices(sub_s_network)
        original_node = sn_subgraph.nodes_of_main_graph[s_node]
        for connecting_edge in vn_decompo.v_edges_master
            if v_subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][v_subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    -dual_costs.flow_conservation[connecting_edge][original_node] , 
                    model[:x][v_node_subgraph, s_node])
            end
            if v_subgraph ∈ keys(vn_decompo.v_nodes_assignment[dst(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[dst(connecting_edge)][v_subgraph]
                add_to_expression!(
                    flow_conservation_cost, 
                    +dual_costs.flow_conservation[connecting_edge][original_node], 
                    model[:x][v_node_subgraph, s_node])
            end
        end
    end


    # departure !
    departure_costs = AffExpr(0.)
    for s_node in vertices(sub_s_network)
        for connecting_edge in vn_decompo.v_edges_master
            if v_subgraph ∈ keys(vn_decompo.v_nodes_assignment[src(connecting_edge)])
                v_node_subgraph = vn_decompo.v_nodes_assignment[src(connecting_edge)][v_subgraph]
                add_to_expression!(
                    departure_costs, 
                    -dual_costs.departure[connecting_edge][sn_subgraph.nodes_of_main_graph[s_node]], 
                    model[:x][v_node_subgraph, s_node])
            end
        end
    end


    @objective(model, Min, 
            -dual_costs.convexity[v_subgraph]
            + placement_cost + routing_cost 
            + flow_conservation_cost 
            + departure_costs);


    return
end




function solve_pricers_sn_decompo(pricer; time_limit = 1000)

    original_s_network_dir = pricer.original_instance.s_network_dir
    sn_subgraph = pricer.sn_subgraph
    vn_subgraph = pricer.vn_subgraph
    v_network = vn_subgraph.graph
    s_network_dir = pricer.sub_instance.s_network_dir
    model = pricer.model

    set_time_limit_sec(model, time_limit)
    optimize!(model)
    status = termination_status(model)

    if status == MOI.FEASIBLE_POINT || status == MOI.OPTIMAL
        # Get the solution
        x_values = value.(model[:x])
        y_values = value.(model[:y])
        true_cost = 0.

        node_placement = []
        for v_node in vertices(v_network)
            for s_node in vertices(s_network_dir)
                if x_values[v_node, s_node] > 0.99
                    real_s_node = sn_subgraph.nodes_of_main_graph[s_node]
                    append!(node_placement, real_s_node)
                    true_cost += original_s_network_dir[real_s_node][:cost]
                end
            end
        end


        edge_routing = Dict()
        for v_edge in edges(v_network)
            if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
                edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
            end
            used_edges = []
            for s_edge in edges(s_network_dir)
                if y_values[v_edge, s_edge] > 0.99
                    real_s_edge = get_edge(original_s_network_dir, sn_subgraph.nodes_of_main_graph[src(s_edge)], sn_subgraph.nodes_of_main_graph[dst(s_edge)])
                    push!(used_edges, real_s_edge)
                    true_cost += original_s_network_dir[src(real_s_edge), dst(real_s_edge)][:cost]
                end
            end
            edge_routing[v_edge] = order_path(original_s_network_dir, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
        end
        mapping = Mapping(v_network, original_s_network_dir, node_placement, edge_routing)
        reduced_cost = objective_value(model)
        return mapping, true_cost, reduced_cost
    else
        return nothing, 9999999, 999999
    end


end



