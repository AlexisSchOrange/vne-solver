using JuMP, CPLEX, Gurobi


# Implementation of : Virtual Network Embedding Based on the Degree and Clustering Coefficient Information


function solve_VNE_DCC(instance)


    v_network = instance.v_network
    s_network = instance.s_network
    s_network = instance.s_network



    # --- STEP 1 --- RANK THE NODES

    # Virtual node ranking
    v_nodes_centrality = closeness_centrality(v_network) # here the article does a bit different (but not really useful)
    most_central_nodes = findmax(v_nodes_centrality)[2]
    distances = desopo_pape_shortest_paths(instance.v_network, most_central_nodes).dists 
    v_nodes_scores = [ v_nodes_centrality[v_node] - distances[v_node] * 10 for v_node in vertices(v_network)]
    v_nodes_ordered = sortperm(v_nodes_scores; rev=true)
    
    # Substrate node ranking
    s_nodes_capacities_around = [s_network[s_node][:cap] * sum(s_network[src(s_edge), dst(s_edge)][:cap] 
                        for s_edge in get_out_edges(s_network, s_node) ∪ get_in_edges(s_network, s_node))
                            for s_node in vertices(s_network)
    ]
    
    s_node_centrality = closeness_centrality(s_network)

    s_nodes_scores = [s_node_centrality[s_node] * s_nodes_capacities_around[s_node] for s_node in vertices(s_network)]
    s_nodes_ordered = sortperm(s_nodes_scores; rev=true)




    # --- STEP 2 --- PLACE THE VIRTUAL NODES
    placement = Dict()
    v_node_placed = []
    s_node_used = []
    
    # init: we place the most central v nodes on the most central s nodes
    placement[v_nodes_ordered[1]] = s_nodes_ordered[1]
    push!(v_node_placed, v_nodes_ordered[1])
    push!(s_node_used, s_nodes_ordered[1])

    println("I shall place $(v_nodes_ordered[1]) on $(s_nodes_ordered[1])")


    for v_node in v_nodes_ordered[2:end]
        #println("Alright lets do $v_node")
        # find the sum of the distance between a s_node and the placements of neighbors of v_node that have been placed
        s_nodes_of_neighbors = []
        for neighbor in neighbors(v_network, v_node)
            if neighbor ∈ v_node_placed
                push!(s_nodes_of_neighbors, placement[neighbor])
            end
        end
        #println("Neighbors: $s_nodes_of_neighbors")

        current_distances_s_nodes = [sum(desopo_pape_shortest_paths(instance.s_network, placement_neighbor).dists[s_node] for placement_neighbor in s_nodes_of_neighbors) for s_node in vertices(s_network)]
        current_s_node_scores = [s_nodes_scores[s_node] - current_distances_s_nodes[s_node] * 10 for s_node in vertices(s_network)]


        #println("Distance score : $current_distances_s_nodes")
        #println("Overall score : $current_s_node_scores")
        current_s_nodes_ranked = sortperm(current_s_node_scores; rev=true)
        keep_on = true
        idx_s_node = 1
        while keep_on
            s_node = current_s_nodes_ranked[idx_s_node]
            #println("I gotta carry on")
            if s_node ∉ s_node_used && s_network[s_node][:cap] >= 1 
                placement[v_node] = s_node
                push!(v_node_placed, v_node)
                push!(s_node_used, s_node)
                keep_on = false
                println("I shall place $v_node on $s_node")
            end
            idx_s_node += 1
        end
    end

    println("Placement of nodes: $placement")

    # STEP 3 --- ROUTE THE VIRTUAL EDGES


    model = Model(CPLEX.Optimizer)
    set_up_problem_placement_restrict(instance, model, placement)
    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        return -999
    end

    return objective_value(model)



end




function set_up_problem_placement_restrict(instance, model, placement_restriction)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);

    

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);




    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
        @constraint(model, sum(x[v_node, s_node] for s_node in placement_restriction[v_node]) == 1)
    end

    # capacity
    for s_node in vertices(s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= s_network[s_node][:cap])
    end


    ## Edges 
    
    # edge capacity (undirected version !)
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( v_network[src(v_edge), dst(v_edge)][:dem] * (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_edge in edges(v_network)) 
            <= 
            instance.s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end

    # Flow departure constraint
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    
end


