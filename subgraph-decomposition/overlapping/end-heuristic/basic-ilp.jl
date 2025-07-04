using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX




function basic_heuristic(instance, vn_decompo, master_problem, time_max)


    model = Model(CPLEX.Optimizer)
    set_up_problem(instance, vn_decompo, model)
    set_time_limit_sec(model, time_max)
    
    for subgraph in vn_decompo.subgraphs
        for column in master_problem.columns[subgraph]
            add_column_ip(model, vn_decompo, instance, subgraph, column)
        end
    end

    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        return -999
    end
    println("Optimal solution : $(objective_value(model))")
    return objective_value(model)
end





function set_up_problem(instance, vn_decompo, model)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    set_silent(model)

    
    ### Variables    
    @variable(model, 0. <= y[
        v_edge in vn_decompo.v_edges_master, 
        s_edge in edges(s_network_dir)] <=1. );
    

    columns = Dict()
    for subgraph in vn_decompo.subgraphs
        columns[subgraph] = []
    end
    
    
    ### Objective
    master_routing_costs = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
        for v_edge in vn_decompo.v_edges_master for s_edge in edges(s_network_dir) ))
    
    @objective(model, Min, master_routing_costs);

    ### Constraints

    # convexity constraints
    @constraint(
        model, 
        mapping_selec[subgraph in vn_decompo.subgraphs],
        0 >= 1
    );


    # overlapping constraint
    overlapping_nodes = keys(vn_decompo.overlapping_nodes)
    @constraint(
        model,
        overlapping[s_node in vertices(s_network), v_node in overlapping_nodes, i_subgraph in 1:(length(vn_decompo.overlapping_nodes[v_node]))],
        0 <= 0
    );


    # capacity of substrate nodes
    @constraint(
        model,
        capacity_s_node[s_node in vertices(s_network)],
        0 <= s_network[s_node][:cap]
    );

    

    # capacity of substrate edges
    # undirected, so both ways !
    @constraint(
        model,
        capacity_s_edge[s_edge in edges(s_network)],
        sum( (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] +  y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))] )
            for v_edge in vn_decompo.v_edges_master)
        + 0
        <= s_network[src(s_edge), dst(s_edge)][:cap]
    );


    # flow conservation constraints
    @constraint(
        model,
        flow_conservation[v_edge in vn_decompo.v_edges_master, s_node in vertices(s_network)],
        0 == 
        sum( y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node))
        - sum( y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
    );


    # Departure constraints
    @constraint(
        model, 
        departure[v_edge in vn_decompo.v_edges_master, s_node in vertices(s_network)],
        0 
        <=
        sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node))
    )
        
    


end



function add_column_ip(model, vn_decompo, instance, subgraph, column)


    lambda = @variable(model, binary=true);
    set_objective_coefficient(model, lambda, column.cost)
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir



    # convexity
    set_normalized_coefficient(model[:mapping_selec][subgraph], lambda, 1)


    # overlapping - this is quite disgusting, but I don't see a way to do it easily...
    for v_node_overlapping in keys(vn_decompo.overlapping_nodes)
        if subgraph in vn_decompo.overlapping_nodes[v_node_overlapping]
            v_node_in_subgraph = vn_decompo.v_nodes_assignment[v_node_overlapping][subgraph]
            s_node_placement = column.mapping.node_placement[v_node_in_subgraph]

            for (i_subgraph_overlapping, subgraph_overlapping) in enumerate(vn_decompo.overlapping_nodes[v_node_overlapping])

                if subgraph_overlapping == subgraph
                    # rhs
                    set_normalized_coefficient(model[:overlapping][s_node_placement, v_node_overlapping, i_subgraph_overlapping], lambda, -1)

                    #lhs
                    if i_subgraph_overlapping == length(vn_decompo.v_nodes_assignment[v_node_overlapping])
                        set_normalized_coefficient(model[:overlapping][s_node_placement, v_node_overlapping, 1], lambda, 1)
                    else
                        set_normalized_coefficient(model[:overlapping][s_node_placement, v_node_overlapping, i_subgraph_overlapping+1], lambda, 1)
                    end
                end
            end
        end
    end


    # capacities on nodes
    for s_node in vertices(s_network)
        usage = 0
        for v_node in vertices(subgraph.graph)
            if column.mapping.node_placement[v_node] == s_node
                original_v_node = subgraph.nodes_of_main_graph[v_node]
                if original_v_node in keys(vn_decompo.overlapping_nodes)
                    usage += 1/length(vn_decompo.overlapping_nodes[original_v_node])
                else
                    usage += 1
                end
            end
        end
        set_normalized_coefficient(model[:capacity_s_node][s_node], lambda, usage)
    end


    # capacities on edges (undirected case!)
    for s_edge in edges(s_network)
        usage = 0
        for v_edge in edges(subgraph.graph)
            s_edge_one = get_edge(s_network_dir, src(s_edge), dst(s_edge))
            if s_edge_one in column.mapping.edge_routing[v_edge].edges
                usage += 1
            end
            s_edge_two = get_edge(s_network_dir, dst(s_edge), src(s_edge))
            if s_edge_two in column.mapping.edge_routing[v_edge].edges
                usage += 1
            end
        end
        set_normalized_coefficient(model[:capacity_s_edge][s_edge], lambda, usage)
    end

    
    # flow conservation constraints (for vnodes that are ends of a cut edge)
    for v_edge in vn_decompo.v_edges_master

        if subgraph in keys(vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]
            if src(v_edge) in keys(vn_decompo.overlapping_nodes)
                coeff = 1/(length(vn_decompo.overlapping_nodes[src(v_edge)]))
            else
                coeff = 1
            end
            set_normalized_coefficient(
                model[:flow_conservation][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                coeff)

        end

        if subgraph in keys(vn_decompo.v_nodes_assignment[dst(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[dst(v_edge)][subgraph]
            if dst(v_edge) in keys(vn_decompo.overlapping_nodes)
                coeff = 1/(length(vn_decompo.overlapping_nodes[dst(v_edge)]))
            else
                coeff = 1
            end
            set_normalized_coefficient(
                model[:flow_conservation][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                -coeff )

        end


    end
    

    # departure constraints (for virtual cut edges)
    for v_edge in vn_decompo.v_edges_master

        if subgraph in keys(vn_decompo.v_nodes_assignment[src(v_edge)])

            v_node_in_subgraph = vn_decompo.v_nodes_assignment[src(v_edge)][subgraph]

            if src(v_edge) in keys(vn_decompo.overlapping_nodes)
                coeff = 1/(length(vn_decompo.overlapping_nodes[src(v_edge)]))
            else
                coeff = 1
            end

            set_normalized_coefficient(
                model[:departure][v_edge, column.mapping.node_placement[v_node_in_subgraph]], 
                lambda, 
                coeff )

        end
    end

end




