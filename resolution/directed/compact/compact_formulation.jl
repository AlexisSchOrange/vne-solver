using Revise, JuMP, CPLEX, Gurobi
includet("../../../utils/import_utils.jl")
includet("solution_compact_fractional.jl")



struct Compact_Formulation
    model
    x_variables
    y_variables
end

function set_up_compact_model(instance, one_to_one = false, departure_cst = false, symmetric = false)
    print("Constructing compact model... ")

    #### Model
    model = Model(CPLEX.Optimizer)

    ### Variables
    x_variables = @variable(model, x[v_network in instance.v_networks, vertices(v_network), vertices(instance.s_network)], binary=true);
    y_variables = @variable(model, y[v_network in instance.v_networks, edges(v_network), edges(instance.s_network)], binary=true);

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_network, v_node, s_node] 
        for v_network in instance.v_networks for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( instance.s_network[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge]
        for v_network in instance.v_networks for v_edge in edges(v_network) for s_edge in edges(instance.s_network) ))
    @objective(model, Min, placement_cost + routing_cost);


    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_network in instance.v_networks
        for v_node in vertices(v_network)
            @constraint(model, sum(x[v_network, v_node, s_node] for s_node in vertices(instance.s_network)) == 1)
        end
    end

    # if one to one : one virtual node per substrate node
    if one_to_one
        for s_node in vertices(instance.s_network)
            for v_network in instance.v_networks
                @constraint(model, sum(x[v_network, v_node, s_node] for v_node in vertices(v_network)) <= 1)
            end
        end
    end



    # node capacity
    for s_node in vertices(instance.s_network)
        @constraint(model, 
            sum( v_network[v_node][:dem] * x[v_network, v_node, s_node] 
                for v_network in instance.v_networks for v_node in vertices(v_network) ) 
            <= 
            instance.s_network[s_node][:cap] )
    end


    ## Edges 
    
    # edge capacity
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( v_network[src(v_edge), dst(v_edge)][:dem] * y[v_network, v_edge, s_edge] 
                for v_network in instance.v_networks for v_edge in edges(v_network)) 
            <= 
            instance.s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(instance.s_network)
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                @constraint(model, 
                    x[v_network, src(v_edge), s_node] - x[v_network, dst(v_edge), s_node] 
                    <=
                    sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) - 
                        sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node))
                )
            end
        end
    end


    ## Additional constraints : Node + Edge
    if one_to_one
        if departure_cst
            for s_node in vertices(instance.s_network)
                for v_network in instance.v_networks
                    for v_node in vertices(v_network)
                        for v_edge in get_out_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_out_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end
                        for v_edge in get_in_edges(v_network, v_node)
                            @constraint(model, sum(y[v_network, v_edge, s_edge] for s_edge in get_in_edges(instance.s_network, s_node)) >= x[v_network, v_node, s_node])
                        end
                    end
                end
            end
        end
    end

    ## Symmetric edges (mostly for undirected)
    if symmetric
        for v_network in instance.v_networks
            for v_edge in edges(v_network)
                for s_edge in edges(instance.s_network)
                    @constraint(model, y[v_network, v_edge, s_edge] == y[v_network, get_edge(v_network, dst(v_edge), src(v_edge)), 
                                                get_edge(instance.s_network, dst(s_edge), src(s_edge))])
                end
            end
        end
    end
    
    println("done.")

    return Compact_Formulation(model, x_variables, y_variables)
end


function get_solution(instance, x_values, y_values)

    mappings = []
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            for s_node in vertices(instance.s_network)
                if x_values[v_network, v_node, s_node] > 0.99
                    append!(node_placement, s_node)
                end
            end
        end

        edge_routing = Dict()
        for v_edge in edges(v_network)
            if node_placement[src(v_edge)] == node_placement[dst(v_edge)]
                edge_routing[v_edge] = Path(src(v_edge), dst(v_edge), [], 0)
            end
            used_edges = []
            for s_edge in edges(instance.s_network)
                if y_values[v_network, v_edge, s_edge] > 0.99
                    push!(used_edges, s_edge)
                end
            end
            edge_routing[v_edge] = order_path(instance.s_network, used_edges, node_placement[src(v_edge)], node_placement[dst(v_edge)]) 
        end
        m = Mapping(v_network, instance.s_network, node_placement, edge_routing)
        push!(mappings, m)
    end

    return mappings
end


function solve_directed_compact_integer(instance, one_to_one = false, departure_cst = false, time_solver = 30, silent = true)

    # Set up the problem
    problem = set_up_compact_model(instance, one_to_one, departure_cst)

    # Solving
    print("Starting solving... ")
    set_time_limit_sec(problem.model, time_solver)
    if silent
        set_silent(problem.model)
    end
    optimize!(problem.model)
    println("done. Solving state: " * string(termination_status(problem.model)) * ", obj value: " * string(objective_value(problem.model)) * ", bound value: " * string(objective_bound(problem.model)))

    # Get the solution
    x_values = value.(problem.model[:x])
    y_values = value.(problem.model[:y])
    mappings = get_solution(instance, x_values, y_values)
    
    return mappings
end



function solve_directed_compact_fractional(instance, one_to_one = false, departure_cst = false, time_solver = 30, silent = true)

    # Set up the problem
    problem = set_up_compact_model(instance, one_to_one, departure_cst)
    relax_integrality(problem.model)

    # Solving
    if silent
        set_silent(problem.model)
    end
    set_time_limit_sec(problem.model, time_solver)
    print("Starting solving... ")
    optimize!(problem.model)

    # Get the solution
    x_values = value.(problem.x_variables)
    y_values = value.(problem.y_variables)
    mappings = []
    for v_network in instance.v_networks
        node_placement = []
        for v_node in vertices(v_network)
            push!(node_placement, [])
            for s_node in vertices(instance.s_network)
                push!(node_placement[v_node], x_values[v_network, v_node, s_node])
            end
        end
        edge_routing = Dict()
        for v_edge in edges(v_network)
            edge_routing[v_edge] = Dict()
            for s_edge in edges(instance.s_network)
                edge_routing[v_edge][s_edge] = y_values[v_network, v_edge, s_edge]
            end
        end
        
        m = MappingCompactFractional(v_network, instance.s_network, node_placement, edge_routing)
        push!(mappings, m)
    end

    return mappings
end

