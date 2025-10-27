
#using Revise, JuMP, CPLEX, Gurobi
using Revise, JuMP, CPLEX
using Graphs, GraphsFlows
using SparseArrays


includet("../utils/import_utils.jl")




# ========== CLASSICAL STUFF
function set_up_problem_ff_plus(instance, model)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);

    

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);




    ###=========== Constraints

    ##---- Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(instance.s_network)) == 1)
    end

    # one to one : one virtual node per substrate node
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end

    # node capacity : NOT USELESS AHHHHHHHHh
    for s_node in vertices(instance.s_network)
        @constraint(model, sum( x[v_node, s_node] for v_node in vertices(v_network)) <= sum(s_network[s_node][:cap]))
    end


    ##--- Edges 
    
    # edge capacity (undirected version !)
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
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

    
    ## Departure constraints    
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    
    
end



function solve_branch_and_cut(instance; time_solver = 100, alpha_acceptance=0.05)
    
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    time_start = time()
    model = Model(CPLEX.Optimizer)
    set_up_problem_ff_plus(instance, model)

    set_time_limit_sec(model, time_solver)
    #set_optimizer_attribute(model, "Threads", 1)
    #set_optimizmer_attribute(model, "CPXPARAM_MIP_Strategy_Search", 1)
    #set_optimizer_attribute(model, "CPX_PARAM_REPEATPRESOLVE", 0)
    #
    #set_optimizer_attribute(model, "CPXPARAM_Parallel", -1)

    # setting things for the branch and cut
    n_s = nv(s_network_dir)
    augmented_network = copy(s_network_dir.graph)
    add_vertex!(augmented_network) # n_s + 1: t(\ebar)
    for i in 1:n_s
        if s_network[i][:cap] > 0
            add_edge!(augmented_network, i, n_s+1)
        end
    end

    nb_cuts_overall = 0
    nb_call = 0
    time_seperation_overall = 0
    time_min_cut = 0
    time_get_values = 0
    
    matrix_flows = zeros(Float32, n_s+1, n_s+1)

    function find_cutset_cuts(cb_data)
        
        nb_call += 1

        if alpha_acceptance > 0.99 # Ahah don't worry (this is ugly but it's just for the tests)
            return
        end

        time_start_separation = time()
        nb_cuts_current = 0
        time_beg_getting_values = time()
        x_values = callback_value.(cb_data, model[:x])
        y_values = callback_value.(cb_data, model[:y])
        time_get_values += time()-time_beg_getting_values

        for v_edge in edges(v_network)

            matrix_flows .= 0.0  # resets in place
            
            for s_edge in edges(s_network_dir)
                matrix_flows[src(s_edge), dst(s_edge)] = y_values[v_edge, s_edge]
            end

            for s_node in vertices(s_network_dir)
                x_values[src(v_edge), s_node] < alpha_acceptance * 1.1 && continue # reduces the number of separations
                x_values[src(v_edge), s_node] > 0.9 && continue # prevent bug in the mincut!

                # Finish flow matrix for that node - remember that no flow between u and t(\ebar) here
                for other_s_node in vertices(s_network_dir)
                    if s_node == other_s_node
                        matrix_flows[other_s_node, n_s+1] = 0.
                    else 
                        matrix_flows[other_s_node, n_s+1] = x_values[dst(v_edge), other_s_node]
                    end
                end
                
                #println("Matrix flow : $matrix_flows")
                # min cut between s_node and n_s+1
                time_beg_cut = time()
                #print("uh")
                #println("Well matrix flow for edge $v_edge and node $s_node:")
                
                (part1, part2, flow) = GraphsFlows.mincut(augmented_network, s_node, n_s+1, matrix_flows, DinicAlgorithm())
                #print("oh")
                time_min_cut += time()-time_beg_cut
                if flow < x_values[src(v_edge), s_node] - alpha_acceptance
                    if (n_s+1) ∈ part1
                        println("WOWOWOW NS+1 IS IN PART1?? $part1 and $part2 wthhh")
                        println("I mean, I'm doing a cut between $s_node and $(n_s+1)... Is that algorithm stupid?")
                        break
                    end
                    if (s_node) ∈ part2
                        println("WOWOWOW $s_node IS IN PART2?? $part1 and $part2 wthhh")
                        println("I mean, I'm doing a cut between $s_node and $(n_s+1)... Is that algorithm stupid?")
                        break
                    end
                    #=
                    for node in part1
                        if s_network[node][:cap] > 0
                            println("Value: $(x_values[dst(v_edge), node])")
                        end
                    end
                    =#
                    cut_s_edges = get_edges_from_S1_to_S2(s_network_dir, part1, part2)
                    cut = @build_constraint(sum(model[:y][v_edge, s_edge] for s_edge in cut_s_edges) + sum( model[:x][dst(v_edge), s_other_node] for s_other_node in part1)
                        >=
                        model[:x][src(v_edge), s_node] + model[:x][dst(v_edge), s_node]
                    )
                    MOI.submit(model, MOI.UserCut(cb_data), cut)

                    nb_cuts_current += 1
                    nb_cuts_overall += 1
                end
            end


        end

        time_seperation_overall += time() - time_start_separation
        #println("Added $nb_cuts_current this time.")
        #println("Added $nb_cuts_overall so far, in $time_overall.")

        return

    end


    MOI.set(model, MOI.UserCutCallback(), find_cutset_cuts) 

    
    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        return ( sol_value= -1,
                lower_bound = -1.,
                gap = -1.,
                node_count = node_count(model),
                time_solving = (time() - time_start),
                time_seperation = time_seperation_overall,
                nb_cuts = nb_cuts_overall
        )
    end




    # Get the solution

    
    

    println("Find the solution $(objective_value(model)) in $(time() - time_start)")
    println("Stats branch and cut: $time_seperation_overall s overall, $time_min_cut on mincut, $time_get_values getting values,  $nb_cuts_overall cuts found")


    return    ( sol_value= objective_value(model),
                lower_bound = objective_bound(model),
                gap = relative_gap(model),
                node_count = node_count(model),
                time_solving = (time() - time_start),
                time_seperation = time_seperation_overall,
                nb_cuts = nb_cuts_overall
    )
end





function get_edges_from_S1_to_S2(graph, S1, S2)

    cut_edges = []
    for edge in edges(graph)
        if src(edge) ∈ S1 && dst(edge) ∈ S2
            push!(cut_edges, edge)
        end
    end

    return cut_edges
end








function solve_simple_cplex(instance; time_solver = 100)
    
    v_network = instance.v_network
    s_network_dir = instance.s_network_dir

    time_start = time()
    model = Model(CPLEX.Optimizer)
    set_up_problem_ff_plus(instance, model)

    set_time_limit_sec(model, time_solver)

    # On
    set_optimizer_attribute(model, "CPXPARAM_Threads", 1)    
    #set_optimizer_attribute(model, "CPXPARAM_MIP_Strategy_Search", 1)

    #set_optimizer_attribute(model, "CPX_PARAM_REPEATPRESOLVE", 0)

    optimize!(model)

    status = primal_status(model)
    if status != MOI.FEASIBLE_POINT
        println("Infeasible or unfinished: $status")
        return ( sol_value= -1,
            lower_bound = -1.,
            gap = -1.,
            node_count = node_count(model),
            time_solving = (time() - time_start)
        )
    end



    # Get the solution



    return ( sol_value= objective_value(model),
        lower_bound = objective_bound(model),
        gap = relative_gap(model),
        node_count = node_count(model),
        time_solving = (time() - time_start),
        nb_call=nb_call
    )
end







# ========== CLASSICAL STUFF


function solve_linear(instance)

    time_beginning = time()

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir

    time_start = time()
    model = Model(CPLEX.Optimizer)
    set_up_problem_ff_plus(instance, model)
    relax_integrality(model)
    set_silent(model)

    alpha = 0.001

    optimize!(model)
    value_ff = round(objective_value(model); digits=3)

    n_s = nv(s_network_dir)
    augmented_network = deepcopy(s_network_dir.graph)
    add_vertex!(augmented_network) # n_s + 1: t(\ebar)
    for i in 1:n_s
        add_edge!(augmented_network, i, n_s+1)
    end

    nb_cuts_overall = 0

    keep_going = true
    iter = 1
    while keep_going

        optimize!(model)

        status = primal_status(model)
        if status != MOI.FEASIBLE_POINT
            println("It's unfeasible!!")
            return (value_ff, 10e9)
        end
        
        println("Iter $iter, value: $(objective_value(model)), nb cuts: $nb_cuts_overall")

        keep_going = false

        x_values = value.(model[:x])
        y_values = value.(model[:y])

        nb_new_cuts = 0
        
        for v_edge in edges(v_network)

            matrix_flows = zeros(n_s+1, n_s+1)
 
            for s_edge in edges(s_network_dir)
                matrix_flows[src(s_edge), dst(s_edge)] = y_values[v_edge, s_edge]
            end

            for s_node in vertices(s_network_dir)
                x_values[src(v_edge), s_node] < 0.001 && continue
                # Finish flow matrix for that node - remember that no flow between u and t(\ebar) here
                for other_s_node in vertices(s_network_dir)
                    if s_node == other_s_node
                        matrix_flows[other_s_node, n_s+1] = 0.
                    else 
                        matrix_flows[other_s_node, n_s+1] = x_values[dst(v_edge), other_s_node]
                    end
                end

                (part1, part2, flow) = GraphsFlows.mincut(augmented_network, s_node, n_s+1, matrix_flows, PushRelabelAlgorithm())

                if flow < x_values[src(v_edge), s_node] - alpha
                    #println("Damn! $part1 to $part2: i got $flow, imma win $(x_values[src(v_edge), s_node] - flow) for $v_edge")
                    cut_s_edges = get_edges_from_S1_to_S2(s_network_dir, part1, part2)
                    @constraint(model, sum(model[:y][v_edge, s_edge] for s_edge in cut_s_edges) + sum( model[:x][dst(v_edge), s_other_node] for s_other_node in part1)
                        >=
                        model[:x][src(v_edge), s_node] + model[:x][dst(v_edge), s_node]
                    )

                    nb_cuts_overall += 1
                    keep_going = true
                    nb_new_cuts += 1
                end
            end


        end

        println("I generated $nb_new_cuts new cuts!")
        iter += 1

    end


    #=
    x_values = value.(model[:x])
    y_values = value.(model[:y])


    node_placement = []
    for v_node in vertices(v_network)
        for s_node in vertices(s_network_dir)
            if x_values[v_node, s_node] > 0.01
                push!(node_placement, s_node)
            end
        end
    end

    for v_edge in edges(v_network)
        for s_edge in edges(s_network_dir)
            if y_values[v_edge, s_edge] > 0.01
                if (s_network_dir[src(s_edge)][:cap] == 0) ||  (s_network_dir[dst(s_edge)][:cap] == 0)
                    println("$v_edge on $s_edge : $(y_values[v_edge, s_edge])...")
                end
                if (s_network_dir[src(s_edge)][:cap] == 0) &&  (s_network_dir[dst(s_edge)][:cap] == 0)
                    println("WTF: $v_edge on $s_edge : $(y_values[v_edge, s_edge])...")
                end
            end
        end
    end
    =#

    value_cutting_plane = round(objective_value(model); digits=3)
    println("\n\n\n FINISHED, obtained $value_cutting_plane in $(time() - time_beginning)s")

    return (value_ff, value_cutting_plane)
end













