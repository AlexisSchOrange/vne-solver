using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf


#general
includet("../../utils/import_utils.jl")

# utils colge
includet("utils/utils-subgraphdecompo.jl")
includet("utils/partition-graph.jl")
includet("utils/checkers.jl")


# init
includet("init/init-paving-simpler.jl")

# pricers
includet("pricers/pricer-subsn.jl")
includet("pricers/pricer-exact.jl")
includet("pricers/mepso-pricer.jl")
includet("pricers/mepso-pricer-subsn.jl")





function crazy_bound(instance)
    
    
    println("Starting...")

    time_master = 0
    nb_iter=0
    lower_bound = 0
    nb_columns=0

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    time_beginning = time()


    # ======= SETTING UP THE DECOMPOSITION ======= #
    nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
    v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)

    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    
    println("Decomposition set: ")
        println("For $v_network, there is $(length(vn_decompo.subgraphs)) subgraphs:")

    for subgraph in vn_decompo.subgraphs
        println("       $(subgraph.graph[][:name]) with $(nv(subgraph.graph)) nodes, $(ne(subgraph.graph)) edges")
    end
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")

    
    vn_subgraphs = vn_decompo.subgraphs

    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")

    # ====== STEP 1 : INIT ======= #
    nb_columns_to_add_init = 200

    println("Paving time...")
    time_0 = time()

    # Get substrate subgraphs - 1
    size_max_v_subgraph = maximum(nv(v_subgraph.graph) for v_subgraph in vn_decompo.subgraphs)
    nb_substrate_subgraphs = floor(Int, nv(s_network) / (size_max_v_subgraph*1.5))

    clusters = partition_graph(s_network.graph, nb_substrate_subgraphs; max_umbalance = 1.3)
    sn_subgraphs = []
    for (i_subgraph, cluster) in enumerate(clusters)
        print("Cluster $i_subgraph has $(length(cluster)) nodes ")
        induced_subg = my_induced_subgraph(s_network, cluster, "sub_sn_$i_subgraph")
        push!(sn_subgraphs,Subgraph(induced_subg, cluster))
    end

    sub_mappings = find_submappings(instance, vn_decompo, sn_subgraphs, solver="mepso", nb_columns=nb_columns_to_add_init)
    println("Mappings gotten! In just $(time() - time_0)")
    for v_subgraph in vn_decompo.subgraphs
        for mapping in sub_mappings[v_subgraph]
            add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
            nb_columns+=1
        end
    end
    print("Submappings added - 1...")



    # ======= ITERATION OF MASTER PROBLEM ======= #
    optimize!(model)
    time_master +=  solve_time(model)

    cg_value = objective_value(model)
    time_overall = time() - time_beginning
    @printf("Iter %2d  RMP value: %10.3f  %5d column    time: %5.2fs  \n",
        nb_iter, cg_value, nb_columns, time_overall
    )


    # ====== STEP 2a) HEURISTIC PAVING PRICING WITH PSO ====== #

    # Get substrate subgraphs 
    println("TIME TO PAVE WITH MY MEPSO PRICER")
    size_max_v_subgraph = maximum(nv(v_subgraph.graph) for v_subgraph in vn_decompo.subgraphs)
    nb_substrate_subgraphs = floor(Int, nv(s_network) / (size_max_v_subgraph*2))

    clusters = partition_graph(s_network.graph, nb_substrate_subgraphs; max_umbalance = 1.3)
    sn_subgraphs = []
    for (i_subgraph, cluster) in enumerate(clusters)
        print("Cluster $i_subgraph has $(length(cluster)) nodes ")
        induced_subg = my_induced_subgraph(s_network, cluster, "sub_sn_$i_subgraph")
        push!(sn_subgraphs,Subgraph(induced_subg, cluster))
    end


    nb_columns_to_put_mepso = 1000
    nb_columns_pricing_mepso = 0
    while nb_columns_pricing_mepso < nb_columns_to_put_mepso

        nb_iter += 1
        #dual_costs = get_empty_duals(instance, vn_decompo)
        dual_costs = get_duals(instance, vn_decompo, master_problem)

        # Associate subvn to a random subsn
        used_sub_s_network = []
        assignment_virtual_substrate_subgraphs = Dict()
        for v_subgraph in vn_subgraphs
            found = false
            while !found
                i_subgraph = rand(1:nb_substrate_subgraphs)
                if i_subgraph ∉ used_sub_s_network
                    push!(used_sub_s_network, i_subgraph)
                    found = true
                    assignment_virtual_substrate_subgraphs[v_subgraph] = sn_subgraphs[i_subgraph]
                end
            end
        end


        for v_subgraph in vn_subgraphs
            s_subgraph = assignment_virtual_substrate_subgraphs[v_subgraph]
            sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)

            sub_mapping, reduced_cost = solve_pricer_mepso_sub_sn(v_subgraph, s_subgraph, dual_costs, vn_decompo, sub_instance, instance, print_things=false)

            #println("Found a columns with reduced cost $reduced_cost")

            if !isnothing(sub_mapping) && reduced_cost < -0.01
                #println("Adding it!")
                #println("Mapping: $sub_mapping")
                #println("Routing: $(sub_mapping.edge_routing)")
                true_cost = get_cost_placement(sub_mapping) + get_cost_routing(sub_mapping)
            
                #println("Cost : $(true_cost)")

                add_column(master_problem, instance, v_subgraph, sub_mapping, true_cost )
                nb_columns += 1
                nb_columns_pricing_mepso += 1
            end
        end


        # ======= ITERATION OF MASTER PROBLEM ======= #
        optimize!(model)
        time_master +=  solve_time(model)

        cg_value = objective_value(model)
        time_overall = time() - time_beginning
        @printf("Iter %2d  RMP value: %10.3f  %5d column    time: %5.2fs  \n",
            nb_iter, cg_value, nb_columns, time_overall
        )

    end



    #= ====== STEP 2 : HEURISTICS PRICING ======= # ITS JUST NOT WORKING BRAH
    nb_columns_to_put_mepso = 200
    nb_columns_pricing_mepso = 0
    while nb_columns_pricing_mepso < nb_columns_to_put_mepso

        nb_iter += 1
        #dual_costs = get_empty_duals(instance, vn_decompo)
        dual_costs = get_duals(instance, vn_decompo, master_problem)


        for v_subgraph in vn_subgraphs

            sub_instance = Instance(v_subgraph.graph, s_network)

            sub_mapping, reduced_cost = solve_pricer_mepso(v_subgraph, dual_costs, vn_decompo, sub_instance; nb_particle=50, time_max=0.5, print_things=false)
            println("Found a columns with reduced cost $reduced_cost")

            if !isnothing(sub_mapping) && reduced_cost < -0.1
                #println("Adding it!")
                #println("Mapping: $sub_mapping")
                #println("Routing: $(sub_mapping.edge_routing)")
                true_cost = get_cost_placement(sub_mapping) + get_cost_routing(sub_mapping)
            
                #println("Cost : $(true_cost)")

                add_column(master_problem, instance, v_subgraph, sub_mapping, true_cost )
                nb_columns += 1
                nb_columns_pricing_mepso += 1
            end
        end


        # ======= ITERATION OF MASTER PROBLEM ======= #
        optimize!(model)
        time_master +=  solve_time(model)

        cg_value = objective_value(model)
        time_overall = time() - time_beginning
        @printf("Iter %2d  RMP value: %10.3f  %5d column    time: %5.2fs  \n",
            nb_iter, cg_value, nb_columns, time_overall
        )



    end
    =#

    # ======= FIRST ITERATION OF MASTER PROBLEM ======= #
    optimize!(model)
    time_master +=  solve_time(model)

    cg_value = objective_value(model)
    time_overall = time() - time_beginning
    @printf("Iter %2d  RMP value: %10.3f  %5d column    time: %5.2fs  \n",
        nb_iter, cg_value, nb_columns, time_overall
    )


    # ====== STEP 2: smaller pricers - paving the network
    println("------- Part 2: Reduced pricers")

    nb_columns_to_put = 1000
    nb_columns_smallpricer = 0

    nb_substrate_subgraph = floor(Int, nv(s_network) / 15)  
    nb_nodes_subgraph = 32
    sn_decompo_clusters = get_sn_decompo(s_network, nb_substrate_subgraph, nb_nodes_subgraph)
    println("We have $nb_substrate_subgraph sub-substrate, with at least $nb_nodes_subgraph capacited nodes")

    pricers_sn_decompo = OrderedDict()
    sub_pricers_last_values = OrderedDict()
    for vn_subgraph in vn_decompo.subgraphs
        pricers_sn_decompo[vn_subgraph] = set_up_pricer_sn_decompo(instance, vn_subgraph, sn_decompo_clusters)
        for pricer in pricers_sn_decompo[vn_subgraph]
            sub_pricers_last_values[pricer] = -99999.
        end
    end

    keep_on = true
    reason = "I don't know"
    while keep_on
        nb_iter += 1

        # ---- pricers part

        dual_costs = get_duals(instance, vn_decompo, master_problem)
        
        for key in keys(sub_pricers_last_values)
            sub_pricers_last_values[key] = sub_pricers_last_values[key]*1.1
        end
        
        sorted_subpb = sort(collect(sub_pricers_last_values), by=x->x[2])
        average_obj = 0
        nb_pricer_to_do = min(5, length(keys(sub_pricers_last_values)))
        for couple in sorted_subpb[1:nb_pricer_to_do]


            pricer_sub_sn = couple[1]


            update_pricer_sn_decompo(vn_decompo, pricer_sub_sn, dual_costs)
            column, true_cost, reduced_cost = solve_pricers_sn_decompo(pricer_sub_sn, time_limit=50)

            if column !== nothing && reduced_cost < -0.001
                add_column(master_problem, instance, pricer_sub_sn.vn_subgraph, column, true_cost)
                nb_columns += 1
                nb_columns_smallpricer += 1
            end

            sub_pricers_last_values[pricer_sub_sn] = reduced_cost

            average_obj += (reduced_cost/nb_pricer_to_do)
        end


        # ---- master problem part

        optimize!(model)
        time_master +=  solve_time(model)

        cg_value = objective_value(model)

        time_overall = time()-time_beginning

        @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs  average reduced cost: %10.3f \n",
                    nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_obj)
    


        
        # ----- useful things

        if nb_columns_smallpricer > nb_columns_to_put
            keep_on = false
        end

    end
    println("\n Step 2 finished, reason: $reason.")
    


    
    # ====== STEP 3: full pricers
    println("\n------- Solving method: Exact pricers")

    pricers_full = Dict()
    for subgraph in vn_decompo.subgraphs
        pricers_full[subgraph] = set_up_pricer(instance, subgraph)
    end
    
    keep_on = true
    reason = "I don't know"
    while keep_on
        nb_iter += 1

        # ---- Pricers things
        dual_costs = get_duals(instance, vn_decompo, master_problem)

        has_found_new_column = false
        sum_pricers_values = 0

        for vn_subgraph in vn_decompo.subgraphs
            pricer = pricers_full[vn_subgraph]

            time_limit_pricer = 500


            sub_mapping, true_cost, reduced_cost = update_solve_pricer(instance, vn_decompo, pricer, dual_costs; time_limit = time_limit_pricer)

            
            if (!isnothing(sub_mapping)) && reduced_cost < -0.0001
                has_found_new_column = true
                add_column(master_problem, instance, vn_subgraph, sub_mapping, true_cost)
                nb_columns += 1
            end
            
            if isnothing(sub_mapping)
                println("Pricer with no solution found, stopping the CG...")
                reason="pricer-unfeasible"
                has_found_new_column = false
                break
            end

            sum_pricers_values += reduced_cost

        end

        current_lower_bound = cg_value + sum_pricers_values
        if current_lower_bound > lower_bound
            lower_bound = current_lower_bound
        end


        # ----- Master problem stuff

        optimize!(model)
        time_master +=  solve_time(model)

        cg_value = objective_value(model)

        time_overall = time()-time_beginning
        average_obj = sum_pricers_values/length(vn_decompo.subgraphs)

        @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs  average reduced cost: %10.3f \n",
            nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_obj)



        if !has_found_new_column
            keep_on = false
            reason="no improving columns"
        end
        keep_on = false


    end

        


    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: idk, time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(cg_value; digits=3))")
    println("====================================================\n")



    
    
    # ======= GETTING A SOLUTION ======= #
    time_cg_heuristic = 600
    value_cg_heuristic, cg_heuristic_solution = basic_heuristic(instance, vn_decompo, master_problem, time_cg_heuristic)


    result = Dict()
    result["solving_time"] = time() - time_beginning
    result["mapping_cost"] = value_cg_heuristic

    return result
end






