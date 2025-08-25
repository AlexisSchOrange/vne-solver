using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf


#general
includet("../../utils/import_utils.jl")

# utils colge
includet("utils/master-problem.jl")
includet("utils/graph-decomposition.jl")
includet("utils/partition-graph.jl")

# pricers
includet("pricers/milp-pricer-subsn-routing.jl")
includet("pricers/local-search-pricer-subsn-routing.jl")
includet("pricers/milp-pricer.jl")


function lower_bound(instance)
    
    
    println("Starting...")  

    time_master = 0
    nb_iter=0
    lower_bound = 0
    nb_columns=0
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    # === SOME PARAMETERS === #
    alpha_colge=0.85
    


    # ======= SETTING UP THE DECOMPOSITION ======= #
    nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
    v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)
    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")





    # === GETTING READY === #
    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")
    empty_dual_costs = get_empty_duals(instance, vn_decompo)

    size_max_v_subgraph = maximum(nv(v_subgraph.graph) for v_subgraph in vn_decompo.subgraphs)
    nb_substrate_nodes_capacited = 0
    for s_node in vertices(s_network)
        if get_attribute_node(s_network, s_node, :cap) >= 1
            nb_substrate_nodes_capacited += 1
        end
    end
    ratio_capacited = nb_substrate_nodes_capacited / nv(s_network)
    nb_substrate_subgraphs = maximum([floor(Int, nv(s_network) / (size_max_v_subgraph*1.5/ratio_capacited)), nb_virtual_subgraph])
    clusters = partition_graph(s_network.graph, nb_substrate_subgraphs; max_umbalance = 1.3)
    sn_subgraphs = []
    for (i_subgraph, cluster) in enumerate(clusters)
        induced_subg = my_induced_subgraph(s_network, cluster, "sub_sn_$i_subgraph")
        push!(sn_subgraphs, Subgraph(induced_subg, cluster))
    end
    println("Substrate network decomposition done:")
    print_stuff_subgraphs(s_network, sn_subgraphs)



    # ====== STEP 1 : PAVING ON STRICT PARTITION ======= #
    nb_columns_to_add_init = 200


    println("Paving time...")
    time_0 = time()


    
    result = find_columns(instance, vn_subgraphs, sn_subgraphs, vn_decompo, empty_dual_costs, solver="local-search", nb_iterations=100, nb_columns=nb_columns_to_add_init)
    sub_mappings = result[:mappings]
    for v_subgraph in vn_decompo.subgraphs
        for mapping in sub_mappings[v_subgraph]
            add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
            nb_columns += 1
        end
    end

    println("\nInitialization dong: \n$nb_columns columns gotten, with local-search, in $(time() - time_0)s\n")




    # ======= ITERATION OF MASTER PROBLEM ======= #
    optimize!(model)
    time_master +=  solve_time(model)

    cg_value = objective_value(model)
    time_overall = time() - time_beginning
    @printf("Iter %2d  RMP value: %10.3f  %5d column    time: %5.2fs  \n",
        nb_iter, cg_value, nb_columns, time_overall
    )


    # ====== STEP 2: PAVING WITH HEURISTICS (local search then milp)
    println("------- Part 2: Reduced pricers")

    # need to complete the substrate graphs!
    nb_nodes_subgraph = 2.5 * size_max_v_subgraph
    sn_clusters_2 = complete_clusters(clusters, s_network, nb_nodes_to_have = nb_nodes_subgraph)
    sn_subgraphs_2 = []
    for (i_subgraph, cluster) in enumerate(sn_clusters_2)
        induced_subg = my_induced_subgraph(s_network, cluster, "sub_sn_$i_subgraph")
        push!(sn_subgraphs_2, Subgraph(induced_subg, cluster))
    end
    println("Adding nodes to the previous substrate decomposition done:")
    print_stuff_subgraphs(s_network, sn_subgraphs_2)


    nb_columns_local_search = 400
    nb_columns_milp = 600


    #!!! For now I won't do "using the most negative reduced costs subproblems first". 
    # Because I'm lazy, and also, because we have no proofs it works.
    # Showing it works would be nice.

    used_dual_costs = empty_dual_costs

    keep_on = true
    reason = "I don't know"
    pricer = "local-search"

    while keep_on
        nb_iter += 1

        # ---- pricers part
        old_dual_costs = used_dual_costs
        current_dual_costs = get_duals(instance, vn_decompo, master_problem)
        used_dual_costs = average_dual_costs(instance, vn_decompo, old_dual_costs, current_dual_costs, alpha=alpha_colge)
        
        result = find_columns(  instance, vn_subgraphs, sn_subgraphs_2, vn_decompo, used_dual_costs; 
                    solver=pricer, nb_iterations=1, nb_columns=10)

        sub_mappings = result[:mappings]
        for v_subgraph in vn_decompo.subgraphs
            for mapping in sub_mappings[v_subgraph]
                add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
                nb_columns += 1
            end
        end



        # ---- master problem part

        optimize!(model)
        time_master +=  solve_time(model)       
        #println("Time last iter master : $(solve_time(model))")
        cg_value = objective_value(model)

        time_overall = time()-time_beginning

        @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs \n",
                    nb_iter, cg_value, lower_bound, nb_columns, time_overall)
    


        
        # ----- useful things

        if nb_columns > nb_columns_local_search && pricer == "local-search"
            pricer="milp"
            println("Switching to milp!")
        end

        if nb_columns > nb_columns_milp
            keep_on = false
            reason="Got enough columns"
        end

    end
    println("\n Step 2 finished, reason: $reason.")
    

    
    
    # ====== STEP 3: full pricers
    println("\n------- Solving method: Exact pricers")

    pricers_full = Dict()
    for subgraph in vn_decompo.subgraphs
        pricers_full[subgraph] = set_up_pricer(instance, subgraph)
    end
    old_cg_value = cg_value
    keep_on = true
    reason = "I don't know"
    while keep_on
        nb_iter += 1

        # ---- Pricers things
        current_dual_costs = get_duals(instance, vn_decompo, master_problem)
        old_dual_costs = used_dual_costs
        used_dual_costs = average_dual_costs(instance, vn_decompo, old_dual_costs, current_dual_costs, alpha=alpha_colge)
        

        has_found_new_column = false
        sum_pricers_values = 0

        for vn_subgraph in vn_decompo.subgraphs
            pricer = pricers_full[vn_subgraph]

            time_limit_pricer = 500


            sub_mapping, true_cost, reduced_cost = update_solve_pricer(instance, vn_decompo, pricer, used_dual_costs; time_limit = time_limit_pricer)

            
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

        current_lower_bound = (1-alpha_colge) * cg_value + (alpha_colge) * old_cg_value + sum_pricers_values
        if current_lower_bound > lower_bound
            lower_bound = current_lower_bound
        end


        # ----- Master problem stuff

        optimize!(model)
        time_master +=  solve_time(model)

        old_cg_value = cg_value
        cg_value = objective_value(model)

        time_overall = time()-time_beginning
        average_obj = sum_pricers_values/length(vn_decompo.subgraphs)

        @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs  average reduced cost: %10.3f \n",
            nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_obj)


        keep_on = true
        if !has_found_new_column
            keep_on = false
            reason="no improving columns"
        end


    end

    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: idk, time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(cg_value; digits=3))")
    println("====================================================\n")



    
    
    # ======= GETTING A SOLUTION ======= #
    

    result = Dict()
    result["solving_time"] = time() - time_beginning

    return result
end






function find_columns(  instance, vn_subgraphs, sn_subgraphs, vn_decompo, dual_costs; 
                        solver="milp", nb_iterations=10, nb_columns=50)


    # Useful things
    v_network = instance.v_network
    s_network = instance.s_network
    iter = 1
    columns_found = 0
    mappings_per_subgraph = Dict()
    for v_subgraph in vn_subgraphs
        mappings_per_subgraph[v_subgraph] = []
    end

    while columns_found < nb_columns && iter <= nb_iterations

        # Associate subvn to a random subsn
        used_sub_s_network = []
        assignment_virtual_substrate_subgraphs = Dict()
        available_s_subgraphs = collect(1:length(sn_subgraphs))
        for v_subgraph in vn_subgraphs
            found = false
            nb_v_nodes = nv(v_subgraph.graph)
            compatible_s_subgraphs = [i_sn_subg 
                    for i_sn_subg in available_s_subgraphs
                    if (nv(sn_subgraphs[i_sn_subg].graph) > 1.1 * nb_v_nodes)
            ]
            i_subgraph=0
            if isempty(compatible_s_subgraphs)
                i_subgraph = rand(available_s_subgraphs)
                filter!(x -> x != i_subgraph, available_s_subgraphs) 
            else
                i_subgraph = rand(compatible_s_subgraphs)
                filter!(x -> x != i_subgraph, available_s_subgraphs) 
            end
            assignment_virtual_substrate_subgraphs[v_subgraph] = sn_subgraphs[i_subgraph]
        end

    


        # FIND THE SUBMAPPING!
        for v_subgraph in vn_subgraphs

            s_subgraph = assignment_virtual_substrate_subgraphs[v_subgraph]
            
            sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)
            
            # Additional cost thing
            additional_costs_routing = []
            for v_node in vertices(v_subgraph.graph)
                current_addition_costs = zeros(nv(s_subgraph.graph))    
                push!(additional_costs_routing, current_addition_costs)
            end
    

            # GETTING THE SUBMAPPING
            if solver == "local-search"
                result = solve_local_search_pricer_subsn_routing(v_subgraph, s_subgraph, sub_instance, instance, vn_decompo, dual_costs, additional_costs_routing, nb_iterations=1000)
                sub_mapping = result[:sub_mapping]
                cost = result[:real_cost]
            elseif solver == "milp"
                result = solve_pricer_milp_routing(v_subgraph, s_subgraph, instance, vn_decompo, additional_costs_routing, dual_costs; time_solver = 60)
                sub_mapping = result[:sub_mapping]
                cost = result[:real_cost]
            else
                println("Pricer unknown???")
                return
            end

            if isnothing(sub_mapping) # invalid submapping!
                #print("A pricer failed. ")
                continue
            end
        
            push!(mappings_per_subgraph[v_subgraph], sub_mapping)
            columns_found += 1
            
        end

        iter += 1
    end

    return (    mappings=mappings_per_subgraph, 
                average_dual_costs=average_dual_costs
    )    


end






function complete_clusters(original_clusters, s_network; nb_nodes_to_have = 30)


    # 2) adding nodes. For now, any adjacent nodes will do.
    i_cluster = 1

    final_clusters = []

    for original_cluster in original_clusters
        cluster = copy(original_cluster)
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

        while nb_nodes_with_capacity < nb_nodes_to_have

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

        sub_s_network = my_induced_subgraph(s_network, cluster, "sub_sn_$i_cluster")
        i_cluster += 1

        push!(final_clusters, cluster)
    end

    return final_clusters

end

