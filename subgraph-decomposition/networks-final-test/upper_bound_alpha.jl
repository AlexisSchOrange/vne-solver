using Printf
using Graphs, MetaGraphsNext
using JuMP, CPLEX


#general
includet("../../utils/import_utils.jl")

# utils colge
includet("utils/master-problem.jl")
includet("utils/graph-decomposition.jl")
includet("utils/partition-graph.jl")

# pricers
includet("pricers/milp-pricer-subsn-routing.jl")
includet("pricers/local-search-pricer-subsn-routing.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")
includet("end-heuristic/local-search-exact.jl")



function solve_price_branch(instance; pricer="milp", nb_virtual_subgraph=0, nb_columns_max=300, time_end_milp=30, alpha_colge=0.5)

    # === SOME USEFUL THINGS === #
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    println("Starting...")
    time_beginning = time()
    nb_columns = 0
    nb_columns_already_found=0

    # === SOME PARAMETERS === #
    nb_columns_init = min(nb_columns_max, 100)




    # ======= SETTING UP THE DECOMPOSITION ======= #

    # ------ Virtual decomposition ------ #
    if nb_virtual_subgraph == 0
        nb_virtual_subgraph = floor(Int, nv(v_network.graph)/10)
    end
    v_node_partitionning = partition_graph(v_network.graph, nb_virtual_subgraph, max_umbalance=1.2)   
    vn_decompo = set_up_decompo(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs
    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")
    
    # ------ Substrate decomposition ------ #
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







    # Starting the algorithm
    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)

    empty_dual_costs = get_empty_duals(instance, vn_decompo)


    # ======= INITIALIZATION  ======= #
    println("STEP 1: INITIALIZATION...")
    time_0=time()

    result = find_columns(instance, vn_subgraphs, sn_subgraphs, vn_decompo, empty_dual_costs, shortest_paths, solver=pricer, nb_iterations=nb_columns_init, nb_columns=nb_columns_init)
    sub_mappings = result[:mappings]
    for v_subgraph in vn_decompo.subgraphs
        for mapping in sub_mappings[v_subgraph]
            add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
            nb_columns += 1
            if !check_if_column_new(master_problem, mapping, v_subgraph)
                nb_columns_already_found += 1
            end
        end
    end

    println("\nInitialization dong: \n$nb_columns columns gotten, with $pricer, in $(time() - time_0)s\n")


    # ======= FIRST ITERATION OF MASTER PROBLEM ======= #
    optimize!(model)
    status = termination_status(model)
    if status != MOI.OPTIMAL
        println("Master problem is infeasible or unfinished: $status")
        return
    end
    rmp_value = objective_value(model)
    time_overall = time() - time_beginning
    nb_iter = 0
    @printf("Iter %2d  RMP value: %10.3f  %5d column    time: %5.2fs  \n", nb_iter, rmp_value, nb_columns, time_overall)


    # ======= COLUMN GENERATION ======= #
    while nb_columns < nb_columns_max

        # --- pricers
        current_dual_costs = get_duals(instance, vn_decompo, master_problem)
        dual_costs_to_use = average_dual_costs(instance, vn_decompo, empty_dual_costs, current_dual_costs, alpha=alpha_colge)
        
        result =  find_columns(instance, vn_subgraphs, sn_subgraphs, vn_decompo, dual_costs_to_use, shortest_paths, solver=pricer, nb_iterations=1, nb_columns=nb_virtual_subgraph)
        sub_mappings = result[:mappings]
        for v_subgraph in vn_decompo.subgraphs
            for mapping in sub_mappings[v_subgraph]
                if !check_if_column_new(master_problem, mapping, v_subgraph)
                    nb_columns_already_found += 1
                end
                add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
                nb_columns += 1
            end
        end


        # ---- master problem part
        optimize!(model)
        rmp_value = objective_value(model)


        # --- printing and utils
        time_overall = time()-time_beginning
        nb_iter += 1
        #average_obj = overall_reduced_costs / (nb_feasible+10e-6)
        #=
        @printf("Iter %2d    RMP value:%10.3f  %5d columns  time: %5.2fs  mean red. cost: %10.3f   feasible %10.3f / %10.3f \n",
            nb_iter, rmp_value, nb_columns, time_overall, average_obj, nb_feasible, nb_virtual_subgraph
        )
        =#
        @printf("Iter %2d    RMP value:%10.3f  %5d columns  time: %5.2fs \n",
            nb_iter, rmp_value, nb_columns, time_overall
        )
    end

    # ======= END HEURISTICS ======= #

    # ---- Restricted master problem heuristic
    value_cg_heuristic, cg_heuristic_solution = basic_heuristic(instance, vn_decompo, master_problem, time_end_milp)

    println("\nSolution obtained: $value_cg_heuristic\n")

    #local_search_changin(instance, cg_heuristic_solution, 300)

    return (mapping_cost=value_cg_heuristic,
        mapping = nothing,
        column_overall = nb_columns,
        columns_repetitive = nb_columns_already_found )
end







function find_columns(  instance, vn_subgraphs, sn_subgraphs, vn_decompo, dual_costs, base_shortest_paths; 
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

        temporary_placement = zeros(Int, nv(v_network))
        
        for v_subgraph in vn_subgraphs
            s_subgraph = assignment_virtual_substrate_subgraphs[v_subgraph]
            cluster = s_subgraph.nodes_of_main_graph
            capacities_score = [ get_attribute_node(s_network, s_node, :cap) * 
                    sum(get_attribute_edge(s_network, get_edge(s_network, s_node, s_neighbor), :cap) for s_neighbor in neighbors(s_network, s_node)) 
                    for s_node in cluster]
            best_node = cluster[findmin(capacities_score)[2]]
            for v_node in v_subgraph.nodes_of_main_graph
                temporary_placement[v_node] = best_node
            end
        end
    


        # FIND THE SUBMAPPING!
        for v_subgraph in vn_subgraphs

            s_subgraph = assignment_virtual_substrate_subgraphs[v_subgraph]
            
            sub_instance = Instance(v_subgraph.graph, s_subgraph.graph)
            
            # Additional cost thing
            additional_costs_routing = []
            for v_node in vertices(v_subgraph.graph)
                current_addition_costs = zeros(nv(s_subgraph.graph))
                original_v_node = v_subgraph.nodes_of_main_graph[v_node]
                
                for v_edge in vn_decompo.v_edges_master
                    if src(v_edge) == original_v_node
                        placement_of_dst_node = temporary_placement[dst(v_edge)]
                        for s_node in vertices(s_subgraph.graph)
                            original_s_node = s_subgraph.nodes_of_main_graph[s_node]
                            current_addition_costs[s_node] += base_shortest_paths.dists[original_s_node, placement_of_dst_node]
                        end
                    end
                    if dst(v_edge) == original_v_node
                        placement_of_dst_node = temporary_placement[src(v_edge)]
                        for s_node in vertices(s_subgraph.graph)
                            original_s_node = s_subgraph.nodes_of_main_graph[s_node]
                            current_addition_costs[s_node] += base_shortest_paths.dists[original_s_node, placement_of_dst_node]
                        end
                    end
                end
                
    
                push!(additional_costs_routing, current_addition_costs)
            end
    

            # GETTING THE SUBMAPPING
            if solver == "local-search"
                result =    solve_local_search_pricer_subsn_routing(v_subgraph, s_subgraph, sub_instance, instance, vn_decompo, dual_costs, additional_costs_routing, nb_iterations=250)
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
            
            # Update temporary placement with the placement chosen
            for v_node in vertices(v_subgraph.graph)
                original_v_node = v_subgraph.nodes_of_main_graph[v_node]
                temporary_placement[original_v_node] = sub_mapping.node_placement[v_node]
            end

        end

        iter += 1
    end

    return (    mappings=mappings_per_subgraph, 
                average_dual_costs=average_dual_costs
    )    


end


