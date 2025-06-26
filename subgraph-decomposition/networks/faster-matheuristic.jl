using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf


#general
includet("../../utils/import_utils.jl")

# utils colge
includet("utils/utils-subgraphdecompo.jl")
includet("utils/partition-vn.jl")
includet("utils/checkers.jl")

# pricers
includet("init/init_uepso.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")
includet("end-heuristic/local-search-exact.jl")



function solve_subgraph_decompo(instance; time_max = 100, v_node_partitionning = [], nb_part = -1, type_pricer="normal")

    println("Starting...")
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network


    # ======= SETTING UP THE DECOMPOSITION ======= #
    if v_node_partitionning == []
        if nb_part<0
            nb_part = floor(Int, nv(v_network.graph)/9)+1
        end
        v_node_partitionning = partition_vn_metis(instance, nb_part)
    end

    vn_decompo = set_up_decompo(instance, v_node_partitionning)

    
    println("Decomposition set: ")
        println("For $v_network, there is $(length(vn_decompo.subgraphs)) subgraphs:")

    for subgraph in vn_decompo.subgraphs
        println("       $(subgraph.graph[][:name]) with $(nv(subgraph.graph)) nodes")
    end
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")

    
    master_problem = set_up_master_problem(instance, vn_decompo)
    model = master_problem.model
    print("Master problem set... ")

    # ====== PAVING THE NETWORK WITH HEURISTIC ======= #

    println("Paving time...")
    time_0 = time()
    nb_column_per_subgraph = 30
    mappings = init_uepso(instance, vn_decompo, nb_column_per_subgraph)
    println("Mappings gotten! In just $(time() - time_0)")
    for v_subgraph in vn_decompo.subgraphs
        for mapping in mappings[v_subgraph]
            add_column(master_problem, instance, v_subgraph, mapping, get_cost_placement(mapping) + get_cost_routing(mapping))
        end
    end


        

    
    # ======= GETTING A SOLUTION ======= #
    time_solution = 20
    val, heur_sol = basic_heuristic(instance, vn_decompo, master_problem, time_solution)
    #local_search(instance, vn_decompo, heur_sol)



    # ======= LOCAL SEARCH TIME ====== #
    time_local_search = 30
    local_search_changin(instance, heur_sol, time_local_search)

    result = Dict()
    result["time"] = time() - time_beginning
    result["heuristic_res"] = heur_sol

    return result
end


