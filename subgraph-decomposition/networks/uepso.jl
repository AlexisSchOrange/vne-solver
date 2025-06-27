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

    heur_sol, cost = solve_UEPSO(instance)
        

    

    # ======= LOCAL SEARCH TIME ====== #
    time_local_search = 30
    local_search_changin(instance, heur_sol, time_local_search)

    result = Dict()
    result["time"] = time() - time_beginning
    result["heuristic_res"] = heur_sol

    return result
end


