using Revise
includet("../../subgraph-decomposition/final/test_stuff/final_test_reduced_costs.jl")



instance = get_instance_from_folder("test_cg/")
# to warm up julia
solve_subgraph_decompo(instance, time_max=10, nb_part=10)




# now the real deal 
solve_subgraph_decompo(instance, time_max=666, nb_part=4)
