

using Revise
includet("../../subgraph-decomposition/final/test_stuff/final_test.jl")



instance = get_instance_from_folder("test_cg/")
# to warm up julia
solve_subgraph_decompo(instance, time_max=60, nb_part=10)


df_stats_integer_sols = DataFrame(Sol2=Int[], Sol5=Int[], Sol10=Int[])

for k in 2:10
    val_2min = solve_subgraph_decompo(instance, time_max=125, nb_part=k)
    val_5min = solve_subgraph_decompo(instance, time_max=333, nb_part=k)
    val_10min= solve_subgraph_decompo(instance, time_max=666, nb_part=k)

    push!(df_stats_integer_sols, (val_2min, val_5min, val_10min)) 

    CSV.write("integer_stats_solving.csv", df_stats)

end




