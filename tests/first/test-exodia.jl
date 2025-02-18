
using DataFrames, CSV
using JuMP, CPLEX, Gurobi


includet("../colge/exodia/exodiavisu.jl")
includet("../compact/compact_undir.jl")







function solve_all_instances(pathvn, pathsn)


    vns = []
    for filename in readdir(pathvn; join=true)
        g, type = read_graph(filename)
        push!(vns, g)
    end

    sns = []
    for filename in readdir(pathsn; join=true)
        g, type = read_graph(filename)
        push!(sns, g)
    end

    overall_df = DataFrame(Gr=String[], Gs=String[], ILP10s=Int[], ILP60s=Int[], ColGe10s=Int[], ColGe60s=Int[])

    for vn in vns
        println("Doing vn $(vn[][:name])...")
        for sn in sns
            println("   for sn $(sn[][:name])...")
            instance = Instance_Undir_VNE_1s(vn, sn)


            val_ilp_10s = ceil(Int, solve_compact(instance, 10, true))
            val_ilp_60s = ceil(Int, solve_compact(instance, 60, true))
            val_colge_10s = ceil(Int, solve_subgraph_decompo(instance, 9, [], 4))
            val_colge_60s = ceil(solve_subgraph_decompo(instance, 55, [], 4))

            push!(overall_df, (instance.v_network[][:name], instance.s_network[][:name], val_ilp_10s, val_ilp_60s, val_colge_10s, val_colge_60s))  # Add a row with values

            CSV.write("experiment_results.csv", overall_df)
        end
    end
    
    println(overall_df)

end



function visu_all(pathvn, pathsn)


    vns = []
    for filename in readdir(pathvn; join=true)
        g, type = read_graph(filename)
        push!(vns, g)
    end

    sns = []
    for filename in readdir(pathsn; join=true)
        g, type = read_graph(filename)
        push!(sns, g)
    end


    for vn in vns
        println("Doing vn $(vn[][:name])...")
        for sn in sns
            println("   for sn $(sn[][:name])...")
            instance = Instance_Undir_VNE_1s(vn, sn)
            solve_subgraph_decompo(instance, 9, [], 4)
        end
    end
    
    println(overall_df)

end
