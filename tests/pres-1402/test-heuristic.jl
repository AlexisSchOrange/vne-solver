
using DataFrames, CSV
using JuMP, CPLEX, Gurobi

includet("../../utils/import_utils.jl")
includet("../../heuristics/VNE-DCC.jl")


function solve_all_instances(path)


    vns = []
    for filename in readdir(path*"vns"; join=true)
        g, type = read_graph(filename)
        push!(vns, g)
    end

    sns = []
    for filename in readdir(path*"sns"; join=true)
        g, type = read_graph(filename)
        push!(sns, g)
    end

    overall_df = DataFrame(Gr=String[], Gs=String[], Algo=String[], Res=Int[])

    for vn in vns
        println("Doing vn $(vn[][:name])...")
        for sn in sns
            println("   for sn $(sn[][:name])...")
            instance = Instance_Undir_VNE_1s(vn, sn)

            
            val_heuristic = ceil(Int, solve_VNE_DCC(instance))
            push!(overall_df, (instance.v_network[][:name], instance.s_network[][:name], "heuristc", val_heuristic)) 
            
            CSV.write("experiment_results.csv", overall_df)
        end
    end
    
    println(overall_df)

end

