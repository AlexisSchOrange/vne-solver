
using DataFrames, CSV
#using JuMP, CPLEX, Gurobi

#includet("../../colge/exodia/exodia2.jl")
includet("../../compact/compact.jl")
includet("../../utils/import_utils.jl")







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

    overall_df = DataFrame(Gr=String[], Gs=String[], Algo=String[], Res5m=Int[], Res20m=Int[])

    for vn in vns
        println("Doing vn $(vn[][:name])...")
        for sn in sns
            println("   for sn $(sn[][:name])...")
            instance = Instance_Undir_VNE_1s(vn, sn)

            
            val_ilp_300s = ceil(Int, solve_compact(instance, time_solver = 15, stay_silent=true, linear=false))
            val_ilp_1200s = ceil(Int, solve_compact(instance, time_solver = 180, stay_silent=true, linear=false))
            push!(overall_df, (instance.v_network[][:name], instance.s_network[][:name], "ILP", val_ilp_300s, val_ilp_1200s)) 

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



function solve_instance(instance) 

    overall_df = DataFrame(Gr=String[], Gs=String[], Algo=String[], Res1s=Int[], Res10s=Int[], Res60s=Int[], Res180s=Int[], Res600s=Int[])


    val_ilp_1s = ceil(Int, solve_compact(instance, 1, true))
    val_ilp_10s = ceil(Int, solve_compact(instance, 10, true))
    val_ilp_60s = ceil(Int, solve_compact(instance, 60, true))
    val_ilp_180s = ceil(Int, solve_compact(instance, 180, true))
    val_ilp_600s = ceil(Int, solve_compact(instance, 600, true))
    push!(overall_df, (instance.v_network[][:name], instance.s_network[][:name], "ILP", val_ilp_1s, val_ilp_10s, val_ilp_60s, val_ilp_180s, val_ilp_600s)) 

    val_cg_1s = ceil(solve_subgraph_decompo(instance, 60, [], 1))
    val_cg_10s = ceil(solve_subgraph_decompo(instance, 60, [], 10))
    val_cg_60s = ceil(solve_subgraph_decompo(instance, 60, [], 60))
    val_cg_180s = ceil(solve_subgraph_decompo(instance, 60, [], 180))
    val_cg_600s = ceil(solve_subgraph_decompo(instance, 60, [], 600))
    push!(overall_df, (instance.v_network[][:name], instance.s_network[][:name], "ILP", val_cg_1s, val_cg_10s, val_cg_60s, val_cg_180s, val_cg_600s)) 


    CSV.write("experiment_results.csv", overall_df)
    
    println(overall_df)


end
