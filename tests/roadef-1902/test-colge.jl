
using DataFrames, CSV
using Revise
#using JuMP, CPLEX, Gurobi

#includet("../../colge/exodia/exodia2.jl")
includet("../../colge/leviathan/leviathan.jl")
includet("../../utils/import_utils.jl")







function solve_all_instances_colge(path)


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

    overall_df = DataFrame(Gr=String[], Gs=String[], Algo=String[], Times=Int[], Res=Int32[], CGval=Float32[], LGBound=Float32[], nbcol=Int32[])

    for vn in vns
        println("Doing vn $(vn[][:name])...")
        for sn in sns
            println("   for sn $(sn[][:name])...")
            instance = Instance_Undir_VNE_1s(vn, sn)
            
            nb_subgraphs=-1
            if vn[][:name] == "C5_K5"
                nb_subgraphs=5
            end
            if vn[][:name] == "W5_AW9"
                nb_subgraphs=5
            end


            for time_limit in [20, 60]
                res, nbsg, CG_bound, LG_bound, nb_col = solve_subgraph_decompo_tests(instance; time_max = time_limit, v_node_partitionning = [], nb_part = nb_subgraphs)
                push!(overall_df, (instance.v_network[][:name], instance.s_network[][:name], "Colge", time_limit, nbsg, res, CG_bound, LG_bound, nb_col)) 
                CSV.write("xp_colge_$(last(path, 4)[1:3]).csv", overall_df)
            end
        end
    end
    
    println(overall_df)

end






function main(ARGS)
    for arg in ARGS
        solve_all_instances_colge(arg)
    end
end

main(ARGS)