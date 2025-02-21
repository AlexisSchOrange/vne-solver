
using DataFrames, CSV
using Revise
#using JuMP, CPLEX, Gurobi

#includet("../../colge/exodia/exodia2.jl")
includet("../../compact/compact.jl")
includet("../../utils/import_utils.jl")







function solve_all_instances_ilp(path)


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

    overall_df = DataFrame(Gr=String[], Gs=String[], Algo=String[], Times=Int[], Res=Int32[], Gap=Float32[])

    for vn in vns
        println("Doing vn $(vn[][:name])...")
        for sn in sns
            println("   for sn $(sn[][:name])...")
            instance = Instance_Undir_VNE_1s(vn, sn)

            for time_limit in [300, 1200, 6000]
                res, gap = solve_compact(instance, time_solver = time_limit, stay_silent=true, linear=false)
                push!(overall_df, (instance.v_network[][:name], instance.s_network[][:name], "ILP", time_limit, res, gap)) 
                CSV.write("xp_ilp_$(last(path, 4)[1:3]).csv", overall_df)
            end
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




function main(ARGS)
    for arg in ARGS
        solve_all_instances_ilp(arg)
    end
end

main(ARGS)