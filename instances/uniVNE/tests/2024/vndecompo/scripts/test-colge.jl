
using DataFrames, CSV
using JuMP, CPLEX, Gurobi
includet("../../../../../utils/import_utils.jl")




struct Resultat
    algo
    vn
    sn
    cg_value
    nb_columns
    time
    parameters
end




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

    overall_df = nothing

    results = []
    for vn in vns
        println("Doing vn $(vn[][:name])...")
        for sn in sns
            println("   for sn $(sn[][:name])...")
            instance = Instance_Undir_VNE_1s(vn, sn)
            result = solve(instance, 100)
            if overall_df === nothing
                overall_df = DataFrame(Dict(n=>[getfield(result, n)] for n in fieldnames(Resultat)))
            else
                current_df = DataFrame(Dict(n=>[getfield(result, n)] for n in fieldnames(Resultat)))
                append!(overall_df, current_df)
            end
            CSV.write("experiment_results.csv", overall_df)
        end
    end
    
    println(overall_df)

end
 


function solve(instance, time_solver = 30)
    
    algo = "local_cuts_both"

    lp_relax = -1.
    int_value = -1.
    lp_relax_end = -1.
    gap = -1.
    nb_node = -1
    time = -1

    nb_vnodes_cuts = 5
    parameters = Dict("nb_vnodes"=>nb_vnodes_cuts)

    # LP
    model_lp = Model(CPLEX.Optimizer)
    set_silent(model_lp)
    set_up_both!(instance, model_lp, nb_vnodes_cuts)
    relax_integrality(model_lp)
    optimize!(model_lp)
    lp_relax = arrondi(objective_value(model_lp))

    # ILP
    model_ip = Model(CPLEX.Optimizer)
    set_silent(model_ip)
    set_time_limit_sec(model_ip, time_solver)
    parameters = set_up_both!(instance, model_ip, nb_vnodes_cuts)
    optimize!(model_ip)
    status_sol = primal_status(model_ip)
    status_solver = termination_status(model_ip)
    if status_sol == MOI.FEASIBLE_POINT
        int_value = arrondi(objective_value(model_ip))
        gap = arrondi(relative_gap(model_ip))
        lp_relax_end = arrondi(int_value / (gap + 1))
    elseif status_solver == MOI.INFEASIBLE
        int_value = -9999999.
        println("       Unfeasible...")
    else
        println("       No solution found...")
    end
    nb_node = node_count(model_ip)
    time = arrondi(solve_time(model_ip))



    return Resultat(
        algo,
        instance.v_network[][:name],
        instance.s_network[][:name],
        lp_relax,
        int_value,
        lp_relax_end,
        gap,
        nb_node,
        time,
        parameters
    )
end