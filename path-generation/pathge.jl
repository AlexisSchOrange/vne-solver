

includet("../utils/import_utils.jl")
includet("utils-pathsge.jl")

using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX



function solve_pathge(instance)

    v_network = instance.v_network

    master_problem = set_up_master_problem(instance)


    # add original paths: those that are just one edge long...
    for s_edge in edges(instance.s_network_dir)
        cost = instance.s_network_dir[src(s_edge), dst(s_edge)][:cost]
        s_path = Path(src(s_edge), dst(s_edge), [s_edge], cost)
        for v_edge in edges(v_network)
            add_column(master_problem, v_edge, s_path)
        end
    end



    # solve the master problem for that shit.


    #=
    pricers = Dict()
    for v_edge in edges(v_network)
        pricers[v_edge] = set_up_pricer(v_network)
    end


    while keep_on

        keep_on = false

        optimize!(master_problem.model)

        print("CG value: $(objective_value(model))")
        for v_edge in edges(g_network)
            column, value = update_and_solve_pricer(pricers[v_edge], dual_values)
            if value < -0.001
                add_column(master_problem, column)
                keep_on = true
            end

        end
    
    end
    =#

    return

end
