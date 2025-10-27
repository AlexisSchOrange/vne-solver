using Graphs, MetaGraphsNext
using JuMP, CPLEX

includet("graph_decomposition_overlapping.jl")
includet("master_problem.jl")
includet("../pricers/pricer-full.jl")


function column_generation(instance, vn_decompo, master_problem; time_max = 900)

    model = master_problem.model

    time_beginning = time()

    #------------ GENERATION DE COLONNES
    nb_columns = 0
    nb_iter = 0

    time_master = 0
    time_subproblems = 0
    time_overall = time()-time_beginning
    
    cg_value = 10e9
    lower_bound =  0

    alpha_smoothing = 0.
    dual_costs = get_empty_duals(instance, vn_decompo, master_problem)
    average_obj = -10000

    print("\n\n==================== Starting CG ====================\n")

    # ====== full pricers
    println("\n------- Solving method: Exact pricers")

    pricers_full = Dict()
    for subgraph in vn_decompo.subgraphs
        pricers_full[subgraph] = set_up_pricer(instance, subgraph)
    end
    keep_on = true
    reason = "I don't know"
    while keep_on
        nb_iter += 1

        optimize!(model)
        time_master +=  solve_time(model)

        status = termination_status(model)
        if status != MOI.OPTIMAL
            println("Infeasible or unfinished: $status")
            return
        end

        if cg_value < 5*10e3 && average_obj > -50.
            alpha_smoothing = 0.85
        end

        old_cg_value = cg_value
        cg_value = (1-alpha_smoothing) * objective_value(model) + alpha_smoothing * old_cg_value

        old_dual_costs = dual_costs
        current_dual_costs = get_duals(instance, vn_decompo, master_problem)
        dual_costs = average_dual_costs(instance, vn_decompo, old_dual_costs, current_dual_costs, alpha=alpha_smoothing)


            

        time_beginning_pricer = time()
        has_found_new_column = false
        sum_pricers_values = 0

        for vn_subgraph in vn_decompo.subgraphs
            pricer = pricers_full[vn_subgraph]

            time_limit_subpb = 100

            column, true_cost, reduced_cost = update_solve_pricer(instance, vn_decompo, pricer, dual_costs; time_limit = time_limit_subpb)

            if reduced_cost < -0.0001
                has_found_new_column = true 
                add_column(master_problem, instance, vn_decompo, vn_subgraph, column, true_cost)
                nb_columns += 1
            end

            sum_pricers_values += reduced_cost
        end

        current_lower_bound = cg_value + sum_pricers_values
        if cg_value < 5*10e3 && current_lower_bound > lower_bound
            lower_bound = current_lower_bound
        end

        time_subproblems += time() - time_beginning_pricer

        average_obj = sum_pricers_values/length(vn_decompo.subgraphs)
        @printf("Iter %2d  CG bound: %10.3f  lower bound: %10.3f  %5d column  time: %5.2fs  average reduced cost: %10.3f \n",
            nb_iter, cg_value, lower_bound, nb_columns, time_overall, average_obj)



        time_overall = time()-time_beginning
        if time_overall < time_max
            keep_on = true
            if !has_found_new_column
                keep_on = false
                reason="no improving columns"
            end
        else
            keep_on = false
            reason="time limit"
        end


    end

        


    print("\n==================== CG finished ====================\nReason: $reason \n")
    println("Time in MP: $(round(time_master; digits=3)) , time in SP: $(round(time_subproblems; digits=3)), time overall: $(round(time_overall; digits=3))")
    println("$nb_iter iters, final value: $(round(cg_value; digits=3))")
    println("====================================================\n")








end