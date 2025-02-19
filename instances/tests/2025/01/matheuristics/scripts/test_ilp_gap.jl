
using DataFrames, CSV
using JuMP, CPLEX, Gurobi
includet("../../../../../../../utils/import_utils.jl")




struct Resultat
    algo
    vn
    n_r
    sn
    n_s
    best_sol
    status
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
 


function solve(instance, time_solver = 100)
    


    # ILP
    model_ip = Model(CPLEX.Optimizer)
    set_silent(model_ip)
    set_time_limit_sec(model_ip, time_solver)
    set_up_basic!(instance, model_ip)
    set_optimizer_attribute(model_ip, "CPXPARAM_MIP_Tolerances_MIPGap", 0.1)
    optimize!(model_ip)
    status_sol = primal_status(model_ip)
    status_solver = termination_status(model_ip)

    status = ""
    best_sol = -1.
    if status_sol == MOI.FEASIBLE_POINT
        gap = arrondi(relative_gap(model_ip))
        if gap > 0.10001
            best_sol = arrondi(objective_value(model_ip))
            println("       Unfinished...")    
            status = "Unfinished"
        else
            best_sol = arrondi(objective_value(model_ip))
            status = "Finished"
        end
    elseif status_solver == MOI.INFEASIBLE
        best_sol = -9999999.
        status = "Unfeasible"
        println("       Unfeasible...")
    else
        best_sol = -9999999.
        status = "No sol found"
        println("       No solution found...")
    end
    time = arrondi(solve_time(model_ip))



    return Resultat(
        "ILP_CPLEX",
        instance.v_network[][:name],
        instance.s_network[][:name],
        nv(instance.v_network),
        nv(instance.s_network),
        best_sol,
        status,
        time,
        Dict()
    )
end


function set_up_basic!(instance, model)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    ### Variables
    @variable(model, x[vertices(v_network), vertices(instance.s_network)], binary=true);
    @variable(model, y[edges(v_network), edges(s_network_dir)], binary=true);

    

    ### Objective
    placement_cost = @expression(model, sum( instance.s_network[s_node][:cost] * v_network[v_node][:dem] * x[v_node, s_node] 
        for v_node in vertices(v_network) for s_node in vertices(instance.s_network) ))
    routing_cost = @expression(model, sum( s_network_dir[src(s_edge), dst(s_edge)][:cost] * v_network[src(v_edge), dst(v_edge)][:dem] * y[v_edge, s_edge]
        for v_edge in edges(v_network) for s_edge in edges(s_network_dir) ))
    @objective(model, Min, placement_cost + routing_cost);




    ### Constraints

    ## Nodes

    # one substrate node per virtual node
    for v_node in vertices(v_network)
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(instance.s_network)) == 1)
    end

    # one to one : one virtual node per substrate node
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(x[v_node, s_node] for v_node in vertices(v_network)) <= 1)
    end

    # node capacity : NOT USELESS AHHHHHHHHh
    for s_node in vertices(instance.s_network)
        @constraint(model, sum(v_network[v_node][:dem] * x[v_node, s_node] for v_node in vertices(v_network)) <= sum(s_network[s_node][:cap]))
    end

    ## Edges 
    
    # edge capacity (undirected version !)
    for s_edge in edges(instance.s_network)
        @constraint(model, 
            sum( v_network[src(v_edge), dst(v_edge)][:dem] * (y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]  )
                for v_edge in edges(v_network)) 
            <= 
            instance.s_network[src(s_edge), dst(s_edge)][:cap] )
    end
    
    # Flow conservation
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, 
                x[src(v_edge), s_node] - x[dst(v_edge), s_node] 
                ==
                sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) - 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node))
            )
        end
    end

    
    # Departure constraints
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end


    #--------------- infeasibility constraints: if there is just not enough capacity on edges yknow...
    nb_var_less = 0
    for v_node in vertices(v_network)
        for s_node in vertices(s_network)
            v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
            necessary_bw = sum(v_network[src(v_edge), dst(v_edge)][:dem] for v_edge in v_edges_incident)

            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident)
            if necessary_bw > available_bw
                nb_var_less += 1
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end
    
    
end



function arrondi(val)
    return (floor(val*1000)/1000)
end