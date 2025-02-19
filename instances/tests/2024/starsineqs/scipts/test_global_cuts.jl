
using DataFrames, CSV
using JuMP, CPLEX, Gurobi
includet("../../../../../../utils/import_utils.jl")




struct Resultat
    vn
    sn
    lp_relax
    int_value
    lp_relax_end
    gap
    nb_node
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
    set_up_global_precise!(instance, model_lp, nb_vnodes_cuts)
    relax_integrality(model_lp)
    optimize!(model_lp)
    lp_relax = arrondi(objective_value(model_lp))

    # ILP
    model_ip = Model(CPLEX.Optimizer)
    set_silent(model_ip)
    set_time_limit_sec(model_ip, time_solver)
    set_up_global_precise!(instance, model_ip, nb_vnodes_cuts)
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
    
    
    
end



function set_up_global_precise!(instance, model, nb_vnodes)

    set_up_basic!(instance, model)


    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    
    # star stuff
    # Get a list of nodes and their degrees
    node_degrees = [(v, degree(v_network, v)) for v in vertices(v_network)]

    # Sort nodes by degree in descending order and take the top five
    dense_v_nodes = sort(node_degrees, by=x -> -x[2])[1:nb_vnodes]


    #---------------- global constraints + precise
    nb_cons_2 = 0
    for (v_node, v_deg) in dense_v_nodes

        trucmuche = @expression(model, v_deg)
        for s_node in vertices(s_network_dir)
            trucmuche += sum(model[:x][v_node, s_node] * (v_deg - degree(instance.s_network, s_node)))
        end

        v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
        
        cons = @constraint(model, 
            trucmuche <= sum(model[:y][v_edge, s_edge] for v_edge in v_edges_incident for s_edge in edges(s_network_dir))
        )
        
        #println("New cons: $cons")
        nb_cons_2 += 1
    end

end


function arrondi(val)
    return (floor(val*1000)/1000)
end