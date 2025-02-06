
using DataFrames, CSV
using JuMP, CPLEX, Gurobi
includet("../../../../../../../utils/import_utils.jl")

using Metis


struct Resultat
    algo
    vn
    n_r
    sn
    n_s
    best_sol
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
            result = solve_partitionning_sn(instance)
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
 





function solve_partitionning_sn(instance)

    v_network = instance.v_network
    s_network = instance.s_network


    time_start = time()


    nb_clusters = Int(floor(nv(s_network) /  nv(v_network)))
    partition = Metis.partition(s_network.graph, nb_clusters)

    parameters = Dict("Nb clusters"=>nb_clusters)

    part = Dict(i=>Vector{Int64}() for i in 1:nb_clusters)
    for i in 1:nv(s_network)
        push!(part[partition[i]], i)
    end


    status = "No solution"
    best_solution = +999999.

    for i_cluster in 1:nb_clusters
        sub_s_network = my_induced_subgraph(s_network, part[i_cluster], "sn_$i_cluster")
        subinstance = Instance_Undir_VNE_1s(v_network, sub_s_network)

        current_model = Model(CPLEX.Optimizer)
        set_up_problem(subinstance, current_model)
        set_time_limit_sec(current_model, 10)
        set_silent(current_model)

        optimize!(current_model)

        status_sol = primal_status(current_model)
    
        if status_sol == MOI.FEASIBLE_POINT
            best_solution = arrondi(objective_value(current_model))
            status = "Solution"
        end
    end

    time_algo = arrondi(time() - time_start)

    return Resultat(
        "partitionning_sn",
        instance.v_network[][:name],
        instance.s_network[][:name],
        nv(instance.v_network),
        nv(instance.s_network),
        best_solution,
        time_algo,
        parameters
    )

end



function set_up_problem(instance, model)

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

    
    ## Departure constraints
    for s_node in vertices(instance.s_network)
        for v_edge in edges(v_network)
            @constraint(model, sum(y[v_edge, s_edge] for s_edge in get_out_edges(s_network_dir, s_node)) 
                >= x[src(v_edge), s_node])
        end
    end
    
    # Outgoing edges cap: pretty stupid but useful
    i = 0
    for v_node in vertices(v_network)
        for s_node in vertices(s_network)
            v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
            necessary_bw = 0 + sum(v_network[src(v_edge), dst(v_edge)][:dem] for v_edge in v_edges_incident; init=0.)

            s_edges_incident = [get_edge(s_network, s_node, neighbor) for neighbor in neighbors(s_network, s_node)]
            available_bw = 0 +sum(s_network[src(s_edge), dst(s_edge)][:cap] for s_edge in s_edges_incident; init=0.)
            if necessary_bw > available_bw
                i+=1
                @constraint(model, model[:x][v_node, s_node] == 0)
            end 
        end
    end
    #println("We get this to delete: $i")
end





function arrondi(val)
    return (floor(val*1000)/1000)
end