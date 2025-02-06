
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
            result = progressive_rounding(instance)
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
 





function progressive_rounding(instance)

    time_start = time()

    v_network = instance.v_network
    s_network = instance.s_network
    parameters = Dict()

    nodes_to_restrict_next = []

    nb_neighbors = []
    nb_neighbors_done = []
    for v_node in vertices(v_network)
        push!(nb_neighbors, length(neighbors(v_network, v_node)))
        push!(nb_neighbors_done, 0)
    end


    restrictions = Dict()
    restricted_nodes = []

    # we simply takes the most dense node and make it start here
    sorted_degree = sort(collect(vertices(v_network)), by=v_node -> -degree(v_network, v_node))
    nodes_to_restrict_next = sorted_degree[1:1]

    current_model = Model(CPLEX.Optimizer)
    set_up_problem(instance, current_model) 
    add_star_constraints(instance, current_model)
    unrelax = relax_integrality(current_model)        
    set_silent(current_model)


    iter=0
    while length(restricted_nodes) != length(vertices(v_network))
        iter+=1
        #println("Current placement restriction: $(restrictions)")
    
        optimize!(current_model)

        if !has_values(current_model)
            println("       Infeasible model !")
            return pls(instance, parameters, time()-time_start)
        end
        #println("Solving time: $(solve_time(current_model)) ")
        
        x_values = value.(current_model[:x])


        for v_node in nodes_to_restrict_next
            selec = []
            vals = [(x_values[v_node, s_node] * 100
                    + sum(x_values[v_node, neighbor] * 50 for neighbor in neighbors(s_network, s_node); init=0.0)) 
                        for s_node in vertices(s_network)]
            sorted = sortperm(vals, rev=true)
            selec = sorted[1:1]
            restrictions[v_node] = selec

            push!(restricted_nodes, v_node)
            @constraint(current_model, current_model[:x][v_node, selec[1]]  == 1)

        end


        # -------- NEXT NODES SELECTION: quite tricky ?
        # idea : 1) most central left, 2) most neighbors 3) less frac 4) all neighbors ?
        if length(restricted_nodes) != length(vertices(v_network))
             # On en prend qu'un, le plus dense parmis les voisins... ?
            possible_nodes = []
            for v_node in restricted_nodes
                for neighbor in neighbors(v_network, v_node)
                    if neighbor ∉ restricted_nodes
                        if neighbor ∉ possible_nodes
                            push!(possible_nodes, neighbor)
                        end
                    end
                end
            end
            ratio = Dict()
            for v_node in possible_nodes
                #ratio[v_node] = -nb_neighbors_done[v_node] / nb_neighbors[v_node] # with less neighbor already fixed first
                #ratio[v_node] = nb_neighbors_done[v_node] / nb_neighbors[v_node] # with most neighbor already fixed first
                ratio[v_node] = degree(v_network, v_node) # with highest degree node
            end
            sorted_keys = sort(collect(keys(ratio)), by=x -> -ratio[x])
            nodes_to_restrict_next = [sorted_keys[1]]

        end

    end

    # model improved
    unrelax()
    optimize!(current_model)
    if !has_values(current_model)
        println("       Infeasible model !")
        return pls(instance, parameters, time()-time_start)
    end

    best_sol = objective_value(current_model)

    time_algo = arrondi(time() - time_start)



    return Resultat(
        "Progressive_rounding",
        instance.v_network[][:name],
        instance.s_network[][:name],
        nv(instance.v_network),
        nv(instance.s_network),
        best_sol,
        time_algo,
        parameters
    )

end



function pls(instance, parameters, time_algo)

    return Resultat(
        "Progressive_rounding",
        instance.v_network[][:name],
        instance.s_network[][:name],
        nv(instance.v_network),
        nv(instance.s_network),
        -9999.,
        time_algo,
        parameters
    )

end


function arrondi(val)
    return (floor(val*1000)/1000)
end




function add_star_constraints(instance, model)

    v_network = instance.v_network
    s_network_dir = instance.s_network_dir
    s_network = instance.s_network

    nb_vnodes = length(vertices(v_network))
    # star stuff
    # Get a list of nodes and their degrees
    node_degrees = [(v, degree(v_network, v)) for v in vertices(v_network)]

    # Sort nodes by degree in descending order and take the top five
    dense_v_nodes = sort(node_degrees, by=x -> -x[2])[1:nb_vnodes]


    #---------------- global constraints + precise
    for (v_node, v_deg) in dense_v_nodes
        trucmuche = @expression(model, v_deg)
        for s_node in vertices(s_network_dir)
            trucmuche += sum(model[:x][v_node, s_node] * (v_deg - degree(instance.s_network, s_node)))
        end

        v_edges_incident = [get_edge(v_network, v_node, neighbor) for neighbor in neighbors(v_network, v_node)]
        
        @constraint(model, 
            trucmuche <= sum(model[:y][v_edge, s_edge] for v_edge in v_edges_incident for s_edge in edges(s_network_dir))
        )
    end




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
        @constraint(model, sum(x[v_node, s_node] for s_node in vertices(s_network)) == 1)
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
