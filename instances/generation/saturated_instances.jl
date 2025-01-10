using Graphs, MetaGraphsNext
using Revise
includet("../../utils/graph.jl")
includet("io.jl")
includet("../../compact/compact_undir.jl")





function get_substrate_shortest_path(g)
    mg = MetaGraph(
        Graph(),
        Int,
        Dict,
        Dict,
        Dict(:name=>"s_" * name, :type=>"substrate")
    )

    for node in vertices(g)
        add_vertex!(mg, node, Dict(:cap=> 1, :cost => 0))
    end

    for edge in edges(g)
        add_edge!(mg, src(edge), dst(edge), Dict(:cap=>1, :cost =>1))
    end

    # Get shortest path and increase by 1 the cap everytime the edge is used in a shortest path
    distmx = ones(Int, size(vertices(g), 1), size(vertices(g), 1))
    top_cap = 1
    for i in vertices(g)
        for j in vertices(g)
            if i!=j
                yen_path = yen_k_shortest_paths(mg, i, j, distmx, 1)
                path = yen_path.paths[1]
                for i_node in 1:length(path)-1
                    mg[path[i_node], path[i_node+1]][:cap] = mg[path[i_node], path[i_node+1]][:cap] +1
                    if mg[path[i_node], path[i_node+1]][:cap] > top_cap
                        top_cap = top_cap + 1
                    end
                end
            end
        end
    end

    cap_max = 5
    for edge in edges(mg)
        mg[src(edge), dst(edge)][:cap] = floor(mg[src(edge), dst(edge)][:cap] * cap_max / top_cap - 0.00001) + 1
    end


    return mg
end



function saturate_edges!(s_network)
    n_s = length(vertices(s_network))
    # Let's put 3 small prism graphs on the sn. we add them until enough nodes are used...
    size_prism = 3
    clg = circular_ladder_graph(size_prism)

    v_network = get_vn_from_base(clg)

    println("The virtual graph look like that:")

    print_graph(v_network)

    println("Original graph:")
    print_graph(s_network)
    nb_node_max_to_be_used = floor(n_s / 4)
    nb_v_nodes = 0
    i_vn = 1
    while nb_v_nodes < nb_node_max_to_be_used
        nb_v_nodes += (size_prism * 2)

        instance = Instance_Undir_VNE_1s(v_network, s_network)
        s_network_dir = instance.s_network_dir

        model = Model(CPLEX.Optimizer)
        set_up_problem_undir_1vn_1t1(instance, model)
        set_time_limit_sec(model, 10.)
        optimize!(model)
        
        #removing capacities
        for v_node in vertices(v_network)
            for s_node in vertices(s_network_dir)
                if value(model[:x][v_node, s_node]) > 0.5
                    s_network[s_node][:cap] = s_network[s_node][:cap] - 1
                end
            end
        end
        for v_edge in edges(v_network)
            for s_edge in edges(s_network_dir)
                if value(model[:y][v_edge, s_edge]) > 0.5
                    s_network[src(s_edge), dst(s_edge)][:cap] = s_network[src(s_edge), dst(s_edge)][:cap] - 1
                end
            end
        end

        println("Graph after placing $(i_vn) graphs:")
        print_graph(s_network)
        i_vn += 1
    end
    

    println("\nFinal graph:")
    print_graph(s_network)


    return s_network


end



function get_vn_from_base(g)

    v_network = MetaGraph(
        Graph(),
        Int,
        Dict,
        Dict,
        Dict(:name=>"vn", :type=>"substrate")
    )

    for node in vertices(g)
        add_vertex!(v_network, node, Dict(:dem=> 1))
    end

    for edge in edges(g)
        add_edge!(v_network, src(edge), dst(edge), Dict(:dem=>1))
    end

    return v_network
end




function affinate_capacities!(s_network)

    for s_node in vertices(s_network)
        if s_network[s_node][:cap] < 0
            s_network[s_node][:cap] = 0
            println("Hey small problem you have a negative cap here $s_node")
        end
    end
    for s_edge in edges(s_network)
        if s_network[src(s_edge), dst(s_edge)][:cap] < 1
            s_network[src(s_edge), dst(s_edge)][:cap] = 1
        end
    end

end



function randomize_costs!(s_network)

    node_cost_max = 0
    node_cost_min = 0
    edge_cost_max = 5
    edge_cost_min = 2

    for node in vertices(s_network) 
        s_network[node][:cost] = rand(node_cost_min:node_cost_max)
    end

    for edge in edges(s_network) 
        s_network[src(edge), dst(edge)][:cost] = rand(edge_cost_min:edge_cost_max)
    end


end

