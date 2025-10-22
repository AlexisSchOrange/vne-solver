
using Revise

using Graphs, MetaGraphsNext
using JuMP, CPLEX
using OrderedCollections
using Printf

#general
includet("../utils/import_utils.jl")

# utils colge
includet("utils/master_problem.jl")
includet("utils/graph_decomposition_overlapping.jl")
includet("utils/column_generation.jl")
includet("../utils/partition-graph.jl")

# init
includet("init/greedy_easy.jl")

# end heuristics
includet("end-heuristic/basic-ilp.jl")



function edge_partition(instance)

    println("Starting...")
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    # ======= SETTING UP THE DECOMPOSITION ======= #

    # AUTOMATIC PARTITION


    v_node_partitionning = []
    for v_edge in edges(v_network)
        push!(v_node_partitionning, [src(v_edge), dst(v_edge)])
    end


    println("Node partitionning: $v_node_partitionning")





    vn_decompo = set_up_decompo_overlapping(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    
    # === COLUMN GENERATION === #

    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")



    # generating first columns. Adapted for overlapping cg...
    for s_edge in edges(s_network_dir)

        for v_subgraph in vn_subgraphs
            node_placement = [src(s_edge), dst(s_edge)]
            cost_mapping = s_network_dir[src(s_edge)][:cost] + s_network_dir[dst(s_edge)][:cost] + s_network_dir[src(s_edge), dst(s_edge)][:cost]
            edge_routing = Dict( collect(edges(v_subgraph.graph))[1] => Path(src(s_edge), dst(s_edge), [s_edge], s_network_dir[src(s_edge), dst(s_edge)][:cost]))
            sub_mapping = Mapping(v_subgraph.graph, s_network_dir, node_placement, edge_routing)
            add_column(master_problem, instance, vn_decompo, v_subgraph, sub_mapping, cost_mapping)
        end

    end


    # column generation!
    column_generation(instance, vn_decompo, master_problem)


    # ======= END HEURISTIC STUFF ======= #

    basic_heuristic(instance, vn_decompo, master_problem, 900)

    return 
end




function start_partitition(instance)


    v_node_partitionning = []
    copy_v_network = copy(instance.v_network.graph)
    
    keep_on = true
    while keep_on

        # Get node with max degree
        node_with_max_degree = 1

        keep_on = false
        if degree(copy_v_network, node_with_max_degree) >= 0
            keep_on = true
        end


        # remove it from the network

        # need to take care of the numbering of nodes in the copy!
    end



end


