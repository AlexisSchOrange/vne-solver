
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



    # generating first columns. For edge...
    #=
    for s_edge in edges(s_network_dir)

        for v_subgraph in vn_subgraphs
            node_placement = [src(s_edge), dst(s_edge)]
            cost_mapping = s_network_dir[src(s_edge)][:cost] + s_network_dir[dst(s_edge)][:cost] + s_network_dir[src(s_edge), dst(s_edge)][:cost]
            edge_routing = Dict( collect(edges(v_subgraph.graph))[1] => Path(src(s_edge), dst(s_edge), [s_edge], s_network_dir[src(s_edge), dst(s_edge)][:cost]))
            sub_mapping = Mapping(v_subgraph.graph, s_network_dir, node_placement, edge_routing)
            add_column(master_problem, instance, vn_decompo, v_subgraph, sub_mapping, cost_mapping)
        end

    end
    =#

    # column generation!
    column_generation(instance, vn_decompo, master_problem)


    # ======= END HEURISTIC STUFF ======= #

    basic_heuristic(instance, vn_decompo, master_problem, 900)

    return 
end




function stars_partition(instance)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    subgraphs = []
    copy_v_network = copy(v_network.graph)
    keep_on = true

    real_indices = collect(1:nv(v_network)) 
    # Graphs.jl is completly stupid when it comes to removing a node. It swaps it with the last node and then removes it. be careful... 

    nodes_in_no_subgraphs = collect(1:nv(v_network))
    v_node_partitionning = []
    while keep_on

        # Get node with max degree
        node_with_max_degree = argmax(degree(copy_v_network))

        keep_on = true
        if degree(copy_v_network, node_with_max_degree) == 0
            break
        end

        #println("most central node: $node_with_max_degree, aka $(real_indices[node_with_max_degree])")

        new_part = [ real_indices[node_with_max_degree] ]

        edges_subg = []
        for neigh in neighbors(copy_v_network, node_with_max_degree)
            #println("Neighbor: $neigh, aka $(real_indices[neigh])")
            push!(new_part, real_indices[neigh])
            push!(edges_subg, get_edge(v_network, real_indices[node_with_max_degree], real_indices[neigh]))
        end

        push!(v_node_partitionning, new_part)
        
        #push!(subgraphs, Dict("nodes"=>new_part, "edges"=>edges_subg))

        for v_node in new_part
            if v_node ∈ nodes_in_no_subgraphs
                filter!(x -> x != v_node, nodes_in_no_subgraphs)
            end
        end

        rem_vertex!(copy_v_network, node_with_max_degree)
        real_indices[node_with_max_degree] = real_indices[length(real_indices)]
        deleteat!(real_indices, length(real_indices))

        
        #println("Real indices: $real_indices")

    end


    #println("Some nodes left? $nodes_in_no_subgraphs")
    #println("Node partitionning: $v_node_partitionning")

    for v_node_left in nodes_in_no_subgraphs
        push!(v_node_partitionning, [v_node_left])
    end

    println("Node partitionning: $v_node_partitionning")





    vn_decompo = set_up_decompo_overlapping(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges: $(vn_decompo.v_edges_master)")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    

    # === COLUMN GENERATION === #

    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")

    # column generation!
    column_generation(instance, vn_decompo, master_problem)

end




# A bit less stars, with a bit less overlapping nodes?
# The center of a star can not be in another star.
function stars_partition_2(instance)

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    copy_v_network = copy(v_network.graph)
    keep_on = true

    real_indices = collect(1:nv(v_network)) 
    # Graphs.jl is completly stupid when it comes to removing a node. It swaps it with the last node and then removes it. be careful... 

    nodes_in_no_subgraphs = collect(1:nv(v_network))
    v_node_partitionning = []
    possible_centers = collect(1:nv(v_network))

    while keep_on

        # Get node with max degree
        nodes_degrees = [degree(copy_v_network, v_node) for v_node in vertices(copy_v_network)]
        node_sorted_degree = sortperm(nodes_degrees, rev=true)
        center = 0
        for node in node_sorted_degree
            if real_indices[node] ∈ possible_centers
                center = node
                break
            end
        end
        #println("The center will be $center")
        if center == 0
            #println("No center found?")
            break
        end
        if degree(copy_v_network, center) == 0
            #println("Degree too low: time to stop.")
            break
        end
        #println("New center! $center")
        #println("most central node: $node_with_max_degree, aka $(real_indices[node_with_max_degree])")

        new_part = [ real_indices[center] ]

        edges_subg = []
        for neigh in neighbors(copy_v_network, center)
            #println("Neighbor: $neigh, aka $(real_indices[neigh])")
            push!(new_part, real_indices[neigh])
            push!(edges_subg, get_edge(v_network, real_indices[center], real_indices[neigh]))
        end

        # Nodes from this star can not be centers anymore.
        for v_node in new_part
            filter!(x -> x != v_node, possible_centers)
        end

        push!(v_node_partitionning, new_part)
        
        for v_node in new_part
            if v_node ∈ nodes_in_no_subgraphs
                filter!(x -> x != v_node, nodes_in_no_subgraphs)
            end
        end

        rem_vertex!(copy_v_network, center)
        real_indices[center] = real_indices[length(real_indices)]
        deleteat!(real_indices, length(real_indices))

        
        #println("Nodes that can be center: $possible_centers")
        #println("Nodes that are still in the graph: $real_indices")
    end


    #println("Some nodes left? $nodes_in_no_subgraphs")
    #println("Node partitionning: $v_node_partitionning")

    for v_node_left in nodes_in_no_subgraphs
        push!(v_node_partitionning, [v_node_left])
    end

    println("Node partitionning: $v_node_partitionning")





    vn_decompo = set_up_decompo_overlapping(instance, v_node_partitionning)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges: $(vn_decompo.v_edges_master)")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    

    # === COLUMN GENERATION === #

    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")

    # column generation!
    column_generation(instance, vn_decompo, master_problem)


end


# Here, looking at real stars:
# Subgraphs are real stars, I don't take edges that are "wheel edges".
# More edges in master problem
function stars_real_2(instance)
    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir

    subgraphs = []
    copy_v_network = copy(v_network.graph)
    keep_on = true

    real_indices = collect(1:nv(v_network)) 
    # Graphs.jl is completly stupid when it comes to removing a node. It swaps it with the last node and then removes it. be careful... 

    nodes_in_no_subgraphs = collect(1:nv(v_network))
    v_node_partitionning = []
    possible_centers = collect(1:nv(v_network))

    while keep_on

        # Get node with max degree
        nodes_degrees = [degree(copy_v_network, v_node) for v_node in vertices(copy_v_network)]
        node_sorted_degree = sortperm(nodes_degrees, rev=true)
        center = 0
        for node in node_sorted_degree
            if real_indices[node] ∈ possible_centers
                center = node
                break
            end
        end
        #println("The center will be $center")
        if center == 0
            #println("No center found?")
            break
        end
        if degree(copy_v_network, center) == 0
            #println("Degree too low: time to stop.")
            break
        end

        #println("most central node: $node_with_max_degree, aka $(real_indices[node_with_max_degree])")

        new_part = [ real_indices[center] ]

        edges_subg = []
        for neigh in neighbors(copy_v_network, center)
            #println("Neighbor: $neigh, aka $(real_indices[neigh])")
            push!(new_part, real_indices[neigh])
            push!(edges_subg, get_edge(v_network, real_indices[center], real_indices[neigh]))
        end

        # Nodes from this star can not be centers anymore.
        for v_node in new_part
            filter!(x -> x != v_node, possible_centers)
        end

        
        push!(subgraphs, Dict("nodes"=>new_part, "edges"=>edges_subg))

        for v_node in new_part
            if v_node ∈ nodes_in_no_subgraphs
                filter!(x -> x != v_node, nodes_in_no_subgraphs)
            end
        end

        rem_vertex!(copy_v_network, center)
        real_indices[center] = real_indices[length(real_indices)]
        deleteat!(real_indices, length(real_indices))

        
        #println("Real indices: $real_indices")

    end


    println("Some nodes left? $nodes_in_no_subgraphs")
    #println("Node partitionning: $v_node_partitionning")




    vn_decompo = set_up_decompo_overlapping_more_info(instance, subgraphs)
    vn_subgraphs = vn_decompo.subgraphs

    println("Virtual network decomposition done:")
    print_stuff_subgraphs(v_network, vn_subgraphs)
    println("   and $(length(vn_decompo.v_edges_master)) cutting edges: $(vn_decompo.v_edges_master)")
    println("   and $(length(vn_decompo.overlapping_nodes)) overlapping nodes : $(vn_decompo.overlapping_nodes)")

    

    # === COLUMN GENERATION === #
    
    # master problem things
    master_problem = set_up_master_problem(instance, vn_decompo)
    print("Master problem set... ")

    # column generation!
    column_generation(instance, vn_decompo, master_problem)

end




# No overlappin nodes!
function strict_stars_partition(instance)

end



# paths !
# overlapping paths!
# maybe triangle and cycles ? and the reminder are paths ?
# It would be nice to use cycle base.