
using Graphs, MetaGraphsNext
includet("../utils/import_utils.jl")
includet("shortest-path-routing.jl")
includet("UEPSO.jl")
includet("../utils/kahip_wrapper.jl")


struct Subgraph
    graph
    nodes_of_main_graph
end



# Based on ODEA, but without overlapping decomposition.
# Only 10% ! It's crazy !
function dea(instance)

    v_network = instance.v_network
    s_network = instance.s_network

    # inialize a mapping randomly
    best_placement = []
    placement_cost = 0
    for v_node in 1:nv(v_network)
        keep_on = true
        while keep_on
            s_node = rand(1:nv(s_network))
            if s_node ∉ best_placement && get_attribute_node(s_network, s_node, :cap)>0
                push!(best_placement, s_node)
                keep_on=false
                placement_cost+= get_attribute_node(s_network, s_node, :cost)
            end
        end
    end
    routing, routing_cost = shortest_path_routing(instance, best_placement)
    overall_cost = placement_cost + routing_cost
    
    best_placement_cost = 999999
    if overall_cost  < best_placement_cost
        best_placement_cost = overall_cost
    end
    # virtual graph partition
    nb_cluster = round(Int, nv(v_network)/8)
    print("Doing a partition into $nb_cluster subgraphs...")
    vn_partition = partition_vn(instance, nb_cluster)

    vn_subgraphs = []
    for (i_subgraph, v_nodes) in enumerate(vn_partition)
        subgraph = Subgraph(my_induced_subgraph(v_network, v_nodes, "subgraph_$i_subgraph"), v_nodes)
        push!(vn_subgraphs, subgraph)
    end


    # let's gongue
    nb_iter = 50
    for i in 1:nb_iter
        print("iter $i... ")
        for v_subgraph in vn_subgraphs
            
            # get forbidden nodes: nodes already used by a virtual node outside of the subgraph
            forbidden_s_nodes = []
            for v_node in vertices(v_network)
                if v_node ∉ v_subgraph.nodes_of_main_graph
                    push!(forbidden_s_nodes, best_placement[v_node])
                end
            end

            subinstance = Instance(v_subgraph.graph, s_network)
            mapping, cost = solve_UEPSO(subinstance, forbidden_s_nodes)
            #print("I've obtained a mapping of cost $cost")
            new_placement = []
            node_assignment = 1
            placement_cost = 0
            for v_node in vertices(v_network)
                if v_node ∈ v_subgraph.nodes_of_main_graph
                    s_node = mapping.node_placement[node_assignment]
                    placement_cost+= get_attribute_node(s_network, s_node, :cost)
                    push!(new_placement, s_node)
                    node_assignment += 1
                else
                    s_node = best_placement[v_node]
                    placement_cost+= get_attribute_node(s_network, s_node, :cost)
                    push!(new_placement, s_node)
                end
            end

            routing, routing_cost = shortest_path_routing(instance, best_placement)
            overall_cost = placement_cost + routing_cost

            if overall_cost < 9999
                println("Well we got a solution of $overall_cost")
            end
            if overall_cost < best_placement_cost
                best_placement = new_placement
                best_placement_cost = overall_cost
                println("We got a new best at $overall_cost !")
            end
        end
    end


end


function partition_vn(instance, nb_clusters)

    graph = instance.v_network.graph

    # 1 : Partitionner
    inbalance = 0.10
    println("$nb_clusters clusters to do, with inbalance $inbalance...")
    partition = partition_kahip(graph, nb_clusters, inbalance)
    clusters = [Vector{Int64}() for i in 1:nb_clusters]
    for s_node in vertices(graph)
        push!(clusters[partition[s_node]], s_node)
    end

    # 2 : Corriger

    for cluster in clusters
        simple_subgraph, vmap = induced_subgraph(graph, cluster)
        if !is_connected(simple_subgraph)
            components = connected_components(simple_subgraph)
            component_sorted = sort(components, by=x->length(x), rev=true)
            for subcluster in component_sorted[2:length(component_sorted)]
                #Let's add all those nodes to a (most) connected subgraph
                nodes_original = [vmap[node] for node in subcluster]
                subgraph_neighbors = zeros(Int, nb_clusters)
                for node in nodes_original
                    for neighbor in neighbors(graph, node)
                        if neighbor ∉ cluster
                            subgraph_neighbors[partition[neighbor]] += 1
                        end
                    end
                end
                most_connected_subgraph = sortperm(subgraph_neighbors, rev=true)
                append!(clusters[most_connected_subgraph[1]], nodes_original)
                #println("Well let's add $nodes_original to cluster $(clusters[most_connected_subgraph[1]])")
                filter!(e->e∉nodes_original, cluster)
            end
        end
    end


    return clusters

end






function solve_UEPSO(instance, forbiden_s_nodes=[])

    v_network = instance.v_network
    s_network = instance.s_network

    # PSO parameters
    nb_particle=25  
    nb_iter=40



    position = []
    velocity = []


    personal_best = []
    personal_best_cost = []

    global_best = [i for i in 1:nv(v_network)]
    global_best_cost = 99999



    # initialization
    for particle in 1:nb_particle

        placement = []
        placement_cost=0
        for i in 1:nv(v_network)
            keep_on = true
            while keep_on
                s_node = rand(1:nv(s_network))
                if s_node ∉ placement && get_attribute_node(s_network, s_node, :cap)>0 && s_node ∉ forbiden_s_nodes
                    push!(placement, s_node)
                    keep_on=false
                    placement_cost+= get_attribute_node(s_network, s_node, :cost)
                end
            end
        end

        routing, routing_cost = shortest_path_routing(instance, placement)

        overall_cost = placement_cost + routing_cost

        push!(position, placement)
        push!(personal_best, position[particle])
        push!(personal_best_cost, overall_cost)

        if overall_cost < global_best_cost
            global_best = position[particle]
            global_best_cost = overall_cost
            #println("We got a new best solution! value $overall_cost")
        end

        push!(velocity, ones(nv(v_network)))
    end



    # iterations
    for iter in 1:nb_iter
        for particle in 1:nb_particle

            velocity[particle] = plus( velocity[particle], 
                                        minus(personal_best[particle], position[particle]), 
                                        minus(global_best, position[particle]))
            position[particle], placement_cost = times(position[particle], velocity[particle], instance, forbiden_s_nodes)

            #println("Current position: $(position[particle])")
            routing, routing_cost = shortest_path_routing(instance, position[particle])

            overall_cost = placement_cost + routing_cost
            if overall_cost < personal_best_cost[particle]
                personal_best[particle] = position[particle]
                personal_best_cost[particle] = overall_cost
            end
            if overall_cost < global_best_cost
                global_best_cost = overall_cost
                global_best = position[particle]
                #println("We got a new best solution! value $overall_cost")
            end
        end
    end

    #println("Final best solution: $global_best")
    routing, routing_cost = shortest_path_routing(instance, global_best)
    final_mapping = Mapping(v_network, s_network, global_best, routing)

    return final_mapping, global_best_cost
end




function minus(pos1, pos2)

    res=[]
    for i in 1:length(pos1)
        if pos1[i] == pos2[i]
            push!(res, 1)
        else
            push!(res, 0)
        end
    end
    return res
end

function plus(vel_inertia, vel_pb, vel_gb)

    p_inertia = 0.1
    p_attraction_personal = 0.2
    p_attraction_global = 0.7

    new_velocity = []
    for i in 1:length(vel_inertia)
        r = rand()
        if r < p_inertia
            push!(new_velocity, vel_inertia[i])
        elseif r < (p_inertia + p_attraction_personal)
            push!(new_velocity, vel_pb[i])
        else
            push!(new_velocity, vel_gb[i])
        end
    end

    return new_velocity

end


function times(position, velocity, instance, forbiden_s_nodes)

    new_placement = []
    placement_cost = 0

    for i in 1:nv(instance.v_network)
        if velocity[i] == 1
            push!(new_placement, position[i])
            placement_cost += get_attribute_node(instance.s_network, position[i], :cost)
        else
            push!(new_placement, -1)
        end
    end

    for i in 1:nv(instance.v_network)
        if new_placement[i] == -1
            keep_on = true
            while keep_on
                s_node = rand(1:nv(instance.s_network))
                if s_node ∉ new_placement && get_attribute_node(instance.s_network, s_node, :cap)>0 && s_node ∉ forbiden_s_nodes
                    new_placement[i]=s_node
                    keep_on=false
                    placement_cost+= get_attribute_node(instance.s_network, s_node, :cost)
                end
            end
        end
    end

    return new_placement, placement_cost
end




