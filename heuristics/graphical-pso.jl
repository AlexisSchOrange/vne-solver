using Graphs, MetaGraphsNext
using Statistics
using  NetworkLayout


includet("../utils/import_utils.jl")





function graphical_pso(instance)




    function shortest_path_routing(v_node_placement, get_routing=false)

        capacities_edges_copy = copy(capacities_edges)
        s_network_copy_dir_copy_ofgraph = nothing
        is_still_original_s_network = true

        edge_routing = Dict{Edge, Path}()
        overall_routing_costs = 0

        for v_edge in edges(v_network)

            s_src = v_node_placement[src(v_edge)]
            s_dst = v_node_placement[dst(v_edge)]

            nodes_of_path = Int[]
            cost_of_routing_current_edge = 0

            # Get the shortest path
            if is_still_original_s_network
                nodes_of_path = [s_dst]
                v = s_dst
                while v != s_src
                    u = shortest_paths.parents[s_src, v]
                    push!(nodes_of_path, u)
                    v = shortest_paths.parents[s_src, v]
                end
                reverse!(nodes_of_path)
                cost_of_routing_current_edge = shortest_paths.dists[s_src, s_dst]
            else # A smarter thing to do here, would be to use the basic path, and if it's using a removed edge, then compute the astar path. But it does not cost so much time so.
                edges_of_path = a_star(s_network_copy_dir_copy_ofgraph, s_src, s_dst, distmx)
                if edges_of_path == [] # No paths found: the graph is full!
                    return Dict(), 10e9
                end
                push!(nodes_of_path, src(edges_of_path[1]))
                for edge in edges_of_path
                    push!(nodes_of_path, dst(edge))
                    cost_of_routing_current_edge += distmx[src(edge), dst(edge)]
                end
            end
            
            overall_routing_costs += cost_of_routing_current_edge

            # Removing the capacities. If no more capacities, removing the edge.
            for i_node in 1:length(nodes_of_path)-1
                current_src = nodes_of_path[i_node]
                current_dst = nodes_of_path[i_node+1]
                capacities_edges_copy[(current_src), (current_dst)] -= 1
                capacities_edges_copy[(current_dst), (current_src)] -= 1

                if capacities_edges_copy[(current_src), (current_dst)] == 0
                    if is_still_original_s_network
                        s_network_copy_dir_copy_ofgraph = deepcopy(instance.s_network_dir.graph)
                        is_still_original_s_network = false
                    end
                    rem_edge!(s_network_copy_dir_copy_ofgraph, current_src, current_dst)
                    rem_edge!(s_network_copy_dir_copy_ofgraph, current_dst, current_src)
                end
            end

            # Sometime it is necessary to get the actual routing, but not all the time. 
            # Since it takes some time to retrieve the edges, due to poor design choices, I only do it when necessary.
            if get_routing
                edges_of_path = Edge[]
                for i_node in 1:length(nodes_of_path)-1
                    current_src = nodes_of_path[i_node]
                    current_dst = nodes_of_path[i_node+1]
                    push!(edges_of_path, get_edge(s_network_dir, current_src, current_dst))
                end
                path = Path(s_src, s_dst, edges_of_path, cost_of_routing_current_edge)
                edge_routing[v_edge] = path
            end
        end

        return edge_routing, overall_routing_costs
    end



    function construct_placement(pos)


        # get closest substrate node to the pos:
        distances = [sqrt((pos[1]-x)^2 + (pos[2]-y)^2) for (x,y) in coords_s_nodes]
        central_s_node = argmin(distances)

        # Useful lists!
        placement = zeros(Int, nv(v_network))
        placement[most_central_v_node] = central_s_node

        already_placed_v_nodes = Int[]
        push!(already_placed_v_nodes, most_central_v_node)


        next_v_nodes = Int[]
        for v_node in already_placed_v_nodes
            for v_neighbor in neighbors_vn_memory[v_node]
                if v_neighbor ∉ next_v_nodes && v_neighbor ∉ already_placed_v_nodes
                    push!(next_v_nodes, v_neighbor)
                end
            end
        end

        while !isempty(next_v_nodes)

            # Take a node of the list
            v_node = popfirst!(next_v_nodes)

            # Get neighbors already placed
            placement_neighbors = Int[]
            for v_neighbor in neighbors_vn_memory[v_node]
                if v_neighbor ∈ already_placed_v_nodes
                    push!(placement_neighbors, placement[v_neighbor])
                end
            end

            # Choose 10 random nodes (or less if the network is not that big)
            number_s_nodes = min(floor(nv(s_network)/5), 10)
            some_s_nodes = Int[]
            while length(some_s_nodes) < number_s_nodes
                s_node = rand(1:nv(s_network))
                if (node_capacities[s_node] >= 1) && (s_node ∉ placement) && (s_node ∉ some_s_nodes)
                    push!(some_s_nodes, s_node)
                end
            end

            # Rank them according distance to already placed nodes and capacity
            distances = [ sum(shortest_paths.dists[s_src, s_node] for s_src in placement_neighbors) for s_node in some_s_nodes]
            distances_norm = (distances .- minimum(distances)) ./ (maximum(distances) - minimum(distances) + 1e-9)

            #capacities = [ s_node_scores[s_node] for s_node in some_s_nodes]
            #capacities_norm = (capacities .- minimum(capacities)) ./ (maximum(capacities) - minimum(capacities) + 1e-9)

            costs = [ node_costs[s_node] for s_node in some_s_nodes]
            costs_norm = (costs .- minimum(costs)) ./ (maximum(costs) - minimum(costs) + 1e-9)

            final_scores = 1. .* distances_norm .+ 0.5 .* costs_norm #.+ 0.5 .* ( 1 .-capacities_norm)

            selected_idx = argmin(final_scores)
            s_node_selected = some_s_nodes[selected_idx]


            # Finish the work
            placement[v_node] = s_node_selected
            push!(already_placed_v_nodes, v_node)
            for v_neighbor in neighbors_vn_memory[v_node]
                if v_neighbor ∉ already_placed_v_nodes && v_neighbor ∉ next_v_nodes 
                    push!(next_v_nodes, v_neighbor)
                end
            end

        end

        placement_cost = 0
        for v_node in vertices(v_network)
            placement_cost += node_costs[placement[v_node]] 
        end

        return placement, placement_cost


    end




    println("FINEEEEE, let's do the graphical thing.")
    time_beginning = time()

    
    time_beginning = time()

    v_network = instance.v_network
    s_network = instance.s_network
    s_network_dir = instance.s_network_dir


    #---- Make sure there are enough capacited nodes
    nodes_with_caps = 0
    for s_node in vertices(s_network)
        if s_network[s_node][:cap] >= 1
            nodes_with_caps += 1
        end
    end
    if nodes_with_caps < nv(v_network)
        println("What the hell? Not enough capacited nodes...")
        return 10e9
    end

        
    #---- Usefull things
    node_capacities = [get_attribute_node(s_network, s_node, :cap) for s_node in vertices(s_network)]
    node_costs = [get_attribute_node(s_network, s_node, :cost) for s_node in vertices(s_network)]
    capacities_edges = Dict{Tuple{Int, Int}, Int}()
    for s_edge in edges(s_network)
        cap = s_network[src(s_edge), dst(s_edge)][:cap]
        capacities_edges[(src(s_edge),dst(s_edge))] = cap
        capacities_edges[(dst(s_edge),src(s_edge))] = cap
    end

    neighbors_sn_memory = []
    for s_node in vertices(s_network)
        push!( neighbors_sn_memory, neighbors(s_network, s_node))
    end

    neighbors_vn_memory = []
    for v_node in vertices(v_network)
        push!( neighbors_vn_memory, neighbors(v_network, v_node))
    end

    # Greedy score based on capacities for nodes: not for now.
    #=
    s_node_scores = [ node_capacities[s_node]  * 
                        sum(capacities_edges[(s_node,s_neighbor)] for s_neighbor in neighbors_sn_memory[s_node]) 
                    for s_node in vertices(s_network)]
    =#

    # shortest path of the substrate network
    distmx = zeros(Int, nv(s_network), nv(s_network))
    for s_edge in edges(s_network_dir)
        distmx[src(s_edge), dst(s_edge)] = get_attribute_edge(s_network_dir, s_edge, :cost)
    end
    shortest_paths = floyd_warshall_shortest_paths(s_network_dir, distmx)


    # Getting the most central virtual node
    most_central_v_node = argmin(closeness_centrality(v_network))


    # Algos for layout: stress, spring, 
    coords_s_nodes = stress(instance.s_network; iterations=500)
    coords_s_nodes = normalize_coords(coords_s_nodes)

    pos_particles = []
    vel_particles = []
    best_pos_particles = []
    best_cost_particles = []

    history_pos = []

    best_pos_overall = (0., 0.)
    best_cost_overall = 10e9

    nb_particle = 25
    nb_iter = 5000

    # initialization
    for particle in 1:nb_particle

        # put a random position
        pos = (rand()*2-1, rand()*2-1)
        vel = ((rand()*2-1)*0.15, (rand()*2-1)*0.15)

        node_placement, placement_cost = construct_placement(pos)
        edge_routing, routing_cost = shortest_path_routing(node_placement)

        total_cost = placement_cost + routing_cost
        if total_cost < best_cost_overall
            best_pos_overall = pos
            best_cost_overall = total_cost
        end

        push!(pos_particles, pos)
        push!(vel_particles, vel)    
        
        push!(best_pos_particles, pos)
        push!(best_cost_particles, total_cost)

        push!(history_pos, [pos])
    end

    println("Best after init: $best_cost_overall")

    for iter in 1:nb_iter

        for particle in 1:nb_particle

            current_vel = vel_particles[particle]
            current_pos = pos_particles[particle]
            
            new_vel =   0.7 .* current_vel .+ 
                        1.5 * rand() .* (best_pos_particles[particle] .- current_pos) .+ 
                        1.5 * rand() .* (best_pos_overall .- current_pos);
            
            new_vel = clamp.(new_vel, -0.10, 0.10)

            new_pos = new_vel .+ current_pos
            new_pos = clamp.(new_pos, -1, 1.)

            node_placement, placement_cost = construct_placement(new_pos)
            edge_routing, routing_cost = shortest_path_routing(node_placement)
    
            total_cost = placement_cost + routing_cost
            if total_cost < best_cost_particles[particle]
                best_pos_particles[particle] = new_pos
                best_cost_particles[particle] = total_cost
            end
            if total_cost <= best_cost_overall
                best_pos_overall = new_pos
                best_cost_overall = total_cost
                println("New best found! $best_cost_overall, at iter $iter")
            end
    
            pos_particles[particle] = new_pos
            vel_particles[particle] = new_vel
            
            push!(history_pos[particle], new_pos)

        end

    end


    print_history_pos(history_pos, nb_particle, nb_iter)

    println("Final best cost: $best_cost_overall, took me $(time() - time_beginning)")

    return best_cost_overall
end



function normalize_coords(coords)
    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)

    return [(2*(x - xmin)/(xmax - xmin) - 1,
             2*(y - ymin)/(ymax - ymin) - 1) for (x,y) in coords]
end


function print_history_pos(history_pos, nb_particle, nb_iter)

    #=
    for particle in 1:nb_particle

        println("For particle $particle:")
        for pos in history_pos[particle]
            println("( $(round(pos[1], digits=3)), $(round(pos[2], digits=3)))")
        end

    end
    =#
    println("End position:")
    for particle in 1:nb_particle
        pos = history_pos[particle][nb_iter]
        println("( $(round(pos[1], digits=3)), $(round(pos[2], digits=3)))")
    end

end

