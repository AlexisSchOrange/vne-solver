
using Graphs, MetaGraphsNext
includet("../utils/import_utils.jl")
includet("optimal_routing.jl")
includet("shortest-path-routing.jl")



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

