
using JuMP, CPLEX, Gurobi
using Polyhedra, CDDLib, XPORTA


includet("../utils/import_utils.jl")



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
    
    
    
    # Simple path constraints, only useful for porta.
    # Note that non-simple path and subtours are possible with the formulation, 
    # but will never appear in practice due to being expensive for nothing.
    for s_node in vertices(instance.s_network)
        for v_node in vertices(v_network)
            for v_edge in get_out_edges(v_network, v_node)
                @constraint(model, 
                    sum(y[v_edge, s_edge] for s_edge in get_in_edges(s_network_dir, s_node)) 
                    <= 1 - x[v_node, s_node] )
            end
        end
    end
    # to remove loops..
    for v_edge in edges(v_network)
        for s_edge in edges(instance.s_network)
            @constraint(model, y[v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))] 
                + y[v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))] 
                <= 1 )
        end
    end

    
    
end



########### GET SOLLSSSSS

function get_sols(instance)

    model = Model(Gurobi.Optimizer)
    set_silent(model)
    s_network_dir = instance.s_network_dir

    set_up_problem(instance, model)
    set_optimizer_attribute(model, "PoolSearchMode", 2)
    set_optimizer_attribute(model, "PoolSolutions", 1000)
    
    print("Model set, starting to solve... ")

    optimize!(model)
    #solution_summary(model)
    println("Done, there are $(result_count(model)) solutions")
    sols = Vector{Vector{Int64}}()
    for i in 1:result_count(model)
        cur_sol = Vector{Int64}()
        for v_node in vertices(instance.v_network)
            for s_node in vertices(s_network_dir)
                push!(cur_sol, round.(value.(model[:x][v_node, s_node]; result = i)))     
            end
        end

        for v_edge in edges(instance.v_network)
            for s_edge in edges(s_network_dir)
                push!(cur_sol, round.(value.(model[:y][v_edge, s_edge]; result = i)))     
            end
        end

        push!(sols, cur_sol)
    end

    name_variables = []
    for v_node in vertices(instance.v_network)
        for s_node in vertices(instance.s_network)
            push!(name_variables, "x_" * string(v_node) * "_" * string(s_node))
        end
    end
    
    for v_edge in edges(instance.v_network)
        for s_edge in edges(s_network_dir)
            push!(name_variables, "y_" * string(src(v_edge)) * string(dst(v_edge)) * "_" * string(src(s_edge)) * string(dst(s_edge)))
        end
    end

    return sols, name_variables
end

function get_sols_simplified(instance)

    model = Model(Gurobi.Optimizer)
    set_silent(model)

    s_network_dir = instance.s_network_dir
    set_up_problem(instance, model)
    set_optimizer_attribute(model, "PoolSearchMode", 2)
    set_optimizer_attribute(model, "PoolSolutions", 10000)
    
    print("Model set, starting to solve... ")

    optimize!(model)
    #solution_summary(model)
    println("Done, there are $(result_count(model)) solutions")
    sols = Vector{Vector{Int64}}()
    for i in 1:result_count(model)
        cur_sol = Vector{Int64}()
        for v_node in vertices(instance.v_network)
            for s_node in vertices(s_network_dir)
                push!(cur_sol, round.(value.(model[:x][v_node, s_node]; result = i)))     
            end
        end

        for v_edge in edges(instance.v_network)
            for s_edge in edges(instance.s_network)
                push!(cur_sol, round.(
                    value.(model[:y][v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))]; result = i) 
                    + value.(model[:y][v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]; result = i) ) )     
            end
        end

        push!(sols, cur_sol)
    end

    name_variables = []
    for v_node in vertices(instance.v_network)
        for s_node in vertices(instance.s_network)
            push!(name_variables, "x_" * string(v_node) * "_" * string(s_node))
        end
    end
    
    for v_edge in edges(instance.v_network)
        for s_edge in edges(instance.s_network)
            push!(name_variables, "y_" * string(src(v_edge)) * string(dst(v_edge)) * "_" * string(src(s_edge)) * string(dst(s_edge)))
        end
    end

    return sols, name_variables



end


function get_sols_simplified_test(instance)

    model = Model(Gurobi.Optimizer)
    set_silent(model)

    s_network_dir = instance.s_network_dir
    set_up_problem(instance, model)
    set_optimizer_attribute(model, "PoolSearchMode", 2)
    set_optimizer_attribute(model, "PoolSolutions", 10000)
    

    # INEQUALTIES TO BREAK SIMETRY ?
    @constraint(model, model[:x][1, 1] == 1)


    print("Model set, starting to solve... ")

    optimize!(model)
    #solution_summary(model)
    println("Done, there are $(result_count(model)) solutions")
    sols = Vector{Vector{Int64}}()
    for i in 1:result_count(model)
        cur_sol = Vector{Int64}()
        for v_node in vertices(instance.v_network)
            for s_node in vertices(s_network_dir)
                push!(cur_sol, round.(value.(model[:x][v_node, s_node]; result = i)))     
            end
        end

        for v_edge in edges(instance.v_network)
            for s_edge in edges(instance.s_network)
                push!(cur_sol, round.(
                    value.(model[:y][v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))]; result = i) 
                    + value.(model[:y][v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]; result = i) ) )     
            end
        end

        push!(sols, cur_sol)
    end

    name_variables = []
    for v_node in vertices(instance.v_network)
        for s_node in vertices(instance.s_network)
            push!(name_variables, "x_" * string(v_node) * "_" * string(s_node))
        end
    end
    
    for v_edge in edges(instance.v_network)
        for s_edge in edges(instance.s_network)
            push!(name_variables, "y_" * string(src(v_edge)) * string(dst(v_edge)) * "_" * string(src(s_edge)) * string(dst(s_edge)))
        end
    end

    return sols, name_variables



end



function get_sols_gigatest(instance)

    model = Model(Gurobi.Optimizer)
    #set_silent(model)

    set_up_problem(instance, model)
    # BREAK THE SYMETRIES ?
    println("What happens with the linear relax ?")
    unrelax = relax_integrality(model)
    optimize!(model)
    unrelax()


    set_optimizer_attribute(model, "PoolSearchMode", 2)
    set_optimizer_attribute(model, "PoolSolutions", 1000)
    
    print("Model set, starting to solve... ")
    
    optimize!(model)
    #solution_summary(model)
    println("Done, there are $(result_count(model)) solutions")
    sols = Vector{Vector{Int64}}()
    for i in 1:result_count(model)
        cur_sol = Vector{Int64}()
        for v_node in vertices(instance.v_network)
            for s_node in vertices(s_network_dir)
                push!(cur_sol, round.(value.(model[:x][v_node, s_node]; result = i)))     
            end
        end

        for v_edge in edges(instance.v_network)
            for s_edge in edges(instance.s_network)
                push!(cur_sol, round.(
                    value.(model[:y][v_edge, get_edge(s_network_dir, src(s_edge), dst(s_edge))]; result = i) 
                    + value.(model[:y][v_edge, get_edge(s_network_dir, dst(s_edge), src(s_edge))]; result = i) ) )     
            end
        end

        push!(sols, cur_sol)
    end

    name_variables = []
    for v_node in vertices(instance.v_network)
        for s_node in vertices(instance.s_network)
            push!(name_variables, "x_" * string(v_node) * "_" * string(s_node))
        end
    end
    
    for v_edge in edges(instance.v_network)
        for s_edge in edges(instance.s_network)
            push!(name_variables, "y_" * string(src(v_edge)) * string(dst(v_edge)) * "_" * string(src(s_edge)) * string(dst(s_edge)))
        end
    end

    return sols, name_variables



end


############### HREP

function get_hrep(sols, names)
    v_rep = vrep(sols);
    poly = polyhedron(v_rep, XPORTA.Library(:float));
    println("There is $(length(sols)) solutions and $(length(names)) variables in your instance.")
    println("Computing the hrep... This might take some time... ")
    h_rep = hrep(poly);
    print("finished.")
    return h_rep
end



function get_dominant_hrep(sols, names)
    rays = []
    for i_var in 1:length(names)
        ray = zeros(Int64, length(names))
        ray[i_var] = 1
        push!(rays, ray)
    end
    v_rep = convexhull(sols...) + conichull(rays...)
    poly = polyhedron(v_rep, CDDLib.Library(:float));
    #print(poly)

    println("There is $(length(sols)) solutions and $(length(names)) variables in your instance.")
    println("Computing the hrep... This might take some time... ")
    h_rep = hrep(poly);

    return(h_rep)

end


# PRINTING

function print_solutions(sols, names)

    for (i_sol, sol) in enumerate(sols)
        println("Solution num $(i_sol)")
        for (i_var, value) in enumerate(sol)
            if value > 0.5
                println("$(names[i_var])")
            end
        end
    end
end



function print_polytope(hr, names_variables)
    println("There are " * string(length(names_variables)) * " variables.\n")
    println("There are " * string(length(hyperplanes(hr))) * " hyperplanes")
    for h in hyperplanes(hr)
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                print("+ ")
                print(floor(Int, h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                print("\t")
            elseif (h.a[i_var] < -0.0001)
                print("- ")
                print(floor(Int, -h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                print("\t")
            else
                print("\t\t")
            end
        end
        print(" = ")
        println(floor(Int,h.β))
    end

    println("\n\n There are " * string(length(halfspaces(hr))) * " halfspaces")

    for h in halfspaces(hr)
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                print("+ ")
                print(floor(Int, h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                print("\t")
            elseif (h.a[i_var] < -0.0001)
                print("- ")
                print(floor(Int, -h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                print("\t")
            else
                print("\t\t")
            end
        end
        print(" ≤ ")
        println(floor(Int,h.β))
    end

end


function print_polytope_simpler(hr, names_variables, print_index=false)
    println("There are " * string(length(names_variables)) * " variables.\n")
    println("There are " * string(length(hyperplanes(hr))) * " hyperplanes")
    for h in hyperplanes(hr)
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                print("+ ")
                print(floor(Int, h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            elseif (h.a[i_var] < -0.0001)
                print("- ")
                print(floor(Int, -h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            end
        end
        print(" = ")
        println(floor(Int,h.β))
    end

    println("\n\n There are " * string(length(halfspaces(hr))) * " halfspaces")

    for h in halfspaces(hr)
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                print("+ ")
                print(floor(Int, h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            elseif (h.a[i_var] < -0.0001)
                print("- ")
                print(floor(Int, -h.a[i_var]))
                print(" ")
                print(names_variables[i_var])
                if print_index
                    print("[" * string(i_var) * "]")
                end
                print("\t")
            end
        end
        print(" ≤ ")
        println(floor(Int,h.β))
    end



end


# todo : order the halfspace by length of nb on the left, so that's its easy to read ?
function print_polytope_better(hr, names_variables)
    println("There are " * string(length(names_variables)) * " variables.\n")
    println("There are " * string(length(hyperplanes(hr))) * " hyperplanes")

    for h in hyperplanes(hr)
        str_right = ""
        str_left = ""
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                str_left *= "+ $(floor(Int, h.a[i_var]))$(names_variables[i_var])  "
            elseif (h.a[i_var] < -0.0001)
                str_right *= "+ $(floor(Int, -h.a[i_var]))$(names_variables[i_var])  "
            end
        end

        if  h.β > 0.00001
            str_right = " $(floor(Int,h.β)) " * str_right
        elseif h.β < - 0.00001
            str_left = " $(floor(Int,-h.β)) " * str_left
        end
        
        if length(str_left) == 0
            str_left = " 0 "
        end
        if length(str_right) == 0
            str_right = " 0 "
        end

        print(str_left)
        print(" = ")
        println(str_right)
    end

    println("\n\n There are " * string(length(halfspaces(hr))) * " halfspaces")
    for h in halfspaces(hr)
        str_right = ""
        str_left = ""
        for i_var in 1:length(names_variables)
            if (h.a[i_var] > 0.0001)
                str_left *= "+ $(floor(Int, h.a[i_var]))$(names_variables[i_var])  "
            elseif (h.a[i_var] < -0.0001)
                str_right *= "+ $(floor(Int, -h.a[i_var]))$(names_variables[i_var])  "
            end
        end

        if  h.β > 0.00001
            str_right = " $(floor(Int,h.β)) " * str_right
        elseif h.β < - 0.00001
            str_left = " $(floor(Int,-h.β)) " * str_left
        end
        
        if length(str_left) == 0
            str_left = " 0 "
        end
        if length(str_right) == 0
            str_right = " 0 "
        end
        print(str_right)
        print(" ≥ ")
        println(str_left)
    end



end



# 13 11 : trying to do star on cycle, but it explodes really fast.
function print_solutions_test(sols, names)

    i_count = 0
    i_var_to_watch = 3

    for (i_sol, sol) in enumerate(sols)
        if sol[i_var_to_watch] == 1     
            i_count += 1
            println("Solution num $(i_sol)")
            for (i_var, value) in enumerate(sol)
                if value > 0.5
                    println("$(names[i_var])")
                end
            end
        end
    end

    println("Y'en a eu $i_count pour $(names[i_var_to_watch])")
end

# trying to remove some symetries.



# AYAAAAA


function print_constraints(hr, instance)

    for h in halfspaces(hr)
        i_var = 1

        print("@constraint(model, ")
        
        for v_node in vertices(instance.v_network)
            for s_node in vertices(instance.s_network)
                if abs((h.a[i_var])) > 0.01
                    print("+ $(floor(Int, h.a[i_var])) * model[:x][$v_node, $s_node] ")
                end
                i_var += 1
            end
        end
        for v_edge in edges(instance.v_network)
            for s_edge in edges(instance.s_network_dir)
                if (abs(h.a[i_var]) > 0.01)
                    print("+ $(floor(Int, h.a[i_var])) * model[:y][get_edge(instance.v_network, $(src(v_edge)), $(dst(v_edge))), get_edge(instance.s_network_dir, $(src(s_edge)), $(dst(s_edge)))] ")
                end
                i_var += 1
            end
        end
        print(" <= +")
        print(floor(Int,h.β))
        println(")")


    end

end


