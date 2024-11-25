
using JuMP, CPLEX, Gurobi
using Polyhedra, CDDLib, XPORTA


includet("../utils/import_utils.jl")
includet("../resolution/undirected/compact_undir.jl")




function get_sols_undir(instance)

    model = Model(Gurobi.Optimizer)
    set_silent(model)

    s_network_dir = generate_dir_sn(instance)
    set_up_problem_undir_1vn_1t1(instance, s_network_dir, model)
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



function get_hrep(sols, names)
    v_rep = vrep(sols);
    poly = polyhedron(v_rep, XPORTA.Library(:float));
    println("There is $(length(sols)) solutions and $(length(names)) variables in your instance.")
    println("Computing the hrep... This might take some time... ")
    h_rep = hrep(poly);
    return poly
    print("finished.")
end


# This is experimental, should work because we are in undirected, so we can remove one of the variable for each substrate edge.
# BUT by doing so, you are not working on the right polytope. 
# But is it a problem ? Because the solution are still valid.
# Finding inequalities from this, might require to adapt them for the original formulation.
function get_sols_undir_simplified(instance)

    model = Model(Gurobi.Optimizer)
    set_silent(model)

    s_network_dir = generate_dir_sn(instance)
    set_up_problem_undir_1vn_1t1(instance, s_network_dir, model)
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



function print_polytope(poly, names_variables)
    hr = hrep(poly);
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











#------------------- TEST stuff


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


# This is experimental, should work because we are in undirected, so we can remove one of the variable for each substrate edge.
#ALSO SOME WEIRD stuff
function get_sols_undir_gigatest(instance)

    model = Model(Gurobi.Optimizer)
    #set_silent(model)

    s_network_dir = generate_dir_sn(instance)
    set_up_problem_undir_1vn_1t1(instance, s_network_dir, model)
    # BREAK THE SYMETRIES ?
    @constraint(model, model[:x][1, 1] +model[:x][1, 3] +model[:x][1, 2]   == 1)
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