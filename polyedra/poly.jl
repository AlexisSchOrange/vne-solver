
using Polyhedra, CDDLib
includet("utils/file_reader.jl")
includet("utils/import_utils.jl")
includet("resolution/directed/compact/compact_formulation.jl")





# takes an instance (with one vr !), returns all the solutions (points) and the names of variables
function get_v_rep_convexhull(instance, time_limit=100, sol_max=1000)

    
    # builds the model and find all optimal points (a bit obscure in cplex, i found this online)
    compact_model = set_up_compact_model(instance)
    model = compact_model.model
    #set_time_limit_sec(model, time_limit)

    set_silent(model)
    optimize!(model)
    #set_optimizer_attribute(model, "CPX_PARAM_DISPLAY", 0)

    set_optimizer_attribute(model, "CPX_PARAM_SOLNPOOLAGAP", 100000000.0)
    set_optimizer_attribute(model, "CPX_PARAM_SOLNPOOLINTENSITY", 4)
    set_optimizer_attribute(model, "CPX_PARAM_POPULATELIM", sol_max)
    
    backend_model = unsafe_backend(model);
    env = backend_model.env;
    lp = backend_model.lp;
    
    CPLEX.CPXpopulate(env, lp);
    
    N_results = CPLEX.CPXgetsolnpoolnumsolns(env, lp)
    if N_results == sol_max
        println("The number of solution max has been reached. Some solutions might not have been compute.")
    end

    # get correct names
    v_network = instance.v_networks[1]
    names_variables_dic = Dict()
    
    for v_node in vertices(v_network)
        for s_node in vertices(instance.s_network)
            names_variables_dic[CPLEX.column(backend_model, compact_model.x_variables[v_network, v_node, s_node].index)] = "x_" * string(v_node) * "_" * string(s_node)
        end
    end
    
    for v_edge in edges(v_network)
        for s_edge in edges(instance.s_network)
            names_variables_dic[CPLEX.column(backend_model, compact_model.y_variables[v_network, v_edge, s_edge].index)] = "y_" * string(src(v_edge)) * string(dst(v_edge)) * "_" * string(src(s_edge)) * string(dst(s_edge))
        end
    end
    
    names_variables = []
    for i_var in 1:length(names_variables_dic)
        push!(names_variables, names_variables_dic[i_var])
    end

    # get all solutions vectors
    all_sols = Vector{Vector{Int64}}()
    for sol in 0:N_results-1
        current_sol = Vector{Int64}()
        for i_variable in 0:length(names_variables)-1
            x = Ref{Cdouble}()  ## Reference to the `Term` variable value
            CPLEX.CPXgetsolnpoolx(env, lp, sol, x, i_variable, i_variable)
            push!(current_sol, convert.(Int64, round.(x[])))
        end
        push!(all_sols, current_sol)
    end

    println("\n\nThere are " * string(N_results) * " solutions")
    return(all_sols, names_variables)
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


function print_polytope_simpler(poly, names_variables, print_index=false)
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



function print_variable(names_variables)
    for i_var in 1:length(names_variables)
        println(string(i_var) * " : " * names_variables[i_var])
    end
end


function eliminate_variables(poly_entry, names_variables, to_elim)
    poly_exit = eliminate(poly_entry, to_elim)
    names_variables_exit = copy(names_variables)
    for i_truc in 1:length(to_elim)
        deleteat!(names_variables_exit, to_elim[length(to_elim) - i_truc + 1])
    end
    return(poly_exit, names_variables_exit)
end


function do_everything(instance)
    all_sols, names_variables = get_v_rep_convexhull(instance);
    v_rep = vrep(all_sols);
    poly = polyhedron(v_rep, CDDLib.Library(:float));
    return poly, names_variables
end

function elimate_and_do_everything(poly_entry, names_variables, to_elim)
    poly_project = eliminate(poly_entry, to_elim)
    names_variables_project = copy(names_variables)
    for i_truc in 1:length(to_elim)
        deleteat!(names_variables_project, to_elim[length(to_elim) - i_truc + 1])
    end
    return poly_project, names_variables_project
end

function do_everything_and_print(instance)
    all_sols, names_variables = get_v_rep_convexhull(instance);
    v_rep = vrep(all_sols);
    poly = polyhedron(v_rep, CDDLib.Library(:float));
    h_rep = hrep(poly);
    print_polytope(h_rep, names_variables);
end