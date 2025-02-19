using JSON
using DataStructures
using Graphs, MetaGraphsNext
includet("../../../../../../utils/graph.jl")

function do_stufff()
    for file in readdir("sn_new")
        mg = read_substrate("sn_new/$file")
        #print_graph(mg)
        write_network_to_json(mg)
    end
end






function write_network_to_json(mg)
    formatted_json_str = "{\n"
    
    # about the graph:
    for k in keys(mg[])
        formatted_json_str *= "\t\"" * string(k) * "\": "
        formatted_json_str *= "\"" * string(mg[][k]) * "\""

        formatted_json_str *= ",\n"
    end

    formatted_json_str *= "\t\"directed\": false,\n"

    # node:
    formatted_json_str *= "\t\"nodes\": [\n"
    for node in vertices(mg)
        formatted_json_str = formatted_json_str * "\t\t{\"id\": " * string(node) * ","
        for k in keys(mg[node])
            formatted_json_str *= " \"" * string(k) * "\": " * string(mg[node][k]) * ","
        end
        formatted_json_str = formatted_json_str[1:end-1] 
        formatted_json_str *= "},\n"
    end
    formatted_json_str = formatted_json_str[1:end-2] 
    formatted_json_str *= "\n\t],"

    # edges
    formatted_json_str *= "\n\t\"edges\": [\n"
    for edge in edges(mg)
        formatted_json_str = formatted_json_str * "\t\t{\"source\": " * string(src(edge)) * ", \"target\": " * string(dst(edge)) * ", "
        for k in keys(mg[src(edge), dst(edge)])
            formatted_json_str *= " \"" * string(k) * "\": " * string(mg[src(edge), dst(edge)][k]) * ","
        end
        formatted_json_str = formatted_json_str[1:end-1] 
        formatted_json_str *= "},\n"
    end
    formatted_json_str = formatted_json_str[1:end-2] 
    formatted_json_str *= "\n\t]"

    formatted_json_str *= "\n}"

    # Write the formatted JSON string to a file
    filename = mg[][:name] * ".json"
    open(filename, "w") do f
        write(f, formatted_json_str)
    end
end






function read_graph(filename::String)
    println("Ok so")
    json_graph = JSON.parsefile(filename)
    if json_graph["type"] == "virtual"
        g = read_virtual(json_graph)
        return g, 1
    end
    if json_graph["type"] == "substrate"
        g = read_substrate(json_graph)
        return g, 2
    end
end

function read_substrate(filename)
    json_graph = JSON.parsefile(filename)

    # testing whether it is directed or undirected. By default: undirected.
    println("Ok ")
    g = MetaGraph(
        Graph(),
        Int,
        Dict,
        Dict,
        Dict(:name => json_graph["name"], :type => "substrate")
    )

    println("Wesh")

    for node in json_graph["nodes"]
        add_vertex!(g, node["id"], Dict{Any, Any}(:cap=> 1, :cost => rand(2:4)))
    end

    for edge in json_graph["edges"]
        add_edge!(g, edge["source"], edge["target"], Dict{Any, Any}(:cap => rand(2:4), :cost => rand(2:4)))
    end

    return g
end

function read_virtual(json_graph)


    # testing whether it is directed or undirected. By default: directed.
    if "directed" in keys(json_graph)
        if json_graph["directed"]
            g = MetaGraph(
                DiGraph(),
                Int,
                Dict,
                Dict,
                Dict(:name => json_graph["name"], :type => "virtual", :directed => json_graph["directed"])
            )
        else
            g = MetaGraph(
                Graph(),
                Int,
                Dict,
                Dict,
                Dict(:name => json_graph["name"], :type => "virtual", :directed => json_graph["directed"])
            )
        end
    else
        g = MetaGraph(
            DiGraph(),
            Int,
            Dict,
            Dict,
            Dict(:name => json_graph["name"], :type => "virtual", :directed => true)
        )
    end

    for node in json_graph["nodes"]
        add_vertex!(g, node["id"], Dict{Any, Any}(:dem => node["dem"]))
    end

    for edge in json_graph["edges"]
        add_edge!(g, edge["source"], edge["target"], Dict{Any, Any}(:dem => edge["dem"]))
    end

    return g
end

