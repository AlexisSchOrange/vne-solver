import Base.parse
import JSON
using Graphs, MetaGraphsNext


function read_graph(filename::String)
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

function read_substrate(json_graph)
    
    # testing whether it is directed or undirected. By default: undirected.
    if "directed" in keys(json_graph)
        if json_graph["directed"]
            g = MetaGraph(
                DiGraph(),
                Int,
                Dict,
                Dict,
                Dict(:name => json_graph["name"], :type => "substrate", :directed => json_graph["directed"])
            )
        else
            g = MetaGraph(
                Graph(),
                Int,
                Dict,
                Dict,
                Dict(:name => json_graph["name"], :type => "substrate", :directed => json_graph["directed"])
            )
        end
    else
        g = MetaGraph(
            DiGraph(),
            Int,
            Dict,
            Dict,
            Dict(:name => json_graph["name"], :type => "substrate", :directed => true)
        )
    end


    for node in json_graph["nodes"]
        add_vertex!(g, node["id"], Dict{Any, Any}(:cap=> node["cap"], :cost => node["cost"]))
    end

    for edge in json_graph["edges"]
        add_edge!(g, edge["source"], edge["target"], Dict{Any, Any}(:cap => edge["cap"], :cost => edge["cost"]))
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



function get_instance_from_folder(folder_path::String)
    virtual_networks = []
    substrate_network = nothing
    for filename in readdir(folder_path; join=true)
        g, type = read_graph(filename)
        if type == 1
            push!(virtual_networks, g)
        else
            substrate_network = g
        end
    end

    instance = InstanceVNE(virtual_networks, substrate_network)

    return(instance)
end

