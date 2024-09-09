using JSON
using DataStructures

function write_s_network_to_json(mg)
    nodes_info = [
    OrderedDict("id" => node, "capacity" => mg[node][:cap], "cost" => mg[node][:cost])
    for node in vertices(mg)
    ]

    edges_info = [
        OrderedDict("source" => src(edge), "target" => dst(edge), "capacity" => mg[src(edge), dst(edge)][:cap], "cost" => mg[src(edge), dst(edge)][:cost])
        for edge in edges(mg)
    ]

    graph_dict = OrderedDict("name" => mg[][:name], "type" => mg[][:type], "nodes" => nodes_info, "edges" => edges_info)

    # Convert the dictionary to a JSON string with 2 spaces of indentation
    json_str = JSON.json(graph_dict, 2)

    # Post-process the JSON string to match the desired formatting
    formatted_json_str = replace(json_str, r"(\s*{\n\s*\"id\": \d+.*\n\s*})" => s -> replace(s, "\n" => " "))

    # Write the formatted JSON string to a file
    filename = "output.json"
    open(filename, "w") do f
        write(f, formatted_json_str)
    end
end

function write_network_to_json(mg)
    formatted_json_str = "{\n"
    
    # about the graph:
    for k in keys(mg[])
        formatted_json_str *= "\t\"" * string(k) * "\": "
        if isa(mg[][k], String)
            formatted_json_str *= "\"" * string(mg[][k]) * "\""
        else
            formatted_json *= string(mg[][k])
        end
        formatted_json_str *= ",\n"
    end

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

function write_v_network_to_json(mg)
    nodes_info = [
    OrderedDict("id" => node, "demand" => mg[node][:dem])
    for node in vertices(mg)
    ]

    edges_info = [
        OrderedDict("source" => src(edge), "target" => dst(edge), "demand" => mg[src(edge), dst(edge)][:dem])
        for edge in edges(mg)
    ]

    graph_dict = OrderedDict("name" => mg[][:name], "type" => "virtual", "nodes" => nodes_info, "edges" => edges_info)

    filename = mg[][:name] * ".json"
    # Write the JSON string to a file
    open(filename,"w") do f
        JSON.print(f, graph_dict, 4)
    end
end



function some_triangle(n)
    for i in 1:n

        mg = MetaDiGraph()
    
        set_prop!(mg, :name, "triangle_" * string(i))
        set_prop!(mg, :type, "virtual")
    
        add_vertices!(mg, 3)
        add_edge!(mg, 1, 2)
        add_edge!(mg, 2, 3)
        add_edge!(mg, 3, 1)
        add_edge!(mg, 3, 2)
    
        for node in vertices(mg)
            set_prop!(mg, node, :dem, rand(1:2))
        end
        for edge in edges(mg)
            set_prop!(mg, edge, :dem, rand(1:2))
        end
        write_graph_to_json(mg)
    end
end