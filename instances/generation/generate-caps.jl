using Graphs, MetaGraphsNext

includet("../../utils/import_utils.jl")
includet("io.jl")



function put_as_substrate_and_save(g, name)
    mg = MetaGraph(
        Graph(),
        Int,
        Dict,
        Dict,
        Dict(:name=>"s_" * name, :type=>"substrate")
    )

    for node in vertices(g)
        cap_to_put = 1
        cost_to_put = rand(2:4)
        add_vertex!(mg, node, Dict(:cap=> cap_to_put, :cost => cost_to_put))
    end

    for edge in edges(g)
        cap_to_put = 5 #rand(2:4)
        cost_to_put = rand(2:4)
        add_edge!(mg, src(edge), dst(edge), Dict(:cap=>cap_to_put, :cost =>cost_to_put))
    end

    write_network_to_json(mg, false)
end

function put_as_virtual_and_save(g, name)
    mg = MetaGraph(
        Graph(),
        Int,
        Dict,
        Dict,
        Dict(:name=>"v_" * name, :type=>"virtual")
    )

    for node in vertices(g)
        add_vertex!(mg, node, Dict(:dem=> 1))
    end

    for edge in edges(g)
        add_edge!(mg, src(edge), dst(edge), Dict(:dem=>1))
    end

    write_network_to_json(mg, false)
end
