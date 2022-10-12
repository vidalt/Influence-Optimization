abstract type Solution end

struct BCSolution <: Solution
    data::DataGlcip
    instance::String
    status::Symbol
    best_bound::Float64 # Upper bound
    objective::Float64 # Cost
    time::Float64 # Time (s)
    alpha::Float64

    app::Dict{String, Any}

    actives::Vector{Int}

    x::Vector{Int}
    y::Vector{Int}
    z::Matrix{Int}
    graph::SimpleDiGraph

    function BCSolution(
        data::DataGlcip,
        instance::String,
        app::Dict{String, Any},
        status::Symbol,
        best_bound::Float64,
        objective::Float64,
        time::Float64,
        x::Vector{Int},
        y::Vector{Int},
        z::Matrix{Int},
    )
        actives = filter(i -> i == 1, x)

        graph = create_graph(x, y, z)

        new(
            data,
            instance,
            status,
            best_bound,
            objective,
            time,
            app["alpha"],
            app,
            actives,
            x,
            y,
            z,
            graph,
        )
    end
end

function create_graph(x, y, z)
    graph = SimpleDiGraph(length(x))

    # Edges
    s = size(z)[1]
    for i = 1:s
        for j = 1:s
            if z[i, j] > 0
                add_edge!(graph, i, j)
            end
        end
    end


    return graph
end

# Checks the feasiblity of a solution
function check_feasibility(sol::Solution)
    return ~is_cyclic(sol.graph) && length(sol.actives) >= sol.alpha
end

# write solution as TikZ figure (.tex)
function draw_solution(
    solution::Solution,
    filename::String = "out/latex/result_graph",
)
    n_incentives = length(unique(solution.y)) - 1

    h(i::Int) = solution.data.hurdle[i]
    d(i::Int, j::Int) = solution.data.influence[j, i]
    incentives_nodes = Dict()

    N = size(solution.z, 2)
    g = DiGraph(N + n_incentives)

    # Edges
    s = size(solution.z)[1]
    for i = 1:s
        for j = 1:s
            if solution.z[i, j] > 0
                # @info "Add edge from $(i+n_incentives) -> $(j+n_incentives)"
                add_edge!(g, i + n_incentives, j + n_incentives)
            end
        end
    end
    influence_labels = Dict(
        (e.src, e.dst) =>
            solution.data.influence[e.src-n_incentives, e.dst-n_incentives]
        for e in edges(g)
    )

    # Incentives
    c = 0
    incentives_labels = Dict()
    main_nodes_styles = Dict()
    for (i, v) in enumerate(solution.y)
        if solution.x[i] > 0
            main_nodes_styles[i+n_incentives] = "draw, fill=green!10"
        else
            main_nodes_styles[i+n_incentives] = "draw, fill=red!10"
        end
        if v != 0
            if !(v in keys(incentives_nodes))
                c += 1
                incentives_nodes[v] = c
            end

            incentives_labels[(incentives_nodes[v], i + n_incentives)] =
                solution.data.incentives[v]
            # @info "Add edge from $(incentives_nodes[v]) -> $(i+n_incentives)"
            add_edge!(g, incentives_nodes[v], i + n_incentives)
        end
    end

    temp = zeros(n_incentives)
    for (k, v) in incentives_nodes
        temp[v] = solution.data.incentives[k]
    end

    incentive_nodes =
        Dict(v => "draw, fill=blue!10" for (k, v) in incentives_nodes)

    labels = vcat(string.(temp), ["$i ($(h(i)))" for i = 1:N])
    node_styles = merge(main_nodes_styles, incentive_nodes)
    edge_labels = merge(influence_labels, Dict())

    t = TikzGraphs.plot(
        g,
        labels = labels,
        node_styles = node_styles,
        edge_labels = edge_labels,
        options = "scale=2",
    )

    TikzPictures.save(TEX(filename), t)
    # TikzPictures.save(PDF("result_graph"), t)
    
    @info("Exported LaTeX solution to $filename")
    @info("\n Use `lualatex $(filename).tex` to generate a PDF with a graph representation of the solution")
end

# # Read solution from file
# function read_solution(filename::String)
#     data = Vector{Float64}(undef, 0)

#     open(filename) do file
#         for line in eachline(file)
#             m = match(r"^\s*(?:#|$)", line)
#             if m === nothing
#                 for peaceofdata in split(line)
#                     push!(data, parse(Float64, peaceofdata))
#                 end
#             end
#         end
#     end

#     open(filename, "r") do f
#         instance = data[1]
#         data_glcip = read_glcip_data(instance)

#         # sol = BCSolution(data_glcip, data[1:7]..., x, y, z)
#     end
# end

# Export solution to a file
function export_solution(sol::Solution, filename::String)
    unixtime = round(Int64, time() * 1000)
    output_filename = "$(filename)_$(unixtime).out"

    open(output_filename, "w") do f
        println(f, sol)
    end

    @info("Exported solution to $output_filename")
end

# Display solution
function Base.show(io::IO, sol::Solution)
    println(io, "Instance = $(sol.instance)")
    println(io, "Status = $(sol.status)")
    println(io, "Best bound = $(sol.best_bound)")
    println(io, "Objective value = $(sol.objective)")
    println(io, "Time (s) = $(sol.time)")
    println(io, "Alpha = $(sol.alpha)")
    println(io, "Length(actives) = $(length(sol.actives))")

    println(io, "Application parameters:")
    for (arg, val) in sol.app
        println(io, "  $arg  =>  $(repr(val))")
    end

    println(io, "\nSolution: ")
    # Variable - X
    for (i, x) in enumerate(sol.x)
        if x > 0
            println(io, "X[$i]= $x")
        end
    end
    println(io)

    # Variable - Z
    for i = 1:size(sol.z, 2)
        for (j, v) in enumerate(sol.z[i, :])
            # @show i, j
            if v > 0
                println(io, "z[$i, $j]= $v")
            end
        end
    end
    println(io)

    W = sol.data.weights
    P = sol.data.incentives

    # Variable - Y
    for (i, p) in enumerate(sol.y)
        if p != 0
            println(io, "y[$i, $p] = 1 | P[$p] = $(P[p]) | W[$i, $p] = $(W[(i,p)])")
        end
    end
end
