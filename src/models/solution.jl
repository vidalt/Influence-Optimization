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

        graph = create_graph(x, z)

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

function create_graph(x, z)
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
    if is_cyclic(sol.graph)
        return false, "Cyclic graph"
    end

    # Notations
    INCENTIVES = sol.data.incentives
    N = i_neighborhood(sol.data)
    d(i::Int, j::Int) = influence(sol.data, j, i)
    h(i::Int) = sol.data.hurdle[i]
    W = sol.data.weights

    # Check for all nodes
    actives = 0
    for (i, x) in enumerate(sol.x)
        influences = []
        for j in N(i)
            if sol.z[j, i] <= x
                influence = d(i,j)
                push!(influences, influence)
            else
                return false, "Influence can only be exerted on arc (j,i) if node j is active"
            end
        end

        if x == 1 # active
            p = sol.y[i]
            incentive = INCENTIVES[p]

            # influences = sum([d(i,j) for j in N(i) if sol.z[i, j] > 0])^sol.app["gamma"]
            diffusion = incentive + sum(influences)^sol.app["gamma"]

            if diffusion >= h(i)
                actives += 1
            else
                return false, "Incentive + influences not suficient to activate node `$i`"
            end
            # println("y[$i, $p] = 1 | P[$p] = $(incentive) | W[$i, $p] = $(W[(i,p)])")
        end
    end

    if actives < sol.alpha * length(sol.data.nodes)
        return false,  "Number of active nodes < alpha |V|"
    end

    test_cost = sum(W[(i,p)] for (i, p) in enumerate(sol.y))

    if sol.objective != test_cost
        return false, "Wrong cost"
    end

    return true, ""
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

# Export solution from a file
function read_solution(filename::String)
    open(filename) do file
        list_results = []
        list_params = []
        command_params = ""

        out_x = Vector{Int}()
        dict_y = Dict()
        dict_W = Dict()
        dict_z = Dict()

        # Extract data from file
        for line in eachline(file)
            
            m = match(r"^\s*(?:#|$)", line) # Comments
            if m === nothing
                # Extract data from results
                regex_results = r"^(?:\w| |\(|\))+ = "
                if match(regex_results, line) !== nothing
                    temp = split(line, "=")
                    value_results = strip(temp[2])
                    push!(list_results, value_results)
                end

                # Extract data from parameters
                regex_params = r"^\s+(?:\w|\%)+ \s* => \s*(\d|\w|\"|\.|\/|-)+$"
                if match(regex_params, line) !== nothing
                    temp = split(line, "=>")
                    value_params = replace(strip(temp[2]), "\"" => "")
                    push!(list_params, value_params)
                end

                # Extract data from command parameters
                regex_command_params = r"Dict{String,Any}\((?:.*)\)"
                match_command_params = match(regex_command_params, line)
                if match_command_params !== nothing
                    command_params = match_command_params.match
                end

                # Extract data from variable `x`
                regex_x = r"^(?:x|X)\[(\d+)\] *= *(\d+)$"
                match_x = match(regex_x, line)
                if match_x !== nothing
                    value_x = parse(Int, match_x[2])
                    push!(out_x, value_x)
                end

                # Extract data from variable `z`
                regex_z = r"^(?:z|Z)\[(\d+(?:, \d+)*)\] *= *(\d+)$"
                match_z = match(regex_z, line)
                if match_z !== nothing
                    z_indices = split(match_z[1], ",")
                    i = parse(Int, strip(z_indices[1]))
                    j = parse(Int, strip(z_indices[2]))
                    z_ij = parse(Int, strip(match_z[2]))

                    dict_z[(i, j)] = z_ij
                end

                # Extract data from variable `y`
                regex_y = r"^(?:y|Y)\[(\d+(?:, \d+)*)\] *= *(\d+)"
                match_y = match(regex_y, line)
                if match_y !== nothing
                    y_indices = split(match_y[1], ",")
                    i = parse(Int, strip(y_indices[1]))
                    # p = parse(Int, strip(y_indices[2]))

                    temp = split(line, " | ")

                    regex_p = r"P\[(\d+)\] *= *((?:\d|\.)+)"
                    match_p = match(regex_p, line)
                    if match_p !== nothing
                        dict_y[i] = parse(Int, strip(match_p[1])) 
                    end

                    regex_W = r"W\[(\d+(?:, \d+)*)\] *= *((?:\d|\.)+)$"
                    match_W = match(regex_W, line)
                    if match_W !== nothing
                        dict_W[i] = parse(Float64, strip(match_W[2])) + 1
                    end
                end
            end
        end

        input_file, status, best_bound, objective, solve_time, alpha, actives = list_results

        upcutoff, latex, output, alpha, instance_type, verbose, gamma, command, filepath = list_params

        # Create list of parameters
        app = Dict(
            "upcutoff"  =>  parse(Float64, upcutoff),
            "latex"  =>  parse(Bool, latex),
            "output"  =>  String(output),
            "alpha"  =>  parse(Float64, alpha),
            "instance"  =>  String(instance_type),
            "verbose"  =>  parse(Bool, verbose),
            "gamma"  =>  parse(Float64, gamma),
            "%COMMAND%"  =>  String(command),
            "filepath"  =>  String(filepath),
            command => eval(Meta.parse(command_params))
        )

        # Read instance data of `filepath`
        data = read_glcip_data(String(filepath), String(instance_type))

        # Create influence matrix  `z_ij` (∀i,j in V)
        out_z = zeros(Int, length(data.nodes), length(data.nodes))
        for (key, value) in dict_z
            i,j = key
            # out_z[j, i] = value
            out_z[i, j] = value
        end

        # Create list of incentives `y_ip` (∀i in V, p in P(i))
        out_y = ones(Int, length(data.nodes))
        for (i, p) in dict_y
            if app["instance"] == "SW"
                out_y[i] = p
            else
                out_y[i] = dict_W[i]
            end
        end

        # Create `solution` structure
        bc_sol = BCSolution(
            data, 
            String(input_file), 
            app, 
            Symbol(status), 
            parse(Float64, best_bound), 
            parse(Float64, objective), 
            parse(Float64, solve_time), 
            out_x, 
            out_y, 
            out_z
        )

        return bc_sol
    end
end

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
            println(io, "x[$i]= $x")
        end
    end
    println(io)

    # Variable - Z
    for i = 1:size(sol.z, 2)
        for (j, v) in enumerate(sol.z[i, :])
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
            if sol.app["instance"] == "SW"
                println(io, "y[$i, $p] = 1 | P[$p] = $(P[p]) | W[$i, $p] = $(W[(i,p)])")
            elseif sol.app["instance"] == "GRZ"
                println(io, "y[$i, $p] = 1 | P[$p] = $(p - 1) | W[$i, $p] = $(p - 1)")
            end
        end
    end
end
