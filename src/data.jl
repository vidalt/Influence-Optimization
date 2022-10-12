# import Base.show, Base.print

struct DataGlcip
    nodes::UnitRange{Int}
    k::Int
    beta::Float64
    dmin::Int
    dmax::Int
    gamma::Float64
    inr::Int

    hurdle::Vector{Int}
    H:: Int
    influence::Matrix{Int}
    weights::Dict{Tuple{Int64,Int64},Float64}
    min_incentives::Array{UnitRange{Int64},1}
    incentives::Vector{Float64}
    graph::SimpleDiGraph

    function DataGlcip(nodes::Float64, k::Float64, beta::Float64, dmin::Float64, dmax::Float64,
            gamma::Float64, inr::Float64, hurdle::Vector{Int}, H::Int,
            influence::Matrix{Int}, weights::Dict{Tuple{Int64,Int64},Float64},
            min_incentives::Array{UnitRange{Int64},1}, incentives::Vector{Float64})
        n = Int(nodes)

        graph = create_graph(n, influence)

        new(1:n, Int(k), beta, Int(dmin), Int(dmax), gamma, Int(inr), hurdle, H,
            influence, weights, min_incentives, incentives, graph)
    end
end

function create_graph(n, influence)
    # Input Graph
    graph = SimpleDiGraph(n)
    for i in 1:n
        for j in 1:n
            if(influence[i,j] > 0)
                add_edge!(graph, i, j)
            end
        end
    end

    return graph
end

function get_incentives(H)
    return [0, 0.25H, 0.5H, 0.75H, H]
end

function read_glcip_data(filename::String, get_incentives::Function)
    data = Vector{Float64}(undef, 0)

    headers = []
    open(filename) do file
        for line in eachline(file)
            m = match(r"^\s*(?:#|$)", line)
            if m === nothing
                splitted_line = split(line)

                if length(splitted_line) >= 7 # Headers Case
                    headers = splitted_line
                else
                    for peaceofdata in splitted_line
                        push!(data, parse(Float64, peaceofdata))
                    end
                end
            end
        end
    end
    headers_data = parse.(Float64, headers)

    n_vertex, n_edges = Int.(data[1:2])

    # Hurdle
    offset = 2
    end_vertices = Int(2*n_vertex+offset)
    hurdle = Int.(data[offset+2:2:end_vertices])

    # Influence values on arcs (i,j)
    influence = zeros(Int, n_vertex, n_vertex)
    for (_, i,j,d) in Iterators.partition(data[end_vertices+1:end],4)
        influence[Int(i)+1, Int(j)+1] = Int(d)
    end

    # Max Hurdle
    if length(headers) == 7
        H = maximum(hurdle)
    elseif length(headers) == 8
        H = parse(Int, headers[8])
    end

    # incentives in P_i are integers obtained by rounding up the frac. values
    P = ceil.(get_incentives(H))

    W_ip = Dict((i,p) => trunc(P[p]^0.9) for i in 1:n_vertex for p in 1:length(P))

    # Preprocessing step to compute the cheapest incentive $p$ that is
    # sufficient to activate a node without receiving influence from its neighbors
    minimum_incentive(i::Int) = minimum(filter(p-> P[p] >= hurdle[i], 1:length(P)))
    cheapest_incentives(i) = 1:length(filter(p->W_ip[i,p] <= W_ip[i, minimum_incentive(i)], 1:length(P)))

    P_i = cheapest_incentives.(1:n_vertex)
    data_glcip = DataGlcip(headers_data[1:7]..., hurdle, H, influence, W_ip, P_i, P)

    return data_glcip
end

i_neighborhood(data::DataGlcip) = i->[j for (j,v) in enumerate(data.influence[:, i]) if v>0]
influence(data::DataGlcip, i::Int, j::Int) = data.influence[i,j]
sum_influence(data::DataGlcip) = i->sum(data.influence[:, i])


# function show(io::IO, d::DataGlcip)
#     println(io, "Generalized Least Cost Influence Propagation dataset.")
#     println(io, "Number of individuals = $(length(d.nodes))")
#
#     println(io, "Thresholds that need to be reached through neighborinf influence: ")
#     for (i, h) in enumerate(d.hurdle)
#         println(io, "\t vertex $i, hurdle = $h")
#     end
# end
