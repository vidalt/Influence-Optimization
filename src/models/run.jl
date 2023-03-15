using JuMP, CPLEX
using Combinatorics
import MathProgBase

include("icc_model.jl")
include("cf_model.jl")

# TODO: caso discreto
function preprocess_incentives_bkp(data::DataGlcip, i::Int)
    # Notations
    N = i_neighborhood(data)
    d(i::Int, j::Int) = influence(data, j, i)
    h(i::Int) = data.hurdle[i] .- 0.5

    # Preprocessing
    preprocessed_sum_influences = Set([0])

    for infl in N(i)
        newer = [s + d(infl, i) for s in preprocessed_sum_influences]
        for s in newer
            if s < h(i)
                push!(preprocessed_sum_influences, s)
            end
        end
    end

    push!(preprocessed_sum_influences, Int(h(i)+0.5))

    return preprocessed_sum_influences
end

function preprocess_incentives(data::DataGlcip, i::Int)
    # Notations
    N = i_neighborhood(data)
    d(i::Int, j::Int) = influence(data, j, i)
    h(i::Int) = data.hurdle[i] .- 0.5
    P(i::Int) = data.min_incentives[i]

    # Preprocessing
    preprocessed_sum_influences = Set([1])

    INCTV(p::Int) = data.incentives[p]

    _INCENTIVES = [INCTV(p) for p in P(i)]
    incentives_idx = Set(p for p in P(i))
    
    for c in combinations(N(i))
        sum_infls = sum(d.(i, c))
        
        if sum_infls < h(i)
            current_idx = first(incentives_idx)
            while !isempty(incentives_idx)
                activation = _INCENTIVES[current_idx] + sum_infls
                if activation >= h(i)
                    push!(preprocessed_sum_influences, current_idx)

                    # Update current index
                    delete!(incentives_idx, current_idx)
                    break
                end
                current_idx += 1
            end

            if length(preprocessed_sum_influences) == length(P(i))
                break
            end
        
        end
    end

    return sort(collect(preprocessed_sum_influences))
end

function build_model(data::DataGlcip, app::Dict{String,Any}, preprocessed_sum_influences::Array{Array{Int,1},1})
    command = app["%COMMAND%"]

    if command == "cf"
        return cf_model(data, app, preprocessed_sum_influences)
    else
        return icc_model(data, app, preprocessed_sum_influences)
    end
end

function solve_glcip(data::DataGlcip, input_file::String, app::Dict{String,Any})
    N = i_neighborhood(data)
    d(i::Int, j::Int) = influence(data, j, i)
    h(i::Int) = data.hurdle[i] .- 0.5
    # P(i::Int) = data.min_incentives[i]
    infls(i::Int) = [d(i,j) for j in N(i)]

    preprocessed_sum_influences = [preprocess_incentives(data, i) for i in data.nodes]
    # P(i::Int) = preprocessed_sum_influences[i]

    (model, variables, P, W) = build_model(data, app, preprocessed_sum_influences)
    status = solve(model)

    # Output
    out_x = zeros(Int, length(data.nodes))
    out_y = zeros(Int, length(data.nodes))
    out_z = zeros(Int, length(data.nodes), length(data.nodes))

    # if status == :Optimal
    # Variable - X
    if haskey(variables, :x)
        x = variables[:x]
        for i in data.nodes
            x_i = round(getvalue(x[i]))
            if x_i > 0
                out_x[i] = x_i
            end
        end
    end

    # Variable - Z
    if haskey(variables, :z)
        z = variables[:z]
        for i in data.nodes
            for j in N(i)
                z_ij = round(getvalue(z[i, j]))
                if (z_ij > 0)
                    out_z[j, i] = z_ij
                end
            end
        end
    end

    # Variable - Y
    if haskey(variables, :y)
        y = variables[:y]
        for i in data.nodes
            # for p in data.min_incentives[i]
            for p in P(i)
                y_ip = getvalue(y[i, p])
                if y_ip !== nothing
                    y_ip_value = round(y_ip)
                    if (y_ip_value > 0)
                        if p > 1
                            if app["instance"] == "SW"
                                out_y[i] = p
                            elseif app["instance"] == "GRZ"
                                out_y[i] = W[(i,p)] + 1
                            end
                        end
                    end
                end
            end
        end 
    end
    # end

    best_bound = getobjbound(model)
    objective = getobjectivevalue(model)
    solve_time = getsolvetime(model)

    bc_sol = BCSolution(data, input_file, app, status, best_bound, objective, solve_time, out_x, out_y, out_z)

    return bc_sol
end
