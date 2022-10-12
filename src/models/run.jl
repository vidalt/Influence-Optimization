using JuMP, CPLEX
import MathProgBase

include("icc_model.jl")
include("cf_model.jl")

function build_model(data::DataGlcip, app)
    command = app["%COMMAND%"]
 
    if command == "cf"
       return cf_model(data, app)
    else
       return icc_model(data, app)
    end
 end

function solve_glcip(data::DataGlcip, input_file, app)
    (model, variables) = build_model(data, app)

    status = solve(model)

    # Output
    N = i_neighborhood(data)
    # W = data.weights
    # INCENTIVES = data.incentives

    if status == :Optimal
        out_x = zeros(Int, length(data.nodes))
        out_y = zeros(Int, length(data.nodes))
        out_z = zeros(Int, length(data.nodes), length(data.nodes))

        # Variable - X
        if haskey(variables, :x)
            x = variables[:x]
            for i in data.nodes

                x_i = round(getvalue(x[i]))
                if x_i > 0
                    out_x[i] = x_i
                    # println("x[$i]: $x_i | h[$i]:$(data.hurdle[i])")
                end
            end
            # println()
        end

        # Variable - Z
        if haskey(variables, :z)
            z = variables[:z]
            for i in data.nodes
                for j in N(i)
                    zij = round(getvalue(z[i, j]))
                    if (zij > 0)
                        out_z[j, i] = zij
                        # println("z[$i, $j]: $(zij) | d($i, $j): ")
                    end
                end
            end
            # println()
        end

        # Variable - Y
        if haskey(variables, :y)
            y = variables[:y]
            for i in data.nodes
                for p in data.min_incentives[i]
                    temp = getvalue(y[i, p])
                    yip = round(temp)
                    if (yip > 0)
                        if p > 1
                            out_y[i] = p
                            # println("y[$i, $p]: $(W[(i,p)]) | $(INCENTIVES[p])")
                        end
                    end
                end
            end 
        end
    end

    best_bound = getobjbound(model)
    objective = getobjectivevalue(model)
    solve_time = getsolvetime(model)

    bc_sol = BCSolution(data, input_file, app, status, best_bound, objective, solve_time, out_x, out_y, out_z)

    return bc_sol
end
