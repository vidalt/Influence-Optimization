# Create a TimerOutput, this is the main type that keeps track of everything.
const to = TimerOutput()

function sliding_window(arr::Vector{Int})
   paths = []
   for (i, j) in zip(arr, arr[2:end])
      push!(paths, (i,j))
   end

   return paths
end

function preprocess_instance(data::DataGlcip)
   N = i_neighborhood(data)
   h(i::Int) = data.hurdle[i]
   d(i::Int, j::Int) = influence(data, j, i)

   W_ip = Dict()
   incentives = Vector{Vector{Int}}()
   for i in data.nodes
      P_i = [1]
      W_ip[(i, 1.0)] = 0

      # Todas as influências são iguais, entao basta pegar a primeira delas
      j = N(i)[1] 
      influence = d(i,j)

      ratio = h(i) / influence
      limit = ceil(ratio)
      for v in 1:limit
         push!(P_i, v+1)
         temp = h(i) - ((limit - v) * influence)
         W_ip[(i, v+1)] = temp < 0 ? h(i) : temp
      end

      push!(incentives, P_i)
   end

   P_new(i::Int) = incentives[i]
   d_new(i::Int, j::Int) = 1
   h_new(i::Int) = incentives[i][end] - 1

   return P_new, W_ip, d_new, h_new
end


function icc_model(data::DataGlcip, app::Dict{String,Any}, preprocessed_sum_influences::Array{Array{Int,1},1})
   # Parameters
   command = app["%COMMAND%"]

   has_step2 = app[command]["step2"]
 
   # Notations
   H = data.H
   INCENTIVES = data.incentives

   N = i_neighborhood(data)

   if app["instance"] == "GRZ"
      P, W, d, h = preprocess_instance(data)
      data.weights = W
   else
      P(i::Int) = data.min_incentives[i]
      P(i::Int) = preprocessed_sum_influences[i]
      W = data.weights
      d(i::Int, j::Int) = influence(data, j, i)
      h(i::Int) = data.hurdle[i] .- 0.5
   end
   
   h_gamma(i::Int) = h(i)^(1 / app["gamma"])
   
   function h_minus_p(i, p)
      temp = h(i) - INCENTIVES[p]

      return max(temp, 0)^(1 / app["gamma"])
   end
 
   function incentive_gamma(i, p)
      return ceil(h_gamma(i)) - ceil(h_minus_p(i, p))
   end

   # Model
   ## Create Model
   model = Model(solver=CplexSolver(CPX_PARAM_THREADS=params_cplex[:threads],
      CPX_PARAM_CUTUP=app["upcutoff"],
      #CPX_PARAM_CUTPASS=-1,
      CPX_PARAM_TILIM=params_cplex[:time_limit],
      CPX_PARAM_MIPDISPLAY=params_cplex[:mip_display],
      CPX_PARAM_MIPINTERVAL=params_cplex[:mip_interval],
      CPX_PARAM_SCRIND=params_cplex[:scrind])
   )

   @variables(model, begin
      x[i in data.nodes], Bin
      z[i in data.nodes, j in N(i)], Bin
      y[i in data.nodes, p in P(i)], Bin
   end)

   @objective(model, Min, sum(W[(i, p)] * y[i, p] for i in data.nodes, p in P(i)))

   @constraints(model, begin
      # diffusion[i in data.nodes], sum(INCENTIVES[p]*y[i,p] for p in P(i)) + sum(d(i,j)*z[i,j] for j in N(i))^bc_params["gamma"] >= h(i) * x[i]
      # diffusion[i in data.nodes], sum(ceil(h_gamma(i)) - ceil((abs(h(i) - INCENTIVES[p]))^(1/bc_params["gamma"]))*y[i,p] for p in P(i)) +
      diffusion[i in data.nodes], sum(incentive_gamma(i, p) * y[i, p] for p in P(i)) +
                                 sum(d(i, j) * z[i, j] for j in N(i)) >= ceil(h_gamma(i)) * x[i]
      incentive[i in data.nodes], sum([y[i, p] for p in P(i)]) == x[i]
      forcing[i in data.nodes, j in N(i); d(j, i) == 0], z[i, j] <= x[j]
      coverage, sum(x) >= ceil(app["alpha"] * (length(data.nodes) - 0.000001))
   end)

   function pre_add_cycles_constraints_up_to_length(limit_length::Int=4)
      histogram_cycles = zeros(Int, limit_length)
   
      @info("BEGIN - add all cycle elimination constraints (1d) for all cycles up to length four")
      
      p = Progress(length(data.nodes); showspeed=true, enabled=app["verbose"])
      for i in data.nodes
         pre_add_cycles_constraints_node_i_up_to_length!(1, limit_length, [i], histogram_cycles, i, i)
         next!(p)
      end
   
      @info("END - add all cycle elimination constraints (1d) for all cycles up to length four")
   
      for (i,v) in enumerate(histogram_cycles)
         @info("Cycles length[$i]: $v")
      end
   end
   
   function pre_add_cycles_constraints_node_i_up_to_length!(cur_length::Int, limit_length::Int, arr::Vector{Int}, histogram_cycles::Vector{Int}, i::Int, j::Int)
      if cur_length <= limit_length
         for j_i in N(j)
            arr_j = [arr; j_i]
            if influence(data, i, j_i) > 0
               path = sliding_window(arr_j)
               
               lhs = append!([z[e...] for e in path], [z[j_i, i]])
               rhs = [x[e[1]] for e in path]
               
               @constraint(model, sum(lhs) <= sum(rhs))
               
               histogram_cycles[cur_length] += 1
            end
      
            cur_length <= limit_length && pre_add_cycles_constraints_node_i_up_to_length!(cur_length+1, limit_length, arr_j, histogram_cycles, i, j_i)
         end
      end
   end

   # Generalized cycle elimination constraints (1d) for all cycles up to length four are added a priori to the model
   pre_add_cycles_constraints_up_to_length(app[command]["pre_add_cycles_up_to_length"])

   stop_cycle_cut = false
   cycle_cut_init_time = time()

   function simple_cycle_elimination_lazy_constraint(cb)
      z_val = zeros(Int, length(data.nodes), length(data.nodes))

      for i in data.nodes
         for j in N(i)
            z_ij = getvalue(z[i, j])
            z_val[i, j] = round(Int, z_ij)
         end
      end

      num_cuts = 0

      NUM_NODES = length(data.nodes)
      graph = create_graph(NUM_NODES, z_val)
      cycles = simplecycles_iter(graph, 1000)

      for cycle in cycles
         i = cycle[1]
         j = cycle[end]

         if influence(data, i, j) > 0
            path = sliding_window(cycle)

            lhs = append!([z[e...] for e in path], [z[j, i]])
            rhs = [x[e[1]] for e in path]

            coefs = vcat([1.0 for i in lhs], [-1.0 for i in rhs])
            vars = vcat(lhs, rhs)
            vars_value = getvalue.(vars)

            if (sum(coefs .* vars_value) > 0 + params_icc[:cut_tolerance])
               num_cuts += 1
               @lazyconstraint(cb, sum(lhs) <= sum(rhs))
               print(".")
            end
         end
      end

      return num_cuts
   end

   function cycle_elimination_cut(cb; is_cut=false)
      ub = min(app["upcutoff"], MathProgBase.cbgetobj(cb))
      gap = round((1.0 - cbgetnodeobjval(cb) / ub) * 100.0; digits=2)
      bestbound=MathProgBase.cbgetbestbound(cb)
   
      # if (is_cut && gap < 0.01)
      if is_cut && stop_cycle_cut
         return 0
      end

      num_cuts = 0
   
      x_val = getvalue.(x[data.nodes])
      z_val = zeros(Float64, length(data.nodes), length(data.nodes))
      w_val = zeros(Float64, length(data.nodes), length(data.nodes))

      for i in data.nodes
         for j in N(i)
            z_ij = getvalue(z[i, j])
            z_val[i, j] = z_ij
   
            w_val[j, i] = x_val[j] - z_ij
         end
      end

      # is_x_val_integral = isinteger.(x_val) |> all
      # is_z_val_integral = isinteger.(z_val) |> all
   
      # if is_x_val_integral && is_z_val_integral
      @info("Pre Floyd Warshall")
      @timeit to "Floyd Warshall" ds = floyd_warshall_shortest_paths(data.graph, w_val)
      show(to)
      println()
   
      has_cut = false
      # @timeit to "Cycles Elimination - Outer" begin
      for k in data.nodes

         # @timeit to "Cycles Elimination - Inner" begin
         for j in N(k)
            # @timeit to "A_star" path = 
            a_star(data.graph, k, j, w_val)
            if (!isempty(path))
               if (ds.dists[k, j] + w_val[j, k] < x_val[k] - params_icc[:cut_tolerance])
                  lhs = append!([z[e.dst, e.src] for e in path], [z[k, j]])
                  rhs = [x[e.dst] for e in path]
   
                  # C = append!([k], [e.dst for e in path])
   
                  coefs = vcat([1.0 for i in lhs], [-1.0 for i in rhs])
   
                  vars = vcat(lhs, rhs)
                  vars_value = getvalue.(vars)
                  
   
                  if (sum(coefs .* vars_value) > 0 + params_icc[:cut_tolerance])
                     # has_cut = true
                     num_cuts += 1
                     if typeof(cb) == Model
                        @constraint(cb, sum(lhs) <= sum(rhs))
                        s = solve(cb)
                     else
                        if is_cut
                           @usercut(cb, sum(lhs) <= sum(rhs))
                        else
                           @lazyconstraint(cb, sum(lhs) <= sum(rhs))
                        end
                     end
                     # @info "Violated generalized cycle."
                     # print(".")

                     current_time = time()
                     delta = current_time - cycle_cut_init_time
                     elapsed_secs = floor(delta)

                     if elapsed_secs % 2 == 0
                        @info "Elapsed secs: $elapsed_secs | Violated generalized cycles: $num_cuts | is_cut: $is_cut | LB: $(floor(bestbound)) | Gap: $gap"
                        # show(to)
                        # println()
                     end
                  end
               end
            end
         end
         # end
      end
      # end

      if is_cut && num_cuts == 0
         stop_cycle_cut = true
      end
   
      return num_cuts
   end
 
   # functions to compute the _y_ and _z_ values
   _x_(i) = getvalue(x[i])
   _y_(i, p) = getvalue(y[i, p])
   _z_(i, j) = getvalue(z[i, j])

   cut_rounds_max_time = get(app[command], "cut-rounds-max-time", 5*60) | params_icc[:cut_rounds_max_time]
   cut_rounds = get(app[command], "cut-rounds", 200) | params_icc[:cut_rounds_init]

   cut_init_time = time()

   function cover_propagation_elimination_cut(cb; is_cut=false, is_x_k=true)
      # @info "Cover Cut | x_k: $is_x_k"
      current_time = time()
      delta = current_time - cut_init_time
      elapsed_secs = floor(delta)
   
      ub = min(app["upcutoff"], MathProgBase.cbgetobj(cb))
      gap = 1.0 - cbgetnodeobjval(cb) / ub
   
      # stop adding cuts after `cut_rounds_max_time` seconds (or after `cut_rounds` cut rounds)
      # if (is_x_k && (elapsed_secs >= cut_rounds_max_time || cut_rounds <= 0)) || (gap < params_icc[:gap_limit])
      if (elapsed_secs >= cut_rounds_max_time || cut_rounds <= 0) || (gap < params_icc[:gap_limit])
         cut_rounds = 0
         return 0
      end

      @show elapsed_secs, cut_rounds, gap
   
      cut_rounds -= 1
   
      num_cuts = 0
   
      V_set = filter(x_k -> getvalue(x[x_k]) > params_icc[:cut_tolerance], data.nodes)
      # @show length(V_set)
      for k in V_set
         x_k = getvalue(x[k])

         # ##
         # for i in data.nodes
         #    temp = [(i, p, getvalue(y[i, p])) for p in P(i) if getvalue(y[i, p]) > 0 && p > 1]
         #    !isempty(temp) && @info temp 
         # end
         # ## 
   
         # Build and solve the separation MIP
         sep = Model(solver=CplexSolver(CPX_PARAM_THREADS=params_cplex[:threads],
            CPX_PARAM_CUTUP=is_x_k ? x_k : 1.0,
            CPX_PARAM_TILIM=2.0,
            CPX_PARAM_MIPDISPLAY=0,
            CPX_PARAM_SCRIND=0))
         @variables(sep, begin
            s[i in data.nodes], Bin
            z1[i in data.nodes, j in N(i)], Bin
            z0[i in data.nodes, j in N(i)], Bin
            y1[i in data.nodes, p in P(i)], Bin
            y0[i in data.nodes, p in P(i)], Bin
         end)
   
         @objective(sep, Min, sum(_y_(i, p) * y1[i, p] for i in data.nodes, p in P(i)) +
                              sum(_z_(i, j) * z1[i, j] for i in data.nodes, j in N(i)))
   
         @constraints(sep, begin
            # cover[i in data.nodes], sum(INCENTIVES[p]*y0[i,p] for p in P(i)) +
            #                         sum(d(i,j)*z0[i,j] for j in N(i)) <= h(i) - 1
            cover[i in data.nodes], sum(incentive_gamma(i, p) * y0[i, p] for p in P(i)) +
                                    sum(d(i, j) * z0[i, j] for j in N(i)) <= ceil(h_gamma(i)) - 1
            forcing_z1[i in data.nodes, j in N(i)], z1[i, j] >= s[i] - s[j] - z0[i, j]
            forcing_y1[i in data.nodes, p in P(i)], y1[i, p] >= s[i] -
                                                               sum(y0[i, q] for q in P(i) if q >= p)
            # forcing_sk, s[k] == 1
         end)
   
         if !is_x_k
            @constraint(sep, separation_x, sum(s) >= floor((1.0 - app["alpha"]) * (length(data.nodes) + 0.000001)) + 1.0)
         else
            @constraint(sep, forcing_sk, s[k] == 1)
   
            if !has_step2
               @constraint(sep, separation_x, sum(s) <= floor((1.0 - app["alpha"]) * (length(data.nodes) + 0.000001)))
            end
   
         end
      #  @info "Build separation model"
   
         status = solve(sep, suppress_warnings=true)

         if status != :Optimal
            continue
         end

         @info "Solve separation model | Node $k | Status: $status"
   
         # check the violation
         lhs_val = getobjectivevalue(sep)
   
         rhs_val = is_x_k ? x_k : 1.0
         # @show lhs_val, rhs_val
         if (lhs_val < rhs_val - params_icc[:cut_tolerance])
            # compute the cut coefficients to be lifted (their coefficients set to zero)
            X = Set([i for i in data.nodes if getvalue(s[i]) > 0.999])
            # @info "Pre lift 1"
            lifted_z = [
               [j for j in N(i) if !(j in X) && getvalue(z0[i, j]) > 0.001]
               for i in data.nodes
            ]
            lifted_y = [
               [p for p in P(i) if sum(getvalue(y0[i, q]) for q in P(i) if q >= p) > 0.001]
               for i in data.nodes
            ]
            # @info "Pre lift 2"
   
            for i in X
               # try to fill the current slack of i with more liftings
               slack = ceil(h_gamma(i))
               # slack = h(i)
               if !isempty(lifted_y[i])
                  slack -= maximum([incentive_gamma(i, p) for p in lifted_y[i]])
               end
               if !isempty(lifted_z[i])
                  slack -= sum(d(i, j) for j in lifted_z[i])
               end
               while slack > 1
                  cand = [j for j in N(i)
                        if !(j in X) && (d(i, j) < slack) && !(j in lifted_z[i])
                  ]
                  if isempty(cand)
                     break
                  else
                     l = cand[argmax([d(i, j) for j in cand])]
                     push!(lifted_z[i], l)
                     slack -= d(i, l)
                  end
               end
               slack += maximum([incentive_gamma(i, p) for p in lifted_y[i]])
               lifted_y[i] = [p for p in P(i) if incentive_gamma(i, p) < slack]
            end
            # @info "Pos lift 1"
            lhs_y = [(i, p) for i in X for p in P(i) if !(p in lifted_y[i])]
            lhs_z = [(i, j) for i in X for j in N(i)
                     if !(j in X) && !(j in lifted_z[i])
            ]
            # @info "Pos lift 2"
            if isempty(lhs_y)
               lhs = sum(z[i, j] for (i, j) in lhs_z)
            elseif isempty(lhs_z)
               lhs = sum(y[i, p] for (i, p) in lhs_y)
            else
               lhs = sum(y[i, p] for (i, p) in lhs_y) + sum(z[i, j] for (i, j) in lhs_z)
            end
   
            rhs = is_x_k ? x[k] : 1.0
            # rhs = x[k]
   
            # @info("Add generalized_propagation_elimination_cut")
            num_cuts += 1
   
            # Save results to file for analysis later
            # writelp(sep, "sep.lp", genericnames=false)
   
            if typeof(cb) == Model
               @constraint(cb, lhs >= rhs)
               status = solve(cb)
            else
               if is_cut
                  @usercut(cb, lhs >= rhs)
               else
                  @lazyconstraint(cb, lhs >= rhs)
               end
            end
            # @info "Added cut to the model | is_x_k? $(is_x_k)"
         end
   
         (!is_x_k) && break
      end
   
      if num_cuts == 0
         cut_rounds = 0 # 10
      end
   
      # if elapsed_secs % 30 == 0
      # @show num_cuts, cut_rounds
      # end
   
      return num_cuts
   end
 
   function cover_propagation_elimination_cut2(cb; is_cut=false, is_x_k=true)
      current_time = time()
      delta = current_time - cut_init_time
      elapsed_secs = floor(delta)
      # @info "Cover Cut 2 | x_k: $is_x_k"

      ub = min(app["upcutoff"], MathProgBase.cbgetobj(cb))
      gap = 1.0 - cbgetnodeobjval(cb) / ub

      # TODO: replicar condição de parar de add cortes similar ao corte anterior
      # stop adding cuts after gap < gap_limit
      # if (gap < params_icc[:gap_limit]) || (is_x_k && gap < params_icc[:gap_limit2]) # elapsed_secs >= cut_rounds_max_time2
      if (elapsed_secs >= cut_rounds_max_time || cut_rounds <= 0) || (gap < params_icc[:gap_limit])
         return 0
      end

      cut_rounds -= 1

      num_cuts = 0
      done = false
      
      V_set = filter(x_k -> getvalue(x[x_k]) > params_icc[:cut_tolerance], data.nodes)
      # @show length(V_set)
      for k in V_set
         x_k = getvalue(x[k])

         # Build and solve the separation MIP
         sep = Model(solver=CplexSolver(CPX_PARAM_THREADS=params_cplex[:threads],
                                       CPX_PARAM_CUTUP= is_x_k ? x_k : 1.0,
                                       CPX_PARAM_TILIM=2.0,
                                       CPX_PARAM_MIPDISPLAY=0,
                                       CPX_PARAM_SCRIND=0))
         @variables(sep, begin
            s[i in data.nodes], Bin
            x1[i in data.nodes], Bin
            x0[i in data.nodes], Bin
            y1[i in data.nodes, p in P(i)], Bin
            y0[i in data.nodes, p in P(i)], Bin
         end)

         @objective(sep, Min, sum(_y_(i,p) * y1[i,p] for i in data.nodes, p in P(i)) +
                              sum(_x_(j) * x1[j] for j in data.nodes))

         @constraints(sep, begin
            cover[i in data.nodes], sum(incentive_gamma(i,p)*y0[i,p] for p in P(i)) +
                                    sum(d(i,j)*x0[j] for j in N(i)) +
                                    sum(d(i,j) for j in N(i))*s[i] <=
                                    sum(d(i,j) for j in N(i)) + ceil(h_gamma(i)) - 1
            forcing_x1[j in data.nodes], x1[j] >= 1 - s[j] - x0[j]
            forcing_y1[i in data.nodes, p in P(i)], y1[i,p] >= s[i] -
                                                   sum(y0[i,q] for q in P(i) if q >= p)
         end)

         # |X| >= floor( (1 - alpha) |V|) + 1
         if !is_x_k
            @constraint(sep, separation_x, sum(s) >= floor((1.0 - app["alpha"]) * (length(data.nodes) + 0.000001)) + 1.0)
         else
            @constraint(sep, separation_x, sum(s) <= floor((1.0 - app["alpha"]) * (length(data.nodes) + 0.000001)))
            @constraint(sep, forcing_sk, s[k] == 1)
         end

         done = true

         status = solve(sep, suppress_warnings=true)
         @info "Solve separation model | Node $k | Status: $status | Cut rounds: $cut_rounds"
         if status != :Optimal
            (!is_x_k) && break

            continue
         end

         # check the violation
         lhs_val = getobjectivevalue(sep)

         # compute the cut coefficients to be lifted (their coefficients set to zero)
         X = Set([i for i in data.nodes if getvalue(s[i]) > 0.999])
         _X_ = [i for i in data.nodes if !(i in X)]

         lifted_x = [j for j in data.nodes if !(j in X) && getvalue(x0[j]) > 0.001]

         lifted_y = [
            [p for p in P(i) if sum(getvalue(y0[i,q]) for q in P(i) if q >= p) > 0.001]
            for i in data.nodes
         ]

         slack = Vector{Float64}(undef, length(X))
         for (k, i) in enumerate(X)
            # try to fill the current slack of i with more liftings
            slack[k] = ceil(h_gamma(i))
            # slack = h(i)

            if !isempty(lifted_y[i])
               slack[k] -= maximum([incentive_gamma(i,p) for p in lifted_y[i]])
            end
            if !isempty(lifted_x)
               slack[k] -= sum(d(i,j) for j in lifted_x)
            end
         end

         while maximum(slack) > 1.0
            cand = [j for j in _X_
               if reduce(*, [d(i,j) for i in X] .< slack) && !(j in lifted_x)
            ]

            if isempty(cand)
               break
            else
               l = cand[argmax([sum(d(i,j) for i in X) for j in cand])]
               push!(lifted_x, l)
               slack .-= [d(i,l) for i in X]
            end
         end

         for (k, i) in enumerate(X)
            slack[k] += maximum([incentive_gamma(i,p) for p in lifted_y[i]])
            lifted_y[i] = [p for p in P(i) if incentive_gamma(i,p) < slack[k]]
         end

         lhs_y = [(i,p) for i in X for p in P(i) if !(p in lifted_y[i])]
         lhs_x = [i for i in _X_ if !(i in lifted_x)]

         if isempty(lhs_y)
            lhs = sum(x[i] for i in lhs_x)
         elseif isempty(lhs_x)
            lhs = sum(y[i,p] for (i,p) in lhs_y)
         else
            lhs = sum(y[i,p] for (i,p) in lhs_y) + sum(x[i] for i in lhs_x)
         end

         rhs = is_x_k ? x_k : 1.0
         
         @show lhs_val, rhs
         if(lhs_val < rhs - params_icc[:cut_tolerance])
            if isempty(lhs_y)
               lhs = sum(x[j] for j in lhs_x)
            elseif isempty(lhs_x)
               lhs = sum(y[i,p] for (i,p) in lhs_y)
            else
               lhs = sum(y[i,p] for (i,p) in lhs_y) + sum(x[j] for j in lhs_x)
            end

            num_cuts += 1

            if typeof(cb) == Model
               if is_x_k
                  @constraint(cb, lhs >= x[k])
               else
                  @constraint(cb, lhs >= 1.0)
               end
               status = solve(cb)
            else
               if is_cut
                  if is_x_k
                     @usercut(cb, lhs >= x[k])
                  else
                     @usercut(cb, lhs >= 1.0)
                  end
               else
                  if is_x_k
                     @lazyconstraint(cb, lhs >= x[k])
                  else
                     @lazyconstraint(cb, lhs >= 1.0)
                  end
               end
            end

         end

         (!is_x_k) && break
      end

      if num_cuts == 0
         cut_rounds = 0 # 10
      end

      return num_cuts
   end
 
   function separate_cut(cb)   
      # cycle_cuts = cycle_elimination_cut(cb, is_cut=true)

      cover_cuts = 0
      if "cover_cut" == app[command]["cuts"]
         # cover_cuts = cover_propagation_elimination_cut(cb, is_cut=true)

         if app[command]["step2"]
            cover_cuts += cover_propagation_elimination_cut(cb, is_cut=true, is_x_k=false)
         end
      end

      cover_cuts2 = 0
      if "cover_cut2" == app[command]["cuts"]
         # @info "cover_cut2"
         # cover_cuts2 = cover_propagation_elimination_cut2(cb, is_cut=true)
         cover_cuts2 += cover_propagation_elimination_cut2(cb, is_cut=true, is_x_k=false)
      end

      if  cover_cuts > 0 || cover_cuts2 > 0
         # @info("Found $cycle_cuts cycle cuts")

         # tell CPLEX to keep calling this routine in this node
         unsafe_store!(cb.userinteraction_p, convert(Cint,2), 1)
      end
   end
 
   # set the callback cut separation function
   addcutcallback(model, separate_cut)
   # addlazycallback(model, cycle_elimination_cut)
   addlazycallback(model, simple_cycle_elimination_lazy_constraint)
   
   
   variables = Dict(
      :x => x,
      :y => y,
      :z => z
   )
     
   if app["export-cplex-log"]
      cplex_output_filename = joinpath(dirname(appfolder), "out", "CPLEX_output.log")
      touch(cplex_output_filename)

      JuMP.build(model)
      CPLEX.set_logfile(model.internalModel.inner.env, cplex_output_filename)

      @info "Export CPLEX log to file: $cplex_output_filename"
   end

   if app["export-cplex-log"]
      lp_filename = joinpath(dirname(appfolder), "out", "glcip.lp")
      writeLP(model, lp_filename, genericnames=false)

      @info "Export .lp to file: $lp_filename"
   end
 
   return model, variables, P, W
end