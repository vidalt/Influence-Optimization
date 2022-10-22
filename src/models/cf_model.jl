function cf_model(data::DataGlcip, app)
   # Parameters
   command = app["%COMMAND%"]
   verbose = Int(app["verbose"]) | params_cplex[:scrind]

   # Notations
   H = data.H
   INCENTIVES = data.incentives

   N = i_neighborhood(data)
   P(i::Int) = data.min_incentives[i]
   h(i::Int) = data.hurdle[i] .- 0.5
   d(i::Int, j::Int) = influence(data, j, i)
   h_gamma(i::Int) = h(i)^(1 / app["gamma"])
   W = data.weights

   # Initialization
   prevobj = -1.0
   countobj = 0

   function h_minus_p(i, p)
      temp = h(i) - INCENTIVES[p]

      return max(temp, 0)^(1 / app["gamma"])
   end

   function incentive_gamma(i, p)
      return ceil(h_gamma(i)) - ceil(h_minus_p(i, p))
   end

   function getNonactivated(p_)
      # start with all individuals in X and no influence
      resistence = [max(0, Int(ceil(h_gamma(i)) - incentive_gamma(i, p_[i]))) for i in data.nodes]
      new_act = [(resistence[i] == 0) for i in data.nodes]

      while count(new_act) > 0
         # propagate all activations
         act = new_act
         new_act = zeros(Bool, length(data.nodes))
         for i in data.nodes
            if resistence[i] > 0
               for j in N(i)
                  if act[j]
                     resistence[i] -= d(i, j)
                     if resistence[i] <= 0
                        resistence[i] = 0
                        new_act[i] = true
                     end
                  end
               end
            end
         end
      end

      nonactivated = Set([i for i in data.nodes if resistence[i] > 0])
      # activated = [i for i in data.nodes if !(i in nonactivated)]

      if length(data.nodes) - length(nonactivated) >= ceil(app["alpha"] * (length(data.nodes) - 0.000001))
         # cost = sum(W[(i, p_[i])] for i in data.nodes)
         # @info("FEASIBLE SOLUTIONs ($(p_), $(activated)) WITH COST $cost")
         return Set()
      end

      return nonactivated
   end

   # Create Model
   model = Model(solver=CplexSolver(CPX_PARAM_THREADS=params_cplex[:threads],
      CPX_PARAM_CUTUP=app["upcutoff"],
      #CPX_PARAM_CUTPASS=-1,
      CPX_PARAM_TILIM=params_cplex[:time_limit],
      CPX_PARAM_MIPDISPLAY=params_cplex[:mip_display],
      #CPX_PARAM_MIPINTERVAL=params_cplex[:mip_interval],
      CPX_PARAM_SCRIND=verbose))

   @variables(model, begin
      y[i in data.nodes, p in P(i)], Bin
   end)

   @objective(model, Min, sum(W[(i, p)] * y[i, p] for i in data.nodes, p in P(i)))

   @constraints(model, begin
      incentive[i in data.nodes], sum([y[i, p] for p in P(i)]) == 1 # x[i]
   end)

   _y_(i, p) = getvalue(y[i, p])

   function cover_propagation_elimination_cut2(cb; is_cut=false)
      # current_time = time()
      # delta = current_time - cut_init_time

      ub = min(app["upcutoff"], MathProgBase.cbgetobj(cb))

      gap = 1.0 - cbgetnodeobjval(cb) / ub
      if is_cut && gap < params_cf[:cut_gap]
         return 0
      end

      if is_cut
         currobj = cbgetnodeobjval(cb)
         if abs(currobj - prevobj) < 1e-4
            countobj += 1
            if countobj >= 100
               return 0
            end
         else
            prevobj = currobj
            countobj = 0
         end
      end

      num_cuts = 0
      status = :Undefined

      is_frac = (count([(_y_(i, p) > 0.001 && _y_(i, p) < 0.999)
                        for i in data.nodes for p in P(i)]) > 0)

      if is_frac
         # Build and solve the separation MIP
         sep = Model(
            solver=CplexSolver(CPX_PARAM_THREADS=params_cplex[:threads],
                                 CPX_PARAM_CUTUP=1.0,
                                 CPX_PARAM_TILIM=(is_cut ? 0.2 : 7200.0), 
                                 CPX_PARAM_MIPDISPLAY=0, 
                                 CPX_PARAM_SCRIND=0
         ))

         @variables(sep, begin
            s0[j in data.nodes], Bin 
            s1[j in data.nodes], Bin 
            y1[i in data.nodes, p in P(i)], Bin
            y0[i in data.nodes, p in P(i)], Bin
         end)

         @objective(sep, Min, sum(_y_(i, p) * y1[i, p] for i in data.nodes, p in P(i)))

         @constraints(sep, begin
            cover[i in data.nodes], sum(incentive_gamma(i, p) * y0[i, p] for p in P(i)) +
                                    sum(d(i, j) * s0[j] for j in N(i)) +
                                    sum(d(i, j) for j in N(i)) * s1[i] <=
                                    sum(d(i, j) for j in N(i)) + ceil(h_gamma(i)) - 1
            forcing_s0[j in data.nodes], s0[j] >= 1 - s1[j]
            forcing_y1[i in data.nodes, p in P(i)], y1[i, p] >= s1[i] -
                                                                  sum(y0[i, q] for q in P(i) if q >= p)
            separation_x, sum(s1) >= floor((1.0 - app["alpha"]) * (length(data.nodes) + 0.000001)) + 1.0 # |X| >= floor( (1 - alpha) |V|) + 1
         end)

         status = solve(sep, suppress_warnings=true)

         # check the violation
         lhs_val = getobjectivevalue(sep)
      else
         X = getNonactivated(
            [p for i in data.nodes for p in P(i) if _y_(i, p) > 0.999])
         lhs_val = isempty(X) ? 1.0 : 0.0
      end

      rhs_val = 1.0

      if (lhs_val < rhs_val - params_icc[:cut_tolerance])
         # compute the cut coefficients to be lifted (their coefficients set to zero)
         if is_frac
            X = Set([i for i in data.nodes if getvalue(s1[i]) > 0.999])
         end
         _X_ = [i for i in data.nodes if !(i in X)]

         #lifted_x = [j for j in data.nodes if !(j in X) && getvalue(x0[j]) < 0.001]
         lifted_x = Int[]

         lifted_y = [
            [p for p in P(i) if ((p == 1) || ((i in X) && (incentive_gamma(i, p) +
                                          (isempty(_X_) ? 0 : sum(d(i, j) for j in _X_))
                                          <= ceil(h_gamma(i)) - 1)))]
            for i in data.nodes
         ]

         lhs_x = [j for j in lifted_x]

         lhs_y = [(i, p) for i in X for p in P(i) if !(p in lifted_y[i])]

         if isempty(lhs_y)
            lhs = sum(x[j] for j in lhs_x)
         elseif isempty(lhs_x)
            lhs = sum(y[i, p] for (i, p) in lhs_y)
         else
            lhs = sum(y[i, p] for (i, p) in lhs_y) + sum(x[j] for j in lhs_x)
         end

         num_cuts += 1

         # Add to the model
         rhs = 1.0
         if typeof(cb) == Model
            @constraint(cb, lhs >= rhs)

            status = solve(cb)
         else
            if is_cut
               @usercut(cb, lhs >= rhs)
            else
               @lazyconstraint(cb, lhs >= 1.0)
            end
         end
      end

      return num_cuts
   end

   function separate_cut(cb)
      cover_cuts2 = 0
      if occursin("cover_cut2", app[command]["cuts"])
         cover_cuts2 = cover_propagation_elimination_cut2(cb, is_cut=true)
      end

      if cover_cuts2 > 0
         # tell CPLEX to keep calling this routine in this node
         unsafe_store!(cb.userinteraction_p, convert(Cint, 2), 1)
      end
   end

   function cut_lazycallback(cb)
      # cycle_elimination_cut(cb)
      cover_propagation_elimination_cut2(cb)
   end

   # set the callback cut separation function
   addcutcallback(model, separate_cut)
   addlazycallback(model, cut_lazycallback)

   # writeLP(model, "glcip.lp", genericnames=false)

   variables = Dict(
      :y => y
   )

   return model, variables
end