function cf_model(data::DataGlcip, app)
    verbose = Int(app["verbose"]) | params_cplex[:scrind]
 
    # Create Model
    model = Model(solver=CplexSolver(CPX_PARAM_THREADS=params_cplex[:threads],
       CPX_PARAM_CUTUP=app["upcutoff"],
       #CPX_PARAM_CUTPASS=-1,
       CPX_PARAM_TILIM=params_cplex[:time_limit],
       CPX_PARAM_MIPDISPLAY=params_cplex[:mip_display],
       #CPX_PARAM_MIPINTERVAL=params_cplex[:mip_interval],
       CPX_PARAM_SCRIND=verbose))
 
    @variables(model, begin
       # 0 <= x[i in data.nodes] <= 1
       # 0 <= z[i in data.nodes, j in N(i)] <= 1
       y[i in data.nodes, p in P(i)], Bin
    end)
 
    @objective(model, Min, sum(W[(i, p)] * y[i, p] for i in data.nodes, p in P(i)))
 
    @constraints(model, begin
       incentive[i in data.nodes], sum([y[i, p] for p in P(i)]) == 1 # x[i]
    end)
 
    return model
 end