const params = Dict(
   :ub => 10000000.0, # upper bound
   :alpha => 0.1, # nodes that need to be activated
   :gamma => 1.0, # activation function (0.9, 1.0, 1.1)
)

const params_icc = Dict(
   :gap_limit => 0.0,
   :gap_limit2 => 0.4,
   :cut_rounds_init => 200,
   :cut_rounds_max_time => 5*60, # seconds

   # :cut_tolerance => 0.1, # cut tol
   :cut_tolerance => 0.01, # cut tol
   :cut_multiplier => 1000000, # new cut multiplier
)

const params_cf = Dict(
   :cut_gap => 0.1 # minimum node gap to add cuts
)

const params_cplex = Dict(
   :time_limit => 72000, # 7200
   :up_cutoff => 10000000.0,
   :mip_display => 4,
   :mip_interval => 200,
   :scrind => 1,
   :threads => 1,
)
