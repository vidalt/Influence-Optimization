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

   :cut_tolerance => 0.1, # cut tol
   :cut_multiplier => 1000000, # new cut multiplier
   :cut_gap => 0.1 # minimum node gap to add cuts

   # hamilton | Dev
   # :cut_tolerance => 0.1, # cut tol
   # :cut_multiplier => 1000000, # new cut multiplier
   # :rounding_cuts => [6], # C+ϵ used in the rouding cut e.g. 6+ϵ
   # :rounding_epsilon => 0.001, # 1+ϵ used in the new cut
   #
   # # :cut => "", # Cut useds
   # :cuts => [], # Cuts used
   #
   # # Covercut1
   # :cut_rounds_init => 2000,
   # :cut_rounds_max_time => 3600, # Seconds
   #
   #  # Covercut2
   #  :cut_rounds_max_time2 => 3600, # Seconds
   #  :cut_rounds2_gap => 5 # (%)
)

const params_cplex = Dict(
   :time_limit => 7200,
   :up_cutoff => 10000000.0,
   :mip_display => 4,
   :mip_interval => 1,
   :scrind => 1,
   :threads => 1,
)
