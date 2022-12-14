function choices(list)
    if length(list) >= 3
        return join(list, ", ", " or ")
    else
        return join(list, " or ")
    end
end

function parse_commandline(args::Array{String,1})
    # usage = "GLCIP executable \n" *
    #             "  On interactive mode, call main([\"arg1\", ..., \"argn\"])",
    s = ArgParseSettings()

    gamma_list = [0.9, 1.0, 1.1]

    # Main commands
    @add_arg_table! s begin
        "--output", "-o"
            help = "Output filename"
        "--verbose", "-v"
            help = "Verbose output"
            action = :store_true
        # "--batch", "-b"
        #     help = "Batch"
        #     action = :store_true
        "--latex", "-l"
            help = "Latex Output"
            action = :store_true
        "--upcutoff", "-u"
            help = "Up cutoff"
            arg_type = Float64
            default = params[:ub]
        "--alpha", "-a"
            help = "α|V| nodes that need to be activated"
            range_tester = (x -> x > 0.0 && x <= 1.0)
            arg_type = Float64
            default = params[:alpha]
        "--gamma", "-g"
            help = "Activation function Γ: " * choices(gamma_list)
            range_tester = (x -> x in gamma_list)
            arg_type = Float64
            default = params[:gamma]
        "filepath"
            help = "Instance file path or path directory"
            default = "/data/socnet-instances-v2/SW-n50-k4-b0.1-d1-10-g0.7-i3"
            required = true
        "icc"
            help = "ICC"
            action = :command
        "icc+"
            help = "ICC+"
            action = :command
        "licc+"
            help = "LICC+"
            action = :command
        "cf"
            help = "CF"
            action = :command
    end

    # Formulations params
    @add_arg_table! s["icc"] begin
        "--cut-rounds", "-r"
            help = "Cover-Cut rounds"
            arg_type = Int
            default = params_icc[:cut_rounds_init]
        "--cut-rounds-max-time"
            help = "Cover-Cut max time (s)"
            arg_type = Int
            default = params_icc[:cut_rounds_max_time]
    end

    @add_arg_table! s["icc+"] begin
        "--cut-rounds", "-r"
            help = "Cover-Cut rounds"
            arg_type = Int
            default = params_icc[:cut_rounds_init]
    end

    # @add_arg_table! s["licc+"] begin
    # end

    # @add_arg_table! s["cf"] begin
    # end

    return parse_args(args, s)
end