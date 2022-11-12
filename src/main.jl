using ArgParse
using ProgressMeter

using LightGraphs
using TikzGraphs, TikzPictures

include("cli.jl")
include("data.jl")
include("models/params.jl")
include("models/run.jl")
include("solution.jl")

appfolder = dirname(@__FILE__)

function run_glcip(app, input_file)
    data = read_glcip_data(input_file, get_incentives)

    command = app["%COMMAND%"]
    if command == "icc"
        @info("<== ICC ==>")

        app[command]["cuts"] = "cover_cut"
        app[command]["step2"] = false
    elseif command == "icc+"
        @info("<== ICC+ ==>")

        app[command]["cuts"] = "cover_cut"
        app[command]["step2"] = true
    elseif command == "licc+"
        @info("<== LICC+ ==>")

        app[command]["cuts"] = "cover_cut2"
        app[command]["step2"] = true
    elseif command == "cf"
        @info("<== CF ==>")

        app[command]["cuts"] = "cover_cut2"
    end

    solution = solve_glcip(data, input_file, app)

    return solution
end

function run(app::Dict{String,Any})
    @info("Application parameters:")
    for (arg, val) in app
        @info("  $arg  =>  $(repr(val))")
    end
    flush(stdout)

    # # DEPRECATED
    # if haskey(app, "batch")
    #     # dir_files = readdir(app["filepath"])
    #     dir_files = filter(x->occursin("SW", x), readdir(app["filepath"]))

    #     # Iterate over files
    #     p = Progress(length(dir_files), 1)
    #     for (i, filename) in enumerate(dir_files)
    #         input_file = joinpath(app["filepath"], filename)
    #         @info(input_file)

    #         sol = run_glcip(app, input_file)

    #         # row = [sol.instance, sol.objective, sol.status, sol.time]
    #         # # TODO: filename as a parameter?
    #         # open("models_track.csv","a") do fp
    #         #     println(fp, join(row, ", "))
    #         # end

    #         update!(p, i)
    #     end
    # else
    # end

    sol = run_glcip(app, app["filepath"])

    # Output 
    println(sol)

    base_filename = "$(basename("$(app["filepath"])"))-a$(sol.alpha)-g$(app["gamma"])-$(app["%COMMAND%"])"

    output_solution_dir = joinpath(appfolder, "..", "out")
    mkpath(output_solution_dir)

    # Export solution output
    output_solution_path = joinpath(output_solution_dir, base_filename)
    export_solution(sol, output_solution_path)

    # Export a graph representation of the solution
    if app["latex"]
        output_latex_dir = joinpath(output_solution_dir, "latex")
        mkpath(output_latex_dir)

        output_latex_dir_gamma = joinpath(output_latex_dir, string(app["gamma"]))
        mkpath(output_latex_dir_gamma)

        output_latex_path = joinpath(output_latex_dir_gamma, base_filename)

        draw_solution(sol, output_latex_path)
    end
end

function main(args)
    app = parse_commandline(args)

    app === nothing && return
    run(app)
end

if isempty(ARGS)
    # main(["--help)

    # main(["data/socnet-instances-v2/SW-n50-k4-b0.1-d1-10-g0.7-i1", "-a", "0.1", "-g", "1.0", "icc", "-r", "2000"])
    # main(["data/socnet-instances-v2/SW-n50-k4-b0.1-d1-10-g0.7-i1", "-l", "cf"])
    # main(["data/socnet-instances-v2/SW-n50-k8-b0.1-d1-10-g0.7-i5", "-a", "0.1", "-u", "59",  "cf"]) # 58
    # main(["data/socnet-instances-v2/SW-n50-k8-b0.1-d1-10-g0.7-i2", "-a", "0.5", "-u", "51", cf"]) # 50
    main(["data/socnet-instances-v2/SW-n50-k8-b0.3-d1-10-g0.7-i1", "-a", "0.1", "-u", "45", "cf"]) # 44
    # main(["data/socnet-instances-v2/SW-n75-k8-b0.1-d1-10-g0.7-i5", "-a", "1.0", "-u", "122", "cf"]) # 

    # arg_string = "data/socnet-instances-v2/SW-n50-k4-b0.1-d1-10-g0.7-i1 -a 0.1 -g 1.0 icc -r 2000"
    # main(split(arg_string))
else
    main(ARGS)
end
