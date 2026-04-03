const PROJECT_ROOT = @__DIR__
const SRC_DIR      = joinpath(PROJECT_ROOT, "src")
const INST_DIR     = joinpath(PROJECT_ROOT, "instances")

include(joinpath(SRC_DIR, "parser.jl"))
include(joinpath(SRC_DIR, "model_mtz.jl"))
include(joinpath(SRC_DIR, "model_sec.jl"))
include(joinpath(SRC_DIR, "heuristic.jl"))
include(joinpath(SRC_DIR, "plot_results.jl"))

using Printf

const DEFAULT_INSTANCES = [
    "instance_6_1.txt",
    "instance_40_1.txt",
    "instance_80_1.txt",
]

function parse_args(args)
    opts = Dict(
        :run_mtz       => true,
        :run_sec       => true,
        :run_heuristic => true,
        :verbose       => false,
        :time_limit    => 120.0,
        :files         => String[],
    )

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--mtz"
            opts[:run_mtz] = true
            opts[:run_sec] = false
            opts[:run_heuristic] = false
        elseif arg == "--sec"
            opts[:run_mtz] = false
            opts[:run_sec] = true
            opts[:run_heuristic] = false
        elseif arg == "--heuristic"
            opts[:run_mtz] = false
            opts[:run_sec] = false
            opts[:run_heuristic] = true
        elseif arg == "--all"
            opts[:run_mtz] = true
            opts[:run_sec] = true
            opts[:run_heuristic] = true
        elseif arg == "--verbose"
            opts[:verbose] = true
        elseif arg == "--time"
            i += 1
            opts[:time_limit] = parse(Float64, args[i])
        elseif startswith(arg, "--")
            println("Unknown option: $arg")
        else
            push!(opts[:files], arg)
        end
        i += 1
    end

    return opts
end

function run_instance(filepath::String, opts::Dict)
    println()
    println("=" ^ 72)
    println("Instance: $filepath")
    println("=" ^ 72)

    # --- Parse ---
    inst = try
        parse_instance(filepath)
    catch e
        println("ERROR parsing instance: $e")
        return
    end
    print_instance_summary(inst)

    d = compute_distances(inst)

    # Validate connectivity briefly
    _, out_n, in_n = build_arc_sets(inst, d)
    println("  Arcs within range R=$(inst.R): $(sum(length.(out_n)))")
    println()

    # Results accumulator
    rows = Any[]

    # --- MTZ ---
    if opts[:run_mtz]
        print("Running MTZ model ... ")
        flush(stdout)
        res = solve_mtz(inst;
                        time_limit = opts[:time_limit],
                        verbose    = opts[:verbose])
        println("done ($(round(res.time_sec, digits=2)) s)")

        push!(rows, (
            method  = "MTZ (ILP)",
            obj     = res.objective,
            status  = res.status,
            time_s  = res.time_sec,
            bb      = res.bb_nodes,
            cuts    = "N/A",
            path    = res.path,
        ))

        if !isempty(res.path)
            print("  Path: ")
            println(join(res.path, " -> "))
            verify_solution(res.path, inst, d)
        end
        println()
    end

    # --- SEC ---
    if opts[:run_sec]
        print("Running SEC model ... ")
        flush(stdout)
        res = solve_sec(inst;
                        time_limit = opts[:time_limit],
                        verbose    = opts[:verbose])
        println("done ($(round(res.time_sec, digits=2)) s)")

        push!(rows, (
            method  = "SEC (ILP)",
            obj     = res.objective,
            status  = res.status,
            time_s  = res.time_sec,
            bb      = res.bb_nodes,
            cuts    = res.sec_cuts,
            path    = res.path,
        ))

        if !isempty(res.path)
            print("  Path: ")
            println(join(res.path, " -> "))
            verify_solution(res.path, inst, d)
        end
        println()
    end

    # --- Heuristic ---
    if opts[:run_heuristic]
        print("Running Heuristic ... ")
        flush(stdout)
        res = solve_heuristic(inst; verbose=opts[:verbose])
        println("done ($(round(res.time_sec, digits=4)) s)")

        push!(rows, (
            method  = "Heuristic",
            obj     = res.objective,
            status  = res.objective < Inf ? "FEASIBLE" : "NO_SOLUTION",
            time_s  = res.time_sec,
            bb      = "N/A",
            cuts    = "N/A",
            path    = res.path,
        ))

        if !isempty(res.path)
            print("  Path: ")
            println(join(res.path, " -> "))
            verify_solution(res.path, inst, d)
        end
        println()
    end

    # --- Summary table ---
    if !isempty(rows)
        println()
        println("Summary table")
        println("-" ^ 72)
        @printf("%-18s  %10s  %18s  %10s  %10s  %6s\n",
                "Method", "Objective", "Status", "Time (s)", "B&B nodes", "SEC cuts")
        println("-" ^ 72)
        for r in rows
            obj_str = isinf(r.obj) ? "---" : string(Int(r.obj))
            @printf("%-18s  %10s  %18s  %10.3f  %10s  %6s\n",
                    r.method,
                    obj_str,
                    r.status[1:min(18,end)],
                    r.time_s,
                    string(r.bb),
                    string(r.cuts))
        end
        println("-" ^ 72)

        # Report optimality gap between heuristic and best ILP
        ilp_objs = [r.obj for r in rows if r.method != "Heuristic" && !isinf(r.obj)]
        heur_rows = [r.obj for r in rows if r.method == "Heuristic" && !isinf(r.obj)]
        if !isempty(ilp_objs) && !isempty(heur_rows)
            best_ilp  = minimum(ilp_objs)
            heur_obj  = heur_rows[1]
            gap_pct   = 100.0 * (heur_obj - best_ilp) / best_ilp
            @printf("  Heuristic gap vs best ILP: %+.2f%%\n", gap_pct)
        end

        # Save plot
        get_path(method) = begin
            r = findfirst(x -> x.method == method, rows)
            r === nothing ? Int[] : rows[r].path
        end
        get_obj(method) = begin
            r = findfirst(x -> x.method == method, rows)
            r === nothing ? Inf : rows[r].obj
        end
        name = splitext(basename(filepath))[1]
        save_plot(inst, d,
                  get_path("MTZ (ILP)"),  get_obj("MTZ (ILP)"),
                  get_path("SEC (ILP)"),  get_obj("SEC (ILP)"),
                  get_path("Heuristic"),  get_obj("Heuristic"),
                  PROJECT_ROOT, name)
    end
end

function main()
    opts = parse_args(ARGS)

    files = opts[:files]

    if isempty(files)
        # Run all default instances
        for fname in DEFAULT_INSTANCES
            fpath = joinpath(INST_DIR, fname)
            if isfile(fpath)
                run_instance(fpath, opts)
            else
                println("Instance file not found: $fpath")
            end
        end
    else
        for fpath in files
            # Allow relative or absolute paths
            resolved = isfile(fpath) ? fpath : joinpath(PROJECT_ROOT, fpath)
            if isfile(resolved)
                run_instance(resolved, opts)
            else
                println("Instance file not found: $fpath")
            end
        end
    end

    println()
    println("Complete.")
end

main()
