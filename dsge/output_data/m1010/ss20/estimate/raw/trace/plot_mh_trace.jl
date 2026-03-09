#!/usr/bin/env julia

using Printf
using Pkg
using HDF5
using DSGE
using Plots

const DEFAULT_INPUT = "../mhsave_vint=250825.h5"
const DEFAULT_OUTPUT = "mh_trace_vint=250825.pdf"
const DEFAULT_SUBSPEC = "ss20"
const DEFAULT_ROWS = 3
const DEFAULT_COLS = 3
const EXPECTED_N_PARAMS = 108

function usage_and_exit(code::Int=0)
    println("Usage: julia plot_mh_trace.jl [--input PATH] [--output PATH] [--subspec STR] [--rows INT] [--cols INT]")
    println("Defaults:")
    println("  --input   $(DEFAULT_INPUT)")
    println("  --output  $(DEFAULT_OUTPUT)")
    println("  --subspec $(DEFAULT_SUBSPEC)")
    println("  --rows    $(DEFAULT_ROWS)")
    println("  --cols    $(DEFAULT_COLS)")
    exit(code)
end

function parse_args(args::Vector{String})
    opts = Dict{String,String}(
        "input"   => DEFAULT_INPUT,
        "output"  => DEFAULT_OUTPUT,
        "subspec" => DEFAULT_SUBSPEC,
        "rows"    => string(DEFAULT_ROWS),
        "cols"    => string(DEFAULT_COLS),
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            usage_and_exit(0)
        end
        if !startswith(arg, "--")
            error("Unexpected positional argument: $arg")
        end
        key = arg[3:end]
        if !haskey(opts, key)
            error("Unknown option: $arg")
        end
        if i == length(args)
            error("Missing value for option: $arg")
        end
        opts[key] = args[i + 1]
        i += 2
    end
    rows = try
        parse(Int, opts["rows"])
    catch
        error("Invalid integer for --rows: $(opts["rows"])")
    end
    cols = try
        parse(Int, opts["cols"])
    catch
        error("Invalid integer for --cols: $(opts["cols"])")
    end
    if rows <= 0 || cols <= 0
        error("--rows and --cols must be positive integers.")
    end
    return opts, rows, cols
end

resolve_path(path::String) = isabspath(path) ? normpath(path) : normpath(joinpath(@__DIR__, path))

function dir_is_writable(path::String)
    probe = joinpath(path, ".write_probe_$(time_ns())")
    try
        open(probe, "w") do io
            write(io, "probe")
        end
        rm(probe; force = true)
        return true
    catch
        return false
    end
end

function dep_field(dep, sym::Symbol, fallback::String="n/a")
    if sym in fieldnames(typeof(dep))
        value = getfield(dep, sym)
        return value === nothing ? fallback : string(value)
    end
    return fallback
end

function pkg_metadata(name::String)
    for (_, dep) in Pkg.dependencies()
        if dep_field(dep, :name, "") == name
            return (
                version = dep_field(dep, :version),
                git_revision = dep_field(dep, :git_revision),
                tree_hash = dep_field(dep, :tree_hash),
            )
        end
    end
    return (version = "n/a", git_revision = "n/a", tree_hash = "n/a")
end

function param_plot_title(p)
    if :tex_label in fieldnames(typeof(p))
        tex = strip(string(getfield(p, :tex_label)))
        if !isempty(tex)
            return "\$$(replace(tex, "\$" => "\\\$"))\$"
        end
    end
    return string(p.key)
end

function trace_subplot(series::AbstractVector{<:Real}, pname::String, n_draws::Int)
    return plot(
        1:n_draws, series;
        linewidth = 0.75,
        color = :steelblue,
        legend = false,
        xlabel = "Draw",
        ylabel = "",
        title = pname,
    )
end

function blank_subplot()
    return plot(legend = false, framestyle = :none, grid = false, axis = false, ticks = false)
end

function write_trace_pdf(input_path::String, output_path::String, subspec::String, rows::Int, cols::Int)
    if !isfile(input_path)
        error("Input file not found: $input_path")
    end
    output_dir = dirname(output_path)
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    if !dir_is_writable(output_dir)
        error("Output directory is not writable: $output_dir")
    end
    if Sys.which("gs") === nothing
        error("Ghostscript executable `gs` not found in PATH.")
    end

    m = Model1010(subspec)
    param_names = [param_plot_title(p) for p in m.parameters]
    n_model_params = length(param_names)
    if n_model_params != EXPECTED_N_PARAMS
        error("Model parameter count mismatch: expected $EXPECTED_N_PARAMS, got $n_model_params")
    end

    mh_raw = h5read(input_path, "mhparams")
    raw_dims = size(mh_raw)

    mh = if size(mh_raw, 1) == n_model_params
        mh_raw
    elseif size(mh_raw, 2) == n_model_params
        permutedims(mh_raw)
    else
        error("Cannot align mhparams dimensions $raw_dims with model parameter count $n_model_params")
    end

    orientation = size(mh_raw, 1) == n_model_params ? "rows_are_parameters" : "columns_are_parameters_transposed"

    n_params, n_draws = size(mh)
    if n_params != EXPECTED_N_PARAMS
        error("Aligned draw matrix has $n_params parameters; expected $EXPECTED_N_PARAMS")
    end
    if n_params != length(param_names)
        error("Parameter name count $(length(param_names)) does not match draw rows $n_params")
    end

    per_page = rows * cols
    n_pages = cld(n_params, per_page)
    page_files = String[]
    page_width = max(1200, 420 * cols)
    page_height = max(900, 260 * rows)

    gr()

    for page in 1:n_pages
        panels = Vector{Any}(undef, per_page)
        base = (page - 1) * per_page
        for slot in 1:per_page
            idx = base + slot
            if idx <= n_params
                panels[slot] = trace_subplot(vec(mh[idx, :]), param_names[idx], n_draws)
            else
                panels[slot] = blank_subplot()
            end
        end
        page_plot = plot(panels...; layout = (rows, cols), size = (page_width, page_height))
        page_file = joinpath(output_dir, @sprintf(".mh_trace_page_%03d.pdf", page))
        savefig(page_plot, page_file)
        push!(page_files, page_file)
    end

    merge_args = String[
        "gs",
        "-dBATCH",
        "-dNOPAUSE",
        "-q",
        "-sDEVICE=pdfwrite",
        "-dCompatibilityLevel=1.4",
        "-sOutputFile=$(output_path)",
    ]
    append!(merge_args, page_files)
    run(Cmd(merge_args))

    if !isfile(output_path)
        error("Ghostscript completed but output PDF was not created: $output_path")
    end

    for f in page_files
        isfile(f) && rm(f)
    end

    dmeta = pkg_metadata("DSGE")
    mcmeta = pkg_metadata("ModelConstructors")

    println("Trace plot generation complete.")
    println("Input file: $input_path")
    println("Raw mhparams dims: $raw_dims")
    println("Orientation: $orientation")
    println("Parameters: $n_params")
    println("Draws: $n_draws")
    println("Pages: $n_pages (layout $(rows)x$(cols), $per_page traces/page)")
    println("Output PDF: $output_path")
    println("DSGE version: $(dmeta.version)")
    println("DSGE git_revision: $(dmeta.git_revision)")
    println("DSGE tree_hash: $(dmeta.tree_hash)")
    println("ModelConstructors version: $(mcmeta.version)")
    println("ModelConstructors git_revision: $(mcmeta.git_revision)")
    println("ModelConstructors tree_hash: $(mcmeta.tree_hash)")
end

function main()
    opts, rows, cols = parse_args(ARGS)
    input_path = resolve_path(opts["input"])
    output_path = resolve_path(opts["output"])
    subspec = opts["subspec"]
    write_trace_pdf(input_path, output_path, subspec, rows, cols)
end

main()
