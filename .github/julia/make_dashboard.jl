using Printf

# --- ANSI colored dots (render fine in GitHub markdown) ---
const RED     = "🔴"
const GREEN   = "🟢"
const YELLOW  = "🟡"
const ORANGE  = "🟠"
const CYAN    = "⚪"
const MAGENTA = "🟣"

# --- Modes ---
modes = ["forward", "vector-forward", "reverse", "vector-reverse"]

RESULTS_DIR = joinpath(@__DIR__, "..", "..", "results")

# ---------------------------
# Load logs
# ---------------------------

results = Dict{String,Dict{String,String}}()

logfiles = filter(f -> endswith(f, ".log"), readdir(RESULTS_DIR; join=true))

for file in logfiles
    for line in eachline(file)

        isempty(strip(line)) && continue

        parts = split(line, ":")

        if length(parts) != 3
            continue
        end

        routine, mode, status = parts

        if !haskey(results, routine)
            results[routine] = Dict(m => "SKIPPED" for m in modes)
        end

        results[routine][mode] = status
    end
end

routines = sort(collect(keys(results)))

# ---------------------------
# Status → dot
# ---------------------------

status_to_dot = Dict(
    "MACHINE_PRECISION" => GREEN,
    "ACCEPTABLE" => GREEN,
    "OUTSIDE_TOLERANCE" => ORANGE,
    "EXECUTION_FAILED" => RED,
    "SKIPPED" => CYAN,
    "TAPENADE_FAILED" => MAGENTA
)

# ---------------------------
# Write dashboard.md
# ---------------------------

open("dashboard.md", "w") do io
    # -------- Legend --------
    println(io, "## Legend")
    println(io, "")
    println(io, "- $(GREEN) MACHINE_PRECISION")
    println(io, "- $(GREEN) ACCEPTABLE")
    println(io, "- $(ORANGE) OUTSIDE_TOLERANCE")
    println(io, "- $(RED) EXECUTION_FAILED")
    println(io, "- $(CYAN) SKIPPED")
    println(io, "- $(MAGENTA) TAPENADE_FAILED")
    println(io, "")
    println(io, "---")
    println(io, "")

    # -------- Summary per mode --------
    println(io, "## Summary by Mode")
    for m in modes
        counts = Dict("MACHINE_PRECISION"=>0, "ACCEPTABLE"=>0, "OUTSIDE_TOLERANCE"=>0,
                      "EXECUTION_FAILED"=>0, "SKIPPED"=>0, "TAPENADE_FAILED"=>0)
        for r in routines
            status = get(results[r], m, "SKIPPED")
            counts[status] = get(counts, status, 0) + 1
        end
        println(io, "- **$m**: " *
            "$(GREEN) $(counts["MACHINE_PRECISION"]) " *
            "$(GREEN) $(counts["ACCEPTABLE"]) " *
            "$(ORANGE) $(counts["OUTSIDE_TOLERANCE"]) " *
            "$(RED) $(counts["EXECUTION_FAILED"]) " *
            "$(CYAN) $(counts["SKIPPED"]) " *
            "$(MAGENTA) $(counts["TAPENADE_FAILED"])")
    end
    println(io, "")
    println(io, "---")
    println(io, "")

    # -------- Table header --------
    println(io, "| Routine | " * join(modes, " | ") * " |")
    println(io, "|---------|" * repeat("------|", length(modes)))

    # -------- Table rows --------
    for r in routines
        line = "| $r "

        for m in modes
            status = get(results[r], m, "SKIPPED")
            dot = get(status_to_dot, status, CYAN)
            line *= "| $dot "
        end

        line *= "|"
        println(io, line)
    end
end

println("Dashboard written to dashboard.md ✅")
