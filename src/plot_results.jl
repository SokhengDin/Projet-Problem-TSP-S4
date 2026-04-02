# plot_results.jl — Save a high-res PNG for each instance showing all three method paths

ENV["GKSwstype"] = "nul"  # suppress GUI window
using Plots
gr()

const REGION_COLORS = [
    :grey, :royalblue, :green3, :red, :darkorange,
    :purple, :deeppink, :saddlebrown, :teal, :gold,
]

function _node_color(r)
    r == 0 && return :grey
    REGION_COLORS[mod1(r, length(REGION_COLORS) - 1) + 1]
end

function _subplot(inst::Instance, d::Matrix{Int}, path::Vector{Int}, title_str::String)
    xs = Float64.(inst.coords[:, 1])
    ys = Float64.(inst.coords[:, 2])

    p = plot(title=title_str, legend=false, aspect_ratio=:equal,
             xlabel="x", ylabel="y", titlefontsize=9, grid=true)

    # All aerodromes coloured by region
    for i in 1:inst.n
        scatter!(p, [xs[i]], [ys[i]], color=_node_color(inst.regions[i]),
                 markersize=6, markerstrokewidth=0.5, markerstrokecolor=:white)
        annotate!(p, xs[i], ys[i] + 2.5, text(string(i), 5, :black))
    end

    # Path arcs — skip zero-length segments to avoid NaN in arrow rendering
    if length(path) >= 2
        for k in 1:length(path)-1
            a, b = path[k], path[k+1]
            # Skip if same node (duplicate stop)
            a == b && continue
            dx = xs[b] - xs[a]
            dy = ys[b] - ys[a]
            # Draw line; add arrow only when segment is long enough
            if sqrt(dx^2 + dy^2) > 0.5
                plot!(p, [xs[a], xs[b]], [ys[a], ys[b]],
                      color=:black, linewidth=1.4,
                      arrow=arrow(:closed, :head, 0.3, 0.2))
            else
                plot!(p, [xs[a], xs[b]], [ys[a], ys[b]],
                      color=:black, linewidth=1.4)
            end
        end
    end

    # Mark s (green star) and t (red diamond) on top
    scatter!(p, [xs[inst.s]], [ys[inst.s]], color=:lime,
             markersize=10, markershape=:star5, markerstrokewidth=1, markerstrokecolor=:darkgreen)
    scatter!(p, [xs[inst.t]], [ys[inst.t]], color=:red,
             markersize=10, markershape=:diamond, markerstrokewidth=1, markerstrokecolor=:darkred)

    return p
end

function save_plot(inst::Instance, d::Matrix{Int},
                   path_mtz::Vector{Int}, obj_mtz,
                   path_sec::Vector{Int}, obj_sec,
                   path_heur::Vector{Int}, obj_heur,
                   outdir::String, name::String)

    obj_s(v) = isinf(v) ? "INF" : string(Int(v))

    p1 = _subplot(inst, d, path_mtz,  "MTZ  obj=$(obj_s(obj_mtz))")
    p2 = _subplot(inst, d, path_sec,  "SEC  obj=$(obj_s(obj_sec))")
    p3 = _subplot(inst, d, path_heur, "Heuristic  obj=$(obj_s(obj_heur))")

    fig = plot(p1, p2, p3, layout=(1, 3), size=(2400, 800), dpi=200,
               plot_title=name, plot_titlefontsize=12)

    out = joinpath(outdir, "$(name)_results.png")
    savefig(fig, out)
    println("  Plot saved: $out")
end
