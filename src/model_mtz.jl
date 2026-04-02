using JuMP, HiGHS

if !isdefined(Main, :Instance)
    include(joinpath(@__DIR__, "parser.jl"))
end

# Solve with MTZ subtour elimination
function solve_mtz(inst::Instance; time_limit::Float64=300.0, verbose::Bool=false)
    n, s, t, vmin, m = inst.n, inst.s, inst.t, inst.vmin, inst.m
    regions = inst.regions

    d = compute_distances(inst)
    arcs, out_neighbors, in_neighbors = build_arc_sets(inst, d)

    if isempty(arcs)
        return (objective=Inf, status="INFEASIBLE", time_sec=0.0, bb_nodes=0, path=Int[])
    end

    model = Model(HiGHS.Optimizer)
    !verbose && set_silent(model)
    set_time_limit_sec(model, time_limit)

    @variable(model, x[arcs], Bin)
    @variable(model, y[1:n], Bin)
    @variable(model, 0 <= u[1:n] <= n)  # position along path

    @objective(model, Min, sum(d[i,j] * x[(i,j)] for (i,j) in arcs))

    # C1: flow at s
    !isempty(out_neighbors[s]) && @constraint(model, sum(x[(s,j)] for j in out_neighbors[s]) == 1)
    !isempty(in_neighbors[s])  && @constraint(model, sum(x[(j,s)] for j in in_neighbors[s])  == 0)

    # C2: flow at t
    !isempty(in_neighbors[t])  && @constraint(model, sum(x[(j,t)] for j in in_neighbors[t])  == 1)
    !isempty(out_neighbors[t]) && @constraint(model, sum(x[(t,j)] for j in out_neighbors[t]) == 0)

    # C3: flow balance
    for i in 1:n
        (i == s || i == t) && continue
        in_sum  = isempty(in_neighbors[i])  ? AffExpr(0) : sum(x[(j,i)] for j in in_neighbors[i])
        out_sum = isempty(out_neighbors[i]) ? AffExpr(0) : sum(x[(i,j)] for j in out_neighbors[i])
        @constraint(model, in_sum  == y[i])
        @constraint(model, out_sum == y[i])
    end

    # C4: force s and t visited
    @constraint(model, y[s] == 1)
    @constraint(model, y[t] == 1)

    # C5: region coverage
    for k in 1:m
        rn = [i for i in 1:n if regions[i] == k]
        !isempty(rn) && @constraint(model, sum(y[i] for i in rn) >= 1)
    end

    # C6: minimum visits
    @constraint(model, sum(y[i] for i in 1:n) >= vmin)

    # C7: MTZ subtour elimination
    for (i,j) in arcs
        (i == s || j == s) && continue
        @constraint(model, u[i] - u[j] + n * x[(i,j)] <= n - 1)
    end

    # C8: link u to y
    @constraint(model, u[s] == 0)
    for i in 1:n
        i != s && @constraint(model, u[i] <= n * y[i])
    end

    elapsed = @elapsed optimize!(model)
    status_str = string(termination_status(model))

    bb_nodes = try Int(MOI.get(model, MOI.NodeCount())) catch; -1 end

    if primal_status(model) == MOI.FEASIBLE_POINT
        obj_val = round(Int, objective_value(model))
        xsol = Dict((i,j) => (value(x[(i,j)]) > 0.5 ? 1 : 0) for (i,j) in arcs)
        path = _reconstruct_path(s, t, n, arcs, xsol)
        return (objective=obj_val, status=status_str, time_sec=elapsed, bb_nodes=bb_nodes, path=path)
    else
        return (objective=Inf, status=status_str, time_sec=elapsed, bb_nodes=bb_nodes, path=Int[])
    end
end

# Follow arcs with xsol==1 from s to t
function _reconstruct_path(s::Int, t::Int, n::Int,
                            arcs::Vector{Tuple{Int,Int}},
                            xsol::Dict{Tuple{Int,Int},Int}) :: Vector{Int}
    next_node = Dict(i => j for (i,j) in arcs if xsol[(i,j)] == 1)
    path = [s]
    current = s
    visited_set = Set([s])
    for _ in 1:n
        !haskey(next_node, current) && break
        nxt = next_node[current]
        if nxt in visited_set && nxt != t
            @warn "Cycle at node $nxt during reconstruction"
            break
        end
        push!(path, nxt); push!(visited_set, nxt)
        current = nxt
        current == t && break
    end
    last(path) != t && @warn "Reconstruction did not reach t=$t"
    return path
end
