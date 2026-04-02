using JuMP, HiGHS

if !isdefined(Main, :Instance)
    include(joinpath(@__DIR__, "parser.jl"))
end

# Solve with iterative SEC cuts
function solve_sec(inst::Instance; time_limit::Float64=300.0, verbose::Bool=false)
    n, s, t, vmin, m = inst.n, inst.s, inst.t, inst.vmin, inst.m
    regions = inst.regions

    d = compute_distances(inst)
    arcs, out_neighbors, in_neighbors = build_arc_sets(inst, d)

    if isempty(arcs)
        return (objective=Inf, status="INFEASIBLE", time_sec=0.0, bb_nodes=0, sec_cuts=0, path=Int[])
    end

    model = Model(HiGHS.Optimizer)
    !verbose && set_silent(model)
    set_time_limit_sec(model, time_limit)

    @variable(model, x[arcs], Bin)
    @variable(model, y[1:n], Bin)

    @objective(model, Min, sum(d[i,j] * x[(i,j)] for (i,j) in arcs))

    # C1
    !isempty(out_neighbors[s]) && @constraint(model, sum(x[(s,j)] for j in out_neighbors[s]) == 1)
    !isempty(in_neighbors[s])  && @constraint(model, sum(x[(j,s)] for j in in_neighbors[s])  == 0)
    # C2
    !isempty(in_neighbors[t])  && @constraint(model, sum(x[(j,t)] for j in in_neighbors[t])  == 1)
    !isempty(out_neighbors[t]) && @constraint(model, sum(x[(t,j)] for j in out_neighbors[t]) == 0)
    # C3
    for i in 1:n
        (i == s || i == t) && continue
        in_sum  = isempty(in_neighbors[i])  ? AffExpr(0) : sum(x[(j,i)] for j in in_neighbors[i])
        out_sum = isempty(out_neighbors[i]) ? AffExpr(0) : sum(x[(i,j)] for j in out_neighbors[i])
        @constraint(model, in_sum  == y[i])
        @constraint(model, out_sum == y[i])
    end
    # C4
    @constraint(model, y[s] == 1)
    @constraint(model, y[t] == 1)
    # C5
    for k in 1:m
        rn = [i for i in 1:n if regions[i] == k]
        !isempty(rn) && @constraint(model, sum(y[i] for i in rn) >= 1)
    end
    # C6
    @constraint(model, sum(y[i] for i in 1:n) >= vmin)

    sec_cuts_added  = 0
    best_obj        = Inf
    best_path       = Int[]
    iteration       = 0
    final_status    = "UNKNOWN"
    cached_bb_nodes = -1
    start_time      = time()

    while true
        iteration += 1
        if time() - start_time >= time_limit
            final_status = "TIME_LIMIT"; break
        end

        optimize!(model)
        status_sym   = termination_status(model)
        final_status = string(status_sym)
        verbose && println("  [SEC] iter=$iteration status=$final_status cuts=$sec_cuts_added")

        cached_bb_nodes = try Int(MOI.get(model, MOI.NodeCount())) catch; cached_bb_nodes end

        if primal_status(model) != MOI.FEASIBLE_POINT
            status_sym == MOI.TIME_LIMIT && (final_status = "TIME_LIMIT")
            verbose && println("  [SEC] No feasible point.")
            break
        end

        xsol        = Dict((i,j) => (value(x[(i,j)]) > 0.5 ? 1 : 0) for (i,j) in arcs)
        ysol        = [value(y[i]) > 0.5 ? 1 : 0 for i in 1:n]
        cur_obj_raw = objective_value(model)

        visited_nodes = [i for i in 1:n if ysol[i] == 1]
        adj = Dict(i => Int[] for i in 1:n)
        for (i,j) in arcs
            xsol[(i,j)] == 1 && push!(adj[i], j)
        end

        reached   = _bfs_reachable(s, adj)
        unreached = Set(v for v in visited_nodes if !(v in reached))

        if isempty(unreached)
            cur_obj  = round(Int, cur_obj_raw)
            cur_path = _reconstruct_path_sec(s, t, n, adj)
            if cur_obj < best_obj
                best_obj = cur_obj; best_path = cur_path
            end
            verbose && println("  [SEC] Valid solution cost=$cur_obj")
            break
        end

        new_cuts = _add_sec_cuts!(model, x, arcs, unreached, adj)
        sec_cuts_added += new_cuts
        new_cuts == 0 && (@warn "No new SEC cuts."; break)
        status_sym == MOI.TIME_LIMIT && (final_status = "TIME_LIMIT"; break)
    end

    return (objective=best_obj, status=final_status,
            time_sec=time()-start_time, bb_nodes=cached_bb_nodes,
            sec_cuts=sec_cuts_added, path=best_path)
end

# BFS from src
function _bfs_reachable(src::Int, adj::Dict{Int,Vector{Int}}) :: Set{Int}
    reached = Set{Int}([src])
    queue   = [src]
    while !isempty(queue)
        curr = popfirst!(queue)
        for nxt in get(adj, curr, Int[])
            !(nxt in reached) && (push!(reached, nxt); push!(queue, nxt))
        end
    end
    return reached
end

# Weakly-connected components of nodes
function _find_components(nodes::Set{Int}, adj::Dict{Int,Vector{Int}}) :: Vector{Vector{Int}}
    undirected = Dict(n => Int[] for n in nodes)
    for (i, nbs) in adj
        i in nodes || continue
        for j in nbs
            j in nodes || continue
            push!(undirected[i], j); push!(undirected[j], i)
        end
    end
    components = Vector{Vector{Int}}()
    unvisited  = copy(nodes)
    while !isempty(unvisited)
        seed = first(unvisited)
        comp = Set{Int}([seed])
        queue = [seed]
        while !isempty(queue)
            curr = popfirst!(queue)
            for nxt in get(undirected, curr, Int[])
                !(nxt in comp) && (push!(comp, nxt); push!(queue, nxt))
            end
        end
        push!(components, collect(comp))
        setdiff!(unvisited, comp)
    end
    return components
end

# Add one SEC cut per subtour component, return cut count
function _add_sec_cuts!(model, x, arcs::Vector{Tuple{Int,Int}},
                        unreached::Set{Int}, adj::Dict{Int,Vector{Int}}) :: Int
    components = _find_components(unreached, adj)
    cuts = 0
    for comp in components
        cs = Set(comp)
        arcs_in = [(i,j) for (i,j) in arcs if i in cs && j in cs]
        isempty(arcs_in) && continue
        @constraint(model, sum(x[(i,j)] for (i,j) in arcs_in) <= length(comp) - 1)
        cuts += 1
    end
    return cuts
end

# Follow used arcs from s to t
function _reconstruct_path_sec(s::Int, t::Int, n::Int, adj::Dict{Int,Vector{Int}}) :: Vector{Int}
    path = [s]; current = s; visited_set = Set([s])
    for _ in 1:n
        nbs = get(adj, current, Int[])
        isempty(nbs) && break
        nxt = first(nbs)
        nxt in visited_set && nxt != t && (@warn "Cycle at $nxt"; break)
        push!(path, nxt); push!(visited_set, nxt)
        current = nxt
        current == t && break
    end
    last(path) != t && @warn "Reconstruction did not reach t=$t"
    return path
end
