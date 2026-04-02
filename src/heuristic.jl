if !isdefined(Main, :Instance)
    include(joinpath(@__DIR__, "parser.jl"))
end

# Entry point: greedy -> 2-opt -> node removal
function solve_heuristic(inst::Instance; verbose::Bool=false)
    start_time = time()
    d = compute_distances(inst)

    path, cost = _greedy_construct(inst, d, verbose)
    if isempty(path)
        return (objective=Inf, time_sec=time()-start_time, path=Int[])
    end

    path, cost = _two_opt!(path, inst, d, verbose)
    path, cost = _node_removal!(path, inst, d, verbose)

    return (objective=cost, time_sec=time()-start_time, path=path)
end

# Greedy nearest-neighbour, prioritising uncovered regions
function _greedy_construct(inst::Instance, d::Matrix{Int}, verbose::Bool)
    n, s, t = inst.n, inst.s, inst.t
    vmin, m  = inst.vmin, inst.m
    regions, R = inst.regions, inst.R

    path    = [s]
    visited = Set([s])
    covered_regions = Set{Int}()
    regions[s] > 0 && push!(covered_regions, regions[s])
    current = s

    while length(covered_regions) < m || length(visited) < vmin
        candidates = [(j, d[current, j]) for j in 1:n
                      if !(j in visited) && j != t && 0 < d[current, j] <= R]

        if isempty(candidates)
            target_nodes = [j for j in 1:n if !(j in visited) && j != t]
            found = false
            if !isempty(target_nodes)
                best_cost = typemax(Int)
                best_sub  = nothing
                for tgt in target_nodes
                    sub = _dijkstra_path(d, R, n, current, tgt)
                    if sub !== nothing
                        c = _path_cost(sub, d)
                        if c < best_cost
                            best_cost = c
                            best_sub  = sub
                        end
                    end
                end
                if best_sub !== nothing && length(best_sub) > 1
                    for node in best_sub[2:end]
                        if !(node in visited)
                            push!(path, node)
                            push!(visited, node)
                            regions[node] > 0 && push!(covered_regions, regions[node])
                        end
                    end
                    current = last(path)
                    found = true
                end
            end
            if !found
                if d[current, t] <= R &&
                   length(covered_regions) >= m && length(visited) >= vmin
                    break
                end
                verbose && println("  [Greedy] Stuck at node $current.")
                return Int[], typemax(Int)
            end
            continue
        end

        uncovered = [(j, dist) for (j, dist) in candidates
                     if regions[j] > 0 && !(regions[j] in covered_regions)]
        chosen = if !isempty(uncovered)
            sort!(uncovered, by=p->p[2]); uncovered[1][1]
        else
            sort!(candidates, by=p->p[2]); candidates[1][1]
        end

        push!(path, chosen)
        push!(visited, chosen)
        regions[chosen] > 0 && push!(covered_regions, regions[chosen])
        current = chosen
    end

    # Top-up to meet vmin
    while length(visited) < vmin
        topup = [(j, d[current, j]) for j in 1:n
                 if !(j in visited) && j != t && 0 < d[current, j] <= R]
        if isempty(topup)
            target_nodes = [j for j in 1:n if !(j in visited) && j != t]
            found = false
            for tgt in target_nodes
                sub = _dijkstra_path(d, R, n, current, tgt)
                if sub !== nothing
                    for node in sub[2:end]
                        if !(node in visited)
                            push!(path, node)
                            push!(visited, node)
                            regions[node] > 0 && push!(covered_regions, regions[node])
                        end
                    end
                    current = last(path)
                    found = true
                    break
                end
            end
            found || break
        else
            sort!(topup, by=p->p[2])
            chosen = topup[1][1]
            push!(path, chosen)
            push!(visited, chosen)
            regions[chosen] > 0 && push!(covered_regions, regions[chosen])
            current = chosen
        end
    end

    length(covered_regions) < m && return Int[], typemax(Int)
    length(visited) < vmin       && return Int[], typemax(Int)

    # Connect to t
    if current != t
        if d[current, t] <= R
            push!(path, t)
        else
            sub = _dijkstra_path(d, R, n, current, t)
            sub === nothing && return Int[], typemax(Int)
            for node in sub[2:end]
                push!(path, node)
            end
        end
    end

    last(path) != t && return Int[], typemax(Int)
    return path, _path_cost(path, d)
end

# 2-opt local search
function _two_opt!(path::Vector{Int}, inst::Instance, d::Matrix{Int}, verbose::Bool)
    R = inst.R
    improved = true
    iters = 0
    while improved
        improved = false
        iters += 1
        L = length(path)
        for i in 2:L-2
            for j in i+1:L-1
                a, b, c, e = path[i-1], path[i], path[j], path[j+1]
                old_cost = d[a,b] + d[c,e]
                new_cost = d[a,c] + d[b,e]
                new_cost >= old_cost && continue
                new_path = vcat(path[1:i-1], reverse(path[i:j]), path[j+1:end])
                feasible = all(d[new_path[k], new_path[k+1]] <= R for k in i-1:j)
                if feasible
                    path = new_path
                    improved = true
                    @goto next_outer
                end
            end
            @label next_outer
        end
    end
    verbose && println("  [2-opt] $iters passes")
    return path, _path_cost(path, d)
end

# Remove non-essential interior nodes (largest saving first)
function _node_removal!(path::Vector{Int}, inst::Instance, d::Matrix{Int}, verbose::Bool)
    vmin, m  = inst.vmin, inst.m
    regions, R = inst.regions, inst.R
    removed = true
    while removed
        removed = false
        L = length(path)
        savings, indices = Float64[], Int[]
        for idx in 2:L-1
            node, prev_node, next_node = path[idx], path[idx-1], path[idx+1]
            d[prev_node, next_node] > R && continue
            length(_covered_regions_except(path, regions, idx)) < m && continue
            length(Set(path)) - (count(==(node), path) == 1 ? 1 : 0) < vmin && continue
            saving = d[prev_node, node] + d[node, next_node] - d[prev_node, next_node]
            saving > 0 && (push!(savings, saving); push!(indices, idx))
        end
        isempty(indices) && break
        best_idx = indices[argmax(savings)]
        verbose && println("  [NodeRemoval] removing $(path[best_idx])")
        deleteat!(path, best_idx)
        removed = true
    end
    return path, _path_cost(path, d)
end

# Total path distance (km)
function _path_cost(path::Vector{Int}, d::Matrix{Int}) :: Int
    sum(d[path[i], path[i+1]] for i in 1:length(path)-1)
end

# Regions covered excluding position skip_idx
function _covered_regions_except(path::Vector{Int}, regions::Vector{Int}, skip_idx::Int) :: Set{Int}
    covered = Set{Int}()
    for (k, node) in enumerate(path)
        k == skip_idx && continue
        r = regions[node]
        r > 0 && push!(covered, r)
    end
    return covered
end

# Dijkstra shortest path from src to dst within range R
function _dijkstra_path(d::Matrix{Int}, R::Int, n::Int, src::Int, dst::Int) :: Union{Vector{Int}, Nothing}
    INF = typemax(Int) / 2
    dist_to = fill(INF, n)
    dist_to[src] = 0
    prev = zeros(Int, n)
    pq = [(0, src)]
    while !isempty(pq)
        sort!(pq, by=p->p[1])
        (cd, u) = popfirst!(pq)
        if u == dst
            path_nodes = [dst]; curr = dst
            while curr != src
                curr = prev[curr]
                pushfirst!(path_nodes, curr)
                curr == 0 && return nothing
            end
            return path_nodes
        end
        cd > dist_to[u] && continue
        for v in 1:n
            v == u && continue
            0 < d[u,v] <= R || continue
            nd = cd + d[u,v]
            if nd < dist_to[v]
                dist_to[v] = nd; prev[v] = u
                push!(pq, (nd, v))
            end
        end
    end
    return nothing
end

# Verify solution feasibility
function verify_solution(path::Vector{Int}, inst::Instance, d::Matrix{Int}) :: Bool
    n, s, t, vmin, m = inst.n, inst.s, inst.t, inst.vmin, inst.m
    regions, R = inst.regions, inst.R
    ok = true
    isempty(path) && (println("  INVALID: empty path"); return false)
    path[1] != s    && (println("  INVALID: starts at $(path[1]), not s=$s"); ok = false)
    last(path) != t && (println("  INVALID: ends at $(last(path)), not t=$t"); ok = false)
    for i in 1:length(path)-1
        dist = d[path[i], path[i+1]]
        dist > R && (println("  INVALID: leg $(path[i])->$(path[i+1]) dist=$dist > R=$R"); ok = false)
    end
    unique_visits = length(Set(path))
    unique_visits < vmin && (println("  INVALID: $unique_visits unique aerodromes < V_min=$vmin"); ok = false)
    covered = Set(regions[v] for v in path if regions[v] > 0)
    for k in 1:m
        !(k in covered) && (println("  INVALID: region $k not covered"); ok = false)
    end
    ok && println("  VALID solution ($(length(path)) stops, $unique_visits unique, cost=$(_path_cost(path,d)))")
    return ok
end
