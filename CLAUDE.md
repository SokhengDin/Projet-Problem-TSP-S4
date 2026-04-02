# CORO Project — Aerodrome Race (Coupe Breitling 100/24)

## Reference Document (`claude.md`)

**ENSIIE 2025–2026 | Prof. Soutil**
**Deadline: 31/03/2026**

---

## 1. Problem Summary

### Input
- `n` aerodromes numbered `1..n`, each with coordinates `(x_i, y_i)` and a region `r_i ∈ {0, 1, ..., m}` (0 = no region)
- Starting aerodrome `s`, arriving aerodrome `t`
- Minimum number of aerodromes to visit: `V_min`
- Number of regions: `m` (each must be visited at least once)
- Maximum flight range (radius): `R`

### Objective
Find a **path** from `s` to `t` visiting at least `V_min` aerodromes, covering all `m` regions, with consecutive aerodromes within distance `R`, **minimizing total Euclidean distance** (rounded to nearest integer).

### Distance
```
d(i, j) = round(sqrt((x_i - x_j)^2 + (y_i - y_j)^2))
```

### Key Differences from Standard TSP
- **Path** (s -> t), not a cycle
- **Subset** selection: we choose which aerodromes to visit (at least V_min)
- **Region covering**: each region 1..m must have at least one visited aerodrome
- **Range constraint**: consecutive stops must be within distance R

---

## 2. Graph Construction

Build a **directed graph** G = (V, A):
- V = {1, ..., n}
- A = {(i, j) : d(i, j) ≤ R, i ≠ j}

Only create arcs between aerodromes reachable within the fuel range R.

---

## 3. ILP Model 1 — MTZ Formulation (Polynomial Constraints)

### Decision Variables
- `x_{ij} ∈ {0, 1}` for (i,j) ∈ A: arc (i,j) is used
- `y_i ∈ {0, 1}` for i ∈ V: aerodrome i is visited
- `u_i ∈ ℝ` for i ∈ V \ {s}: MTZ position variable (order of visit)

### Objective
```
min  Σ_{(i,j) ∈ A}  d(i,j) · x_{ij}
```

### Constraints

#### (C1) Flow conservation — starting aerodrome s
```
Σ_j x_{sj} = 1          (exactly one arc out)
Σ_j x_{js} = 0          (no arc into s — it's the start)
```

#### (C2) Flow conservation — arriving aerodrome t
```
Σ_j x_{jt} = 1          (exactly one arc in)
Σ_j x_{tj} = 0          (no arc out of t — it's the end)
```

#### (C3) Flow conservation — intermediate visited nodes
For all i ∈ V \ {s, t}:
```
Σ_j x_{ji} = y_i        (enter i iff visited)
Σ_j x_{ij} = y_i        (leave i iff visited)
```

#### (C4) Forcing s and t as visited
```
y_s = 1
y_t = 1
```

#### (C5) Region covering
For each region k = 1, ..., m:
```
Σ_{i : r_i = k}  y_i  ≥  1
```

#### (C6) Minimum visits
```
Σ_i y_i  ≥  V_min
```

#### (C7) MTZ subtour elimination
For all (i,j) ∈ A with i ≠ s, j ≠ s:
```
u_i - u_j + n · x_{ij}  ≤  n - 1
```
This ensures: if x_{ij} = 1, then u_j ≥ u_i + 1 (j comes after i in the path).

Bounds: u_i ∈ [1, n] for i ∈ V \ {s}. We can set u_s = 0 (or not define it).

#### (C8) Linking u to y (optional tightening)
```
u_i ≤ n · y_i          for all i ≠ s
```
If node i is not visited (y_i = 0), then u_i = 0.

### MTZ Recap (from CORO Lecture — TSP)

The MTZ formulation introduces **position variables** u_i that track the order in which nodes are visited. The key constraint:

```
x_{ij} = 1  ⟹  u_i - u_j + n ≤ n - 1  ⟹  u_j ≥ u_i + 1
```

This means: if we travel from i to j, then j's position is strictly after i's position. This makes subtours impossible because a subtour would require u values to simultaneously increase and loop back.

**Weakness**: The LP relaxation of MTZ is weaker than SEC (subtour elimination constraints). MTZ gives fractional u values that allow weak bounds. However, it has only O(n²) constraints.

---

## 4. ILP Model 2 — Subtour Elimination Constraints (SEC / Exponential)

### Decision Variables
Same as Model 1 but **without** u_i. Only x_{ij} and y_i.

### Constraints
(C1)–(C6) are identical to Model 1.

Instead of MTZ, we use **subtour elimination constraints**:

#### (SEC) For every subset S ⊆ V \ {s, t} with S ≠ ∅:
```
Σ_{i ∈ S, j ∈ S} x_{ij}  ≤  Σ_{i ∈ S} y_i  -  1
```

Meaning: within any subset S that doesn't contain s or t, the number of arcs used inside S is at most |visited nodes in S| - 1. This prevents disconnected cycles.

**Equivalently** (for a path problem), for any S ⊆ V \ {s, t}:
```
Σ_{i ∈ S, j ∉ S} x_{ij}  ≥  y_k    for any k ∈ S
```
(If k is visited, there must be an arc leaving S.)

### Implementation: Lazy Constraint Callback

Since there are exponentially many subsets S, we add SEC **on the fly**:

```
Algorithm: SEC with Lazy Constraints
─────────────────────────────────────
1. Solve the ILP without any SEC (just C1–C6 + integrality)
2. At each integer-feasible incumbent:
   a. Extract arcs with x_{ij} = 1
   b. Find connected components of the solution
   c. For each connected component S that doesn't contain s:
      - Add: Σ_{i∈S, j∈S} x_{ij} ≤ |S| - 1
   d. If no violated SEC found -> solution is valid
3. Continue solving
```

### Why SEC is Stronger than MTZ

From CORO lectures (LP relaxation quality):
- **SEC LP relaxation** = subtour polytope (tight description for TSP)
- **MTZ LP relaxation** ⊂ subtour polytope (strictly weaker)

This means: z_LP(SEC) ≥ z_LP(MTZ), so SEC gives better bounds and typically fewer B&B nodes. The tradeoff is implementation complexity (need callbacks).

---

## 5. Separation Problem (from CORO Lecture on Cutting Planes)

The **separation problem** for SEC asks:

> Given the current fractional solution x̂, is there a subset S such that the SEC is violated?

For the integer case (lazy constraints), it reduces to **finding connected components** in the solution graph. For fractional solutions (user cuts), it becomes a **minimum cut problem**:

For each pair of nodes (i, t), compute min s-t cut in the support graph with capacities x̂_{ij}. If the min cut value < y_i (or 1 for visited nodes), the corresponding SEC is violated.

This connects to the **max-flow min-cut theorem** and the general theory of separation oracles from CORO.

---

## 6. Theoretical Foundations from CORO

### 6.1 LP Duality (Cours 1)

For the LP relaxation of both models, the dual provides a **lower bound**:

- **Weak duality**: z_LP ≤ z* (for minimization)
- **Strong duality**: z_LP = z_D at optimality
- **Complementary slackness**: x_i > 0 ⟹ dual constraint i is binding

When we solve the LP relaxation of our ILP, the optimal LP value gives a lower bound on the integer optimum. The gap between z_LP and z_ILP is the **integrality gap**.

### 6.2 Branch-and-Bound (Cours 2)

Both models are solved using **Branch-and-Bound**:

```
Algorithm: Branch-and-Bound (best-first)
────────────────────────────────────────
1. Solve LP relaxation at root -> z_LP (lower bound)
2. If integer -> done
3. Maintain priority queue of nodes, incumbent z*
4. Select node with best bound
5. Branch on fractional variable (x_{ij} = 0 or 1)
6. Solve LP at each child node
7. Prune if:
   a. LP infeasible
   b. LP bound ≥ incumbent z*
   c. LP solution is integer (update z* if better)
8. Repeat until queue empty
```

**Branching rule**: Branch on the fractional variable with the largest fractional part (from Prof. Soutil's lecture).

**Three pruning criteria** (from CORO Cours 2):
1. **Optimality**: LP relaxation gives integer solution
2. **Bound**: LP bound ≥ best known incumbent
3. **Infeasibility**: LP relaxation is infeasible

### 6.3 Gomory Cutting Planes (Cours 3)

Not directly used in this project, but the **SEC with lazy constraints** follows the same cutting plane philosophy:

```
Cutting Plane Algorithm:
1. Solve LP relaxation -> x̂
2. If x̂ is integer -> optimal
3. Find violated inequality (separation problem)
4. Add cut, re-solve
5. Repeat
```

The Gomory cut for a fractional basic variable x_i with row Σ_j a̅_{ij} x_j = b̅_i is:

```
Σ_{j ∈ N}  frac(a̅_{ij}) · x_j  ≥  frac(b̅_i)
```

where frac(z) = z - ⌊z⌋.

### 6.4 Lagrangian Relaxation (Cours 4)

An alternative approach to get better bounds: relax the complicating constraints (e.g., region covering) into the objective with Lagrange multipliers λ:

```
L(λ) = min  Σ d_{ij} x_{ij} + Σ_{k=1}^{m} λ_k (1 - Σ_{i:r_i=k} y_i)
       s.t. flow constraints, range constraints, x,y binary
```

The **Lagrangian dual** maximizes L(λ) over λ ≥ 0 using the **subgradient method**:

```
λ_k^{t+1} = max(0, λ_k^t + θ_t · (1 - Σ_{i:r_i=k} y_i^t))
```

with Polyak step size:
```
θ_t = α_t · (z̄ - L(λ^t)) / ‖g^t‖²
```

**Key relationship** (from CORO):
```
z_LP ≤ z_LD ≤ z_ILP
```
The Lagrangian bound z_LD is at least as tight as the LP relaxation.

### 6.5 Column Generation / Dantzig-Wolfe (Cours 5)

Not directly applicable here (no block-diagonal structure), but the pricing subproblem idea connects to the separation problem in SEC:

- **Pricing subproblem** in CG = find column with negative reduced cost
- **Separation problem** in cutting planes = find violated inequality

Both are "find the most useful object to add" to the current formulation.

### 6.6 TSP Lower Bound: 1-Tree Relaxation

From the CORO TSP lecture, the **1-tree** relaxation gives a lower bound:

1. Remove node 1 from graph
2. Compute MST on remaining nodes
3. Add 2 cheapest edges incident to node 1
4. This gives a **1-tree** (spanning tree + one extra edge at node 1)

For our path problem, we can adapt: compute MST on the set of "must-visit" nodes (one per region + s + t), then use the MST cost as a lower bound.

---

## 7. Instance Parser (Julia)

```julia
function parse_instance(filename::String)
    lines = readlines(filename)
    idx = 1
    
    n = parse(Int, strip(lines[idx])); idx += 1
    s = parse(Int, strip(lines[idx])); idx += 1
    t = parse(Int, strip(lines[idx])); idx += 1
    vmin = parse(Int, strip(lines[idx])); idx += 1
    m = parse(Int, strip(lines[idx])); idx += 1
    
    # Region of each aerodrome
    regions = parse.(Int, split(strip(lines[idx])))
    idx += 1
    
    # Radius
    R = parse(Int, strip(lines[idx])); idx += 1
    
    # Coordinates
    coords = zeros(Int, n, 2)
    for i in 1:n
        vals = parse.(Int, split(strip(lines[idx])))
        coords[i, 1] = vals[1]
        coords[i, 2] = vals[2]
        idx += 1
    end
    
    return (n=n, s=s, t=t, vmin=vmin, m=m, regions=regions, R=R, coords=coords)
end

function compute_distances(coords, n)
    d = zeros(Int, n, n)
    for i in 1:n, j in 1:n
        if i != j
            dx = coords[i,1] - coords[j,1]
            dy = coords[i,2] - coords[j,2]
            d[i,j] = round(Int, sqrt(dx^2 + dy^2))
        end
    end
    return d
end
```

---

## 8. Model 1 — MTZ (Julia/JuMP)

```julia
using JuMP, HiGHS  # or Cbc

function solve_mtz(inst)
    n, s, t, vmin, m = inst.n, inst.s, inst.t, inst.vmin, inst.m
    regions, R, coords = inst.regions, inst.R, inst.coords
    
    d = compute_distances(coords, n)
    
    # Build arc set: only arcs within range R
    arcs = [(i,j) for i in 1:n, j in 1:n if i != j && d[i,j] <= R]
    
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    
    # Variables
    @variable(model, x[arcs], Bin)          # arc usage
    @variable(model, y[1:n], Bin)           # node visited
    @variable(model, 0 <= u[1:n] <= n)      # MTZ position
    
    # Objective
    @objective(model, Min, sum(d[i,j] * x[(i,j)] for (i,j) in arcs))
    
    # Arc set as a dict for fast lookup
    out_arcs = Dict(i => [(i,j) for (i,j) in arcs if true] for i in 1:n)
    # Better: precompute
    out_neighbors = [Int[] for _ in 1:n]
    in_neighbors  = [Int[] for _ in 1:n]
    for (i,j) in arcs
        push!(out_neighbors[i], j)
        push!(in_neighbors[j], i)
    end
    
    # (C1) Flow from s
    @constraint(model, sum(x[(s,j)] for j in out_neighbors[s]) == 1)
    @constraint(model, sum(x[(j,s)] for j in in_neighbors[s]) == 0)
    
    # (C2) Flow into t
    @constraint(model, sum(x[(j,t)] for j in in_neighbors[t]) == 1)
    @constraint(model, sum(x[(t,j)] for j in out_neighbors[t]) == 0)
    
    # (C3) Flow conservation for intermediate nodes
    for i in 1:n
        if i == s || i == t
            continue
        end
        @constraint(model, sum(x[(j,i)] for j in in_neighbors[i]) == y[i])
        @constraint(model, sum(x[(i,j)] for j in out_neighbors[i]) == y[i])
    end
    
    # (C4) s and t must be visited
    @constraint(model, y[s] == 1)
    @constraint(model, y[t] == 1)
    
    # (C5) Region covering
    for k in 1:m
        region_nodes = [i for i in 1:n if regions[i] == k]
        if !isempty(region_nodes)
            @constraint(model, sum(y[i] for i in region_nodes) >= 1)
        end
    end
    
    # (C6) Minimum visits
    @constraint(model, sum(y[i] for i in 1:n) >= vmin)
    
    # (C7) MTZ subtour elimination
    for (i,j) in arcs
        if i != s && j != s
            @constraint(model, u[i] - u[j] + n * x[(i,j)] <= n - 1)
        end
    end
    
    # (C8) Link u to y
    for i in 1:n
        if i != s
            @constraint(model, u[i] <= n * y[i])
        end
    end
    @constraint(model, u[s] == 0)
    
    optimize!(model)
    
    return model, x, y, u
end
```

---

## 9. Model 2 — SEC with Lazy Constraints (Julia/JuMP)

```julia
using JuMP, HiGHS
import Graphs

function solve_sec(inst)
    n, s, t, vmin, m = inst.n, inst.s, inst.t, inst.vmin, inst.m
    regions, R, coords = inst.regions, inst.R, inst.coords
    
    d = compute_distances(coords, n)
    arcs = [(i,j) for i in 1:n, j in 1:n if i != j && d[i,j] <= R]
    
    out_neighbors = [Int[] for _ in 1:n]
    in_neighbors  = [Int[] for _ in 1:n]
    for (i,j) in arcs
        push!(out_neighbors[i], j)
        push!(in_neighbors[j], i)
    end
    
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    
    @variable(model, x[arcs], Bin)
    @variable(model, y[1:n], Bin)
    
    @objective(model, Min, sum(d[i,j] * x[(i,j)] for (i,j) in arcs))
    
    # Same constraints C1–C6 as MTZ model (no u variables, no C7/C8)
    @constraint(model, sum(x[(s,j)] for j in out_neighbors[s]) == 1)
    @constraint(model, sum(x[(j,s)] for j in in_neighbors[s]) == 0)
    @constraint(model, sum(x[(j,t)] for j in in_neighbors[t]) == 1)
    @constraint(model, sum(x[(t,j)] for j in out_neighbors[t]) == 0)
    
    for i in 1:n
        if i == s || i == t; continue; end
        @constraint(model, sum(x[(j,i)] for j in in_neighbors[i]) == y[i])
        @constraint(model, sum(x[(i,j)] for j in out_neighbors[i]) == y[i])
    end
    
    @constraint(model, y[s] == 1)
    @constraint(model, y[t] == 1)
    
    for k in 1:m
        region_nodes = [i for i in 1:n if regions[i] == k]
        if !isempty(region_nodes)
            @constraint(model, sum(y[i] for i in region_nodes) >= 1)
        end
    end
    
    @constraint(model, sum(y[i] for i in 1:n) >= vmin)
    
    # Subtour elimination via iterative approach
    # (JuMP lazy callbacks require solver support; 
    #  alternative: iterative solve-and-cut loop)
    
    while true
        optimize!(model)
        
        if termination_status(model) != MOI.OPTIMAL
            println("No optimal solution found: ", termination_status(model))
            break
        end
        
        # Extract solution
        xsol = Dict((i,j) => round(Int, value(x[(i,j)])) for (i,j) in arcs)
        ysol = round.(Int, value.(y))
        
        # Build solution graph and find connected components
        # among visited nodes
        visited = [i for i in 1:n if ysol[i] == 1]
        
        # Find connected components via BFS/DFS on the arc solution
        adj = Dict(i => Int[] for i in visited)
        for (i,j) in arcs
            if xsol[(i,j)] == 1
                push!(adj[i], j)
            end
        end
        
        # BFS from s
        reached = Set{Int}()
        queue = [s]
        push!(reached, s)
        while !isempty(queue)
            curr = popfirst!(queue)
            for next in get(adj, curr, Int[])
                if next ∉ reached
                    push!(reached, next)
                    push!(queue, next)
                end
            end
        end
        
        # Check if all visited nodes are reachable from s
        unreached_visited = setdiff(Set(visited), reached)
        
        if isempty(unreached_visited)
            println("Valid path found! No subtours.")
            break
        end
        
        # Add SEC for each disconnected component
        # Find connected components among unreached visited nodes
        remaining = collect(unreached_visited)
        while !isempty(remaining)
            # BFS to find one component
            comp = Set{Int}()
            bfs_queue = [first(remaining)]
            push!(comp, first(remaining))
            while !isempty(bfs_queue)
                curr = popfirst!(bfs_queue)
                for next in get(adj, curr, Int[])
                    if next ∈ unreached_visited && next ∉ comp
                        push!(comp, next)
                        push!(bfs_queue, next)
                    end
                end
                # Also check reverse arcs for undirected connectivity
                for prev in visited
                    if haskey(xsol, (prev, curr)) && xsol[(prev, curr)] == 1
                        if prev ∈ unreached_visited && prev ∉ comp
                            push!(comp, prev)
                            push!(bfs_queue, prev)
                        end
                    end
                end
            end
            
            S = collect(comp)
            # SEC: Σ_{i,j ∈ S} x_{ij} ≤ |S| - 1
            arcs_in_S = [(i,j) for (i,j) in arcs if i ∈ comp && j ∈ comp]
            if !isempty(arcs_in_S)
                @constraint(model, 
                    sum(x[(i,j)] for (i,j) in arcs_in_S) <= length(S) - 1)
            end
            
            remaining = setdiff(remaining, comp) |> collect
        end
        
        println("Added SEC constraints, re-solving...")
    end
    
    return model, x, y
end
```

---

## 10. Heuristic Approach

### Strategy: Greedy Construction + 2-opt Local Search

#### Phase 1: Greedy Construction

```
Algorithm: Greedy Path Construction
────────────────────────────────────
1. Start at s, mark all regions as uncovered
2. While not all regions covered OR |visited| < V_min:
   a. Among reachable neighbors (within R):
      - Prioritize nodes in uncovered regions
      - Among those, pick the nearest one
      - If no uncovered-region node reachable, pick nearest overall
   b. Move to chosen node, mark its region as covered
3. Once all regions covered and V_min reached:
   a. Find shortest path to t via reachable nodes (Dijkstra on range graph)
   b. Append this path
```

#### Phase 2: 2-opt Improvement

```
Algorithm: 2-opt for Path Improvement
──────────────────────────────────────
1. Let path = [s, v1, v2, ..., vk, t]
2. For all i < j in the path (excluding s and t endpoints):
   a. Reverse the segment [v_i, ..., v_j]
   b. Check feasibility (all consecutive distances ≤ R)
   c. If feasible and shorter -> accept
3. Repeat until no improvement found
```

#### Phase 3: Node Insertion/Removal

```
- Try removing each non-essential node (keep region coverage + V_min)
- Try inserting skipped nodes if they create a shortcut
```

### Alternative: Simulated Annealing (from MESIM course)

```
Algorithm: SA for Aerodrome Path
────────────────────────────────
1. Start with greedy solution, cost C_current
2. T ← T_0 (initial temperature)
3. Repeat:
   a. Generate neighbor by one of:
      - Swap two intermediate nodes
      - Insert a random reachable node
      - Remove a non-essential node
      - 2-opt reversal
   b. Check feasibility (range, regions, V_min)
   c. ΔC = C_new - C_current
   d. Accept if ΔC < 0 or with probability exp(-ΔC / T)
   e. T ← α · T  (cooling, α ≈ 0.995)
4. Until frozen
```

---

## 11. Julia Heuristic Code

```julia
function greedy_heuristic(inst)
    n, s, t, vmin, m = inst.n, inst.s, inst.t, inst.vmin, inst.m
    regions, R, coords = inst.regions, inst.R, inst.coords
    d = compute_distances(coords, n)
    
    # Track state
    path = [s]
    visited = Set([s])
    covered_regions = Set{Int}()
    if regions[s] > 0
        push!(covered_regions, regions[s])
    end
    
    current = s
    
    # Phase 1: Greedy with region priority
    while length(covered_regions) < m || length(visited) < vmin
        # Find reachable unvisited neighbors
        candidates = [(j, d[current, j]) for j in 1:n 
                       if j ∉ visited && d[current, j] <= R && d[current, j] > 0]
        
        if isempty(candidates)
            println("Warning: stuck at node $current, no reachable unvisited nodes")
            break
        end
        
        # Prioritize uncovered regions
        uncovered_cands = [(j, dist) for (j, dist) in candidates 
                           if regions[j] > 0 && regions[j] ∉ covered_regions]
        
        if !isempty(uncovered_cands)
            sort!(uncovered_cands, by=x->x[2])
            next = uncovered_cands[1][1]
        else
            sort!(candidates, by=x->x[2])
            next = candidates[1][1]
        end
        
        push!(path, next)
        push!(visited, next)
        if regions[next] > 0
            push!(covered_regions, regions[next])
        end
        current = next
    end
    
    # Phase 2: Go to t
    if current != t
        if d[current, t] <= R
            push!(path, t)
        else
            # Dijkstra to find shortest path to t within range
            shortest = dijkstra_path(d, R, n, current, t, visited)
            if shortest !== nothing
                append!(path, shortest[2:end])
            else
                println("Warning: cannot reach t from current position")
            end
        end
    end
    
    # Compute total distance
    total_dist = sum(d[path[i], path[i+1]] for i in 1:length(path)-1)
    
    return path, total_dist
end

function dijkstra_path(d, R, n, src, dst, already_visited)
    # Dijkstra on range graph from src to dst
    dist_to = fill(typemax(Int), n)
    dist_to[src] = 0
    prev = zeros(Int, n)
    pq = [(0, src)]  # (distance, node)
    
    while !isempty(pq)
        sort!(pq, by=x->x[1])
        (cd, u) = popfirst!(pq)
        
        if u == dst
            # Reconstruct path
            path = [dst]
            curr = dst
            while curr != src
                curr = prev[curr]
                pushfirst!(path, curr)
            end
            return path
        end
        
        if cd > dist_to[u]
            continue
        end
        
        for v in 1:n
            if v != u && d[u,v] <= R && d[u,v] > 0
                nd = cd + d[u,v]
                if nd < dist_to[v]
                    dist_to[v] = nd
                    prev[v] = u
                    push!(pq, (nd, v))
                end
            end
        end
    end
    
    return nothing  # unreachable
end
```

---

## 12. Verification: Instance `instance_6_1.txt`

```
n=6, s=1, t=5, V_min=4, m=2
regions = [1, 2, 1, 2, 0, 1]
R = 6
coords: (0,10), (1,3), (1,7), (5,3), (2,8), (6,6)
```

**Distance matrix** (rounded Euclidean):
```
d(1,2) = round(sqrt(1+49)) = round(7.07) = 7
d(1,3) = round(sqrt(1+9))  = round(3.16) = 3
d(1,5) = round(sqrt(4+4))  = round(2.83) = 3
d(2,3) = round(sqrt(0+16)) = 4
d(2,4) = round(sqrt(16+0)) = 4
d(2,5) = round(sqrt(1+25)) = round(5.10) = 5
d(3,4) = round(sqrt(16+16))= round(5.66) = 6
d(3,5) = round(sqrt(1+1))  = round(1.41) = 1
d(4,5) = round(sqrt(9+25)) = round(5.83) = 6
d(4,6) = round(sqrt(1+9))  = round(3.16) = 3
d(5,6) = round(sqrt(16+4)) = round(4.47) = 4
```

With R=6, arcs exist for d ≤ 6. So d(1,2)=7 > R -> no arc (1,2).

**Optimal route**: 1 -> 3 -> 2 -> 5, distance = 3 + 4 + 5 = 12. ✓

This visits: aerodrome 1 (region 1), 3 (region 1), 2 (region 2), 5 (no region) -> both regions covered, 4 aerodromes visited ≥ V_min=4. ✓

---

## 13. Report Checklist

For each instance, report:

| Metric | Model 1 (MTZ) | Model 2 (SEC) | Heuristic |
|--------|---------------|---------------|-----------|
| Objective value (distance) | | | |
| Optimality status | Optimal / Gap% | Optimal / Gap% | No guarantee |
| Computation time (s) | | | |
| Number of B&B nodes | | | |
| Number of lazy cuts (SEC) | N/A | | N/A |

### Expected Observations
- **MTZ**: Faster to code, weaker LP relaxation, more B&B nodes
- **SEC**: Stronger LP relaxation, fewer nodes, but needs callback/iterative implementation
- **Heuristic**: Scales to large instances, fast, but no optimality guarantee
- Both ILP models should find the same optimal solution
- For large instances (80+ aerodromes), ILP may time out -> heuristic is essential

---

## 14. Useful Julia Packages

```julia
using JuMP       # Modeling
using HiGHS      # LP/MIP solver (free, fast)
using Cbc        # Alternative MIP solver (free)
using GLPK       # Alternative (free)
using Printf     # Formatting output
using Plots      # Visualization (optional)
```

---

## 15. References

- CORO Cours 1: LP Duality, Complementary Slackness
- CORO Cours 2: Branch-and-Bound, ILP
- CORO Cours 3: Gomory Cutting Planes, Separation Problem
- CORO Cours 4: Lagrangian Relaxation, Subgradient Method
- CORO Cours 5: Column Generation, Dantzig-Wolfe Decomposition
- Miller, Tucker, Zemlin (1960): "Integer Programming Formulation of TSP"
- Dantzig, Fulkerson, Johnson (1954): "Solution of a Large-Scale TSP" (SEC formulation)