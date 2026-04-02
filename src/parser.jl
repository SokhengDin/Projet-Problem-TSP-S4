struct Instance
    n       :: Int
    s       :: Int
    t       :: Int
    vmin    :: Int
    m       :: Int
    regions :: Vector{Int}
    R       :: Int   # range in km
    coords  :: Matrix{Int}
end

# Parse instance file
function parse_instance(filename::String) :: Instance
    lines = readlines(filename)
    lines = [strip(l) for l in lines if !isempty(strip(l)) && !startswith(strip(l), '#')]

    idx = 1
    n    = parse(Int, lines[idx]); idx += 1
    s    = parse(Int, lines[idx]); idx += 1
    t    = parse(Int, lines[idx]); idx += 1
    vmin = parse(Int, lines[idx]); idx += 1
    m    = parse(Int, lines[idx]); idx += 1
    regions = parse.(Int, split(lines[idx])); idx += 1
    @assert length(regions) == n
    R = parse(Int, lines[idx]); idx += 1

    coords = zeros(Int, n, 2)
    for i in 1:n
        vals = parse.(Int, split(lines[idx])); idx += 1
        coords[i, 1] = vals[1]; coords[i, 2] = vals[2]
    end

    @assert 1 <= s <= n && 1 <= t <= n && s != t
    @assert vmin >= 2 && m >= 1 && R > 0
    return Instance(n, s, t, vmin, m, regions, R, coords)
end

# Rounded Euclidean distance matrix (n x n)
function compute_distances(inst::Instance) :: Matrix{Int}
    n = inst.n
    d = zeros(Int, n, n)
    for i in 1:n, j in 1:n
        if i != j
            dx = inst.coords[i,1] - inst.coords[j,1]
            dy = inst.coords[i,2] - inst.coords[j,2]
            d[i,j] = round(Int, sqrt(dx^2 + dy^2))
        end
    end
    return d
end

# Arcs with d[i,j] <= R
function build_arc_sets(inst::Instance, d::Matrix{Int})
    n, R = inst.n, inst.R
    arcs = Tuple{Int,Int}[]
    out_neighbors = [Int[] for _ in 1:n]
    in_neighbors  = [Int[] for _ in 1:n]
    for i in 1:n, j in 1:n
        if i != j && d[i,j] <= R
            push!(arcs, (i,j))
            push!(out_neighbors[i], j)
            push!(in_neighbors[j],  i)
        end
    end
    return arcs, out_neighbors, in_neighbors
end

# Print instance summary
function print_instance_summary(inst::Instance)
    println("  n=$(inst.n)  s=$(inst.s)  t=$(inst.t)  V_min=$(inst.vmin)  m=$(inst.m)  R=$(inst.R)")
end
