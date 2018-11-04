using Random, Base.Iterators, Plots

function isvalid_row(posn, soln)
    all(skipmissing(.!xor.(soln, posn))) # XNOR gate applied to (posn, soln) elementwise. if any bool values differ, returns false. returns true if all missing
end

function isvalid_board(soln, specs)
    n_rows, n_cols = size(soln)
    retval = try
        for i=1:n_rows
            spec = specs[('r', i)]
            exp_spec = get_expected_spec(soln[i,:])
            for (s, es) in zip(spec, exp_spec)
                if s != es
                    return false
                end
            end
        end
        for j=1:n_cols
            spec = specs[('c', j)]
            exp_spec = get_expected_spec(soln[:,j])
            for (s, es) in zip(spec, exp_spec)
                if s != es
                    return false
                end
            end
        end
        return true
    catch
        return false
    end
    return retval
end

function iscomplete_board(soln)
    return all(soln .!== missing)
end

function get_expected_spec(row_soln)
    spec = []
    c = 0
    for cell in row_soln
        if ismissing(cell) || cell == false
            c > 0 && push!(spec, c)
            c = 0
        else
            c += 1
        end
    end
    c > 0 && push!(spec, c)
    return spec
end

function construct_row(spec, posn, n)
    row = zeros(n)
    for (s,p) in zip(spec, posn)
        row[p:p+s-1] .= 1
    end
    return row
end

function get_fill_freqs(current_soln, current_spec)  # i think there's a bug in here about how it handles the last spec in the row
    n_cells = length(current_soln)
    if length(current_spec) == 0 || current_spec[1] == 0
        isvalid_row(falses(n_cells), current_soln) && return (1, zeros(n_cells))
        return (0, zeros(n_cells))
    end
    current_run = current_spec[1]
    limiter = sum(current_spec[2:end]) + length(current_spec[2:end])
    n_aggd_solutions = 0
    aggd_freqs = zeros(n_cells)
    for i in 1:n_cells-limiter-current_run+1
        current_cell = i + current_run
        this_freqs = zeros(n_cells)
        this_freqs[i:current_cell-1] .= 1
        current_cell = min(current_cell, n_cells)  # this is to avoid BoundsError
        this_posn = map(x->convert(Bool, x), this_freqs[1:current_cell]) # note: includes blankspace after run if run is not last in spec
        if isvalid_row(this_posn, current_soln[1:current_cell])
            n_valid_solutions, fill_freqs = get_fill_freqs(current_soln[current_cell+1:end], current_spec[2:end])
            this_freqs .*= n_valid_solutions
            this_freqs[current_cell+1:end] = fill_freqs
            n_aggd_solutions += n_valid_solutions
            aggd_freqs .+= this_freqs
        end
    end
    return n_aggd_solutions, aggd_freqs
end

function get_specs(fname)
    """Reads .NIN files."""
    n_rows, n_cols, specs = open(fname) do f
        lines = readlines(f)
        ns = split(lines[1])
        n_cols = parse(Int, ns[1]); n_rows = parse(Int, ns[2])  # note: n_cols first
        rowspecs = Dict{Tuple{Char, Int}, Array{Int}}(tuple('r', i) => map(x->parse(Int, x), split(lines[i+1]))        for i in 1:n_rows)
        colspecs = Dict{Tuple{Char, Int}, Array{Int}}(tuple('c', i) => map(x->parse(Int, x), split(lines[i+n_rows+1])) for i in 1:n_cols)
        specs = merge(rowspecs, colspecs)
        return n_rows, n_cols, specs
    end
end

function solve(fname)
    n_rows, n_cols, specs = get_specs(fname)
    solution = Array{Union{Missing,Bool}}(missing, n_rows, n_cols)
    probs = init_probs(n_rows, n_cols, specs)
    priorities = init_priorities(n_rows, n_cols, specs)
    q = init_queue(priorities)
    solution, probs, priorities = deduce(solution, probs, q, priorities, specs)
    return solution
end

function update_row(row_soln, row_spec)
    n_positions, fill_frequencies = get_fill_freqs(row_soln, row_spec)
    previously_undetermined_cells = (1:length(row_soln))[ismissing.(row_soln)]
    fills = (1:length(row_soln))[fill_frequencies.==n_positions]; row_soln[fills] .= true
    blanks = (1:length(row_soln))[fill_frequencies.==0]; row_soln[blanks] .= false
    modified_cells = intersect(previously_undetermined_cells, union(fills, blanks))
    fill_probabilities = fill_frequencies ./ n_positions
    return row_soln, modified_cells, fill_probabilities, n_positions
end

function deduce(solution, probs, q, priorities, specs, depth=0)
    old_solution = deepcopy(solution)
    old_probs = deepcopy(probs)
    old_priorities = deepcopy(priorities)
    score = 0
    #q = init_queue(priorities)
    while length(q) > 0
        next = popfirst!(q)
        pprint(solution)
        k = next[2]; ind = k[2]
        if k[1] == 'r'
            row_before = solution[ind,:]
        elseif k[1] == 'c'
            row_before = solution[:,ind]
        end
        n_cells = length(row_before)
        row_soln = deepcopy(row_before)
        n_posns, row_states = get_fill_freqs(row_before, specs[k])
        previously_undetermined_cells = (1:n_cells)[ismissing.(row_before)]
        row_probs = n_posns == 0 ? zeros(n_cells) : row_states ./ n_posns
        fills  = (1:n_cells)[row_probs.==1]; row_soln[fills]  .= true
        blanks = (1:n_cells)[row_probs.==0]; row_soln[blanks] .= false
        modified_cells = intersect(previously_undetermined_cells, union(fills, blanks))
        if !isvalid_row(row_before, row_soln)
            return old_solution, old_probs, old_priorities
        elseif iscomplete_board(solution)
            return solution, probs, priorities
        end
        update_soln!(solution, k, row_soln)
        update_probs!(probs, k, row_probs)
        update_priorities!(priorities, k, n_posns)
        update_queue!(q, k, modified_cells, priorities)
    end
    if !iscomplete_board(solution)
        trial_solution = similar(solution)
        normalize_probs!(probs, solution)
        #selected_cell = find_min_state_count(probs, priorities, solution)
        #selected_cell = find_max_state_count(probs, priorities, solution)
        #selected_cell = find_wgtd_lowest_entropy(probs, priorities, solution)
        #selected_cell = find_wgtd_highest_entropy(probs, priorities, solution)
        #selected_cell = find_most_decisive(probs)
        #selected_cell = find_wgtd_most_decisive(probs, priorities)
        selected_cell = find_lowest_entropy(probs, solution)
        #selected_cell = find_highest_entropy(probs, solution)
        #selected_cell = uniform_random(probs)
        for (p, (i,j)) in selected_cell
            trial_solution = deepcopy(solution)
            trial_solution[i, j] = convert(Bool, round(p))
            ks = [tuple('r', i), tuple('c', j)]
            q = sort([tuple(priorities[k], k) for k in ks])
            trial_solution, probs, priorities = deduce(trial_solution, probs, q, priorities, specs, depth+1)
            iscomplete_board(trial_solution) && isvalid_board(trial_solution, specs) && break
        end
        solution = trial_solution
    end
    return solution, probs, priorities
end

function entropy(probs)
    p = probs[:, :, 1] .* probs[:, :, 2]
    return (probs[:,:,1].+probs[:,:,2])/2, map(x->(x==0. || x==1.) ? 0. : -x*log(x), p)
end

function delta(probs)
    return abs.(probs[:, :, 1] .- probs[:, :, 2])
end

function state_count(probs, priorities)
    n_rows, n_cols, two = size(probs)
    row_state_counts = [priorities[('r', i)] for i in 1:n_rows]
    col_state_counts = [priorities[('c', j)] for j in 1:n_cols]
    return row_state_counts * col_state_counts'
end

function uniform_random(probs)
    n_rows, n_cols, two = size(probs)
    i = rand(1:n_rows); j = rand(1:n_cols)
    loc = CartesianIndex{2}((i,j))
    p = probs[:,:,1] .* probs[:,:,2]
    return [tuple(  p[loc], convert(Tuple, loc)),
            tuple(1-p[loc], convert(Tuple, loc))]
end

function find_lowest_entropy(probs, soln)
    p, ent = entropy(probs)
    ent[.!ismissing.(soln)] .= 1.
    loc = findmin(ent)[2]
    return [tuple(  p[loc], convert(Tuple, loc)),
            tuple(1-p[loc], convert(Tuple, loc))]  # only need these two. if cell `loc` is not filled, then it's blank
end

function find_highest_entropy(probs, soln)
    p, ent = entropy(probs)
    ent[.!ismissing.(soln)] .= 0.
    loc = findmax(ent)[2]
    return [tuple(  p[loc], convert(Tuple, loc)),
            tuple(1-p[loc], convert(Tuple, loc))]  # only need these two. if cell `loc` is not filled, then it's blank
end

function find_most_decisive(probs, soln)
    d = delta(probs)
    d[.!ismissing.(soln)] .= 0.
    loc = findmax(d)[2]
    return [tuple(0, convert(Tuple, loc)),
            tuple(1, convert(Tuple, loc))]
end

function find_most_sure(probs)
    p = (probs[:, :, 1] .+ probs[:, :, 2]) ./ 2
    d = delta(probs)
    d[probs[:,:,1].==0.] .= 1.
    d[probs[:,:,1].==1.] .= 1.
    d[probs[:,:,2].==0.] .= 1.
    d[probs[:,:,2].==1.] .= 1.
    loc = findmin(d)[2]
    return [tuple(  p[loc], convert(Tuple, loc)),
            tuple(1-p[loc], convert(Tuple, loc))]
end

function find_min_state_count(probs, priorities, soln)
    p = probs[:, :, 1] .* probs[:, :, 2]
    state_counts = state_count(probs, priorities) .* p
    state_counts[.!ismissing.(soln)] .= 1e24
    loc = findmin(state_counts)[2]
    return [tuple(  p[loc], convert(Tuple, loc)),
            tuple(1-p[loc], convert(Tuple, loc))]
end

function find_max_state_count(probs, priorities, soln)
    p = probs[:, :, 1] .* probs[:, :, 2]
    state_counts = state_count(probs, priorities) .* p
    state_counts[.!ismissing.(soln)] .= 0
    loc = findmax(state_counts)[2]
    return [tuple(  p[loc], convert(Tuple, loc)),
            tuple(1-p[loc], convert(Tuple, loc))]
end

function find_wgtd_lowest_entropy(probs, priorities, soln)
    p, ent = entropy(probs)
    state_counts = state_count(probs, priorities) .* p
    state_counts .*= (log(2)/2) / maximum(state_counts)
    wgtd_ent = ent .- state_counts
    wgtd_ent[.!ismissing.(soln)] .= 1.1 * log(2)/2
    loc = findmin(wgtd_ent)[2]
    return [tuple(  p[loc], convert(Tuple, loc)),
            tuple(1-p[loc], convert(Tuple, loc))]
end

function find_wgtd_highest_entropy(probs, priorities, soln)
    p, ent = entropy(probs)
    state_counts = state_count(probs, priorities) .* p
    state_counts .*= (log(2)/2) / maximum(state_counts)
    wgtd_ent = ent .+ state_counts
    wgtd_ent[.!ismissing.(soln)] .= 0.
    loc = findmax(wgtd_ent)[2]
    return [tuple(  p[loc], convert(Tuple, loc)),
            tuple(1-p[loc], convert(Tuple, loc))]
end

function find_wgtd_most_decisive(probs, priorities)
    d = delta(probs)
    p, ent = entropy(probs)
    state_counts = state_count(probs, priorities)
    wgtd_d = d ./ state_counts
    loc = findmax(wgtd_d)[2]
    return [tuple(  p[loc], convert(Tuple, loc)),
            tuple(1-p[loc], convert(Tuple, loc))]
end

function init_priorities(n_rows, n_cols, specs)
    priorities = Dict{Tuple{Char, Int}, Int}()
    for (k,spec) in specs
        if k[1] == 'r'
            L = n_cols
        elseif k[1] == 'c'
            L = n_rows
        end
        if L-(sum(spec)+length(spec)-1) < maximum(spec)  # in this condition, there are guaranteed filled blocks in this row
            priorities[k] = -sum([max(s-(L-(sum(spec)+length(spec))), 0) for s in spec])
        else  # no guaranteed blocks
            priorities[k] = prod([L-(sum(spec)+length(spec)-1)-i for i in 0:length(spec)-1]) #binomial(L-(sum(spec)+length(spec)), L-(sum(spec)+length(spec))-2*length(spec))
        end
    end
    return priorities
end

function update_soln!(soln, k, row)
    if k[1] == 'r'
        soln[k[2],:] = row
    elseif k[1] == 'c'
        soln[:,k[2]] = row
    end
end

function init_probs(n_rows, n_cols, specs)
    probs = zeros(n_rows, n_cols, 2)
    for (k,spec) in specs
        if k[1] == 'r'
            probs[k[2], :, 1] .= sum(spec) / n_cols
        elseif k[1] == 'c'
            probs[:, k[2], 2] .= sum(spec) / n_rows
        end
    end
    return probs
end

function update_probs!(probs, k, row_probs)
    if k[1] == 'r'
        probs[k[2], :, 1] .= row_probs
    elseif k[1] == 'c'
        probs[:, k[2], 2] .= row_probs
    end
end

function normalize_probs!(probs, solution)
    n_rows, n_cols = size(solution)
    for i=1:n_rows
        for j=1:n_cols
            if !ismissing(solution[i,j])
                if solution[i,j]
                    probs[i,j,:] .= 1.
                else
                    probs[i,j,:] .= 0.
                end
            end
        end
    end
end

function init_queue(priorities)
    q = sort([tuple(v,k) for (k,v) in pairs(priorities)])
    return q
end

function update_queue!(q, k, modified, priorities)
    if k[1]=='r'
        newrc = 'c'
    elseif k[1]=='c'
        newrc = 'r'
    end
    for m in modified
        newk = tuple(newrc, m)
        filter!(x->x[2]!=newk, q)                # remove old entry
        push!(q, tuple(priorities[newk], newk))  # add new entry with (hopefully) updated priority
    end
    sort!(q)
end

function update_priorities!(priorities, k, n_posns)
    priorities[k] = n_posns
end

function pprint(soln)
    out = ""
    n_rows, n_cols = size(soln)
    out *= "+" * "-"^(2*n_cols) * "+\n"
    for i in 1:n_rows
        out *= "|"
        for j in 1:n_cols
            if ismissing(soln[i,j])
                out *= "\u2591"^2
            elseif soln[i,j]
                out *= "\u2588"^2
            else
                out *= " "^2
            end
        end
        out *= "|\n"
    end
    out *= "+" * "-"^(2*n_cols) * "+\n"
    print(out)
end

function pplot(soln, probs)

end

function init_plots()
    gr(reuse=true)
    gui()
end

for arg in ARGS
    soln = solve(arg)
end
