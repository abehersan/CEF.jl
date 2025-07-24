function simulated_annealing(objective, initial_solution; 
                             temp_init=1.0, temp_min=1e-3, alpha=0.9, 
                             max_iter=1000, neighbor_fn=x -> x .+ 0.01 .* randn(length(x)))
    # objective: function to minimize
    # initial_solution: starting point
    # temp_init: initial temperature
    # temp_min: minimum temperature
    # alpha: cooling rate
    # max_iter: maximum number of iterations
    # neighbor_fn: function to generate a neighbor solution

    current = initial_solution
    best = current
    current_score = objective(current)
    best_score = current_score
    temp = temp_init

    for iter in 1:max_iter
        if temp < temp_min
            break
        end

        # Generate neighbor
        neighbor = neighbor_fn(current)
        neighbor_score = objective(neighbor)

        Δ = neighbor_score - current_score
        if Δ < 0 || rand() < exp(-Δ / temp)
            current = neighbor
            current_score = neighbor_score
            if current_score < best_score
                best = current
                best_score = current_score
            end
        end

        temp *= alpha
    end

    return best, best_score
end