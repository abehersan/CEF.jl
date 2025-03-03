@doc raw"""
    blm_dframe(blm_dict::Dict{String, <:Real})::DataFrame

Given a dictionary of Stevens coefficients of the form Blm -> Value, return
a DataFrame with equivalent information.
"""
function blm_dframe(blm_dict::Dict{String, <:Real})::DataFrame
    l = Int[]
    m = Int[]
    B = Float64[]
    for k in keys(blm_dict)
        ll, mm = parse_blm(k)
        BB = blm_dict[k]
        push!(l, ll)
        push!(m, mm)
        push!(B, BB)
    end
    return DataFrame("B"=>B, "l"=>l, "m"=>m)
end


function parse_blm(b::String)::Tuple{Int, Int}
    if b[3] == 'm'
        l = parse(Int, b[2])
        m = -1*parse(Int, b[end])
    else
        l = parse(Int, b[2])
        m = parse(Int, b[end])
    end
    return (l, m)
end


@doc raw"""
    get_alm!(bfactors::DataFrame, single_ion::mag_ion)::DataFrame

Factorization of the Stevens B_lm parameters defined as
B_lm = A_lm * <r^l> * theta_l,
where <r^l> is the expectation value of the radial wavefunction of the
4f electron density (tabulated) and theta_l are the Stevens geometrical factors.
"""
function get_alm!(bfactors::DataFrame, single_ion::mag_ion)
    alpha, beta, gamma = single_ion.stevens_factors
    r2, r4, r6 = single_ion.rad_wavefunction
    alm = zeros(Float64, nrow(bfactors))
    for (i, r) in enumerate(eachrow(bfactors))
        if r.l == 2
            alm[i] = r.B / (alpha * r2)
        elseif r.l == 4
            alm[i] = r.B / (beta * r4)
        elseif r.l == 6
            alm[i] = r.B / (gamma * r6)
        else
            alm[i] = r.B
            err_message =
            "Given Blm parameter has invalid l and/or m.\n"*
            "l=$l and m=$m were parsed and are not supported.\n"*
            "Writing unscaled Stevens parameter."
            @error err_message
        end
    end
    bfactors[!, :A] .= alm
    return nothing
end