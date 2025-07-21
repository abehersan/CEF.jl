function hc_units(units::Symbol)
    if isequal(units, :SI)
        return Rg
    else
        @error "Units not supported in heat capacity calculations. Use SI, (J/K/mol)"
    end
end


function mag_heatcap(Ep::Vector{Float64}, T::Real)::Float64
    np = population_factor(Ep, T)
    heatcap = sum( ( Ep ./ (kB*T) ) .^2 .* np) - sum( Ep ./ (kB*T) .* np )^2
    return heatcap
end


function mag_entropy(HC::Vector{Float64}, T::Vector{Float64})::Vector{Float64}
    S = similar(HC)
    @views @inbounds for i in eachindex(T)
        S[i] = trapz(T[1:i], HC[1:i] ./ T[1:i])
    end
    return S
end


function cef_entropy!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; B::Vector{<:Real}=[0.0,0.0,0.0], units::Symbol=:SI, method::Symbol=:O)::Nothing
    unit_factor=hc_units(units)
    E=eigvals(cef_hamiltonian(ion,cefparams;B=B,method=method))
    E .-= minimum(E)
    dfcalc[:,:HC_CALC].=[mag_heatcap(E,tt) for tt in dfcalc[:,:T]]
    dfcalc[:,:HC_CALC]*=unit_factor
    dfcalc[:,:SM_CALC]=mag_entropy(dfcalc[:,:HC_CALC],dfcalc[:,:T])
    return nothing
end


function cef_entropy!(lfield::local_env, dfcalc::DataFrame; B::Vector{<:Real}=[0.0,0.0,0.0], units::Symbol=:SI, method::Symbol=:O)::Nothing
    if isempty(lfield.cefparams)
        calc_cefparams!(lfield)
    end
    cef_entropy!(lfield.ion,lfield.cefparams,dfcalc;B,units,method)
    return nothing
end


function cef_entropy_speclevels!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; B::Vector{<:Real}=[0.0,0.0,0.0], levels::UnitRange=1:4, units::Symbol=:SI, method::Symbol=:O)::Nothing
    # only levels specified contribute (2J+1 levels total)
    unit_factor=hc_units(units)
    E=eigvals(cef_hamiltonian(ion,cefparams;B=B,method=method))
    E .-= minimum(E)
    E=E[levels]
    dfcalc[:,:HC_CALC]*=unit_factor
    dfcalc[:,:SM_CALC]=mag_entropy(dfcalc[:,:HC_CALC],dfcalc[:,:T])
    return nothing
end


function cef_entropy_speclevels!(lfield::local_env, dfcalc::DataFrame; B::Vector{<:Real}=[0.0,0.0,0.0], levels::UnitRange=1:4, units::Symbol=:SI, method::Symbol=:O)::Nothing
    if isempty(lfield.cefparams)
        calc_cefparams!(lfield)
    end
    cef_entropy_speclevels!(lfield.ion,lfield.cefparams,dfcalc;B,levels,units,method)
    return nothing
end