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


function cef_entropy!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; units::Symbol=:SI, method::Symbol=:EO)::Nothing
    convfac = hc_units(units)
    @eachrow! dfcalc begin
        @newcol :HC_CALC::Vector{Float64}
        B=[:Bx,:By,:Bz]
        E=eigvals(cef_hamiltonian(ion,cefparams;B=B,method=method))
        E .-= minimum(E)
        :HC_CALC=round(mag_heatcap(E,:T)*convfac,digits=SDIG)
    end
    dfcalc[:,:SM_CALC]=round.(mag_entropy(dfcalc[:,:HC_CALC],dfcalc[:,:T]),digits=SDIG)
    return nothing
end


function cef_entropy_speclevels!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; levels::UnitRange=1:4, units::Symbol=:SI, method::Symbol=:EO)::Nothing
    # only levels specified contribute (2J+1 levels total)
    convfac = hc_units(units)
    @eachrow! dfcalc begin
        @newcol :HC_CALC::Vector{Float64}
        B=[:Bx,:By,:Bz]
        E=eigvals(cef_hamiltonian(ion,cefparams;B=B,method=method))
        E .-= minimum(E)
        E=E[levels]
        :HC_CALC=round(mag_heatcap(E,:T)*convfac,digits=SDIG)
    end
    dfcalc[:,:SM_CALC]=round.(mag_entropy(dfcalc[:,:HC_CALC],dfcalc[:,:T]),digits=SDIG)
    return nothing
end