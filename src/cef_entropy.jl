function hc_units(units::Symbol)
    if isequal(units, :SI)
        return Rg
    else
        @error "Units not supported in heat capacity calculations. Use SI, [J/K/mol]"
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


@doc raw"""
    cef_entropy!(cefsys::cef_system, dfcalc::DataFrame; units::String="SI")::Nothing

Calculate the temperature-dependence of the magnetic entropy and specific
heat capacity of a magnetic system consisting of non-interacting magnetic ions.
The details of the CEF system are given in `cefsys`, which is a `cef_system` struct.

The calculation is performed on a `DataFrame` `dfcalc`.
`dfcalc` must have columns `[:T, :Bx, :By, :Bz]`.

The entropy and heat capacity are calculated in `:SI` units (J/mol/K).
"""
function cef_entropy!(cefsys::cef_system, dfcalc::DataFrame; units::Symbol=:SI, method::Symbol=:EO)::Nothing
    convfac = hc_units(units)
    @eachrow! dfcalc begin
        @newcol :HC_CALC::Vector{Float64}
        B=[:Bx,:By,:Bz]
        E=eigvals(cef_hamiltonian(cefsys.ion,cefsys.cefparams;B=B,method=method))
        E .-= minimum(E)
        :HC_CALC=mag_heatcap(E,:T)*convfac
    end
    dfcalc[:,:SM_CALC] = mag_entropy(dfcalc[:,:HC_CALC],dfcalc[:,:T])
    return nothing
end


@doc raw"""
    cef_entropy_speclevels!(single_ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; levels::UnitRange=1:4, units::String="SI")::Nothing

Calculate the temperature-dependence of the magnetic entropy and specific
heat capacity of a magnetic system consisting of magnetic ions of type `single_ion`.
The CEF parameters are given in the `cefparams` `DataFrame`.
The calculation is performed on the `DataFrame` `dfcalc`.

`dfcalc` must have column `[:T, :Bx, :By, :Bz]`.

`levels` is a step range that specifies the indices of the contributing crystal
field energy levels.

As per the default it is `1:4` which means that the first
4 levels contribute to the entropy of the system.

The entropy and heat capacity are calculated in `SI` units [J/mol/K].
"""
function cef_entropy_speclevels!(cefsys::cef_system, dfcalc::DataFrame;
        levels::UnitRange=1:4, units::Symbol=:SI, method::Symbol=:EO)::Nothing
    # only levels specified contribute (2J+1 levels total)
    convfac = hc_units(units)
    @eachrow! dfcalc begin
        @newcol :HC_CALC::Vector{Float64}
        B=[:Bx,:By,:Bz]
        E=eigvals(cef_hamiltonian(cefsys.ion,cefsys.cefparams;B=B,method=method))
        E .-= minimum(E)
        E=E[levels]
        :HC_CALC=mag_heatcap(E,:T)*convfac
    end
    dfcalc[:,:SM_CALC] = mag_entropy(dfcalc[:,:HC_CALC],dfcalc[:,:T])
    return nothing
end