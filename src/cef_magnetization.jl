function mag_units(units::Symbol)::Float64
    if isequal(units, :SI)
        return 5.5849397            # NA * muB  ( J/T/mol )
    elseif isequal(units, :CGS)
        return 5.5849397*1000.0     # NA * muB  ( emu/mol )
    elseif isequal(units, :ATOMIC)
        return 1.0                  # units of Bohr magneton per mol
    else
        @error "Units $units not understood. Use one of either :SI, :CGS or :ATOMIC"
    end
end


function cef_magnetization_crystal!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; T::Real=1.0, units::Symbol=:ATOMIC, method::Symbol=:O, mode::Function=real)
    unit_factor=mag_units(units)
    @eachrow! dfcalc begin
        @newcol :M_CALC::Vector{Float64}
        extfield = [:Bx,:By,:Bz]
        E, V = eigen(cef_hamiltonian(ion,cefparams;B=extfield,method=method))
        E .-= minimum(E)
        :M_CALC=thermal_average(Ep=E,Vp=V,op=ion.Jx,T=T,mode=mode)+
                thermal_average(Ep=E,Vp=V,op=ion.Jy,T=T,mode=mode)+
                thermal_average(Ep=E,Vp=V,op=ion.Jz,T=T,mode=mode)
    end
    dfcalc[:,:M_CALC]*=(ion.gj*unit_factor)
    return nothing
end


function cef_magnetization_crystal!(lfield::local_env, dfcalc::DataFrame; T::Real=1.0, units::Symbol=:ATOMIC, method::Symbol=:O, mode::Function=real)
    if isempty(lfield.cefparams)
        calc_cefparams!(lfield)
    end
    cef_magnetization_crystal!(lfield.ion,lfield.cefparams,dfcalc;T,units,method,mode)
    return nothing
end


function cef_magnetization_powder!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; T::Real=1.0, units::Symbol=:ATOMIC, method::Symbol=:O, mode::Function=real)
    unit_factor = mag_units(units)
    @eachrow! dfcalc begin
        @newcol :M_CALC::Vector{Float64}
        E, V = eigen(cef_hamiltonian(ion,cefparams; B=[:B,0.0,0.0],method=method))
        E .-= minimum(E)
        MX = thermal_average(Ep=E,Vp=V,op=ion.Jx,T=T,mode=mode)

        E, V = eigen(cef_hamiltonian(ion,cefparams; B=[0.0,:B,0.0],method=method))
        E .-= minimum(E)
        MY = thermal_average(Ep=E,Vp=V,op=ion.Jy,T=T,mode=mode)

        E, V = eigen(cef_hamiltonian(ion,cefparams; B=[0.0,0.0,:B],method=method))
        E .-= minimum(E)
        MZ = thermal_average(Ep=E,Vp=V,op=ion.Jz,T=T,mode=mode)

        :M_CALC=sqrt(MX^2 + MY^2 + MZ^2)
    end
    dfcalc[:,:M_CALC]*=(ion.gj*unit_factor)
    return nothing
end


function cef_magnetization_powder!(lfield::local_env, dfcalc::DataFrame; T::Real=1.0, units::Symbol=:ATOMIC, method::Symbol=:O, mode::Function=real)
    if isempty(lfield.cefparams)
        calc_cefparams!(lfield)
    end
    cef_magnetization_powder!(lfield.ion,lfield.cefparams,dfcalc;T,units,method,mode)
    return nothing
end