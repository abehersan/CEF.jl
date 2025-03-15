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


function cef_magneticmoment_crystal!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame, T::Real; units::Symbol=:ATOMIC, method::Symbol=:EO, mode::Function=real)
    unit_factor=mag_units(units)
    spinops=[ion.Jx,ion.Jy,ion.Jz]
    @eachrow! dfcalc begin
        @newcol :M_CALC::Vector{Float64}
        extfield = [:Bx,:By,:Bz]
        E, V = eigen(cef_hamiltonian(ion,cefparams;B=extfield,method=method))
        E .-= minimum(E)
        :M_CALC=sum([
                thermal_average(Ep=E,Vp=V,op=spinops[1],T=T,mode=mode),
                thermal_average(Ep=E,Vp=V,op=spinops[2],T=T,mode=mode),
                thermal_average(Ep=E,Vp=V,op=spinops[3],T=T,mode=mode)
            ])*ion.gj*unit_factor
    end
    return nothing
end


function cef_magneticmoment_powder!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame, T::Real; units::Symbol=:ATOMIC, method::Symbol=:EO, mode::Function=real)
    unit_factor = mag_units(units)
    spinops = [ion.Jx,ion.Jy,ion.Jz]
    @eachrow! dfcalc begin
        @newcol :M_CALC::Vector{Float64}
        E, V = eigen(cef_hamiltonian(ion,cefparams; B=[:B,0.0,0.0],method=method))
        E .-= minimum(E)
        MX = thermal_average(Ep=E,Vp=V,op=spinops[1],T=T,mode=mode)

        E, V = eigen(cef_hamiltonian(ion,cefparams; B=[0.0,:B,0.0],method=method))
        E .-= minimum(E)
        MY = thermal_average(Ep=E,Vp=V,op=spinops[2],T=T,mode=mode)

        E, V = eigen(cef_hamiltonian(ion,cefparams; B=[0.0,0.0,:B],method=method))
        E .-= minimum(E)
        MZ = thermal_average(Ep=E,Vp=V,op=spinops[3],T=T,mode=mode)

        :M_CALC=round(((MX + MY + MZ) / 3.0) * unit_factor,digits=SDIG)
    end
    return nothing
end