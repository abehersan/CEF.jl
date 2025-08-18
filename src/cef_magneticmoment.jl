function mag_units(units::Symbol)::Float64
    if isequal(units,:SI)
        return 5.5849397            # NA * muB  ( J/T/mol )
    elseif isequal(units,:CGS)
        return 5.5849397*1000.0     # NA * muB  ( emu/mol )
    elseif isequal(units,:ATOMIC)
        return 1.0                  # units of Bohr magneton per mol
    else
        @error "Units $units not understood. Use one of either :SI, :CGS or :ATOMIC"
    end
end


function cef_magneticmoment_crystal!(ion::mag_ion,cefparams::DataFrame,dfcalc::DataFrame;NJ::Real=0.0,T::Real=1.0,units::Symbol=:ATOMIC,method::Symbol=:O,mode::Function=real)
    unit_factor=mag_units(units)
    @eachrow! dfcalc begin
        @newcol :M_CALC::Vector{Float64}
        extfield=[:Bx,:By,:Bz]
        E,V=eigen(cef_hamiltonian(ion,cefparams;B=extfield,method=method))
        E .-= minimum(E)
        mux=thermal_average(Ep=E,Vp=V,op=ion.Jx,T=T,mode=mode)
        muy=thermal_average(Ep=E,Vp=V,op=ion.Jy,T=T,mode=mode)
        muz=thermal_average(Ep=E,Vp=V,op=ion.Jz,T=T,mode=mode)
        if iszero(NJ)
            mutot=mux+muy+muz
        else
            beff=extfield .+ (NJ/(ion.gj*muB)^2)*[mux,muy,muz] # molecular field
            E,V=eigen(cef_hamiltonian(ion,cefparams;B=beff,method=method))
            E .-= minimum(E)
            muxp=thermal_average(Ep=E,Vp=V,op=ion.Jx,T=T,mode=mode)
            muyp=thermal_average(Ep=E,Vp=V,op=ion.Jy,T=T,mode=mode)
            muzp=thermal_average(Ep=E,Vp=V,op=ion.Jz,T=T,mode=mode)
            mutot=muxp+muyp+muzp
        end
        :M_CALC=mutot
    end
    dfcalc[:,:M_CALC]*=(ion.gj*unit_factor)
    return nothing
end


function cef_magneticmoment_crystal!(lfield::local_env,dfcalc::DataFrame;NJ::Real=0.0,T::Real=1.0,units::Symbol=:ATOMIC,method::Symbol=:O,mode::Function=real)
    if isempty(lfield.cefparams)
        calc_cefparams!(lfield)
    end
    cef_magneticmoment_crystal!(lfield.ion,lfield.cefparams,dfcalc;NJ,T,units,method,mode)
    return nothing
end


function cef_magneticmoment_powder!(ion::mag_ion,cefparams::DataFrame,dfcalc::DataFrame;NJ::Real=0.0,T::Real=1.0,units::Symbol=:ATOMIC,method::Symbol=:O,mode::Function=real)
    unit_factor = mag_units(units)
    @eachrow! dfcalc begin
        @newcol :M_CALC::Vector{Float64}
        E,V=eigen(cef_hamiltonian(ion,cefparams; B=[:B,0.0,0.0],method=method))
        E .-= minimum(E)
        mux=thermal_average(Ep=E,Vp=V,op=ion.Jx,T=T,mode=mode)

        E,V=eigen(cef_hamiltonian(ion,cefparams; B=[0.0,:B,0.0],method=method))
        E .-= minimum(E)
        muy=thermal_average(Ep=E,Vp=V,op=ion.Jy,T=T,mode=mode)

        E,V=eigen(cef_hamiltonian(ion,cefparams; B=[0.0,0.0,:B],method=method))
        E .-= minimum(E)
        muz=thermal_average(Ep=E,Vp=V,op=ion.Jz,T=T,mode=mode)

        if iszero(NJ)
            :M_CALC=sqrt(mux^2 + muy^2 + muz^2)
        else
            beff=(NJ/(ion.gj*muB)^2)*[mux,muy,muz] # molecular field
            E,V=eigen(cef_hamiltonian(ion,cefparams; B=[:B,0.0,0.0] .+ beff,method=method))
            E .-= minimum(E)
            muxp=thermal_average(Ep=E,Vp=V,op=ion.Jx,T=T,mode=mode)

            E,V=eigen(cef_hamiltonian(ion,cefparams; B=[0.0,:B,0.0] .+ beff,method=method))
            E .-= minimum(E)
            muyp=thermal_average(Ep=E,Vp=V,op=ion.Jy,T=T,mode=mode)

            E,V=eigen(cef_hamiltonian(ion,cefparams; B=[0.0,0.0,:B] .+ beff,method=method))
            E .-= minimum(E)
            muzp=thermal_average(Ep=E,Vp=V,op=ion.Jz,T=T,mode=mode)

            :M_CALC=sqrt(muxp^2 + muyp^2 + muzp^2)
        end
    end
    dfcalc[:,:M_CALC]*=(ion.gj*unit_factor)
    return nothing
end


function cef_magneticmoment_powder!(lfield::local_env,dfcalc::DataFrame;NJ::Real=0.0,T::Real=1.0,units::Symbol=:ATOMIC, method::Symbol=:O, mode::Function=real)
    if isempty(lfield.cefparams)
        calc_cefparams!(lfield)
    end
    cef_magneticmoment_powder!(lfield.ion,lfield.cefparams,dfcalc;NJ,T,units,method,mode)
    return nothing
end