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


function calc_magmom(g::MVEC{3}, spinops::Vector{Matrix{ComplexF64}}, Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real, mode::Function)::Float64
    spin_expval = zeros(Float64, 3)
    @inbounds for i in eachindex(spin_expval)
        if iszero(spinops[i])
            continue
        end
        spin_expval[i] = thermal_average(Ep=Ep, Vp=Vp, op=spinops[i], T=T, mode=mode)
    end
    return g[1]*spin_expval[1]+g[2]*spin_expval[2]+g[3]*spin_expval[3]
end


function cef_magneticmoment_crystal!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; units::Symbol=:ATOMIC, method::Symbol=:EO, mode::Function=real)
    unit_factor = mag_units(units)
    spinops = [ion.Jx,ion.Jy,ion.Jz]
    @eachrow! dfcalc begin
        @newcol :Mp_CALC::Vector{Float64}
        @newcol :Mx_CALC::Vector{Float64}
        @newcol :My_CALC::Vector{Float64}
        @newcol :Mz_CALC::Vector{Float64}
        extfield = [:Bx, :By, :Bz]
        E, V = eigen(cef_hamiltonian(ion,cefparams;B=extfield,method=method))
        E .-= minimum(E)
        if iszero(extfield)
            spin_proj = spinops
        else
            spin_proj = spinops .* normalize(extfield)
        end
        spin_expval = zeros(Float64, 3)
        @inbounds for i in eachindex(spin_expval)
            if iszero(spinops[i])
                continue
            end
            spin_expval[i] = thermal_average(Ep=E,Vp=V,op=spinops[i],T=:T,mode=mode)
        end
        :Mp_CALC=round(norm(ion.g .* spin_expval)*unit_factor,digits=SDIG)
        :Mx_CALC=round(thermal_average(Ep=E,Vp=V,op=spinops[1],T=:T,mode=mode)*ion.g[1]*unit_factor,digits=SDIG)
        :My_CALC=round(thermal_average(Ep=E,Vp=V,op=spinops[2],T=:T,mode=mode)*ion.g[2]*unit_factor,digits=SDIG)
        :Mz_CALC=round(thermal_average(Ep=E,Vp=V,op=spinops[3],T=:T,mode=mode)*ion.g[3]*unit_factor,digits=SDIG)
    end
    return nothing
end


function cef_magneticmoment_powder!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; units::Symbol=:ATOMIC, method::Symbol=:EO, mode::Function=real)
    unit_factor = mag_units(units)
    spinops = [ion.Jx,ion.Jy,ion.Jz]
    @eachrow! dfcalc begin
        @newcol :M_CALC::Vector{Float64}
        E, V = eigen(cef_hamiltonian(ion,cefparams; B=[:B,0.0,0.0],method=method))
        E .-= minimum(E)
        MX = calc_magmom(ion.g,spinops .* [1.0, 0.0, 0.0],E,V,:T,mode)

        E, V = eigen(cef_hamiltonian(ion,cefparams; B=[0.0,:B,0.0],method=method))
        E .-= minimum(E)
        MY = calc_magmom(ion.g,spinops .* [0.0, 1.0, 0.0],E,V,:T,mode)

        E, V = eigen(cef_hamiltonian(ion,cefparams; B=[0.0,0.0,:B],method=method))
        E .-= minimum(E)
        MZ = calc_magmom(ion.g,spinops .* [0.0, 0.0, 1.0],E,V,:T,mode)

        :M_CALC=round(((MX + MY + MZ) / 3.0) * unit_factor,digits=SDIG)
    end
    return nothing
end