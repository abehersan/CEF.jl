function chi_units(units::Symbol)::Float64
    if isequal(units, :SI)
        return 4.062426*1e-7    # N_A * muB(erg/G) * muB(meV/G)
    elseif isequal(units, :CGS)
        return 0.03232776       # N_A * muB(J/T) * muB(meV/T) * mu0
    elseif isequal(units, :ATOMIC)
        return 1.0              # muB / magnetic ion
    else
        @error "Units $units not understood. Use one of either :SI, :CGS or :ATOMIC"
    end
end


function calc_chialphaalpha(; op_alpha::Matrix{ComplexF64}, Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real, mode::Function=real)::Float64
    chi_alphaalpha::Float64 = 0.0
    np = population_factor(Ep, T)
    @views @inbounds for p in eachindex(Ep)
        ep = Ep[p]
        @views @inbounds for pp in eachindex(Ep)
            epp = Ep[pp]
            if isapprox(ep, epp; atol=PREC)
                if isapprox(np[p], 0.0, atol=PREC)
                    continue
                else
                    m_element_alpha = dot(Vp[:,p], op_alpha, Vp[:,pp])
                    m_element = m_element_alpha*conj(m_element_alpha)
                    chi_alphaalpha += m_element*np[p]/(kB*T)
                end
            else
                if isapprox(np[p], 0.0, atol=PREC) && isapprox(np[pp], 0.0, atol=PREC)
                    continue
                else
                    m_element_alpha = dot(Vp[:,p], op_alpha, Vp[:,pp])
                    m_element = m_element_alpha*conj(m_element_alpha)
                    pop_diff = np[p] - np[pp]
                    chi_alphaalpha += (m_element*pop_diff)/(epp-ep)
                end
            end
        end
    end
    t_avg_alpha = thermal_average(Ep=Ep,Vp=Vp,op=op_alpha,T=T,mode=mode)
    chi_alphaalpha -= (t_avg_alpha^2)/(kB*T)
    return chi_alphaalpha
end


function calc_chi(g::MVEC{3}; spin_proj::Vector{Matrix{ComplexF64}}, Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, T::Real, mode::Function)::Float64
    chi_vec = zeros(Float64, 3)
    @inbounds for i in eachindex(spin_proj)
        if iszero(spin_proj[i])
            continue
        else
            chi_vec[i]=calc_chialphaalpha(op_alpha=spin_proj[i],Ep=Ep,Vp=Vp,T=T,mode=mode)
        end
    end
    return norm((g .^2) .* chi_vec)
end


function cef_susceptibility_crystal!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; units::Symbol=:CGS, method::Symbol=:EO, mode::Function=real)::Nothing
    unit_factor = chi_units(units)
    spin_ops = [ion.Jx,ion.Jy,ion.Jz]
    extfield = [mean(dfcalc.Bx), mean(dfcalc.By), mean(dfcalc.Bz)]
    E, V = eigen(cef_hamiltonian(ion,cefparams;B=extfield,method=method))
    E .-= minimum(E)
    @eachrow! dfcalc begin
        @newcol :CHIx_CALC::Vector{Float64}
        @newcol :CHIy_CALC::Vector{Float64}
        @newcol :CHIz_CALC::Vector{Float64}
        :CHIx_CALC=round(calc_chi(ion.g,spin_proj=spin_ops.*[1.0,0.0,0.0],Ep=E,Vp=V,T=:T,mode=mode)*unit_factor,digits=SDIG)
        :CHIy_CALC=round(calc_chi(ion.g,spin_proj=spin_ops.*[0.0,1.0,0.0],Ep=E,Vp=V,T=:T,mode=mode)*unit_factor,digits=SDIG)
        :CHIz_CALC=round(calc_chi(ion.g,spin_proj=spin_ops.*[0.0,0.0,1.0],Ep=E,Vp=V,T=:T,mode=mode)*unit_factor,digits=SDIG)
    end
    return nothing
end


function cef_susceptibility_powder!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; units::Symbol=:CGS, method::Symbol=:EO, mode::Function=real)::Nothing
    unit_factor = chi_units(units)
    spin_ops = [ion.Jx,ion.Jy,ion.Jz]
    E, V = eigen(cef_hamiltonian(ion,bfactors;method=method))
    E .-= minimum(E)
    @eachrow! dfcalc begin
        @newcol :CHI_CALC::Vector{Float64}
        chix = calc_chi(ion.g,spin_proj=spin_ops.*[1.0,0.0,0.0],Ep=E,Vp=V,T=:T,mode=mode)*unit_factor
        chiy = calc_chi(ion.g,spin_proj=spin_ops.*[0.0,1.0,0.0],Ep=E,Vp=V,T=:T,mode=mode)*unit_factor
        chiz = calc_chi(ion.g,spin_proj=spin_ops.*[0.0,0.0,1.0],Ep=E,Vp=V,T=:T,mode=mode)*unit_factor
        :CHI_CALC=round(((chix+chiy+chiz)/3),digits=SDIG)
    end
    return nothing
end