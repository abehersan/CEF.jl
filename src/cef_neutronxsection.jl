function dipolar_formfactor(ion::mag_ion, Q::Real)::Float64
    A_j0, a_j0, B_j0, b_j0, C_j0, c_j0, D_j0 = ion.ff_coeff_j0
    A_j2, a_j2, B_j2, b_j2, C_j2, c_j2, D_j2 = ion.ff_coeff_j2
    s = Q / 4pi
    ff_j0 = A_j0 * exp(-a_j0*s^2) + B_j0 * exp(-b_j0*s^2) + C_j0 * exp(-c_j0*s^2) + D_j0
    ff_j2 = A_j2*s^2 * exp(-a_j2*s^2) + B_j2*s^2 * exp(-b_j2*s^2) + C_j2*s^2 * exp(-c_j2*s^2) + D_j2*s^2
    return ff_j0 + ion.C2 * ff_j2
end


function calc_polmatrix(Qcart::Vector{Float64})::Matrix{Float64}
    polmat = Matrix{Float64}(undef, (3, 3))
    Q = normalize(Qcart)
    for I in CartesianIndices(polmat)
        a, b = Tuple(I)
        polmat[I] = (isequal(a, b)*1.0 - Q[a]*Q[b])
    end
    return polmat
end


function calc_transitions(ion::mag_ion, i::Int64, Vp::Matrix{ComplexF64})::Vector{VEC{9}}
    SIGMAS = VEC{9}[]
    @views @inbounds for j in 1:size(Vp, 1)
        sxx = real( dot( Vp[:, i], ion.Jx, Vp[:, j] ) * dot( Vp[:, j], ion.Jx, Vp[:, i] ) )
        sxy = real( dot( Vp[:, i], ion.Jy, Vp[:, j] ) * dot( Vp[:, j], ion.Jx, Vp[:, i] ) )
        sxz = real( dot( Vp[:, i], ion.Jz, Vp[:, j] ) * dot( Vp[:, j], ion.Jx, Vp[:, i] ) )
        syx = real( dot( Vp[:, i], ion.Jx, Vp[:, j] ) * dot( Vp[:, j], ion.Jy, Vp[:, i] ) )
        syy = real( dot( Vp[:, i], ion.Jy, Vp[:, j] ) * dot( Vp[:, j], ion.Jy, Vp[:, i] ) )
        syz = real( dot( Vp[:, i], ion.Jz, Vp[:, j] ) * dot( Vp[:, j], ion.Jy, Vp[:, i] ) )
        szx = real( dot( Vp[:, i], ion.Jx, Vp[:, j] ) * dot( Vp[:, j], ion.Jz, Vp[:, i] ) )
        szy = real( dot( Vp[:, i], ion.Jy, Vp[:, j] ) * dot( Vp[:, j], ion.Jz, Vp[:, i] ) )
        szz = real( dot( Vp[:, i], ion.Jz, Vp[:, j] ) * dot( Vp[:, j], ion.Jz, Vp[:, i] ) )
        push!(SIGMAS, [sxx, sxy, sxz, syx, syy, syz, szx, szy, szz]*ion.gj^2)
    end
    return SIGMAS
end


function calc_neutronspectrum_xtal(ion::mag_ion, Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, Qcart::Vector{<:Real}, T::Real)::Vector{VEC{2}}
    np = population_factor(Ep, T)
    ffactor = dipolar_formfactor(ion, norm(Qcart))
    polfactors = reshape(calc_polmatrix(Qcart), 9)
    NXS = VEC{2}[]
    @inbounds for i in eachindex(Ep)
        if isapprox(np[i], 0.0, atol=PREC)
            continue
        end
        sigmas = calc_transitions(ion, i, Vp)
        @inbounds for j in eachindex(Ep)
            dE = Ep[j] - Ep[i]
            NINT = dot(polfactors, sigmas[j])*np[i]*abs2(ffactor)
            if dE < 0.0     # detailed balance
                NINT *= exp( -abs(dE)/(kB*T) )
            end
            push!(NXS, [dE, NINT])
        end
    end
    return NXS
end


function calc_neutronspectrum_powd(ion::mag_ion, Ep::Vector{Float64}, Vp::Matrix{ComplexF64}, Qlen::Real, T::Real)::Vector{VEC{2}}
    np = population_factor(Ep, T)
    ffactor = dipolar_formfactor(ion, Qlen)
    NXS = VEC{2}[]
    @inbounds for i in eachindex(Ep)
        if isapprox(np[i], 0.0, atol=PREC)
            continue
        end
        sigmas = calc_transitions(ion, i, Vp)
        @inbounds for j in eachindex(Ep)
            dE = Ep[j] - Ep[i]
            NINT = sum(sigmas[j])*np[i]*abs2(ffactor)*2.0/3.0
            if dE < 0.0     # detailed balance
                NINT *= exp( -abs(dE)/(kB*T) )
            end
            push!(NXS, [dE, NINT])
        end
    end
    return NXS
end


function simulate_Escan(NXS::Vector{VEC{2}}, Es::AbstractVector, R::Function=TAS_resfunc)::Vector{Float64}
    Is = zeros(length(Es))
    @inbounds for i in eachindex(Es)
        @inbounds for j in eachindex(NXS)
            E, I = NXS[j]
            Is[i] += I * R(Es[i], E)
        end
    end
    return Is
end


function cef_neutronxsection_crystal!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; Q::Vector{<:Real}, T::Real=1.0, B::Vector{<:Real}=[0.0,0.0,0.0], resfunc::Function=TAS_resfunc, method::Symbol=:EO)::Nothing
    E, V = eigen(cef_hamiltonian(ion,cefparams,B=B,method=method))
    E .-= minimum(E)
    NINT = calc_neutronspectrum_xtal(ion,E,V,Q,T)
    EN = dfcalc.EN
    II = round.(simulate_Escan(NINT, EN, resfunc),digits=SDIG)
    dfcalc[!, :I_CALC] = II
    return nothing
end


function cef_neutronxsection_powder!(ion::mag_ion, cefparams::DataFrame, dfcalc::DataFrame; Q::Real, T::Real=1.0, resfunc::Function=TAS_resfunc, method::Symbol=:EO)::Nothing
    E, V = eigen(cef_hamiltonian(ion,cefparams;method=method))
    E .-= minimum(E)
    NINT = calc_neutronspectrum_powd(ion,E,V,Q,T)
    EN = dfcalc.EN
    II = round.(simulate_Escan(NINT, EN, resfunc),digits=SDIG)
    dfcalc[!, :I_CALC] = II
    return nothing
end


function TAS_resfunc(E::Float64, Epeak::Float64, width::Function=x->0.2/(2*sqrt(2*log(2))))::Float64
    return gaussian(x=E, mu=Epeak, A=1.0, sigma=width(E))
end


function gaussian(; x::Real, A::Real, mu::Real, sigma::Real)::Float64
    return (A/(sqrt(2pi)*sigma))*exp(-(x-mu)^2/(2*sigma^2))
end


function lorentz(; x::Real, A::Real, mu::Real, gamma::Real)::Float64
    return (A/pi)*(gamma/((x-mu)^2+gamma^2))
end