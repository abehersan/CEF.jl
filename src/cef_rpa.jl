function calc_chi0(spinops::Vector{Matrix{ComplexF64}},hw::Float64,E::Vector{Float64},V::Matrix{ComplexF64},T::Float64;epsilon::Real=PREC,mode::Function=real)::CMAT{3}
    n=population_factor(E,T)
    chi0=CMAT{3}(zeros(3,3))
    for I in CartesianIndices(chi0)
        opa=spinops[I[1]]
        opb=spinops[I[2]]
        chi0val=ComplexF64(0.0)
        for p in eachindex(E)
            ep=E[p]
            for pp in eachindex(E)
                epp=E[pp]
                if isapprox(ep,epp,atol=PREC)
                    if isapprox(n[p],0.0,atol=PREC)
                        continue
                    else
                        chi0val+=(1im*epsilon/(hw+1im*epsilon))*
                            dot(V[:,p],opa,V[:,pp])*dot(V[:,pp],opb,V[:,p])*n[p]/(kB*T)-
                            thermal_average(;Ep=E,Vp=V,op=opa,T=T,mode=mode)*
                            thermal_average(;Ep=E,Vp=V,op=opa,T=T,mode=mode)*n[p]/(kB*T)
                    end
                else
                    if isapprox(n[p],0.0,atol=PREC) && isapprox(n[pp],0.0,atol=PREC)
                        continue
                    else
                        chi0val+=dot(V[:,p],opb,V[:,pp])*dot(V[:,pp],opa,V[:,p])*(n[p]-n[pp])/
                            (epp-ep-hw-1im*epsilon)
                    end
                end
            end
        end
        chi0[I]=chi0val
    end
    return chi0
end