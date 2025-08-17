function calc_chi0tensor(spinops,esample,E,V;T=2.0,maxt=-1,epsilon=0.001)
    np=population_factor(E,T)
    if maxt > 0
        E=E[1:maxt]
        V=V[:,1:maxt]
    end
    chi0s=Matrix{ComplexF64}[]
    @views @inbounds for i in eachindex(esample)
        chi0=zeros(ComplexF64,(3,3))
        @views @inbounds for (p,ep) in enumerate(E), (pp,epp) in enumerate(E)
            if isapprox(ep,epp,atol=epsilon)
                continue
            elseif iszero(np[p]) & iszero(np[pp])
                continue
            else
                mxx=(dot(V[:,p],spinops[1],V[:,pp])*dot(V[:,pp],spinops[1],V[:,p])*(np[p]-np[pp]))/(epp-ep-esample[i]-1im*epsilon)
                mxy=(dot(V[:,p],spinops[1],V[:,pp])*dot(V[:,pp],spinops[2],V[:,p])*(np[p]-np[pp]))/(epp-ep-esample[i]-1im*epsilon)
                mxz=(dot(V[:,p],spinops[1],V[:,pp])*dot(V[:,pp],spinops[3],V[:,p])*(np[p]-np[pp]))/(epp-ep-esample[i]-1im*epsilon)
                myx=(dot(V[:,p],spinops[2],V[:,pp])*dot(V[:,pp],spinops[1],V[:,p])*(np[p]-np[pp]))/(epp-ep-esample[i]-1im*epsilon)
                myy=(dot(V[:,p],spinops[2],V[:,pp])*dot(V[:,pp],spinops[2],V[:,p])*(np[p]-np[pp]))/(epp-ep-esample[i]-1im*epsilon)
                myz=(dot(V[:,p],spinops[2],V[:,pp])*dot(V[:,pp],spinops[3],V[:,p])*(np[p]-np[pp]))/(epp-ep-esample[i]-1im*epsilon)
                mzx=(dot(V[:,p],spinops[3],V[:,pp])*dot(V[:,pp],spinops[1],V[:,p])*(np[p]-np[pp]))/(epp-ep-esample[i]-1im*epsilon)
                mzy=(dot(V[:,p],spinops[3],V[:,pp])*dot(V[:,pp],spinops[2],V[:,p])*(np[p]-np[pp]))/(epp-ep-esample[i]-1im*epsilon)
                mzz=(dot(V[:,p],spinops[3],V[:,pp])*dot(V[:,pp],spinops[3],V[:,p])*(np[p]-np[pp]))/(epp-ep-esample[i]-1im*epsilon)
                chi0[1,1]+=mxx
                chi0[1,2]+=mxy
                chi0[1,3]+=mxz
                chi0[2,1]+=myx
                chi0[2,2]+=myy
                chi0[2,3]+=myz
                chi0[3,1]+=mzx
                chi0[3,2]+=mzy
                chi0[3,3]+=mzz
            end
        end
        push!(chi0s,chi0)
    end
    return chi0s
end


function calc_rpapeaks(jq,chi0s)
    rpas=Matrix{Float64}[]
    @views @inbounds for i in eachindex(chi0s)
        push!(rpas,imag( inv( diagm(ones(ComplexF64,3)) .- jq .*chi0s[i] ) .*chi0s[i] ) )
    end
    return rpas
end


function rpa_unpolarized_neutronxsection(esample,rpas;T=2.0,ion,polmat,qcart)
    @assert length(esample)==length(rpas)
    IQE=zeros(length(esample))
    @views @inbounds for i in eachindex(IQE)
        neutronrpa=sum(polmat .*rpas[i])*( 1.0/(1.0-exp(-(esample[i])/(kB*T))) )
        if esample[i]<0.0
            neutronrpa*=exp( -abs(esample[i])/(kB*T) )
        end
        if isnan(neutronrpa) || isinf(neutronrpa)
            neutronrpa=0.0
        end
        IQE[i]+=neutronrpa
    end
    IQE*=CEF.CC*(1.0/pi)*(ion.gj*dipolar_formfactor(ion,sqrt(dot(qcart,qcart))))^2
    return IQE
end

function convgauss(E,epeaks,peakamps;sigmas::Function)
    IE=zeros(length(E))
    for i in eachindex(E)
        for j in eachindex(epeaks)
            IE[i]+=gauss(x=E[i],mu=epeaks[j],A=1.0,sigma=sigmas(E[i]))*peakamps[j]
        end
    end
    return IE
end

function convlorentz(E,epeaks,peakamps;gammas::Function)
    IE=zeros(length(E))
    for i in eachindex(E)
        for j in eachindex(epeaks)
            IE[i]+=lorentz(x=E[i],mu=epeaks[j],A=1.0,gamma=gammas(E[i]))*peakamps[j]
        end
    end
    return IE
end