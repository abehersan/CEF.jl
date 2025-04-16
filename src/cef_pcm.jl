Base.@kwdef mutable struct local_env
    ion::mag_ion
    lparams::VEC{6}
    lvecs::MAT3
    cartesian_pointcs::Vector{VEC{4}}
    spherical_pointcs::Vector{VEC{4}}
    cefparams::DataFrame
end


function Base.show(io::IO, ::MIME"text/plain", pcenv::local_env)
    display(pcenv.ion)
    println()
    printstyled(io,"Lattice parameters [a, b, c, α, β, γ]: $(pcenv.lparams)\n")
    printstyled(io,"Lattice vector a: $(pcenv.lvecs[:,1])\n")
    printstyled(io,"Lattice vector b: $(pcenv.lvecs[:,2])\n")
    printstyled(io,"Lattice vector c: $(pcenv.lvecs[:,3])\n")
    println()
    printstyled(io, "Cartesian ligand field (Å, Ze):\n\tx\t\ty\t\tz\t\tZ\n")
    for pc in pcenv.cartesian_pointcs
        x,y,z,Z=pc
        LBL=@sprintf("\t%+3.5f\t%+3.5f\t%+3.5f\t%+3.5f\n",x,y,z,Z)
        printstyled(io,LBL)
    end
    println()
    printstyled(io, "Spherical ligand field (Å, deg., Ze):\n\tr\t\tθ\t\tϕ\t\tZ\n")
    for pc in pcenv.spherical_pointcs
        rr,th,ph,Z=pc
        LBL=@sprintf("\t%+3.5f\t%+3.5f\t%+3.5f\t%+3.5f\n",rr,th*180/pi,ph*180/pi,Z)
        printstyled(io,LBL)
    end
    return nothing
end


function ligand_field(ion,lparams,pointcharges)::local_env
    return local_env(
        ion=ion,
        lparams=VEC{6}(lparams),
        lvecs=lattice_vectors(lparams...),
        cartesian_pointcs=[VEC{4}(pc) for pc in pointcharges],
        spherical_pointcs=[[to_spherical(pc[1:3])...,pc[4]] for pc in pointcharges],
        cefparams=DataFrame()
    )
end


function lattice_vectors(a,b,c,alpha,beta,gamma)::MAT3
    @assert all(0 < x < 180 for x in (alpha, beta, gamma))
    sgamma=sind(gamma)
    cgamma=cosd(gamma)
    cbeta=cosd(beta)
    calpha=cosd(alpha)
    v1=[a,0.0,0.0]
    v2=[b*cgamma,b*sgamma,0]
    v3x=c*cbeta
    v3y=(c/sgamma)*(calpha-cbeta*cgamma)
    v3z=(c/sgamma)*sqrt(sgamma^2-calpha^2-cbeta^2+2*calpha*cbeta*cgamma)
    v3=[v3x,v3y,v3z]
    latvecs = MAT3(hcat(v1, v2, v3))
    return latvecs
end


function to_spherical(cartvec)::VEC3
    x,y,z=cartvec
    r=sqrt(x^2+y^2+z^2)
    theta=if isapprox(r,0,atol=PREC)
        0.0
    else
        acos(clamp(z/r,-1,1))
    end
    phi=atan(y,x)
    return VEC3(r,theta,phi)
end


function tesseral_harmonics(l::Int64,m::Int64,x::Real,y::Real,z::Real,r::Real)::Real
    if isequal(l,2)
        if isequal(m,-2)
            Z2m2=sqrt(15/pi)*(1/4)*((2*x*y)/r^2)
            return Z2m2
        elseif isequal(m,-1)
            Z2m1=sqrt(15/pi)*(1/2)*((y*z)/r^2)
            return Z2m1
        elseif isequal(m,0)
            Z20=sqrt(5/pi)*(1/4)*((3*z^2-r^2)/r^2)
            return Z20
        elseif isequal(m,+1)
            Z21=sqrt(15/pi)*(1/2)*((x*z)/r^2)
            return Z21
        elseif isequal(m,+2)
            Z22=sqrt(15/pi)*(1/4)*((x^2-y^2)/r^2)
            return Z22
        end
    end
    if isequal(l,4)
        if isequal(m,-4)
            Z4m4=sqrt(35/pi)*(3/16)*(4*(x^3*y-x*y^3)/r^4)
            return Z4m4
        elseif isequal(m,-3)
            Z4m3=sqrt(70/pi)*(3/8)*(z*(3*x^2*y-y^3)/r^4)
            return Z4m3
        elseif isequal(m,-2)
            Z4m2=sqrt(5/pi)*(3/8)*(2*x*y*(7*z^2-r^2)/r^4)
            return Z4m2
        elseif isequal(m,-1)
            Z4m1=sqrt(5/(2*pi))*(3/4)*(y*z*(7*z^2-3*r^2)/r^4)
            return Z4m1
        elseif isequal(m,0)
            Z40=sqrt(1/pi)*(3/16)*((35*z^4-30*z^2*r^2+3*r^4)/r^4)
            return Z40
        elseif isequal(m,1)
            Z41=sqrt(5/(2*pi))*(3/4)*(x*z*(7*z^2-3*r^2)/r^4)
            return Z41
        elseif isequal(m,2)
            Z42=sqrt(5/pi)*(3/8)*((x^2-y^2)*(7*z^2-r^2)/r^4)
            return Z42
        elseif isequal(m,3)
            Z43=sqrt(70/pi)*(3/8)*(z*(x^3-3*x*y^2)/r^4)
            return Z43
        elseif isequal(m,4)
            Z44=sqrt(35/pi)*(3/16)*((x^4-6*x^2*y^2+y^4)/r^4)
            return Z44
        end
    end
    if isequal(l,6)
        if isequal(m,-6)
            Z6m6=sqrt(26/(231*pi))*(231/64)*((6*x^5*y-20*x^3*y^3+6*x*y^5)/r^6)
            return Z6m6
        elseif isequal(m,-5)
            Z6m5=sqrt(9009/(512*pi))*(z*(5*x^4*y-10*x^2*y^3+y^5)/r^6)
            return Z6m5
        elseif isequal(m,-4)
            Z6m4=sqrt(13/(7*pi))*(21/32)*(4*(x^3*y-x*y^3)*(11*z^2-r^2)/r^6)
            return Z6m4
        elseif isequal(m,-3)
            Z6m3=sqrt(2730/pi)*(1/32)*((3*x^2*y-y^3)*(11*z^3-3*z*r^2)/r^6)
            return Z6m3
        elseif isequal(m,-2)
            Z6m2=sqrt(2730/pi)*(1/64)*(2*x*y*(33*z^4-13*z^2*r^2+5*r^4)/r^6)
            return Z6m2
        elseif isequal(m,-1)
            Z6m1=sqrt(273/(4*pi))*(1/8)*(x*z*(33*z^4-30*z^2*r^2+5*r^4)/r^6)
            return Z6m1
        elseif isequal(m,0)
            Z60=sqrt(13/pi)*(1/32)*((231*z^6-315*z^4*r^2+105*z^2*r^4-5*r^6)/r^6)
            return Z60
        elseif isequal(m,1)
            Z61=sqrt(273/(4*pi))*(1/8)*(x*z*(33*z^4-30*z^2*r^2+5*r^4)/r^6)
            return Z61
        elseif isequal(m,2)
            Z62=sqrt(2730/pi)*(1/64)*((x^2-y^2)*(33*z^4-18*z^2*r^2+r^4)/r^6)
            return Z62
        elseif isequal(m,3)
            Z63=sqrt(2730/pi)*(1/32)*((x^3-3*x*y^2)*(11*z^3-3*z*r^2)/r^6)
            return Z63
        elseif isequal(m,4)
            Z64=sqrt(13/(7*pi))*(21/32)*((x^4-6*x^2*y^2+y^4)*(11*z^2-r^2)/r^6)
            return Z64
        elseif isequal(m,5)
            Z65=sqrt(9009/(512*pi))*(z*(x^5-10*x^3*y^2+5*x*y^4)/r^6)
            return Z65
        elseif isequal(m,6)
            Z66=sqrt(26/(231*pi))*(231/4)*((x^6-15*x^4*y^2+15*x^2*y^4-y^6)/r^6)
            return Z66
        end
    end
    @error("Values of l=$(l) and/or m=$(m) invalid. l must be one of [2,4,6] and m takes values between -l and l.")
end


function calc_cefparams!(lfield::local_env;shielded::Bool=true)
    cefparams=DataFrame(B=Float64[],l=Int[],m=Int[])
    if shielded
        radwav=lfield.ion.rad_wavefunction_shielded
    else
        radwav=lfield.ion.rad_wavefunction_unshielded
    end
    sfactors=lfield.ion.stevens_factors
    for l in [2,4,6]
        unit_factor=1.44e4*(0.5292^l)
        if isequal(l,2)
            rl=radwav[1]
            al=sfactors[1]
        elseif isequal(l,4)
            rl=radwav[2]
            al=sfactors[2]
        elseif isequal(l,6)
            rl=radwav[3]
            al=sfactors[3]
        end
        for m in -l:1:l
            Alm=0.0
            for pc in lfield.cartesian_pointcs
                x,y,z,Z=pc
                rxtal=[x,y,z]
                rcart=lfield.lvecs*rxtal
                x,y,z=rcart
                R=sqrt(x^2+y^2+z^2)
                Zlm=tesseral_harmonics(l,m,x,y,z,R)
                Alm+=(Z*Zlm)/(R^(l+1))
            end
            Alm*=(4pi)/(2*l+1)
            if shielded
                if isequal(l,2)
                    sig=lfield.ion.shielding_factors[1]
                elseif isequal(l,4)
                    sig=lfield.ion.shielding_factors[2]
                elseif isequal(l,6)
                    sig=lfield.ion.shielding_factors[3]
                end
                al=(1-sig)*al
            end
            Blm=Alm*rl*al*unit_factor
            if Blm < CUTOFF
                continue
            end
            append!(cefparams,DataFrame(:B=>Blm,:l=>l,:m=>m))
        end
    end
    lfield.cefparams=cefparams
    return nothing
end