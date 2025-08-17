Base.@kwdef mutable struct local_env
    ion::mag_ion
    lparams::VEC{6}
    dlattvecs::MAT3
    rlattvecs::MAT3
    cartesian_pointcs::Vector{VEC{4}}
    spherical_pointcs::Vector{VEC{4}}
    cefparams::DataFrame
end


function Base.show(io::IO, ::MIME"text/plain", pcenv::local_env)
    display(pcenv.ion)
    println()
    printstyled(io,"Lattice parameters [a, b, c, α, β, γ]: $(pcenv.lparams)\n")
    println()
    printstyled(io,"Lattice vector a (Å): \t\t$(pcenv.dlattvecs[:,1])\n")
    printstyled(io,"Lattice vector b (Å): \t\t$(pcenv.dlattvecs[:,2])\n")
    printstyled(io,"Lattice vector c (Å): \t\t$(pcenv.dlattvecs[:,3])\n")
    printstyled(io,"Lattice vector a* (1/Å): \t$(pcenv.rlattvecs[:,1])\n")
    printstyled(io,"Lattice vector b* (1/Å): \t$(pcenv.rlattvecs[:,2])\n")
    printstyled(io,"Lattice vector c* (1/Å): \t$(pcenv.rlattvecs[:,3])\n")
    println()
    println("Point charge coordinates")
    printstyled(io, "Cartesian (Å, Ze):\n\tx\t\ty\t\tz\t\tZ\n")
    for pc in pcenv.cartesian_pointcs
        x,y,z,Z=pc
        LBL=@sprintf("\t%+3.5f\t%+3.5f\t%+3.5f\t%+3.5f\n",x,y,z,Z)
        printstyled(io,LBL)
    end
    println()
    printstyled(io, "Spherical (Å, deg., Ze):\n\tr\t\tθ\t\tϕ\t\tZ\n")
    for pc in pcenv.spherical_pointcs
        rr,th,ph,Z=pc
        LBL=@sprintf("\t%+3.5f\t%+3.5f\t%+3.5f\t%+3.5f\n",rr,th*180/pi,ph*180/pi,Z)
        printstyled(io,LBL)
    end
    return nothing
end


function make_pcm(ion,lparams,pointcharges;coords::Symbol=:cartesian)::local_env
    dlattvecs,rlattvecs=lattice_vectors(lparams...)
    cartesian_pointcs=VEC{4}[]
    spherical_pointcs=VEC{4}[]
    if isequal(coords,:cartesian)
        for pc in pointcharges
            pccart=dlattvecs*pc[1:3]
            pcsphe=to_spherical(pccart)
            push!(cartesian_pointcs,VEC{4}(pccart...,pc[end]))
            push!(spherical_pointcs,VEC{4}(pcsphe...,pc[end]))
        end
    elseif isequal(coords,:spherical)
        for pc in pointcharges
            pcsphe=pc[1:3]
            pccart=to_cartesian(pcsphe)
            push!(cartesian_pointcs,VEC{4}(pccart...,pc[end]))
            push!(spherical_pointcs,VEC{4}(pcsphe...,pc[end]))
        end
    else
        @error "Coordinates $units not understood. Use one of either :spherical or :cartesian"
    end
    lenv=local_env(
        ion=ion,
        lparams=VEC{6}(lparams),
        dlattvecs=dlattvecs,
        rlattvecs=rlattvecs,
        cartesian_pointcs=cartesian_pointcs,
        spherical_pointcs=spherical_pointcs,
        cefparams=DataFrame()
    )
    return lenv
end


function lattice_vectors(a,b,c,alpha,beta,gamma)
    @assert all(0 < x < 180 for x in (alpha, beta, gamma))
    sgamma=sind(gamma)
    cgamma=cosd(gamma)
    cbeta=cosd(beta)
    calpha=cosd(alpha)
    a1=[a,0.0,0.0]
    a2=[b*cgamma,b*sgamma,0]
    a3x=c*cbeta
    a3y=(c/sgamma)*(calpha-cbeta*cgamma)
    a3z=(c/sgamma)*sqrt(sgamma^2-calpha^2-cbeta^2+2*calpha*cbeta*cgamma)
    a3=[a3x,a3y,a3z]
    dlattvecs=MAT3(hcat(a1,a2,a3))

    VV=dot(a1,cross(a2,a3))
    b1=(2pi/VV)*cross(a2,a3)
    b2=(2pi/VV)*cross(a3,a1)
    b3=(2pi/VV)*cross(a1,a2)
    rlattvecs=MAT3(hcat(b1,b2,b3))

    return (dlattvecs,rlattvecs)
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


function to_cartesian(sphvec)::VEC3
    r,th,ph=sphvec
    x=r*sind(th)*cosd(ph)
    y=r*sind(th)*sind(ph)
    z=r*cosd(th)
    return VEC3(x,y,z)
end


function plot_pcm(pcm::local_env; path="./pcm.gp")
    open(path, "w") do io
        println(io, "set terminal wxt size 720,720 enhanced title 'CEF.jl point charges'")
        println(io, "set xlabel 'x (Å)'")
        println(io, "set ylabel 'y (Å)'")
        println(io, "set zlabel 'z (Å)'")
        println(io, "set view equal xyz")
        println(io, "splot \\")
        println(io, "    '-' using 1:2:3 with points pointtype 7 pointsize 1.5 lc rgb 'purple' title 'RE ion', \\")
        println(io, "    '-' using 1:2:3 with points pointtype 7 pointsize 1.5 lc rgb 'orange' title 'PCs', \\")
        println(io, "    '-' using 1:2:3:4:5:6 with vectors head filled lc rgb 'red' lw 2 title 'a', \\")
        println(io, "    '-' using 1:2:3:4:5:6 with vectors head filled lc rgb 'green' lw 2 title 'b', \\")
        println(io, "    '-' using 1:2:3:4:5:6 with vectors head filled lc rgb 'blue' lw 2 title 'c'")

        # central ion data block
        println(io, "0.0 0.0 0.0")
        println(io, "e")

        # charges data block
        for pc in pcm.cartesian_pointcs
            println(io, join(pc[1:3], " "))
        end
        println(io, "e")

        # basis vectors data block (from origin)
        aa=pcm.dlattvecs[:,1]
        bb=pcm.dlattvecs[:,2]
        cc=pcm.dlattvecs[:,3]
        println(io, "0 0 0 ", join(aa, " ")*"\ne")
        println(io, "0 0 0 ", join(bb, " ")*"\ne")
        println(io, "0 0 0 ", join(cc, " ")*"\ne")

        println(io, "pause -1")
    end
    println("Gnuplot .gp file generated: $(path)")
    return nothing
end


function tesseral_harmonics(l::Int64,m::Int64,x::Real,y::Real,z::Real,r::Real)::Real
    if isequal(l,0)
        T00=+sqrt(1/(4*pi))
        return T00
    end
    if isequal(l,1)
        if isequal(m,-1)
            T1m1=+sqrt(3/(4*pi))*(y/r)
            return T1m1
        elseif isequal(m,0)
            T10=+sqrt(3/(4*pi))*(z/r)
            return T10
        elseif isequal(m,+1)
            T11=-sqrt(3/(4*pi))*(x/r)
            return T11
        end
    end
    if isequal(l,2)
        if isequal(m,-2)
            T2m2=+sqrt(15/pi)*(1/4)*((2*x*y)/r^2)
            return T2m2
        elseif isequal(m,-1)
            T2m1=+sqrt(15/pi)*(1/2)*((y*z)/r^2)
            return T2m1
        elseif isequal(m,0)
            T20=+sqrt(5/pi)*(1/4)*((3*z^2-r^2)/r^2)
            return T20
        elseif isequal(m,+1)
            T21=-sqrt(15/pi)*(1/2)*((x*z)/r^2)
            return T21
        elseif isequal(m,+2)
            T22=+sqrt(15/pi)*(1/4)*((x^2-y^2)/r^2)
            return T22
        end
    end
    if isequal(l,3)
        if isequal(m,-3)
            T3m3=+sqrt(35/(32*pi))*((3*x^2*y-y^3)/r^3)
            return T3m3
        elseif isequal(m,-2)
            T3m2=+sqrt(105/(16*pi))*((2*x*y*z)/r^3)
            return T3m2
        elseif isequal(m,-1)
            T3m1=+sqrt(21/(32*pi))*(y*(5*z^2-r^2)/r^3)
            return T3m1
        elseif isequal(m,0)
            T30=+sqrt(7/(16*pi))*(z*(5*z^2-3*r^2)/r^3)
            return T30
        elseif isequal(m,+1)
            T31=-sqrt(21/(32*pi))*(x*(5*z^2-r^2)/r^3)
            return T31
        elseif isequal(m,+2)
            T32=+sqrt(105/(16*pi))*(z*(x^2-y^2)/r^3)
            return T32
        elseif isequal(m,+3)
            T33=-sqrt(35/(32*pi))*((x^3-3*x*y^2)/r^3)
            return T33
        end
    end
    if isequal(l,4)
        if isequal(m,-4)
            T4m4=+sqrt(35/pi)*(3/16)*(4*(x^3*y-x*y^3)/r^4)
            return T4m4
        elseif isequal(m,-3)
            T4m3=+sqrt(70/pi)*(3/8)*(z*(3*x^2*y-y^3)/r^4)
            return T4m3
        elseif isequal(m,-2)
            T4m2=+sqrt(5/pi)*(3/8)*(2*x*y*(7*z^2-r^2)/r^4)
            return T4m2
        elseif isequal(m,-1)
            T4m1=+sqrt(5/(2*pi))*(3/4)*(y*z*(7*z^2-3*r^2)/r^4)
            return T4m1
        elseif isequal(m,0)
            T40=+sqrt(1/pi)*(3/16)*((35*z^4-30*z^2*r^2+3*r^4)/r^4)
            return T40
        elseif isequal(m,1)
            T41=+sqrt(5/(2*pi))*(3/4)*(x*z*(7*z^2-3*r^2)/r^4)
            return T41
        elseif isequal(m,2)
            T42=+sqrt(5/pi)*(3/8)*((x^2-y^2)*(7*z^2-r^2)/r^4)
            return T42
        elseif isequal(m,3)
            T43=+sqrt(70/pi)*(3/8)*(z*(x^3-3*x*y^2)/r^4)
            return T43
        elseif isequal(m,4)
            T44=+sqrt(35/pi)*(3/16)*((x^4-6*x^2*y^2+y^4)/r^4)
            return T44
        end
    end
    if isequal(l,5)
        if isequal(m,-5)
            T5m5=sqrt(693/(512*pi))*((5*x^4*y-10*x^2*y^3+y^5)/r^5)
            return T5m5
        elseif isequal(m,-4)
            T5m4=sqrt(3465/(256*pi))*(4*z*(x^3*y-x*y^3)/r^5)
            return T5m4
        elseif isequal(m,-3)
            T5m3=sqrt(385/(512*pi))*((3*x^2*y-y^3)*(9*z^2-r^2)/r^5)
            return T5m3
        elseif isequal(m,-2)
            T5m2=sqrt(1155/(64*pi))*(2*x*y*(3*z^3-z*r^2)/r^5)
            return T5m2
        elseif isequal(m,-1)
            T5m1=sqrt(165/(256*pi))*(y*(21*z^4-14*z^2*r^2+r^4)/r^5)
            return T5m1
        elseif isequal(m,0)
            T50=sqrt(11/(256*pi))*((63*z^5-70*z^3*r^2+15*z*r^4)/r^5)
            return T50
        elseif isequal(m,+1)
            T51=-sqrt(165/(256*pi))*(x*(21*z^4-14*z^2*r^2+r^4)/r^5)
            return T51
        elseif isequal(m,+2)
            T52=sqrt(1155/(64*pi))*((x^2-y^2)*(3*z^3-z*r^2)/r^5)
            return T52
        elseif isequal(m,+3)
            T53=-sqrt(385/(512*pi))*((x^3-3*x*y^2)*(9*z^2-r^2)/r^5)
            return T53
        elseif isequal(m,+4)
            T54=sqrt(3465/(256*pi))*((x^4-6*x^2*y^2+y^4)/r^5)
            return T54
        elseif isequal(m,+5)
            T55=-sqrt(693/(512*pi))*((x^5-10*x^3*y^2+5*x*y^4)/r^5)
            return T55
        end
    end
    if isequal(l,6)
        if isequal(m,-6)
            T6m6=sqrt(26/(231*pi))*(231/64)*((6*x^5*y-20*x^3*y^3+6*x*y^5)/r^6)
            return T6m6
        elseif isequal(m,-5)
            T6m5=sqrt(9009/(512*pi))*(z*(5*x^4*y-10*x^2*y^3+y^5)/r^6)
            return T6m5
        elseif isequal(m,-4)
            T6m4=sqrt(13/(7*pi))*(21/32)*(4*(x^3*y-x*y^3)*(11*z^2-r^2)/r^6)
            return T6m4
        elseif isequal(m,-3)
            T6m3=sqrt(2730/pi)*(1/32)*((3*x^2*y-y^3)*(11*z^3-3*z*r^2)/r^6)
            return T6m3
        elseif isequal(m,-2)
            T6m2=sqrt(2730/pi)*(1/64)*(2*x*y*(33*z^4-13*z^2*r^2+5*r^4)/r^6)
            return T6m2
        elseif isequal(m,-1)
            T6m1=sqrt(273/(4*pi))*(1/8)*(x*z*(33*z^4-30*z^2*r^2+5*r^4)/r^6)
            return T6m1
        elseif isequal(m,0)
            T60=sqrt(13/pi)*(1/32)*((231*z^6-315*z^4*r^2+105*z^2*r^4-5*r^6)/r^6)
            return T60
        elseif isequal(m,1)
            T61=-sqrt(273/(4*pi))*(1/8)*(x*z*(33*z^4-30*z^2*r^2+5*r^4)/r^6)
            return T61
        elseif isequal(m,2)
            T62=sqrt(2730/pi)*(1/64)*((x^2-y^2)*(33*z^4-18*z^2*r^2+r^4)/r^6)
            return T62
        elseif isequal(m,3)
            T63=-sqrt(2730/pi)*(1/32)*((x^3-3*x*y^2)*(11*z^3-3*z*r^2)/r^6)
            return T63
        elseif isequal(m,4)
            T64=sqrt(13/(7*pi))*(21/32)*((x^4-6*x^2*y^2+y^4)*(11*z^2-r^2)/r^6)
            return T64
        elseif isequal(m,5)
            T65=-sqrt(9009/(512*pi))*(z*(x^5-10*x^3*y^2+5*x*y^4)/r^6)
            return T65
        elseif isequal(m,6)
            T66=sqrt(26/(231*pi))*(231/4)*((x^6-15*x^4*y^2+15*x^2*y^4-y^6)/r^6)
            return T66
        end
    end
    @error("Values of l=$(l) and/or m=$(m) invalid. l must be one of [1,2,3,4,5,6] and m takes values between -l and l.")
end


function calc_cefparams!(pcm::local_env)
    cefparams=DataFrame(B=Float64[],l=Int[],m=Int[])
    radwav=pcm.ion.rad_wavefunction
    sfactors=pcm.ion.stevens_factors
    ahc=1.43996e4
    a0=0.52917721067
    for l in [2,4,6]
        if isequal(l,2)
            rl=radwav[1]
            al=sfactors[1]
            sig=pcm.ion.shielding_factors[1]
        elseif isequal(l,4)
            rl=radwav[2]
            al=sfactors[2]
            sig=pcm.ion.shielding_factors[2]
        elseif isequal(l,6)
            rl=radwav[3]
            al=sfactors[3]
            sig=pcm.ion.shielding_factors[3]
        end
        for m in -l:1:l
            Alm=0.0
            for pc in pcm.cartesian_pointcs
                x,y,z,Z=pc
                R=sqrt(x^2+y^2+z^2)
                Tlm=tesseral_harmonics(l,m,x,y,z,R)
                Alm+=((4pi)/(2*l+1))*(Z*Tlm)/(R^(l+1))
            end
            Blm=al*(1-sig)*rl*Alm*ahc*a0^l
            if iszero(Blm)
                continue
            end
            append!(cefparams,DataFrame(:B=>Blm,:l=>l,:m=>m))
        end
    end
    pcm.cefparams=cefparams
    return nothing
end