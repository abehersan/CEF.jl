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
    printstyled(io, "Cartesian ligand field (Å, -Ze):\n\tx\t\ty\t\tz\t\tZ\n")
    for pc in pcenv.cartesian_pointcs
        x,y,z,Z=pc
        LBL=@sprintf("\t%+3.5f\t%+3.5f\t%+3.5f\t%+3.5f\n",x,y,z,Z)
        printstyled(io,LBL)
    end
    println()
    printstyled(io, "Spherical ligand field (Å, deg., -Ze):\n\tr\t\tθ\t\tϕ\t\tZ\n")
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
    cgamma=sind(gamma)
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


function calc_stevensparams()
    # basis vectors from lattice constants
    # ligand positions relative to the lanthanide and effective charges
    # spherical coordinates trafo
    # calculate Blm using formula
    # rotation for compact wavefunctions? most likely post processing function
    # return Blm dframe
end