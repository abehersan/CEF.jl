function print_cef_diagonalization(ion::mag_ion, cefparams::DataFrame; B::Vector{<:Real}=zeros(Float64, 3), method::Symbol=:EO)::Nothing
    cef_matrix = cef_hamiltonian(ion,cefparams;B=B,method=method)
    @assert ishermitian(cef_matrix)
    E = eigvals(cef_matrix)
    E .-= minimum(E)
    printstyled("CEF matrix diagonalization results.\n\n", color=:underline, bold=true)
    display(ion)
    LBL=@sprintf("%1.7f",ion.gj)
    println("Landé g-factor: $(LBL)\n")
    println("External magnetic field in (Tesla) [Bx, By, Bz]: $B\n")
    println("CEF energy levels in (meV) and in (K):")
    for i in eachindex(E)
        EmeV = round(E[i], digits=SDIG)
        EK = round(E[i]/meV_per_K, digits=SDIG)
        println(@sprintf("%.7f\t%.7f", EmeV, EK))
    end
    return nothing
end


function print_cef_diagonalization(lfield::local_env; shielded::Bool=true, B::Vector{<:Real}=zeros(Float64, 3), method::Symbol=:EO)::Nothing
    calc_cefparams!(lfield;shielded)
    print_cef_diagonalization(lfield.ion,lfield.cefparams;B,method)
    return nothing
end


function cef_hamiltonian(ion::mag_ion, cefparams::DataFrame; B::Vector{<:Real}=zeros(Float64, 3), method::Symbol=:EO)::HERMITIANC64
    if iszero(B)
        return H_cef(ion, cefparams, method)
    else
        return H_cef(ion, cefparams, method) + H_zeeman(ion, B)
    end
end


function cef_hamiltonian(ion::mag_ion, D::Real, E::Real; B::Vector{<:Real}=zeros(Float64, 3), method::Symbol=:EO)::HERMITIANC64
    m_dim = Int(2*ion.J+1)
    h_cef = zeros(ComplexF64, (m_dim, m_dim))
    h_cef += D * ion.Jz^2
    h_cef += E * (ion.Jp^2 + ion.Jm^2)/2
    if iszero(B)
        return Hermitian(h_cef)
    else
        return Hermitian(h_cef) + H_zeeman(ion, B)
    end
end


function H_zeeman(ion::mag_ion, extfield::Vector{<:Real})::HERMITIANC64
    Bx, By, Bz = extfield
    gjBxJx = (ion.gj*Bx)*ion.Jx
    gjByJy = (ion.gj*By)*ion.Jy
    gjBzJz = (ion.gj*Bz)*ion.Jz
    zeeman_matrix = (-1.0 * muB) * (gjBxJx .+ gjByJy .+ gjBzJz)
    return Hermitian(zeeman_matrix)
end


function H_cef(ion::mag_ion, cefparams::DataFrame, method::Symbol)::HERMITIANC64
    m_dim = Int(2*ion.J+1)
    cef_matrix = zeros(ComplexF64, (m_dim, m_dim))
    if isequal(method,:EO)
        mode=stevens_EO
    elseif isequal(method,:O)
        mode=stevens_O
    else
        @error "Method: $(method) not one of `:EO` or `:O`."
    end
    @eachrow! cefparams begin
        cef_matrix += :B * mode(ion, :l, :m)
    end
    return Hermitian(cef_matrix)
end


function ryabov_clm(l::Int, m::Int)::Float64
    lmax = 13
    if l > lmax
        @error "Invalid l, l<lmax, where l: $l, lmax: $lmax"
    elseif !(m in -l:1:l)
        @error "Invalid m, m in {-l, l}, where m: $m, l: $l"
    end
    # Flm coefficients calculated by Stoll and implemented in EasySpin
    # see: https://github.com/StollLab/EasySpin/blob/main/easyspin/stev.m
    F = SMatrix{13, 13, Int}([
            1           0           0           0           0       0       0       0       0       0   0   0   0;
            2           1           0           0           0       0       0       0       0       0   0   0   0;
            4           2           1           0           0       0       0       0       0       0   0   0   0;
            24          6           6           1           0       0       0       0       0       0   0   0   0;
            48          24          8           4           1       0       0       0       0       0   0   0   0;
            480         240         240         10          10      1       0       0       0       0   0   0   0;
            2880        1440        360         60          12      6       1       0       0       0   0   0   0;
            40320       5040        1680        168         168     14      14      1       0       0   0   0   0;
            80640       40320       40320       6720        672     336     16      8       1       0   0   0   0;
            1451520     725700      725700      60480       60480   864     288     18      18      1   0   0   0;
            14515200    7257600     1209600     604800      86400   2880    360     180     20      10  1   0   0;
            319334400   79833600    79833600    13305600    2661120 23760   7920    1320    1320    22  22  1   0;
            1916006400  958003200   958003200   31933440    3991680 1995840 31680   15840   1584    264 24  12  1;
        ]
    )
    Flm = F[l+1, abs(m)+1]

    if Bool(mod(l, 2)) # odd l
        alpha = 1.0
    else # even l
        if Bool(mod(m, 2)) # odd m
            alpha = 1.0/2.0
        else # even m
            alpha = 1.0
        end
    end
    clm = alpha/(Flm) # Ryabov (1999) Eq.(22) with Nlm set to 1
    return clm
end


function stevens_EO(ion::mag_ion, l::Int, m::Int)::HERMITIANC64
    T = ion.Jp^l # T^l_l
    for _ in l-1:-1:abs(m) # Eqn (1 and 2) of Ryabov (1999), see Stoll EasySpin
        T = ion.Jm*T - T*ion.Jm
    end
    # Construction of cosine and sine tesseral operators, Ryabov, Eq.(21)
    clm = ryabov_clm(l, m)
    if m >= 0
        Op = clm/2 * (T + adjoint(T))
    else
        Op = clm/2im * (T - adjoint(T))
    end
    return Hermitian(Op)
end


function stevens_O(ion::mag_ion, l::Int, m::Int)::HERMITIANC64
    J=ion.J
    m_dim = Int(2*J+1)
    Jz = ion.Jz
    Jp = ion.Jp
    Jm = ion.Jm
    X = Diagonal(fill(J*(J+1), m_dim))
    if l == 1
        if m == 1
            O11 = 1.0/2.0 * (Jp + Jm)
            return Hermitian(O11)
        elseif m == -1
            O1m1 = -1.0im/2.0 * (Jp - Jm)
            return Hermitian(O1m1)
        elseif m == 0
            O10 = Jz
            return Hermitian(O10)
        else
            err_message =
            "Given values of l=$l and m=$m not supported.\n"*
            "Only l=[2, 4, 6] and m = -l:1:l  are currently supported."
            @error err_message
        end
    elseif l == 2
        if m == -2
            O2m2 = -1.0im/2.0 * (Jp^2 - Jm^2)
            return Hermitian(O2m2)
        elseif m == -1
            O2m1 = -1.0im/4.0 * (Jz*(Jp - Jm) + (Jp - Jm)*Jz)
            return Hermitian(O2m1)
        elseif m == 0
            O20 = 3.0*Jz^2 - X
            return Hermitian(O20)
        elseif m == 1
            O21 = 1.0/4.0 * (Jz*(Jp + Jm) + (Jp + Jm)*Jz)
            return Hermitian(O21)
        elseif m == 2
            O22 = 1.0/2.0 * (Jp^2 + Jm^2)
            return Hermitian(O22)
        else
            err_message =
            "Given values of l=$l and m=$m not supported.\n"*
            "Only l=[2, 4, 6] and m = -l:1:l  are currently supported."
            @error err_message
        end
    elseif l == 4
        if m == -4
            O4m4 = -1.0im/2.0 * (Jp^4 - Jm^4)
            return Hermitian(O4m4)
        elseif m == -3
            O4m3 = -1.0im/4.0 * ((Jp^3 - Jm^3)*Jz + Jz*(Jp^3 - Jm^3))
            return Hermitian(O4m3)
        elseif m == -2
            O4m2 = -1.0im/4.0 * ((Jp^2 - Jm^2)*(7.0*Jz^2 - X - 5.0I) + (7.0*Jz^2 - X - 5.0I)*(Jp^2 - Jm^2))
            return Hermitian(O4m2)
        elseif m == -1
            O4m1 = -1.0im/4.0 * ((Jp - Jm)*(7.0*Jz^3 - (3.0*X + 1.0I)*Jz) + (7.0*Jz^3 - (3.0*X + 1.0I)*Jz)*(Jp - Jm))
            return Hermitian(O4m1)
        elseif m == 0
            O40 = 35.0*Jz^4 - (30.0*X - 25.0I)*Jz^2 + 3.0*X^2 - 6.0*X
            return Hermitian(O40)
        elseif m == 1
            O41 = 1.0/4.0 * ((Jp + Jm)*(7.0*Jz^3 - (3.0*X + 1.0I)*Jz) + (7.0*Jz^3 - (3.0*X + 1.0I)*Jz)*(Jp + Jm))
            return Hermitian(O41)
        elseif m == 2
            O42 = 1.0/4.0 * ((Jp^2 + Jm^2)*(7.0*Jz^2 - X - 5.0I) + (7.0*Jz^2 - X - 5.0I)*(Jp^2 + Jm^2))
            return Hermitian(O42)
        elseif m == 3
            O43 = 1.0/4.0 * (Jz*(Jp^3 + Jm^3) + (Jp^3 + Jm^3)*Jz)
            return Hermitian(O43)
        elseif m == 4
            O44 = 1.0/2.0 * (Jp^4 + Jm^4)
            return Hermitian(O44)
        else
            err_message =
            "Given values of l=$l and m=$m not supported.\n"*
            "Only l=[2, 4, 6] and m = -l:1:l  are currently supported."
            @error err_message
        end
    elseif l == 6
        if m == -6
            O6m6 = -1.0im/2.0 * (Jp^6 - Jm^6)
            return Hermitian(O6m6)
        elseif m == -5
            O6m5 = -1.0im/4.0 * ((Jp^5 - Jm^5)*Jz + Jz*(Jp^5 - Jm^5))
            return Hermitian(O6m5)
        elseif m == -4
            O6m4 = -1.0im/4.0 * ((Jp^4 - Jm^4)*(11.0*Jz^2 - X - 38.0I) + (11.0*Jz^2 - X - 38.0I)*(Jp^4 - Jm^4))
            return Hermitian(O6m4)
        elseif m == -3
            O6m3 = -1.0im/4.0 * ((Jp^3 - Jm^3)*(11.0*Jz^3 - (3.0*X + 59.0I)*Jz) + (11.0*Jz^3 - (3.0*X + 59.0I)*Jz)*(Jp^3 - Jm^3))
            return Hermitian(O6m3)
        elseif m == -2
            O6m2 = -1.0im/4.0 * ((Jp^2 - Jm^2)*(33.0*Jz^4 - (18.0*X + 123.0I)*Jz^2 + X^2 + 10.0*X + 102I) + (33.0*Jz^4 - (18.0*X + 123.0I)*Jz^2 + X^2 + 10.0*X + 102.0I)*(Jp^2 - Jm^2))
            return Hermitian(O6m2)
        elseif m == -1
            O6m1 = -1.0im/4.0 * ((Jp - Jm)*(33.0*Jz^5 - (30.0*X - 15.0I)*Jz^3 + (5.0*X^2 - 10.0*X + 12.0I)*Jz) + (33.0*Jz^5 - (30.0*X - 15.0I)*Jz^3 + (5.0*X^2 - 10.0*X + 12.0I)*Jz)*(Jp - Jm))
            return Hermitian(O6m1)
        elseif m == 0
            O60 = 231.0*Jz^6 - (315.0*X - 735.0I)*Jz^4 + (105.0*X^2 - 525.0*X + 294.0I)*Jz^2 - 5.0*X^3 + 40.0*X^2 - 60.0*X
            return Hermitian(O60)
        elseif m == 1
            O61 = 1.0/4.0 * ((Jp + Jm)*(33.0*Jz^5 - (30.0X - 15.0I)*Jz^3 + (5.0*X^2 - 10.0*X + 12.0I)*Jz) + (33.0*Jz^5 - (30.0*X - 15.0I)*Jz^3 + (5.0*X^2 - 10.0*X + 12.0I)*Jz)*(Jp + Jm))
            return Hermitian(O61)
        elseif m == 2
            O62 = 1.0/4.0 *((Jp^2 + Jm^2)*(33.0*Jz^4 - (18*X + 123.0I)*Jz^2 + X^2 + 10.0*X + 102.0I) + (33.0*Jz^4 - (18.0*X + 123.0I)*Jz^2 + X^2 + 10.0*X + 102.0I)*(Jp^2 + Jm^2))
            return Hermitian(O62)
        elseif m == 3
            O63 = 1.0/4.0 * ((Jp^3 + Jm^3)*(11.0*Jz^3 - (3.0*X + 59.0I)*Jz) + (11.0*Jz^3 - (3.0*X + 59.0I)*Jz)*(Jp^3 + Jm^3))
            return Hermitian(O63)
        elseif m == 4
            O64 = 1.0/4.0 * ((Jp^4 + Jm^4)*(11.0*Jz^2 - X - 38.0I) + (11.0*Jz^2 - X - 38.0I)*(Jp^4 + Jm^4))
            return Hermitian(O64)
        elseif m == 5
            O65 = 1.0/4.0 * ((Jp^5 + Jm^5)*Jz + Jz*(Jp^5 + Jm^5))
            return Hermitian(O65)
        elseif m == 6
            O66 = 1.0/2.0 * (Jp^6 + Jm^6)
            return Hermitian(O66)
        else
            err_message =
            "Given values of l=$l and m=$m not supported.\n"*
            "Only l=[2, 4, 6] and m = -l:1:l  are currently supported."
            @error err_message
        end
    else
        err_message =
        "Given values of l=$l and m=$m not supported.\n"*
        "Only l=[2, 4, 6] and m = -l:1:l  are currently supported."
        @error err_message
    end
end