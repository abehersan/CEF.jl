const re_gs=Dict(
    # ion J gj <r2> <r4> <r6> s2 s4 s6 alpha beta gamma
    # J. Jensen and A. R. Mackintosh, Rare earth magnetism
    # S. Edvardsson and M. Klintenberg, Journal of Alloys and Compounds 275-277, 230 (1998)
    "Ce3+"=>[ 5/2,  6/7,  1.456, 5.437, 42.26, 0.510, 0.0132, -0.0294, -5.7140e-2, 63.490e-4, 00.000000],
    "Pr3+"=>[ 4/1,  4/5,  1.327, 4.537, 32.65, 0.514, 0.0150, -0.0302, -2.1010e-2, -7.346e-4, 60.990e-6],
    "Nd3+"=>[ 9/2,  8/11, 1.222, 3.875, 26.12, 0.515, 0.0164, -0.0307, -0.6428e-2, -2.911e-4, -37.99e-6],
    "Pm3+"=>[ 4/1,  3/5,  1.135, 3.366, 21.46, 0.512, 0.0175, -0.0309, 0.77140e-2, 4.0760e-4, 60.780e-6],
    "Sm3+"=>[ 5/2,  2/7,  1.061, 2.964, 17.99, 0.507, 0.0184, -0.0309, 4.12700e-2, 25.010e-4, 00.000000],
    "Tb3+"=>[ 6/1,  3/2,  0.893, 2.163, 11.75, 0.486, 0.0193, -0.0300, -1.0101e-2, 1.2240e-4, -1.121e-6],
    "Dy3+"=>[ 15/2, 4/3,  0.849, 1.977, 10.44, 0.477, 0.0193, -0.0295, -0.6349e-2, -0.592e-4, 1.0350e-6],
    "Ho3+"=>[ 8/1,  5/4,  0.810, 1.816, 9.345, 0.469, 0.0192, -0.0289, -0.2222e-2, -0.333e-4, -1.294e-6],
    "Er3+"=>[ 15/2, 6/5,  0.773, 1.677, 8.431, 0.460, 0.0190, -0.0283, 0.25400e-2, 0.4440e-4, 2.0700e-6],
    "Tm3+"=>[ 6/1,  7/6,  0.740, 1.555, 7.659, 0.450, 0.0188, -0.0277, 1.01010e-2, 1.6320e-4, -5.606e-6],
    "Yb3+"=>[ 7/2,  8/7,  0.710, 1.448, 7.003, 0.441, 0.0185, -0.0270, 3.17500e-2, -17.32e-4, 148.00e-6]
)

const re_ffparams=Dict(
    # j0 coeffs j2 coeffs
    # https://mcphase.github.io/webpage/manual/node148.html data from J. Brown
    "Ce3+"=>[ 0.2291,  18.180, 0.7897,  5.8070, -0.0191, 0.0000,  0.0000, 2.1284,  8.9174, 1.1229, 2.8371, 0.01108, 0.0000, 0.0000 ],
    "Pr3+"=>[ 0.0504, 24.9989, 0.2572, 12.0377,  0.7142, 5.0039, -0.0219, 0.8734, 18.9876, 1.5594, 6.0872, 0.81420, 2.4150, 0.0111 ],
    "Nd3+"=>[ 0.0540, 25.0293, 0.3101, 12.1020,  0.6575, 4.7223, -0.0216, 0.6751, 18.3421, 1.6272, 7.2600, 0.96440, 2.6016, 0.0150 ],
    "Pm3+"=>[ 0.0000,  0.0000, 0.0000,  0.0000,  0.0000, 0.0000,  0.0000, 0.0000,  0.0000, 0.0000, 0.0000, 0.00000, 0.0000, 0.0000 ],
    "Sm3+"=>[ 0.0288, 25.2068, 0.2973, 11.8311,  0.6954, 4.2117, -0.0213, 0.4707, 18.4301, 1.4261, 7.0336, 0.95740, 2.4387, 0.0182 ],
    "Tb3+"=>[ 0.0177, 25.5095, 0.2921, 10.5769,  0.7133, 3.5122, -0.0231, 0.2892, 18.4973, 1.1678, 6.7972, 0.94370, 2.2573, 0.0232 ],
    "Dy3+"=>[ 0.1157, 15.0732, 0.3270,  6.7991,  0.5821, 3.0202, -0.0249, 0.2523, 18.5172, 1.0914, 6.7362, 0.93450, 2.2082, 0.0250 ],
    "Ho3+"=>[ 0.0566, 18.3176, 0.3365,  7.6880,  0.6317, 2.9427, -0.0248, 0.2188, 18.5157, 1.0240, 6.7070, 0.92510, 2.1614, 0.0268 ],
    "Er3+"=>[ 0.0586, 17.9802, 0.3540,  7.0964,  0.6126, 2.7482, -0.0251, 0.1710, 18.5337, 0.9879, 6.6246, 0.90440, 2.1004, 0.0278 ],
    "Tm3+"=>[ 0.0581, 15.0922, 0.2787,  7.8015,  0.6854, 2.7931, -0.0224, 0.1760, 18.5417, 0.9105, 6.5787, 0.89700, 2.0622, 0.0294 ],
    "Yb3+"=>[ 0.0416, 16.0949, 0.2849,  7.8341,  0.6961, 2.6725, -0.0229, 0.1570, 18.5553, 0.8484, 6.5403, 0.88800, 2.0367, 0.0318 ]
)

function single_ion(ion::String)
    try
        iongs=re_gs[ion]
        ionff=re_ffparams[ion]
        J, gj, r2, r4, r6, s2, s4, s6, alpha, beta, gamma=iongs
        A_j0, a_j0, B_j0, b_j0, C_j0, c_j0, D_j0, A_j2, a_j2, B_j2, b_j2, C_j2, c_j2, D_j2 = ionff
        C2=(2-gj)/gj
        return mag_ion(
                ion, J,
                spin_operators(J, "x"),
                spin_operators(J, "y"),
                spin_operators(J, "z"),
                spin_operators(J, "+"),
                spin_operators(J, "-"),
                gj,
                [alpha, beta, gamma],
                [r2, r4, r6],
                [s2, s4, s6],
                C2,
                [A_j0, a_j0, B_j0, b_j0, C_j0, c_j0, D_j0],
                [A_j2, a_j2, B_j2, b_j2, C_j2, c_j2, D_j2],
            )
    catch y
        mag_ions = collect(keys(re_gs))
        err_message =
        "$y\n"*
        "Given ion $ion not supported.\n"*
        "Available magnetic ions: $mag_ions"
        @error err_message
    end
end


Base.@kwdef mutable struct mag_ion
    ion::String
    J::Float64
    Jx::Matrix{ComplexF64}
    Jy::Matrix{ComplexF64}
    Jz::Matrix{ComplexF64}
    Jp::Matrix{ComplexF64}
    Jm::Matrix{ComplexF64}
    gj::Float64
    stevens_factors::VEC{3}
    rad_wavefunction::VEC{3}
    shielding_factors::VEC{3}
    C2::Float64
    ff_coeff_j0::VEC{7}
    ff_coeff_j2::VEC{7}
end


function Base.show(io::IO, ::MIME"text/plain", ion::mag_ion)
    printstyled(io, "Magnetic ion: $(ion.ion)\n")
    print(io, "Quantum number J: $(ion.J).\nHamiltonian matrix dimension: $(Int(2*ion.J+1))×$(Int(2*ion.J+1)).")
    return nothing
end


function spin_operators(J::Float64, a::String)::Matrix{ComplexF64}
    mJ = -J:1:J
    if isequal(a, "x")
        jp_eigval = @. sqrt(J*(J+1)-mJ*(mJ+1))
        jm_eigval = @. sqrt(J*(J+1)-mJ*(mJ-1))
        Jp = diagm(1=>jp_eigval[1:end-1])
        Jm = diagm(-1=>jm_eigval[2:end])
        Jx = (Jp + Jm)/2.0
        return Jx
    elseif isequal(a, "y")
        jp_eigval = @. sqrt(J*(J+1)-mJ*(mJ+1))
        jm_eigval = @. sqrt(J*(J+1)-mJ*(mJ-1))
        Jp = diagm(1=>jp_eigval[1:end-1])
        Jm = diagm(-1=>jm_eigval[2:end])
        Jy = (Jp - Jm)/2.0im
        return Jy
    elseif isequal(a, "z")
        Jz = diagm(mJ)
        return Jz
    elseif isequal(a, "+")
        jp_eigval = @. sqrt(J*(J+1)-mJ*(mJ+1))
        Jp = diagm(1=>jp_eigval[1:end-1])
        return Jp
    elseif isequal(a, "-")
        jm_eigval = @. sqrt(J*(J+1)-mJ*(mJ-1))
        Jm = diagm(-1=>jm_eigval[2:end])
        return Jm
    else
        @error "String $(a) not understood. Choose one of either {x, y, z, +, -}"
    end
end