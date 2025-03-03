const re_hundsrules = Dict(
    "Ce3+"=>[
        5/2, 6/7,                                           # J, gJ
        -5.7143/1e2, 63.4921/1e4, 0,                        # <r2>, <r4>, <r6>
        0.3666, 0.3108, 0.5119,                             # Stevens geometrical factors
        0.2291, 18.18, 0.7897, 5.807, -0.0191, 0.0, 0.0,    # for j0
        2.1284, 8.9174, 1.1229, 2.8371, 0.01108, 0.0, 0.0   # for j2
        ],
    "Pr3+"=>[
        4/1, 4/5,
        -2.101/1e2, -7.3462/1e4, 60.994/1e6,
        0.3350, 0.2614, 0.4030,
        0.0504, 24.9989, 0.2572, 12.0377, 0.7142, 5.0039, -0.0219,
        0.8734, 18.9876, 1.5594, 6.0872, 0.8142, 2.4150, 0.0111
        ],
    "Nd3+"=>[
        9/2, 8/11,
        -0.6428/1e2, -2.9111/1e4, -37.9980/1e6,
        0.3120, 0.2282, 0.3300,
        0.0540, 25.0293, 0.3101, 12.1020, 0.6575, 4.7223, -0.0216,
        0.6751, 18.3421, 1.6272, 7.2600, 0.9644, 2.6016, 0.0150
        ],
    "Pm3+"=>[
        4/1, 3/5,
        0.7713/1e2, 4.0755/1e4, 60.7807/1e6,
        0.2899, 0.1991, 0.2755,
        0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0
        ],
    "Sm3+"=>[
        5/2, 2/7,
        4.127/1e2, 25.0120/1e4, 0,
        0.2728, 0.1772, 0.2317,
        0.0288, 25.2068, 0.2973, 11.8311, 0.6954, 4.2117, -0.0213,
        0.4707, 18.4301, 1.4261, 7.0336, 0.9574, 2.4387, 0.0182
        ],
    "Tb3+"=>[
        6/1, 3/2,
        -1.0101/1e2, 1.2244/1e4, -1.1212/1e6,
        0.2302, 0.1295, 0.1505,
        0.0177, 25.5095, 0.2921, 10.5769, 0.7133, 3.5122, -0.0231,
        0.2892, 18.4973, 1.1678, 6.7972, 0.9437, 2.2573, 0.0232
        ],
    "Dy3+"=>[
        15/2, 4/3,
        -0.6349/1e2, -0.592/1e4, 1.035/1e6,
        0.2188, 0.1180, 0.1328,
        0.1157, 15.0732, 0.3270, 6.7991, 0.5821, 3.0202, -0.0249,
        0.2523, 18.5172, 1.0914, 6.7362, 0.9345, 2.2082, 0.0250
        ],
    "Ho3+"=>[
        8/1, 5/4,
        -0.2222/1e2, -0.333/1e4, -1.2937/1e6,
        0.2085, 0.1081, 0.1181,
        0.0566, 18.3176, 0.3365, 7.6880, 0.6317, 2.9427, -0.0248,
        0.2188, 18.5157, 1.0240, 6.7070, 0.9251, 2.1614, 0.0268
        ],
    "Er3+"=>[
        15/2, 6/5,
        0.254/1e2, 0.444/1e4, 2.0699/1e6,
        0.1991, 0.0996, 0.1058,
        0.0586, 17.9802, 0.3540, 7.0964, 0.6126, 2.7482, -0.0251,
        0.1710, 18.5337, 0.9879, 6.6246, 0.9044, 2.1004, 0.0278
        ],
    "Tm3+"=>[
        6/1, 7/6,
        1.0101/1e2, 1.632/1e4, -5.606/1e6,
        0.1905, 0.0921, 0.0953,
        0.0581, 15.0922, 0.2787, 7.8015, 0.6854, 2.7931, -0.0224,
        0.1760, 18.5417, 0.9105, 6.5787, 0.8970, 2.0622, 0.0294
        ],
    "Yb3+"=>[
        7/2, 8/7,
        3.175/1e2, -17.3160/1e4, 148.0/1e6,
        0.1826, 0.0854, 0.0863,
        0.0416, 16.0949, 0.2849, 7.8341, 0.6961, 2.6725, -0.0229,
        0.1570, 18.5553, 0.8484, 6.5403, 0.8880, 2.0367, 0.0318
        ],
)


@doc raw"""
    single_ion(ion::String)

Given a magnetic ion name, generate a `mag_ion` struct that contains
the following data:

- Total angular momentum quantum number J, isotropic Landé-factor gJ
- Stevens geometric factors θ_l, α, β, γ
- Expectation value of radial wave functions <r^l>, l={2, 4, 6}
- Magnetic form factor coefficients

Currently only trivalent rare-earth ions are supported.
See `custom_spin` if other ions are required.

Coefficients taken from Rotter Bauer (2010),
the McPhase online manual: https://mcphase.github.io/webpage/manual/node129.html
and Brown's coefficients https://www.ill.eu/sites/ccsl/ffacts/ffactnode1.html
"""
function single_ion(ion::String)
    try
        J, g,
        alpha, beta, gamma,
        r2, r4, r6,
        A_j0, a_j0, B_j0, b_j0, C_j0, c_j0, D_j0,
        A_j2, a_j2, B_j2, b_j2, C_j2, c_j2, D_j2 = re_hundsrules[ion]
        C2=(2-g)/g
        return mag_ion(
                ion, J,
                spin_operators(J, "x"),
                spin_operators(J, "y"),
                spin_operators(J, "+"),
                spin_operators(J, "-"),
                spin_operators(J, "z"),
                [g, g, g],
                [alpha, beta, gamma],
                [r2, r4, r6],
                C2,
                [A_j0, a_j0, B_j0, b_j0, C_j0, c_j0, D_j0],
                [A_j2, a_j2, B_j2, b_j2, C_j2, c_j2, D_j2],
            )
    catch y
        mag_ions = collect(keys(re_hundsrules))
        err_message =
        "$y\n"*
        "Given ion $ion not supported.\n"*
        "Available magnetic ions: $mag_ions"
        @error err_message
    end
end


function custom_ion(S::Real, g::Real, C2::Real)::mag_ion
    return mag_ion(
        ion="S", J=S,
        Jx=spin_operators(S,"x"),
        Jy=spin_operators(S,"y"),
        Jz=spin_operators(S,"z"),
        Jp=spin_operators(S,"+"),
        Jm=spin_operators(S,"-"),
        g=[g,g,g],
        stevens_factors=zeros(3),
        rad_wavefunction=zeros(3),
        C2=C2,
        ff_coeff_j0=zeros(7),
        ff_coeff_j2=zeros(7),
    )
end


Base.@kwdef mutable struct mag_ion
    ion::String
    J::Float64
    Jx::Matrix{ComplexF64}
    Jy::Matrix{ComplexF64}
    Jp::Matrix{ComplexF64}
    Jm::Matrix{ComplexF64}
    Jz::Matrix{ComplexF64}
    g::MVEC{3}
    stevens_factors::VEC{3}
    rad_wavefunction::VEC{3}
    C2::Float64
    ff_coeff_j0::VEC{7}
    ff_coeff_j2::VEC{7}
end


function Base.show(io::IO, ::MIME"text/plain", ion::mag_ion)
    printstyled(io, "Magnetic ion: $(ion.ion)\n")
    if isequal(ion.ion, "S")
        @warn "Custom ion! Geometric factors, radial wavefunction coefficients and form-factor coefficients are zero. Manually initialize them."
    end
    print(io, "Quantum number J: $(ion.J).\nHilbert space dimension: $(Int(2*ion.J+1)).")
    return nothing
end


@doc raw"""
    spin_operators(J::Float64, a::String)::Matrix{ComplexF64}

Explicit matrices for spin operators of given total angular momentum
quantum number `J`.
The string `a` is one of either `x`, `y`, `z` for the Cartesian spin-operators.
If `a` is either `+` or `-` explicit matrices for the ladder operators
are returned.
"""
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