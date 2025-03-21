module CEF

using DataFrames
using DataFramesMeta
using LinearAlgebra
using OffsetArrays
using StaticArrays
using Statistics
using Printf
using Trapz

const PREC::Float64 = 1.0e-7    # for degeneracy calculations
const SDIG::Int64 = 7           # for numerical cutoffs

# static vector/array definitions
const VEC3 = SVector{3, Float64}
const MAT3 = MMatrix{3, 3, Float64, 9}
const VEC{N} = SVector{N, Float64}
const MVEC{N} = MVector{N, Float64}
const CVEC{N} = SVector{N, ComplexF64}
const CMAT{N} = MMatrix{N, N, ComplexF64}
const HERMITIANC64 = Hermitian{ComplexF64, Matrix{ComplexF64}}

# physical constants
const meV_per_K = 0.086173332621451774  # divide E in meV by meV_per_K to get E in K
const mu0 = 1.25663706212e-6            # [N/A]
const muB = 0.057883738013331           # [meV/T]
const kB = 0.08617333262                # [meV/K]
const NA = 6.02214076e23                # Avogadro constant, [1/mol]
const Rg = 8.314462618                  # Ideal gas constant, [J/mol/K]

include("./single_ion.jl")
export single_ion
export custom_ion
export mag_ion
export spin_operators
export re_hundsrules

include("./cef_utils.jl")
export blm_dframe
export get_alm!

include("./cef_system.jl")
export print_cef_diagonalization
export cef_hamiltonian
export stevens_EO, stevens_O

include("./thermodynamical_quantities.jl")
export partition_function
export population_factor
export thermal_average

include("./cef_entropy.jl")
export cef_entropy!
export cef_entropy_speclevels!

include("./cef_magnetization.jl")
export cef_magneticmoment_crystal!
export cef_magneticmoment_powder!

include("./cef_susceptibility.jl")
export cef_susceptibility_crystal!
export cef_susceptibility_powder!

include("./cef_neutronxsection.jl")
export cef_neutronxsection_crystal!
export cef_neutronxsection_powder!

include("./cef_rotation.jl")
export get_euler_angles
export rotate_blm

include("./cef_rpa.jl")
export calc_chi0

include("./cef_pcm.jl")
export local_env
export ligand_field
export lattice_vectors
export tesseral_harmonics
export calc_cefparams!

end