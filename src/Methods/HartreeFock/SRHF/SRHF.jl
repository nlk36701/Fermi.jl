using Fermi.DIIS

using TensorOperations
using LinearAlgebra
using Formatting
import Base: show

export SRHF

"""
    Fermi.HartreeFock.RHFAlgorithm

Abstract type for RHF implementations.
"""
abstract type SRHFAlgorithm end


"""
    Fermi.HartreeFock.actual_SRHF

Wave function object for Restricted Hartree-Fock methods

# High Level Interface 
Run a RHF computation and return the RHF object:
```
julia> @energy rhf
```
Equivalent to
```
julia> Fermi.HartreeFock.RHF()
```
Computes RHF using information from `Fermi.Options.Current`

# Fields

| Name   |   Description     |
|--------|---------------------|
| `molecule` |   Molecule object |
| `energy`   |   RHF Energy      |
| `ndocc`    | Number of doubly occupied spatial orbitals |
| `nvir`  | Number of virtual spatial orbitals |
| `orbitals` |    RHF Orbitals object      |
| `e_conv`   | ΔE from the last iteration  |
| `d_conv`   |  Orbitals RMS change from the last iteration|

# Relevant options 

These options can be set with `@set <option> <value>`

| Option         | What it does                      | Type      | choices [default]     |
|----------------|-----------------------------------|-----------|-----------------------|
| `rhf_alg`      | Picks RHF algorithm               | `Int`     | [1]                   |
| `scf_max_rms`  | RMS density convergence criterion | `Float64` | [10^-9]               |
| `scf_max_iter` | Max number of iterations          | `Int`     | [50]                  |
| `scf_e_conv`   | Energy convergence criterion      | `Float64` | [10^-10]              |
| `basis`        | What basis set to use             | `String`  | ["sto-3g"]            |
| `df`           | Whether to use density fitting    | `Bool`    | `true` [`false`]      |
| `jkfit`        | What aux. basis set to use for JK | `String`  | ["auto"]              |
| `diis`         | Whether to use DIIS               | `Bool`    | [`true`] `false`      |
| `oda`          | Whether to use ODA                | `Bool`    | [`true`] `false`      |
| `oda_cutoff`   | When to turn ODA off (RMS)        | `Float64` | [1E-1]                |
| `oda_shutoff`  | When to turn ODA off (iter)       | `Int`     | [20]                  |
| `scf_guess`    | Which guess density to use        | `String`  | "core" ["gwh"]        |

# Struct tree

**RHF** <: AbstractHFWavefunction <: AbstractWavefunction
"""
struct SRHF <: AbstractHFWavefunction
    molecule::Molecule
    energy::Float64
    ndocc::Int
    nvir::Int
    orbitals::SRHFOrbitals
    e_conv::Float64
    d_conv::Float64
end

function SRHF(x...)
    if !any(i-> i isa SRHFAlgorithm, x)
        SRHF(x..., get_srhf_alg())
    else
        # Print the type of arguments given for a better feedback
        args = "("
        for a in x[1:end-1]
            args *= "$(typeof(a)), "
        end
        args = args[1:end-2]*")"
        throw(FermiException("invalid arguments for SRHF method: $args"))
    end
end

"""
    Fermi.HartreeFock.get_rhf_alg()

Returns a singleton type corresponding to a RHF implementation.
"""
function get_srhf_alg(N::Int = Options.get("srhf_alg"))
    try 
        return get_srhf_alg(Val(N))
    catch MethodError
        throw(FermiException("implementation number $N not available for RHF."))
    end
end

# Actual RHF routine is in here
# For each implementation a singleton type must be create
struct SRHFa <: SRHFAlgorithm end
include("SRHFa.jl")
include("SRHFHelper.jl")
# And a number is assigned to the implementation
get_srhf_alg(x::Val{1}) = SRHFa()

#struct SDirect <: SRHFAlgorithm end
#include("SDirect.jl")
#include("SDirectHelper.jl")
#get_srhf_alg(x::Val{2}) = SDirect()


# Gradient methods
#include("Gradients/RHFgrad.jl")

### MISCELLANEOUS
# Pretty printing
function string_repr(X::RHF)
    out = ""
    out = out*" ⇒ Fermi Restricted Hartree--Fock Wave function\n"
    out = out*" ⋅ Basis:                  $(X.orbitals.basis)\n"
    out = out*" ⋅ Energy:                 $(X.energy)\n"
    out = out*" ⋅ Occ. Spatial Orbitals:  $(X.ndocc)\n"
    out = out*" ⋅ Vir. Spatial Orbitals:  $(X.nvir)\n"
    out = out*"Convergence: " 
    out = out*"ΔE => $(format("{:1.2e}",abs(X.e_conv)))"
    out = out*" Dᵣₘₛ => $(format("{:1.2e}",abs(X.d_conv)))"
    return out
end

function show(io::IO, ::MIME"text/plain", X::RHF)
    print(io, string_repr(X))
end

