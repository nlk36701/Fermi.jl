# Auxiliar functions for RHF computations
"""
    function RHF_core_guess(ints::IntegralHelper)

Compute an initial orbital coefficient matrix using the core guess.
"""
function SRHF_core_guess(ints::IntegralHelper)
    output("Using Core Guess")
    S = ints["S"]
    Λ = S^(-1/2)
    F = ints["T"] + ints["V"]
    Ft = Λ*F*Λ'

    # Get orbital energies and transformed coefficients
    _, Ct = LinearAlgebra.eigen(Symmetric(Ft), sortby=x->x)

    # Reverse transformation to get MO coefficients
    C = Λ*Ct

    return C, Λ
end

"""
    function RHF_gwh_guess(ints::IntegralHelper)

Compute an initial orbital coefficient matrix using the GWH guess.
"""
function SRHF_gwh_guess(ints::IntegralHelper)

    # Form GWH guess
    output("Using GWH Guess")
    molecule = ints.molecule
    S = ints["S"]
    d, U = LinearAlgebra.eigen(Symmetric(S), sortby = x->1/abs(x))
    Λ = Hermitian(S)^(-1/2) |> Array
    idxs = [abs(d[i]) > 1E-7 for i = eachindex(d)]

    H = ints["T"] + ints["V"]
    ndocc = molecule.Nα
    nvir = size(S,1) - ndocc
    F = similar(S)

    for i = 1:(ndocc+nvir)
        F[i,i] = H[i,i]
        for j = i+1:ndocc+nvir
            F[i,j] = 0.875*S[i,j]*(H[i,i] + H[j,j])
            F[j,i] = F[i,j]
        end
    end
    Ft = Λ'*F*Λ

    # Get orbital energies and transformed coefficients
    _, Ct = LinearAlgebra.eigen(Symmetric(Ft), sortby=x->x)

    # Reverse transformation to get MO coefficients
    C = Λ*Ct

    return C, Λ         
end

"""
    RHFEnergy(D::Array{Float64}, H::Array{Float64}, F::Array{Float64})

Compute RHF energy given a density matrix `D`, Core Hamiltonian `H` and Fock matrix `F`.
"""
function SRHFEnergy(D::Array{Float64}, H::Array{Float64}, F::Array{Float64})
    return sum(D .* (H .+ F))
end

"""
    build_fock!(F::Array{Float64}, H::Array{Float64}, D::Array{Float64}, ERI::Array{Float64,4})

Build a Fock matrix into `F` using the Core Hamiltonian `H`, density matrix `D` and two-electron repulsion integral `ERI`.
"""
function build_fock!(F::Array{Float64}, H::Array{Float64}, D::Array{Float64}, ints::IntegralHelper{Float64,Chonky,AtomicOrbitals})
    ERI = ints["ERI"]
    F .= H
    @tensor F[m,n] += 2*D[r,s]*ERI[m,n,r,s]
    @tensor F[m,n] -= D[r,s]*ERI[m,r,n,s]
end

"""
    build_fock!(F::Array{Float64,2}, H::Array{Float64,2}, D::Array{Float64,2}, ints::IntegralHelper{Float64,<:AbstractDFERI,AtomicOrbitals})

Build a Fock matrix into `F` using the Core Hamiltonian `H`, density matrix `D` and two-electron repulsion integral `ERI`
approximated by density fitting.
"""
function build_fock!(F::Array{Float64,2}, H::Array{Float64,2}, D::Array{Float64,2}, ints::IntegralHelper{Float64,<:AbstractDFERI,AtomicOrbitals})
    b = ints["ERI"]
    F .= H
    @tensoropt F[m,n] += 2*D[r,s]*b[Q,m,n]*b[Q,r,s]
    @tensoropt F[m,n] -= D[r,s]*b[Q,m,r]*b[Q,n,s]
end

function newindex(i, irreplength)
    sum = 0
    for (a, Γ) in enumerate(irreplength)
        if i <= Γ + sum
            return i - sum, a
        end
    sum += Γ
    end
end 
#build fock in subspace of irreps instead of massive fock
function smallfock!(otherF, H, D, symsparse, symdex, nbas, so_irrep, irreplength)
    D = D
    otherF = smallFock(irreplength)
    Ft = smallFock(irreplength)
    otherF .= deepcopy(H)
    for z = eachindex(symdex)
        indx = symdex[z]
        ν = symsparse[z]
        i,j,k,l = indx[1], indx[2], indx[3], indx[4] #.+ 1
        I, iΓ = newindex(i, irreplength)
        J, jΓ = newindex(j, irreplength)
        K, kΓ = newindex(k, irreplength)
        L, lΓ = newindex(l, irreplength)
        

        #ir, jr, kr, lr = so_irrep[i], so_irrep[j], so_irrep[k], so_irrep[l]
        ij = Fermi.index2(i,j) 
        kl = Fermi.index2(k,l) 
        
        # Logical auxiliar: γpq (whether p and q are different) Xpq = δpq + 1
        γij = i !== j
        γkl = k !== l
        γab = ij !== kl

        Xik = i === k ? 2.0 : 1.0
        Xjk = j === k ? 2.0 : 1.0
        Xil = i === l ? 2.0 : 1.0
        Xjl = j === l ? 2.0 : 1.0
        if γij && γkl && γab
            #J
            if kΓ == lΓ && iΓ == jΓ
                Ft[iΓ][I,J] += 4.0*D[kΓ][K,L]*ν
                Ft[kΓ][K,L] += 4.0*D[iΓ][I,J]*ν
            end
            # K
            if jΓ == lΓ && kΓ == iΓ 
                Ft[iΓ][I,K] -= Xik*D[jΓ][J,L]*ν
                Ft[jΓ][J,L] -= Xjl*D[iΓ][I,K]*ν
            end
            if iΓ == lΓ && jΓ == kΓ
                Ft[jΓ][J,K] -= Xjk*D[iΓ][I,L]*ν
                Ft[iΓ][I,L] -= Xil*D[jΓ][J,K]*ν
            end
        elseif γkl && γab
            # J
            if kΓ == lΓ && iΓ == jΓ 
                Ft[iΓ][I,J] += 4.0*D[kΓ][K,L]*ν
                Ft[kΓ][K,L] += 2.0*D[iΓ][I,J]*ν
            end
            # K
            if jΓ == lΓ && kΓ == iΓ
                Ft[iΓ][I,K] -= Xik*D[jΓ][J,L]*ν
            end
            if jΓ == kΓ && iΓ == lΓ
                Ft[iΓ][I,L] -= Xil*D[jΓ][J,K]*ν
            end
        elseif γij && γab
            # J
            if kΓ == lΓ && iΓ == jΓ 
                Ft[iΓ][I,J] += 2.0*D[kΓ][K,L]*ν
                Ft[kΓ][K,L] += 4.0*D[iΓ][I,J]*ν
            end
            # K
            if jΓ == lΓ && iΓ == kΓ 
                Ft[iΓ][I,K] -= Xik*D[jΓ][J,L]*ν
            end
            if iΓ == lΓ && jΓ == kΓ
                Ft[jΓ][J,K] -= Xjk*D[iΓ][I,L]*ν
            end

        elseif γij && γkl

            # Only possible if i = k and j = l
            # and i < j ⇒ i < l

            # J
            if kΓ == lΓ && iΓ == jΓ && jΓ == kΓ
                Ft[iΓ][I,J] += 4.0*D[kΓ][K,L]*ν
                Ft[iΓ][I,L] -= D[jΓ][J,K]*ν
            end

            # K
            #not adding a statement here since j = l and i = k....
            if jΓ == lΓ && iΓ == kΓ
                Ft[iΓ][I,K] -= D[jΓ][J,L]*ν
                Ft[jΓ][J,L] -= D[iΓ][I,K]*ν
            end
        elseif γab
            # J
            if kΓ == lΓ && iΓ == jΓ
                Ft[iΓ][I,J] += 2.0*D[kΓ][K,L]*ν
                Ft[kΓ][K,L] += 2.0*D[iΓ][I,J]*ν
            end
            # K
            if jΓ == lΓ && iΓ == kΓ
                Ft[iΓ][I,K] -= Xik*D[jΓ][J,L]*ν
            end
        else
            if kΓ == lΓ && iΓ == jΓ
                Ft[iΓ][I,J] += 2.0*D[kΓ][K,L]*ν
            end
            if jΓ == lΓ && iΓ == kΓ
                Ft[iΓ][I,K] -= D[jΓ][J,L]*ν
            end
        end
    end
    for (Γ, l) in enumerate(irreplength)
        if l != 0
            for i in 1:l 
                otherF[Γ][i,i] += Ft[Γ][i,i] 
                for j = (i+1):l
                    otherF[Γ][i,j] += Ft[Γ][i,j] + Ft[Γ][j,i]
                    otherF[Γ][j,i] = otherF[Γ][i,j]
                end
            end
        end

    end
    otherF 
end