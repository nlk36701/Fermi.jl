using GaussianBasis
using Printf

function SRHF(Alg::A) where A <: SRHFAlgorithm
    ints = IntegralHelper{Float64}(eri_type = Chonky())
    SRHF(ints, Alg)
end

#future will include options to symmetrize molecule to highest possible point group
#before molecule is used. 
function block_SRHF(mol, basis_string)
    Fermi.HartreeFock.srhf_header()
    output("Ascertaining Molecular Point group...")
    ts = @elapsed begin
    mole, symtext = Fermi.Molecules.Symmetry.CharacterTables.symtext_from_mol(mol)
    end
    output("Done in {:10.5f} s", ts)
    println("The full molecular point-group is $(symtext.pg)")
    println("Running in $(symtext.pg) symmetry")
    display(symtext.ctab)
    bigboi = atom_to_string(mole)
    moles = Molecule(molstring = bigboi)
    #println("mole $mole") 
    #call the integral helper with the molecule, chonky integral, and basis string
    ints =  IntegralHelper{Float64}(molecule=moles, eri_type = Chonky(), basis = basis_string)
    
    output("Collecting necessary integrals...")
    t = @elapsed begin
        ints["S"]
        ints["T"]
        ints["V"]
        ints["ERI"]
    end
    output("Done in {:10.5f} s", t)

    #generate basis set object for SALC creation 
    bset = BasisSet(basis_string, mole)
    output("Generating SALCs...")
    
    #generate the SALCS
    #salcs, bigg, so_irrep, stevie_boy, rotated = GaussianBasis.SALCs.ProjectionOp(mole, bset, symtext)
    #bigg = Matrix(I, size(bigg)[1], size(bigg)[1])
    bigg = Matrix(I, 7,7)
    #oei
    TT = newaotoso(ints["T"], bigg)
    VV = newaotoso(ints["V"], bigg)
    SS = newaotoso(ints["S"], bigg)
    SERI = newaotoso_2(ints["ERI"], bigg, bset)
    Vnuc = Molecules.nuclear_repulsion(ints.molecule.atoms) 
    ndocc = try
        Int((ints.molecule.Nα + ints.molecule.Nβ)/2)
    catch InexactError
        throw(FermiException("Invalid number of electrons $(molecule.Nα + molecule.Nβ) for RHF method."))
    end
    #core guess, gwh not available yet
    println("Core guess...")
    
     
    
    molecule = ints.molecule
    output(Fermi.string_repr(molecule))
    #maxit = Options.get("scf_max_iter") 
    maxit = 100
    Etol  = Options.get("scf_e_conv")
    Dtol  = Options.get("scf_max_rms")
    Dtol = 1e-6
    C, A, H = newanti_symmetrize(SS, TT, VV)
    
    nvir = size(C,2) - ndocc
    nao = size(C,1)
    Vnuc = Molecules.nuclear_repulsion(molecule.atoms)

    output("Nuclear repulsion: {:15.10f}", Vnuc)
    output(" Number of AOs:                        {:5.0d}", nao)
    output(" Number of Doubly Occupied Orbitals:   {:5.0d}", ndocc)
    output(" Number of Virtual Spatial Orbitals:   {:5.0d}", nvir)
    
    Co = C[:, 1:ndocc]
    @tensor D[u,v] := Co[u,m]*Co[v,m]
    println("HERES THE D")
    println(D)
    D_old = deepcopy(D)
    eps = zeros(Float64,ndocc+nvir)

    # Build the inital Fock Matrix and diagonalize
    
    seri = permutedims(SERI, (1, 3, 2, 4))
    ite = 1
    E = 0.0
    ΔE = 1.0
    Drms = 1.0
    Drms = 1.0
    diis = false
    damp = 0.0 
    converged = false

    @tensoropt J[a,b] := 2 * D[c, d] * SERI[a, b, c, d]
    #@tensoropt K[a,b] := D[c, d] * seri[a, c, b, d]
    @tensoropt K[a,b] := D[c, d] * SERI[a, c, b, d]
    F = J - K + H
    println("The fock")
    println(F)
    F̃ = deepcopy(F)
    D̃ = deepcopy(D)
    N = length(D) # Number of elements in D (For RMS computation)
    output(" Guess Energy {:20.14f}", RHFEnergy(D,TT+VV,F))

    output("\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8} {:>8}", "E[RHF]", "ΔE", "Dᵣₘₛ", "t", "DIIS", "damp")
    output(repeat("-",80))
    t = @elapsed while ite ≤ maxit
        t_iter = @elapsed begin
            # Produce Ft
            Ft = A'*F*A

            # Get orbital energies and transformed coefficients
            eps, Ct = LinearAlgebra.eigen(Symmetric(Ft), sortby=x->x)
            println("Et $eps")
            # Reverse transformation to get MO coefficients
            C = A*Ct

            # Produce new Density Matrix
            Co = C[:,1:ndocc]
            @tensor D[u,v] = Co[u,m]*Co[v,m]

            # Build the Fock Matrix
            @tensoropt J[a,b] := 2 * D[c, d] * SERI[a, b, c, d]
            @tensoropt K[a,b] := D[c, d] * SERI[a, c, b, d]
            F = J - K + H
            #println("New F")
            #println(F)
            evals, evecs = eigen(F)
            Eelec = newRHFEnergy(D, TT + VV, F)

            # Compute Energy
            Enew = Eelec + Vnuc
            ΔD = D - D_old
            Drms = √(sum(ΔD.^2) / N)

            # Compute Energy Change
            ΔE = Enew - E
            E = Enew
            D_old .= D
        end
        output("    {:<3} {:>15.10f} {:>11.3e} {:>11.3e} {:>8.4f} {:>8}    {:5.2f}", ite, E, ΔE, Drms, t_iter, diis, damp)
        ite += 1

        if (abs(ΔE) < Etol) & (Drms < Dtol) & (ite > 5)
            converged = true
            MP2 = true
            if MP2
                println("HELL YEAH MP2 TIME!")
                #MP2time(evals, C, SERI, ndocc, nao)
                MP2time(eps, Ct, SERI, ndocc, nao)
                #MP2time(evals, C, ints["ERI"], ndocc, nao)
            end 
            break
        end
    end
end
function MP2time(e, C, SERI, ndocc, nao)
    println("These are the evals")
    println(e)
    #first, transform ao -> mo
    G = aotomo(SERI, C, nao)
    E = 0.0
    for i = 1:ndocc
        for j = 1:ndocc
            for a = ndocc + 1:nao
                for b = ndocc + 1:nao
                    #println("$i $j $a $b")
                    E += G[i,j,a,b] * ( 2*G[i,j,a,b] - G[i,j,b,a] ) / (e[i] + e[j] - e[a] - e[b])
                end
            end
        end
    end
    println("MP2 correlation energy: $E")
end
function aotomo(SERI, C, nao)
    G = zeros(nao, nao, nao, nao)
    @tensoropt G[i,j,k,l] =  SERI[μ, ν, ρ, σ]*C[μ, i]*C[ν, j]*C[ρ, k]*C[σ, l]
    return G
end


function newaotoso(OEI, bigg)
        @tensor soei[p,q] := OEI[u,v] * bigg[u,p] * bigg[v,q]
    return soei
end

#function that symmetry adapts TEIs
function newaotoso_2(eri, bigg, bset)
    SERI = zeros(bset.nbas, bset.nbas, bset.nbas, bset.nbas)
    @tensoropt SERI[i,j,k,l] =  eri[μ, ν, ρ, σ]*bigg[μ, i]*bigg[ν, j]*bigg[ρ, k]*bigg[σ, l]
    return SERI
end
function newanti_symmetrize(S, T, V)
    A = real(S^(-1/2))
    F = T + V
    Ft = A * F * A'
    #energies, ct = LinearAlgebra.eigen((ft), sortby=x->x)
    energies, Ct = LinearAlgebra.eigen(Symmetric(Ft), sortby=x->x)
    println("energies $energies")
    C = A*Ct
    H = T + V
    return C, A, H
end
function newRHFEnergy(D, H, F)
    return sum(D .* (H .+ F))
end