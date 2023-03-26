using GaussianBasis
using Printf

function SRHF(Alg::A) where A <: SRHFAlgorithm
    ints = IntegralHelper{Float64}(eri_type = Chonky())
    SRHF(ints, Alg)
end

function SRHF(mol::Molecule, Alg::A) where A <: SRHFAlgorithm
    SRHF(IntegralHelper{Float64}(molecule=mol, eri_type = Chonky()), Alg)
end
function bd_from_full(b, irreplength)
    dabtda = []
    c = 1
    for i in irreplength
        start = c
        stop = c + i - 1
        if i == 0
            push!(dabtda, [])
        else
            smallboi = zeros(i, i)
            smallboi = b[start:stop, start:stop]
            push!(dabtda, smallboi)
            c += i
        end
    end
    return dabtda
end

function full_mat(b)
    fshape = 0
    for (i, block) in enumerate(b)
        fshape += size(block)[1]
    end
    fullboi = zeros(fshape,fshape)
    c = 1
    for (i, block) in enumerate(b)
        s = size(block)[1]
        start = c
        stop = c + s -1

        fullboi[start:stop, start:stop] = block
        c += s
    end
    return fullboi
end

#Function that symmetry adapts OEIs
function aotoso(OEI, salcs)
    B = []
    for (z, salc) in enumerate(salcs)
        if length(salcs[z].lcao) == 0
            push!(B, [])
        else
            @tensor soei[p,q] := OEI[u,v] * salcs[z].lcao[u,p] * salcs[z].lcao[v,q]
            push!(B, soei)
        end
    end
    return B
end

#function that symmetry adapts TEIs
function aotoso_2(eri, bigg, bset)
    SERI = zeros(bset.nbas, bset.nbas, bset.nbas, bset.nbas)
    @tensoropt SERI[i,j,k,l] =  eri[μ, ν, ρ, σ]*bigg[μ, i]*bigg[ν, j]*bigg[ρ, k]*bigg[σ, l]
    return SERI
end

function firststaotoso_2(eri, bigg, bset)
    firsthalfSERI = zeros(bset.nbas, bset.nbas, bset.nbas, bset.nbas)
    @tensoropt firsthalfSERI[p,q,k,l] =  eri[p, q, r, s]*bigg[r, k]*bigg[s, l]
    return firsthalfSERI
end

function secondaotoso_2(firsthalfSERI, bigg, bset)
    @tensoropt SERI[i,j,k,l] :=  firsthalfSERI[p, q, k, l]*bigg[p, i]*bigg[q, j]
    return SERI
end

#function that enforces the overlap diagonal is the identity
#see Susi's comments in PSI4 code, it isn't necessary, which is
#weird
function normalize_antisymmetrize(S)
    A = []
    normlists = []
    for (z, soei) in enumerate(S)
        if length(soei) == 0
            #apend empty array per convention
            push!(A, [])
        else
            normlist = []
            for i in 1:size(soei)[1]
                norm1 = 1/sqrt(soei[i,i])
                push!(normlist, norm1)
                for j in 1:size(soei)[1]
                    norm2 = 1/sqrt(soei[j,j])
                    soei[i,j] = soei[i,j] * norm1 * norm2
                end
            end
            eigval = eigvals(soei)
            eigvec = eigvecs(soei)
            Us = copy(eigvec)
            for k in 1:length(eigval)
                Us[:,k] = eigvec[:,k] * 1.0/sqrt(eigval[k])
            end
            for k in 1:length(eigval)
                Us[k,:] = Us[k,:] * normlist[k]
            end
            anti = Us'eigvec
            push!(A, anti)
        end
    end
    return A
end
#creates indices for sparse Fock build
#8-fold permutational symmetry baked in
function unidx(bset)
    d = []
    N = bset.nbas

    for i = 1:N 
        for j = i:N 
            for k = 1:N 
                for l = k:N 
                    if idx2(i,j) < idx2(k,l)
                        continue
                    end
                    push!(d, [i, j, k, l])
                end
            end
        end
    end
    return d
end

#assistant function to above function,
#evaluated at fock build
function idx2(i, j)
    if i < j
        return (j * (j + 1)) >> 1 + i
    else
        return (i * (i + 1)) >> 1 + j
    end
end

function index4(i,j,k,l)
    return index2(index2(i,j), index2(k,l))
end

#reduction coefficient for TSIR 
#if result/group order !=0, contains A1

function containsTSIR(i, j, k, l, symtext, so_irrep)
    ct = symtext.ctab
    class_map = symtext.class_map
    ir = so_irrep[i] 
    jr = so_irrep[j] 
    kr = so_irrep[k] 
    lr = so_irrep[l]
    group_order = sum(ct.class_orders)
    #println("Grup order $group_order")
    irc = ct.characters[ir,:]
    jrc = ct.characters[jr,:]
    krc = ct.characters[kr,:]
    lrc = ct.characters[lr,:]
    DirectProduct = ct.class_orders .* irc .* jrc .* lrc .* krc
    if sum(DirectProduct) > 0
        return(i, j, k, l)
    else
        #println("Does not contain TSIR $i $j $k $l)")
    end
end

function containsTSIR_new(i, j, k, l, symtext)
    ct = symtext.ctab
    ir = ct.characters[i,:]
    jr = ct.characters[j,:]
    kr = ct.characters[k,:]
    lr = ct.characters[l,:]
    DirectProduct = ct.class_orders .*ir .*jr .* lr .* kr
    if sum(DirectProduct) > 0
        return i, j, k, l
    end
end
#creates dictionary of irreps Gamma and their associated
#number of occupied orbitals
function countmemb1(y)
    d = Dict{Int, Int}()
    for val in y
        if isnan(val)
            continue
        end
        if val in keys(d)
            d[val] += 1
        else
            d[val] = 1
        end
    end
    return d
end


#Create antisymmetrizer, initiate Fock
function anti_symmetrize(S, Tint, V, ndocc, symtext)
     
    #AAA = normalize_antisymmetrize(S)
    #println("AAA $AAA")
    AA = []
    C = []
    CT = []
    Energies = []
    cont = []
    doccirrep = []
    F = []
    H = []
    for (i,s) in enumerate(S)
        if length(s) == 0
            push!(AA, [])
            push!(C, [])
            push!(CT, [])
            push!(Energies, [])
            push!(F, [])
            push!(H, [])
        else
            a = real(s^(-1/2))
            #println(a)
            f = Tint[i] + V[i]
            ft = a * f * a'
            #println(ft)
            energies, ct = LinearAlgebra.eigen((ft), sortby=x->x)
            #energies, ct = LinearAlgebra.eigen((ft))
            #println("energies $energies")
            energies, ct = LinearAlgebra.eigen(Symmetric(ft), sortby=x->x)
            #println("energies $energies")
            #println("orbitals $ct")
            #println("real energies $realenergies")
            #realenergies, ct = LinearAlgebra.eigen(Hermitian(ft))
            #println("real energies $realenergies")
            c = a*ct
            push!(AA, a)
            push!(C, c)
            push!(CT, ct)
            push!(H, Tint[i] + V[i])
            doccirrep = vcat(doccirrep, fill(i, size(s)[1]))
            cont = vcat(cont, energies...)
            push!(Energies, energies)

        end
        #println("Gottem")
    end
    #println("What the hell is going on????")
    order = sortperm(cont)
    doccirrep = doccirrep[order]
    doccirrep = doccirrep[1:ndocc]
    doccdict = countmemb1(doccirrep)
    return C, AA, doccdict, H
end

#move to Molecules. This function only exists because Fermi's molecule object can't handle scientific notation in string format
function atom_to_string(mol)
    bigboi = " "
    for atom in mol
        symbol = Fermi.PhysicalConstants.num_atom[string(atom.Z)]
        coords = atom.xyz
        piece = string(symbol) * " " * Printf.@sprintf("%.20f",coords[1]) * " " * Printf.@sprintf("%.20f", coords[2]) * " " * Printf.@sprintf("%.20f", coords[3]) * "\n"
        bigboi *= piece 
    end
    return bigboi
end


#Initialize empty Fock matrix for sparse build
function smallFock(irreplength)
    F = []
    for x in irreplength
        if x == 0
            push!(F, [])
        else
            push!(F, zeros(x, x))
        end
    end
    return F
end


#function direct_SERI(bset, bigg, SERI)
#    println("hi")
#    #println(GaussianBasis.ERI_2e4c(bset, 1, 1, 1, 4)
#    
#    #integral = 0.0
#    #for p in 1:size(bigg)[1]
#    #    for q in 1:size(bigg)[1]
#    #        for r in 1:size(bigg)[1]
#    #            for s in 1:size(bigg)[1]
#
#    #                integral = GaussianBasis.ERI_2e4c(bset, p, q, r, s)
#
#end
#function actual_SRHF(ints::IntegralHelper{Float64}, Alg::A) where A <: SRHFAlgorithm

#Symmetry adapted RHF code as a Sub-Method of Fermi.Hartree--Fock
#arguments are currently a molecule object and a basis set string.
#THIS MODULE WILL ROTATE YOUR MOLECULE, AND NOT ROTATE IT BACK

#future will include options to symmetrize molecule to highest possible point group
#before molecule is used. 
function actual_SRHF(mol, basis_string)
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
    salcs, bigg, so_irrep, stevie_boy, rotated = GaussianBasis.SALCs.ProjectionOp(mole, bset, symtext)
    
    irreplength = []
    for (i, s) in enumerate(salcs)
        push!(irreplength, size(salcs[i].lcao)[2])
    end
    #generate symmetry-adapted one-electron ints
    Soverlap   = aotoso(ints["S"], salcs)
    Skinetic   = aotoso(ints["T"], salcs)
    Spotential = aotoso(ints["V"], salcs)

    #generate symmetry-adapted two-electron ints
    eri = ints["ERI"]
    #@time SERI = aotoso_2(eri, bigg, bset)
    println("First half transform") 
    @time firsthalfSERI = firststaotoso_2(eri, bigg, bset)
    
    println("Second half transform") 
    @time SERI = secondaotoso_2(firsthalfSERI, bigg, bset)
    
    #BS_SERI = blocksparse_SERI(eri, bigg, irreplength, symtext)
    #println("Length!")
    #println(length(BS_SERI)) 
    #compute nuclear repulsion
    Vnuc = Molecules.nuclear_repulsion(ints.molecule.atoms) 
    
    #grab unique integrals based on 8-fold permutational symmetry
    
    #uniqueindex = unidx(bset)
    #println("unique index")
    #println(uniqueindex[5]) 
    #println(uniqueindex[5][1]) 
    #println(uniqueindex[5][2]) 
    #println(uniqueindex[5][3]) 
    #println(uniqueindex[5][4]) 
    
    #need to filter out non-TSIR ERIs from sparse index
    #symdex takes uniqueindex and further filters the integrals
    
    #dont use for now, just do regular sparse
    #symdex = []
    #for indx in uniqueindex
    #    i, j, k, l = indx[1], indx[2], indx[3], indx[4]
    #    index = containsTSIR(i,j,k,l,symtext, so_irrep)
    #    if index != nothing
    #        push!(symdex, index)
    #    end
    #end
    
    #creates sparse ERI based on eight-fold permutational symmetry and TSIR
    #symsparse contains the integrals, symdex are the idxs to loop over
    
    #symsparse = chonk_to_sparse(SERI, symdex)
    
    #symmetry allowed sparse eri built. Now built antisymmetrizer, core guess
    ndocc = try
        Int((ints.molecule.Nα + ints.molecule.Nβ)/2)
    catch InexactError
        throw(FermiException("Invalid number of electrons $(molecule.Nα + molecule.Nβ) for RHF method."))
    end
    #core guess, gwh not available yet
    println("Core guess...")
    C, AA, doccdict, smallH = anti_symmetrize(Soverlap, Skinetic, Spotential, ndocc, symtext)
    ### = imported from SRHF function 
    
    ### 
    molecule = ints.molecule
    output(Fermi.string_repr(molecule))
    #maxit = Options.get("scf_max_iter") 
    maxit = 100
    Etol  = Options.get("scf_e_conv")
    Dtol  = Options.get("scf_max_rms")
    Dtol = 1e-6
    do_diis = Options.get("diis")
    oda = Options.get("oda")
    oda_cutoff = Options.get("oda_cutoff")
    oda_shutoff = Options.get("oda_shutoff")

    # Variables that will get updated iteration-to-iteration
    ite = 1
    E = 0.0
    ΔE = 1.0
    Drms = 1.0
    Drms = 1.0
    diis = false
    damp = 0.0 
    converged = false
     
    # Build a diis_manager, if needed

    #for symmetry, initialize a diis manager for each irrep?
    #dm object in diis_set for irreps that have a subspace of rank 0 is just [] placeholder
    if do_diis
        diis_set = []
        for Gamma in irreplength
            if Gamma == 0
                push!(diis_set, [])
            else
                DM = Fermi.DIIS.DIISManager{Float64,Float64}(size=Options.get("ndiis"))
                #println("DM")
                #println(DM)
                push!(diis_set, DM)
                #diis_start = Options.get("diis_start")
                #println("diis_start")
                #println(diis_start)
            end
        end
    end
    diis_start = Options.get("diis_start")
    if do_diis
        DM = Fermi.DIIS.DIISManager{Float64,Float64}(size=Options.get("ndiis"))
        diis_start = Options.get("diis_start")
    end
    
    #compute nvir and nao based on irrep

    nvir = 0
    nao = 0
    for c in C
        if length(c) == 0
            continue
        else
            nvir += size(c,2) 
            nao += size(c,1)
        end
    end
    nvir -= ndocc
    Vnuc = Molecules.nuclear_repulsion(molecule.atoms)

    output("Nuclear repulsion: {:15.10f}", Vnuc)
    output(" Number of AOs:                        {:5.0d}", nao)
    output(" Number of Doubly Occupied Orbitals:   {:5.0d}", ndocc)
    output(" Number of Virtual Spatial Orbitals:   {:5.0d}", nvir)
    ### 
    #construct the initial density matrix  
    smallD = BuildD(C, doccdict) 
    smallD_old = deepcopy(smallD)
    
    otherF = smallFock(irreplength)
    smallF = chonkyfock(otherF, smallH, smallD, SERI, irreplength, symtext, so_irrep)
    #println(smallF)
    #smallF = smallfock!(otherF, smallH, smallD, symsparse, symdex, bset.nbas, so_irrep, irreplength)
    #blocksparsefock(otherF, smallH, smallD, BS_SERI, irreplength, symtext, so_irrep)
    
    F̃ = deepcopy(smallF)
    D̃ = deepcopy(smallD)
    N = irreplength # Number of elements in D (For RMS computation)
    doccirrep = nothing
    sorted_e = nothing
    output(" Guess Energy {:20.14f}", rhf_energy_new(smallF, smallD, smallH))
    output("\n Iter.   {:>15} {:>10} {:>10} {:>8} {:>8} {:>8}", "E[RHF]", "ΔE", "Dᵣₘₛ", "t", "DIIS", "damp")
    output(repeat("-",80))
    t = @elapsed while ite ≤ maxit
        t_iter = @elapsed begin
            # Produce Ft
            # Get orbital energies and transformed coefficients
            # Reverse transformation to get MO coefficients
            Ft, Ct, Et = backt(smallF, AA)
            #println("Ft $Ft")
            otherF = smallFock(irreplength)
            # Produce new Density Matrix
            doccdict, doccirrep, sorted_e = orderenergy(Ct, Et, ndocc) 
            D = BuildD(Ct, doccdict)
            
            smallD = D 
            # Build the Fock Matrix
            
            #smallF = smallfock!(otherF, smallH, smallD, symsparse, symdex, bset.nbas, so_irrep, irreplength)
            
            smallF = chonkyfock(otherF, smallH, smallD, SERI, irreplength, symtext, so_irrep)
            #smallF = blocksparsefock(otherF, smallH, smallD, BS_SERI, irreplength, symtext, so_irrep)
            
            ΓEelec = rhf_energy_new(smallF, smallD, smallH) 
            # Compute Energy
            ΓEnew = ΓEelec + Vnuc
            ## Store vectors for DIIS
            ###
            aaa = full_mat(AA)
            SSS = full_mat(Soverlap)
            FFF = full_mat(smallF)
            DDD = full_mat(smallD)

            if do_diis
                err = transpose(aaa)*(FFF*DDD*SSS - SSS*DDD*FFF)*aaa
                push!(DM, FFF, err)
            end
            #######################################################
            
            #No dm object in diis_set for irreps that have a subspace of rank 1 or less?
            #if do_diis
            #    for (i, dm) in enumerate(diis_set)
            #        if irreplength[i] == 0
            #            continue
            #            #println("Again, no dice!")
            #        else
            #            Gamma_err =ΓDice(AA, smallF, smallD, Soverlap)
            #            #println("Γ Error")
            #            #println(Gamma_err)
            #            push!(dm, smallF[i], Gamma_err[i])
            #            #println("dm $i")
            #            #println(dm)
            #        end
            #    end
            #end
            #do_diis = false
            if do_diis && ite > diis_start
                diis = true
                F = Fermi.DIIS.extrapolate(DM)
                smallF = bd_from_full(F, irreplength) 
            end
            #if do_diis && ite > diis_start
            #    newF = []
            #    diis = true
            #    for (x, i) in enumerate(irreplength)
            #        #println("This is i $i and x $x")
            #        if i == 0
            #            push!(newF, [])
            #        elseif i == 1
            #            push!(newF, smallF[x])
            #        else
            #            #println("Is it here?")
            #            smallfock = Fermi.DIIS.extrapolate(diis_set[x])
            #            push!(newF, smallfock)
            #        end
            #    end
            #    #println("newF")
            #    #println(newF)
            #    #smallF = newF
            #end
            
            ###
            # Branch for ODA vs DIIS convergence aids
            #diis = false
            #damp = 0.0
            # Use ODA damping?
            #######################################################
            #if oda && Drms > oda_cutoff && ite < oda_shutoff
            #    diis = false
            #    dD = D - D̃
            #    s = tr(F̃ * dD)
            #    c = tr((F - F̃) * (dD))
            #    if c <= -s/(2*c)
            #        λ = 1.0
            #    else
            #        λ = -s/(2*c)
            #    end
            #    F̃ .= (1-λ)*F̃ + λ*F
            #    D̃ .= (1-λ)*D̃ + λ*D
            #    damp = 1-λ
            #    F .= F̃
            
            ## Or Use DIIS?
            ######################################################

            # Compute the Density RMS
            ΓΔD = smallD - smallD_old
            ΓDrms = 0
            for (i, d) in enumerate(ΓΔD)
                if length(d) != 0 
                    ΓDrms += √(sum(d.^2) / irreplength[i])
                end
            end
            # Compute Energy Change
            ΔE = ΓEnew - E
            E = ΓEnew
            smallD_old .= smallD
        end
        output("    {:<3} {:>15.10f} {:>11.3e} {:>11.3e} {:>8.4f} {:>8}    {:5.2f}", ite, E, ΔE, ΓDrms, t_iter, diis, damp)
        ite += 1

        #if (abs(ΔE) < Etol) & (Drms < Dtol) & (ite > 5)
        if (abs(ΔE) < Etol) & (ΓDrms < Dtol) & (ite > 5)
            converged = true
            break
        end
    end

    output(repeat("-",80))
    output("    SRHF done in {:>5.4f}s", t)
    output("    @Final SRHF Energy     {:>20.12f} Eₕ", E)
    output("\n   • Orbitals Summary",)
    output("\n {:>12}  {:>8} {:>15}   {:>10}", "Orbital(s)", "Irrep", "Energy", "Occupancy")
    i = 1
    while i < length(sorted_e) + 1
        if salcs[doccirrep[i]].irrep == "E"
            output(" {:>12} {:>8}  {:> 15.10f}   {:>7}", (i,  i + 1), salcs[doccirrep[i]].irrep, sorted_e[i], (i ≤ ndocc ? "↿⇂↿⇂" : ""))
            i += 2 
        else 
            output(" {:>10} {:>10}  {:> 15.10f}   {:>6}", i, salcs[doccirrep[i]].irrep, sorted_e[i], (i ≤ ndocc ? "↿⇂" : ""))
            i += 1 
        end
    end
    output("")
    if converged
        output("   ✔  SCF Equations converged 😄")
    else
        output("❗ SCF Equations did not converge in {:>5} iterations ❗", maxit)
    end
    output(repeat("-",80))

    #Orbitals = SRHFOrbitals(molecule, ints.basis, eps, E, C)

    #return SRHF(molecule, E, ndocc, nvir, Orbitals, ΔE, Drms)
end

function ERI_indicies(irreplength)
    vec = []
    start = 1
    total = 0
    for (a, i) in enumerate(irreplength)
        if i == 0
            push!(vec, nothing)
            continue
        elseif a > 1
            push!(vec, total+1:total+i)
        else
            push!(vec, start:i)
        end
        total += i
        start = i
    end
    return vec 
end

function blocksparse_SERI(ERI, bigg, irreplength, symtext)
    BS_SERI = []
    tsir = 0
    nontsir = 0
    total = 0
    ind = ERI_indicies(irreplength)
    for (i, p) in enumerate(irreplength)
        if p != 0
            for (j, q) in enumerate(irreplength)
                if q != 0
                    for (k, r) in enumerate(irreplength)
                        if r != 0
                            for (l, s) in enumerate(irreplength)
                                if s != 0
                                    doit = containsTSIR_new(i,j,k,l,symtext)
                                    if doit != nothing
                                        tsir += 1
                                        if p == q && r == s 
                                            ibigg = bigg[:, ind[i]]
                                            jbigg = bigg[:, ind[j]]
                                            kbigg = bigg[:, ind[k]]
                                            lbigg = bigg[:, ind[l]]
                                            @tensoropt SERI[I,J,K,L] :=  ERI[μ, ν, ρ, σ]*ibigg[μ, I]*jbigg[ν, J]*kbigg[ρ, K]*lbigg[σ, L]
                                            #println(SERI)
                                            total += length(SERI)
                                            push!(BS_SERI, [SERI, [i, j, k, l]])
                                        end
                                    else
                                        nontsir += 1
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    println("Total length")
    println(total)
    println(nontsir)
    return BS_SERI
end

function blocksparsefock(otherF, smallH, smallD, BS_SERI, irreplength, symtext, so_irrep)
    otherF += smallH
    #loop over nonzero integral blocks
    for ay = 1:length(BS_SERI)
        index = BS_SERI[ay][2]
        seri = BS_SERI[ay][1]
        i, j, k, l = index[1], index[2], index[3], index[4]
        #ind = ERI_indicies(irreplength)
        println("$i $j $k $l")
        IKJL = permutedims(seri, (1, 3, 2, 4)) 
        println("bad")
        println(IKJL)
        #D = smallD[k]
        #@tensoropt first[a,b] := 2 * D[c, d] * seri[a, b, c, d]
        #@tensoropt second[a,b] := D[c,d] * IKJL[a, c, b, d] 
        #otherF[i] +=  first #- second

        #if i == j && k == l
        #    D = smallD[k]
        #    @tensoropt first[a,b] := 2 * D[c, d] * seri[a, b, c, d]
        #    otherF[i] += first
        #    #println("otherF")
        #    #println(otherF) 
        #elseif i == k && j == l
        #    D = smallD[j]
        #    #println("D, wrong")
        #    #println(D)
        #    #seri = permutedims(seri, (1, 3, 2, 4))
        #    @tensoropt second[a,b] := D[c,d]* seri[a, c, b, d] 
        #    otherF[i] -= second
        #    #println("second")
        #    #println(second)
        #end
    end
    println("WRONG")
    println(otherF)
    return otherF 
end
function chonkyfock(smallF, smallH, smallD, eri, irreplength, symtext, so_irrep)
    #loop over the indices and check for TSIR in direct product
    #continue if != TSIR
    count = 0
    closeone = 0
    smallF = smallFock(irreplength) .+ smallH
    ind = ERI_indicies(irreplength)
    for h in 1:length(irreplength)
        if ind[h] != nothing
            for g in 1:length(irreplength)
                if ind[g] != nothing
                    I = ind[h]
                    J = ind[h]
                    K = ind[g]
                    L = ind[g]
                    IJKL = eri[I, J, K, L]
                    #println(IJKL)
                    IKJL = eri[I, K, J, L]
                    #println("good")
                    #println("$I $J $K $L")
                    #println(IKJL)
                    #println(IJKL)
                    D = smallD[g]
                    @tensoropt small[a,b] := 2 * D[c, d] * IJKL[a, b, c, d]
                    @tensoropt big[a,b] := D[c, d] * IKJL[a, c, b, d]
                    smallF[h] += small - big
                    #println("smallF")
                    #println(smallF)
                    #println("big")
                    #println(big)
                end
            end
        end
    end
    return smallF
    #return nothing
end

function ΓDice(AA, smallF, smallD, Soverlap)
    err = []
    for g in 1:length(AA)
        if length(AA[g]) == 0
            push!(err, 0)
        else
            error = transpose(AA[g])*(smallF[g]*smallD[g]*Soverlap[g] - Soverlap[g]*smallD[g]*smallF[g])*(AA[g])
            push!(err, error)
        end
    end
    return err
end

function orderenergy(Ct, Et, ndocc)
    doccirrep = []
    cont = []
    for (i, e) in enumerate(Et)
        cont = vcat(cont, e...)
        doccirrep = vcat(doccirrep, fill(i, size(e)))
    end
    order = sortperm(cont)
    sorted_e = cont[order]

    doccirrep = doccirrep[order]
    doccirrep_occ = doccirrep[1:ndocc]
    doccdict = countmemb1(doccirrep_occ)
    return doccdict, doccirrep, sorted_e 
end
function backt(F, AA)
    Ft = []
    Ct = []
    Et = []
    for (i, f) in enumerate(F)
        if size(f)[1] == 0
            push!(Ft, [])
            push!(Ct, [])
            push!(Et, [])
        else
            a = AA[i]
            ft = a'*f*a
            eps, ct = eigen(ft)
            ct = a*ct
            push!(Ft, ft)
            push!(Ct, ct)
            push!(Et, eps)
        end
    end
    return Ft, Ct, Et
end

function rhf_energy_new(F, D, H)
    total = 0
    for (i, f) in enumerate(F)
        if length(f) != 0
            total += sum(D[i] .* (H[i] .+ F[i]))
        end
    end
    return total
end
#
function BuildD(C, doccdict)
    D = []
    for (i, c) in enumerate(C)
        #check to see if there are any irreps that do not contain a set
        #of orbitals. Not the case in cc-pvdz
        if length(c) == 0
            push!(D, [])
        else
            if haskey(doccdict, i)
                val = doccdict[i]
                co = c[:, 1:val]
                @tensor d[u,v] := co[u,m]*co[v,m]
                push!(D, d)
            else
                d = zeros(size(c)[1],size(c)[1])
                push!(D, d)
                continue
            end
        end
    end
    return D
end
function chonk_to_sparse(SERI, symdex)
    symsparse = []
    for i in symdex
        idx = [i...]
        push!(symsparse, SERI[i...])
    end
    return symsparse
end
