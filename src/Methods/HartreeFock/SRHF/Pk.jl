using TensorOperations
using LinearAlgebra
using Formatting
import Base: show


function fast_SRHF(mol, basis_string)
    Fermi.HartreeFock.srhf_header()
    output("Ascertaining Molecular Point group...")
    ts = @elapsed begin
    mole, symtext = Fermi.Molecules.Symmetry.CharacterTables.symtext_from_mol(mol)
    end
    output("Done in {:10.5f} s", ts)
    println("The full molecular point-group is $(symtext.pg)")
    println("Running in $(symtext.pg) because this program is neat-o")
    
    bigboi = atom_to_string(mole)
    moles = Molecule(molstring = bigboi)
    println("mole $mole") 
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
    println("First half transform") 
    @time firsthalfSERI = firststaotoso_2(eri, bigg, bset)
    
    println("Second half transform") 
    @time SERI = secondaotoso_2(firsthalfSERI, bigg, bset)
    
    Vnuc = Molecules.nuclear_repulsion(ints.molecule.atoms)
     
    pk_symoffset, pk_pairs = compute_offset(irreplength) 
    
    println("pk_symoffset $pk_symoffset")
    #println("pk_pairs $pk_pairs")
    #println("sum pk_offset")
    pk_size = sum(pk_pairs)
    
    #initialize PK array 
    pk = Vector{Float64}(undef, index2(pk_size, pk_size))
    println("pk $pk")

    #create vector of integral indices using 8-fold perm-sym
    uniqueindex = unidx(bset)
    ind = ERI_indicies(irreplength)
    
    for (id, idx) in enumerate(uniqueindex)
        i = idx[1]
        j = idx[2]
        k = idx[3]
        l = idx[4]
        ii, is = newindex(i, irreplength)
        jj, js = newindex(j, irreplength)
        kk, ks = newindex(k, irreplength)
        ll, ls = newindex(l, irreplength)
        value = SERI[id]
        #J
        if is == js && ks == ls
            bra = index2(ii, jj) + pk_symoffset[is]
            ket = index2(kk, ll) + pk_symoffset[ks]
            #pk_symoffset corrects for the symmetry offset in the pk vectors
            braket = index2(bra, ket)
            pk[braket] += value
            #println("braket $braket")
            # K/2, 2nd sort
            if ii != jj && kk != ll
                if is == ls && js == ks
                    bra = index2(ii, ll) + pk_symoffset[is]
                    ket = index2(jj, kk) + pk_symoffset[js]
                    braket = index2(bra, ket)
                    pk[braket] += value
                    if ii == ll || jj == kk
                        pk[braket] -= 0.5 * value
                    else
                        pk[braket] -= 0.25 * value
                    end
                end
            end
        end
        # K/2, 1st sort
        if is == ks && js == ls
            bra = index2(ii, kk) + pk_symoffset[is]
            ket = index2(jj, ll) + pk_symoffset[js]
            braket = index2(bra, ket)
            if ii == kk || jj == ll
                pk[braket] -= 0.5 * value
            else
                pk[braket] -= 0.25 * value
            end
        end

        #println("value $value")
    end
    for i = 1:pk_size
        pk[index2(i, i)] *= 0.5
    end
    println("PK NOW")
    println(pk)
end

#compute number of indices in upper triangle of each pair that transforms as a1

function compute_offset(irreplength)
    pk_symoffset = []
    pk_pairs = []
    symoffset = 0
    for irr in irreplength
        push!(pk_symoffset, symoffset)
        pairs = Int64(irr * (irr + 1) / 2)
        symoffset += pairs
        push!(pk_pairs, pairs)
    end
    return pk_symoffset, pk_pairs
end

#assistant function to above function,
#evaluated at fock build
function index2(i, j)
    if i < j
        return (j * (j + 1)) >> 1 + i
    else
        return (i * (i + 1)) >> 1 + j
    end
end

function indx4(i,j,k,l)
    return index2(index2(i,j), index(k,l))
end
