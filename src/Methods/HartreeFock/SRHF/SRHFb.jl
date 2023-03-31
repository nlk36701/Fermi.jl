
function degen_SRHF(mol, basis_string)
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
    #for (ir, irrep) in enumerate(symtext.ctab.class_orders)
    for (ir, irrep) in enumerate(symtext.ctab.irreps)
        println("ir $ir $irrep")
    end
    irreplength = []
    for (i, s) in enumerate(salcs)
        push!(irreplength, size(salcs[i].lcao)[2])
    end
end