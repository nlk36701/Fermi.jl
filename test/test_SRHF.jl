@reset 
using Molecules
using GaussianBasis
path = joinpath(@__DIR__, "xyz/water.xyz")
mol = Molecules.parse_file(path)
#mol, symtext = Molecules.Symmetry.CharacterTables.symtext_from_mol(mol)
println(mol)
@set scf_guess core
@set basis "sto-3g"
@molecule mol
wfn = @energy mol => srhf