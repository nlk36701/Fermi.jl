"""
    Fermi.PhysicalConstants

Module for storing physical constants and conversion factors used in computations.

# Functions:

    atomic_number   Given an element symbol, return the atomic number.
"""
module PhysicalConstants

export atomic_number, num_core_electrons

"""
    Fermi.PhysicalConstants.AvogadroNumber

Float64 object with the Avogadro Number in mol^-1
Source: CODATA 2018 [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?na)
"""
const AvogadroNumber = 6.02214076e23 

"""
    Fermi.PhysicalConstants.kcal_to_kj

Float64 object with the convesion factor from kcal to kJ
Source: [NIST SP 1038 - 2006 page 12](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication1038.pdf)
"""
const kcal_to_kj = 4.184

"""
    Fermi.PhysicalConstants.bohr_to_angstrom

Float64 object with the conversion factor from Bohr to Angstrom.
Source: CODATA 2018 [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0)
"""
# Note. Psi4 constant is
#const bohr_to_angstrom = 0.529177210670000
const bohr_to_angstrom = 0.529177210903

"""
    Fermi.PhysicalConstants.hartree_to_kjmol

Float64 object with the conversion factor from Hartree energy to kJ/mol
Source: CODATA 2018 [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?hr)
"""
const hartree_to_kjmol = 4.3597447222071e-18 * AvogadroNumber / 1000 # (E_H in J) * (Na in 1/mol) / (1000 J/kJ) 

"""
    Fermi.PhysicalConstants.hartree_to_kcalmol

Float64 object with the conversion factor from Hartree energy to kcal/mol.
Source: CODATA 2018 [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0)
"""
const hartree_to_kcalmol = hartree_to_kjmol / kcal_to_kj

"""
    Fermi.PhysicalConstants.hartree_to_cm

Float64 object with the conversion factor from Hartree energy to cm^-1 (wavenumber)
Source: CODATA 2018 [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?hrminv)
"""
const hartree_to_cm = 2.1947463136320e7 / 100 # Value in 1/meter * (1/100) meter/cm

"""
    Fermi.PhysicalConstants.hartree_to_eV

Float64 object with the conversion factor from Hartree energy to eV.
Source: CODATA 2018 [NIST](https://physics.nist.gov/cgi-bin/cuu/Value?hrev)
"""
const hartree_to_eV = 27.21138505

"""
    Fermi.PhysicalConstants.atomic_number(atom::String)

Given a string with the element symbol, return its atomic number
"""
function atomic_number(atom::String)
    
    if haskey(atom_num, atom)
        return atom_num[atom]
    else
        throw(KeyError("atomic symbol $atom not defined."))
    end
end

function num_core_electrons(atom::String)
    
    Z = atomic_number(atom)

    for k in keys(core_elec_perrow)
        if Z in k
            return core_elec_perrow[k]
        end
    end
end

const atom_num = Dict(
    "H" => 1,
    "He" => 2,
    "Li" => 3,
    "Be" => 4,
    "B" => 5,
    "C" => 6,
    "N" => 7,
    "O" => 8,
    "F" => 9,
    "Ne" => 10,
    "Na" => 11,
    "Mg" => 12,
    "Al" => 13,
    "Si" => 14,
    "P" => 15,
    "S" => 16,
    "Cl" => 17,
    "Ar" => 18,
    "K" => 19,
    "Ca" => 20,
    "Sc" => 21,
    "Ti" => 22,
    "V" => 23,
    "Cr" => 24,
    "Mn" => 25,
    "Fe" => 26,
    "Co" => 27,
    "Ni" => 28,
    "Cu" => 29,
    "Zn" => 30,
    "Ga" => 31,
    "Ge" => 32,
    "As" => 33,
    "Se" => 34,
    "Br" => 35,
    "Kr" => 36,
    "Rb" => 37,
    "Sr" => 38,
    "Y" => 39,
    "Zr" => 40,
    "Nb" => 41,
    "Mo" => 42,
    "Tc" => 43,
    "Ru" => 44,
    "Rh" => 45,
    "Pd" => 46,
    "Ag" => 47,
    "Cd" => 48,
    "In" => 49,
    "Sn" => 50,
    "Sb" => 51,
    "Te" => 52,
    "I" => 53,
    "Xe" => 54,
    "Cs" => 55,
    "Ba" => 56,
    "La" => 57,
    "Ce" => 58,
    "Pr" => 59,
    "Nd" => 60,
    "Pm" => 61,
    "Sm" => 62,
    "Eu" => 63,
    "Gd" => 64,
    "Tb" => 65,
    "Dy" => 66,
    "Ho" => 67,
    "Er" => 68,
    "Tm" => 69,
    "Yb" => 70,
    "Lu" => 71,
    "Hf" => 72,
    "Ta" => 73,
    "W" => 74,
    "Re" => 75,
    "Os" => 76,
    "Ir" => 77,
    "Pt" => 78,
    "Au" => 79,
    "Hg" => 80,
    "Tl" => 81,
    "Pb" => 82,
    "Bi" => 83,
    "Po" => 84,
    "At" => 85,
    "Rn" => 86,
    "Fr" => 87,
    "Ra" => 88,
    "Ac" => 89,
    "Th" => 90,
    "Pa" => 91,
    "U" => 92,
    "Np" => 93,
    "Pu" => 94,
    "Am" => 95,
    "Cm" => 96,
    "Bk" => 97,
    "Cf" => 98,
    "Es" => 99,
    "Fm" => 100,
    "Md" => 101,
    "No" => 102,
    "Lr" => 103,
    "Rf" => 104,
    "Db" => 105,
    "Sg" => 106,
    "Bh" => 107,
    "Hs" => 108,
    "Mt" => 109,
    "Ds" => 110,
    "Rg" => 111,
    "Cn" => 112,
    "Nh" => 113,
    "Fl" => 114,
    "Mc" => 115,
    "Lv" => 116,
    "Ts" => 117,
    "Og" => 118
)

const num_atom = Dict(
"1"  => "H",
"2"  => "He",
"3"  => "Li",
"4"  => "Be",
"5"  => "B",
"6"  => "C",
"7"  => "N",
"8"  => "O",
"9"  => "F",
"10" => "Ne",
"11" => "Na",
"12" => "Mg",
"13" => "Al",
"14" => "Si",
"15" => "P",
"16" => "S",
"17" => "Cl",
"18" => "Ar",
"19" => "K",
"20" => "Ca",
"21" => "Sc",
"22" => "Ti",
"23" => "V",
"24" => "Cr",
"25" => "Mn",
"26" => "Fe",
"27" => "Co",
"28" => "Ni",
"29" => "Cu",
"30" => "Zn",
"31" => "Ga",
"32" => "Ge",
"33" => "As",
"34" => "Se",
"35" => "Br",
"36" => "Kr",
"37" => "Rb",
"38" => "Sr",
"39" => "Y",
"40" => "Zr",
"41" => "Nb",
"42" => "Mo",
"43" => "Tc",
"44" => "Ru",
"45" => "Rh",
"46" => "Pd",
"47" => "Ag",
"48" => "Cd",
"49" => "In",
"50" => "Sn",
"51" => "Sb",
"52" => "Te",
"53" => "I",
"54" => "Xe",
"55" => "Cs",
"56" => "Ba",
"57" => "La",
"58" => "Ce",
"59" => "Pr",
"60" => "Nd",
"61" => "Pm",
"62" => "Sm",
"63" => "Eu",
"64" => "Gd",
"65" => "Tb",
"66" => "Dy",
"67" => "Ho",
"68" => "Er",
"69" => "Tm",
"70" => "Yb",
"71" => "Lu",
"72" => "Hf",
"73" => "Ta",
"74" => "W",
"75" => "Re",
"76" => "Os",
"77" => "Ir",
"78" => "Pt",
"79" => "Au",
"80" => "Hg",
"81" => "Tl",
"82" => "Pb",
"83" => "Bi",
"84" => "Po",
"85" => "At",
"86" => "Rn",
"87" => "Fr",
"88" => "Ra",
"89" => "Ac",
"90" => "Th",
"91" => "Pa",
"92" => "U",
"93" => "Np",
"94" => "Pu",
"95" => "Am",
"96" => "Cm",
"97" => "Bk",
"98" => "Cf",
"99" => "Es",
"100" => "Fm",
"101" => "Md",
"102" => "No",
"103" => "Lr",
"104" => "Rf",
"105" => "Db",
"106" => "Sg",
"107" => "Bh",
"108" => "Hs",
"109" => "Mt",
"110" => "Ds",
"111" => "Rg",
"112" => "Cn",
"113" => "Nh",
"114" => "Fl",
"115" => "Mc",
"116" => "Lv",
"117" => "Ts",
"118" => "Og"
)


const core_elec_perrow = Dict(
    1:2    => 0,
    3:10   => 2,
    11:18  => 10,
    19:36  => 18,
    37:54  => 36,
    55:86  => 54,
    87:118 => 86
)
end #module
