"""Reference data for homonuclear diatomics (H2 ... F2).

Experimental total energy of a molecule is built from
    E_exp = 2 * E_atom(exact, non-relativistic) - D_e
where D_e is the experimental well depth (spectroscopic dissociation energy).
All energies in Hartree, distances in Angstrom.
"""

EV_TO_HARTREE = 1.0 / 27.211386

# Z = atomic number of ONE atom
MOLECULES = {
    #        Z   spin(2S)  Re_exp(A)  De_exp(eV)  E_atom_exact(Eh)  term symbol
    "H2": dict(Z=1, spin=0, re=0.7414, de=4.747, e_atom=-0.500000, term="1Σg+"),
    "Li2": dict(Z=3, spin=0, re=2.6729, de=1.056, e_atom=-7.478060, term="1Σg+"),
    "B2": dict(Z=5, spin=2, re=1.5900, de=3.090, e_atom=-24.653910, term="3Σg-"),
    "C2": dict(Z=6, spin=0, re=1.2425, de=6.320, e_atom=-37.845000, term="1Σg+"),
    "N2": dict(Z=7, spin=0, re=1.0977, de=9.905, e_atom=-54.589300, term="1Σg+"),
    "O2": dict(Z=8, spin=2, re=1.2075, de=5.213, e_atom=-75.067300, term="3Σg-"),
    "F2": dict(Z=9, spin=0, re=1.4119, de=1.660, e_atom=-99.733900, term="1Σg+"),
}

def exp_energy(name):
    m = MOLECULES[name]
    return 2 * m["e_atom"] - m["de"] * EV_TO_HARTREE