"""Hartree-Fock energies of diatomics with PySCF."""
import numpy as np
from pyscf import gto, scf
from scipy.optimize import minimize_scalar
from data import MOLECULES


def hf_energy(name, r, basis="sto-3g"):
    """HF total energy (Hartree) of molecule `name` at bond length r (Angstrom).
    Closed-shell -> RHF, triplets (B2, O2) -> UHF with 2S = 2."""
    spin = MOLECULES[name]["spin"]
    atom = name.rstrip("0123456789")  # 'H2' -> 'H', 'Li2' -> 'Li'
    mol = gto.M(atom=f"{atom} 0 0 0; {atom} 0 0 {r}", basis=basis,
                spin=spin, unit="Angstrom", verbose=0)
    mf = scf.RHF(mol) if spin == 0 else scf.UHF(mol)
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:                       # fall back to 2nd-order solver
        mf = mf.newton()
        mf.kernel()
    return mf.e_tot


def scan(name, basis="sto-3g", rmin=0.7, rmax=2.6, n=30):
    """Potential energy curve; r range given as multiples of experimental Re."""
    re = MOLECULES[name]["re"]
    rs = np.linspace(rmin * re, rmax * re, n)
    es = np.array([hf_energy(name, r, basis) for r in rs])
    return rs, es


def equilibrium(name, basis="sto-3g"):
    """Optimise the bond length -> (r_eq, E_min)."""
    re = MOLECULES[name]["re"]
    res = minimize_scalar(lambda r: hf_energy(name, r, basis),
                          bounds=(0.8 * re, 1.3 * re), method="bounded",
                          options={"xatol": 1e-4})
    return res.x, res.fun