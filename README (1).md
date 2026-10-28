# Hartree–Fock (HF) Solver for Homonuclear Diatomic Molecules

An interactive Python app that uses **PySCF** to estimate the ground-state energies and
equilibrium bond lengths of **H₂, Li₂, B₂, C₂, N₂, O₂ and F₂** at the Hartree–Fock level,
and compares them with experiment. The interface is built with **Streamlit**.

## Features

**Tab 1 – E vs r**
- Select one or more molecules (default: H₂) and plot their HF potential-energy curves E(r) on the same axes.
- Marks the HF equilibrium bond length (★ and dashed line) and the experimental bond length (dotted line); both values appear in the legend.
- Optional experimental ground-state energy line (dash-dot).
- Choose the energy axis: absolute (Eh) or relative to each HF minimum, and the bond-length axis: r (Å) or r / r_e(exp).
- Summary table: HF and experimental bond lengths, HF and experimental energies, ΔE and percentage absolute error.

**Tab 2 – Ground-state energy vs Z**
- HF and experimental ground-state energies plotted against the atomic number Z of the atom (H = 1 … F = 9).
- Bar chart of the percentage absolute error, plus a summary table.

## Method

- Basis set: **STO-3G** (fixed).
- Closed-shell ground states (H₂, Li₂, C₂, N₂, F₂) use **RHF**; the triplet ground states (B₂, O₂) use **UHF** with 2S = 2.
- E(r) is scanned from 0.7 to 2.6 × the experimental r_e; the equilibrium bond length is found by bounded minimisation (SciPy) of the HF energy.
- Experimental total energy of a molecule:

      E_exp = 2 · E_atom(exact, non-relativistic) − D_e

  using exact atomic energies and experimental dissociation energies D_e (stored in `data.py`).
- Percentage absolute error = |E_exp − E_HF| / |E_exp| × 100.

## Files

| File | Purpose |
|------|---------|
| `app.py` (or `App.py`) | Streamlit interface and plots |
| `hf_solver.py` | PySCF calculations: energy at a given r, E(r) scan, equilibrium search |
| `data.py` | Experimental reference data (r_e, D_e, atomic energies, spin states) |
| `requirements.txt` | Python dependencies |

## Installation and usage

    pip install -r requirements.txt
    streamlit run app.py

Use the exact file name on your machine (e.g. `streamlit run App.py`). If `streamlit` is not
found, use `python3 -m streamlit run App.py`. The app opens in your browser at
`http://localhost:8501`. Running the file with `python App.py` will not start the interface.

PySCF works best on macOS and Linux; on Windows use WSL or a conda environment.

## Notes and limitations

- STO-3G is a minimal basis, so HF energies lie well above experiment; the gap includes both basis-set error and the missing electron-correlation energy.
- Restricted HF dissociates to the wrong limit for closed-shell molecules, so the energy rises at large r instead of levelling off.
- The experimental D_e and atomic-energy values were entered by hand; check them against a source such as NIST before quoting them (B₂ and F₂ are the least certain).
- F₂ is predicted to be barely bound by HF, and B₂/C₂ can be difficult to converge (a second-order solver is used as a fallback).
