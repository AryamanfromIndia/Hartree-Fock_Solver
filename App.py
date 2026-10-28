import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

from data import MOLECULES, exp_energy
from hf_solver import scan, equilibrium

st.set_page_config(page_title="HF diatomics", layout="wide")
st.title("Hartree–Fock (HF) energies of homonuclear diatomics (PySCF)")
st.caption("**HF = Hartree–Fock**, **Exp. = experimental**, **Eh = Hartree** (energy unit), **r_e = equilibrium bond length**")

basis = "sto-3g"
st.sidebar.markdown("**Basis set:** STO-3G")
npts = st.sidebar.slider("Points in E(r) scan", 15, 50, 30)
st.sidebar.caption("Closed-shell molecules use RHF; B2 and O2 (triplets) use UHF.")


@st.cache_data(show_spinner=False)
def cached_scan(name, basis, n):
    return scan(name, basis, n=n)

@st.cache_data(show_spinner=False)
def cached_eq(name, basis):
    return equilibrium(name, basis)


tab1, tab2 = st.tabs(["E vs r (one or more molecules)", "Ground-state energy vs Z"])

# ---------------------------------------------------------------- tab 1
with tab1:
    names = st.multiselect("Molecules (select one or more)", list(MOLECULES), default=["H2"])

    c1, c2, c3 = st.columns(3)
    mode = c1.radio("Energy axis",
                    ["Absolute energy (Eh)", "Relative to each HF minimum (Eh)"],
                    help="Use 'relative' when comparing molecules - their absolute energies differ hugely.")
    xmode = c2.radio("Bond-length axis", ["r (Å)", "r / r_e(exp)"])
    show_exp = c3.checkbox("Show experimental energy / bond length", value=True)
    relative = mode.startswith("Relative")
    ymax = None
    if relative:
        ymax = st.slider("Max energy above minimum shown (Eh)", 0.1, 3.0, 1.0, 0.1)

    if not names:
        st.info("Select at least one molecule.")
    else:
        fig, ax = plt.subplots(figsize=(9, 5.5))
        rows = []
        for i, name in enumerate(names):
            m = MOLECULES[name]
            col = f"C{i % 10}"
            with st.spinner(f"Running Hartree–Fock for {name}..."):
                rs, es = cached_scan(name, basis, npts)
                r_eq, e_eq = cached_eq(name, basis)
            e_exp = exp_energy(name)

            scale = 1.0 if xmode == "r (Å)" else 1.0 / m["re"]
            shift = e_eq if relative else 0.0

            ax.plot(rs * scale, es - shift, "o-", ms=3, color=col, label=f"{name} Hartree–Fock (HF)")
            ax.axvline(r_eq * scale, color=col, ls="--", lw=1,
                       label=f"{name} HF $r_e$ = {r_eq:.3f} Å")            # HF bond length
            ax.plot(r_eq * scale, e_eq - shift, "*", color=col, ms=14, mec="k")
            if show_exp:
                ax.axvline(m["re"] * scale, color=col, ls=":", lw=1.5,
                           label=f"{name} Exp. $r_e$ = {m['re']:.3f} Å")     # exp. bond length
                if not relative:
                    ax.axhline(e_exp, color=col, ls="-.", lw=1)            # exp. energy
            rows.append({"Molecule": name, "State": m["term"],
                         "r_e HF (Å)": round(r_eq, 4), "r_e exp (Å)": m["re"],
                         "E HF (Eh)": round(e_eq, 5), "E exp (Eh)": round(e_exp, 5),
                         "ΔE = HF − exp (Eh)": round(e_eq - e_exp, 5),
                         "% abs. error": round(100 * abs(e_eq - e_exp) / abs(e_exp), 3)})

        ax.set_xlabel("Bond length " + xmode)
        ax.set_ylabel("E − E_min (Hartree)" if relative else "Energy (Hartree)")
        if relative:
            ax.set_ylim(-0.05 * ymax, ymax)
        ax.set_title("HF potential-energy curves  (★ HF minimum, dashed = HF $r_e$"
                     + (", dotted = exp. $r_e$" + ("" if relative else ", dash-dot = exp. E") if show_exp else "") + ")",
                     fontsize=10)
        ax.legend()
        ax.grid(alpha=.3)
        st.pyplot(fig)

        st.table(rows)
        st.caption("Restricted HF dissociates to the wrong limit (energy rises at large r for "
                   "closed-shell molecules) - a well-known static-correlation failure.")

# ---------------------------------------------------------------- tab 2
with tab2:
    if st.button("Compute all 7 molecules (may take a minute)") or st.session_state.get("done"):
        st.session_state["done"] = True
        allnames = list(MOLECULES)
        Z = [MOLECULES[n]["Z"] for n in allnames]
        prog = st.progress(0.0)
        hf, ex, re_hf = [], [], []
        for i, n in enumerate(allnames):
            r, e = cached_eq(n, basis)
            re_hf.append(r); hf.append(e); ex.append(exp_energy(n))
            prog.progress((i + 1) / len(allnames))
        hf, ex = np.array(hf), np.array(ex)

        fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8))
        a1.plot(Z, hf, "o-", label=f"Hartree–Fock (HF) / {basis.upper()}")
        a1.plot(Z, ex, "s--", label="Experiment")
        for z, n, e in zip(Z, allnames, hf):
            a1.annotate(n, (z, e), textcoords="offset points", xytext=(5, 6))
        a1.set_xlabel("Atomic number Z of the atom")
        a1.set_ylabel("Ground-state energy (Hartree)")
        a1.set_title("Total energy vs Z")
        a1.legend(); a1.grid(alpha=.3)

        pct = 100 * np.abs(ex - hf) / np.abs(ex)
        bars = a2.bar(Z, pct, width=.7, color="C2")
        for z, p in zip(Z, pct):
            a2.annotate(f"{p:.2f}%", (z, p), ha="center", xytext=(0, 3), textcoords="offset points", fontsize=8)
        a2.set_xlabel("Atomic number Z of the atom")
        a2.set_ylabel("|E_exp − E_HF| / |E_exp|  (%)")
        a2.set_title("Percentage absolute error of HF energy")
        a2.grid(alpha=.3)
        st.pyplot(fig)

        st.table({
            "Molecule": allnames,
            "Z (atom)": Z,
            "r_e HF (Å)": [round(x, 4) for x in re_hf],
            "r_e exp (Å)": [MOLECULES[n]["re"] for n in allnames],
            "E HF (Eh)": [round(x, 5) for x in hf],
            "E exp (Eh)": [round(x, 5) for x in ex],
            "|ΔE| (Eh)": [round(abs(a - b), 5) for a, b in zip(hf, ex)],
            "% abs. error": [round(float(p), 3) for p in pct],
        })