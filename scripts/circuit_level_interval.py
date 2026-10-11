#!/usr/bin/env python
"""Syndrome-interval analysis, end to end, for the distances you ask for.

Two backends (see src/qec_timing/circuit_level/ for what each noise model can and cannot express):

  --backend selene   Selene's native depolarizing model on the compiled circuit. For each d:
    1. measure the effective noise model of one extraction round -- charged 1q ops per data
       qubit and per Z-ancilla, and the Rzz counts when p_2q > 0 (notebooks 04s and 04c);
    2. run the memory at fixed R, sweeping the idle interval tau, decoded with PyMatching
       weighted by the step-1 counts (notebook 07's design).

  --backend stim     circuit-level Stim with emulator-style channels. For each d:
    1. build the memory circuit with a time-based idle channel (depolarizing + linear and
       quadratic dephasing), in the Z or X basis;
    2. sweep tau at fixed R, decoded with PyMatching on the exact detector error model.

Both then score each tau by the paper's figure of merit pL / T, with T = R * tau, alongside
the decay rate Gamma = -ln(1 - 2 pL) / 2T -- equal to pL / T while pL << 1, but not rewarding
a failed code as pL -> 0.5 (derivation in notebooks/07, "The objective").

Outputs (<tag> starts with the backend name)
  results/circuit_level/cache/effective_noise_model/<code>_d<d>.json   selene step 1, one file per distance
  results/circuit_level/cache/stim_circuits/<tag>_d<d>_tau<t>.stim     stim step 1, the circuit at the first tau
  results/circuit_level/cache/interval_analysis/<tag>.csv              every (d, tau) point
  results/circuit_level/plots/pL_per_T_<tag>.png                      the plot
Raw runs are cached in results/circuit_level/cache/<experiment>/ via runcache; Selene runs use the same
names the notebooks use, so anything a notebook already ran is reused.

Examples
  python scripts/circuit_level_interval.py --distances 3 5
  python scripts/circuit_level_interval.py --distances 5 --p2q 0.005 --k 3
  python scripts/circuit_level_interval.py --backend stim --basis X --distances 3 5 7 \\
      --p2q 0.005 --dephase-linear 0.001
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time

HERE = pathlib.Path(__file__).resolve()
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))

import matplotlib                                          # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                            # noqa: E402
import numpy as np                                         # noqa: E402
import pandas as pd                                        # noqa: E402

from qec_timing.circuit_level import runcache             # noqa: E402
from qec_timing.circuit_level import surface_memory as sm  # noqa: E402

PLOTS = REPO / "results" / "circuit_level" / "plots"
TAUS = [1, 2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 32, 48, 64, 96]   # 07's grid
PALETTE = ['#2a78d6', '#eb6834', '#1baf7a', '#8a5cd6', '#c9a227', '#d64f8a']
INK = '#0b0b0b'


def parse():
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--backend", choices=("selene", "stim"), default="selene")
    ap.add_argument("--distances", type=int, nargs="+", required=True,
                    help="odd code distances, e.g. 3 5 7")
    ap.add_argument("--code", choices=sorted(sm.SPEC), default="rotated")
    ap.add_argument("--rounds", type=int, default=8,
                    help="syndrome rounds R, held fixed while tau is swept (must be >= d)")
    ap.add_argument("--taus", type=int, nargs="+", default=TAUS,
                    help="idle intervals to sweep, in tau units")
    ap.add_argument("--p1q", type=float, default=0.001, help="1q depolarizing rate")
    ap.add_argument("--p2q", type=float, default=None,
                    help="2q depolarizing rate (default: 0 for selene, 0.005 for stim, where "
                         "the CX is the only noisy part of extraction)")
    ap.add_argument("--pmeas", type=float, default=0.002, help="readout flip rate")
    ap.add_argument("--shots", type=int, default=None,
                    help="shots per replica (default: 6000 selene, 100000 stim)")
    ap.add_argument("--replicas", type=int, default=2,
                    help="independent seeds per point, pooled")
    sel = ap.add_argument_group("selene only")
    sel.add_argument("--k", type=int, default=1,
                     help="idle ops per tau unit; lambda scales with k")
    st = ap.add_argument_group("stim only")
    st.add_argument("--basis", choices=("Z", "X"), default="Z",
                    help="memory basis; Z-type dephasing is only visible in X")
    st.add_argument("--idle-depol", type=float, default=None,
                    help="depolarizing per tau unit of idling (default: p1q, i.e. one "
                         "idle op per tau unit, as in the selene arm)")
    st.add_argument("--dephase-linear", type=float, default=0.0,
                    help="linear dephasing: Z with probability rate * tau per idle period")
    st.add_argument("--dephase-quad", type=float, default=0.0,
                    help="quadratic dephasing frequency f: Z with probability sin^2(f tau/2)")
    st.add_argument("--pinit", type=float, default=0.0, help="reset flip rate")
    a = ap.parse_args()

    bad = [d for d in a.distances if d < 3 or d % 2 == 0]
    if bad:
        ap.error(f"distances must be odd and >= 3: {bad}")
    if a.rounds < max(a.distances):
        ap.error(f"--rounds {a.rounds} is below d={max(a.distances)}; the paper's "
                 f"figure of merit needs R >= d")
    if a.backend == "selene":
        stim_only = [f for f, v, dflt in (("--basis", a.basis, "Z"),
                                          ("--idle-depol", a.idle_depol, None),
                                          ("--dephase-linear", a.dephase_linear, 0.0),
                                          ("--dephase-quad", a.dephase_quad, 0.0),
                                          ("--pinit", a.pinit, 0.0)) if v != dflt]
        if stim_only:
            ap.error(f"{', '.join(stim_only)}: requires --backend stim; Selene's error model has "
                     f"no duration-based or biased channel and this memory is Z-basis only")
        a.p2q = 0.0 if a.p2q is None else a.p2q
        a.shots = a.shots or 6000
    else:
        if a.k != 1:
            ap.error("--k is selene only; in stim the idle channel is set by "
                     "--idle-depol / --dephase-linear / --dephase-quad")
        a.p2q = 0.005 if a.p2q is None else a.p2q
        a.idle_depol = a.p1q if a.idle_depol is None else a.idle_depol
        a.shots = a.shots or 100_000
        if a.basis == "Z" and (a.dephase_linear or a.dephase_quad):
            print("WARNING: dephasing is Z-type and invisible to a Z-basis memory; "
                  "use --basis X to see it")
    return a


def score(d, tau, rounds, pL, s):
    T = rounds * tau
    sat = pL > 0.5 - 3 * s          # a coin flip: the code has failed
    return {'d': d, 'tau': tau, 'T': T, 'pL': pL, 'sigma': s,
            'pL_per_T': pL / T, 'pL_per_T_sigma': s / T,
            'Gamma': np.nan if sat else -np.log(1 - 2 * pL) / (2 * T),
            'Gamma_sigma': np.nan if sat else s / ((1 - 2 * pL) * T),
            'saturated': bool(sat)}


# ---- selene ------------------------------------------------------------------------------
def selene_logical_error(G, a, tau):
    """Pooled (pL, sigma) over replicas, each replica cached separately."""
    F = N = 0
    for r in range(a.replicas):
        params = dict(code=a.code, distance=G.d, tau=tau, rounds=a.rounds,
                      shots=a.shots, shot_offset=r * a.shots, p_1q=a.p1q, p_2q=a.p2q,
                      p_meas=a.pmeas, decoder='pymatching', priors='measured-per-qubit')
        # Keys match 07_memory_error_rate at its defaults, so its runs are reused; the
        # extras appear only when they change what is simulated or how it is decoded.
        if a.k != 1:
            params['idle_ops_per_tau'] = a.k
        if a.p2q:
            params['priors'] = 'measured-per-qubit+p2q'
        hit = runcache.load('07_memory_error_rate', params)
        if hit is None:
            fails = G.failures(tau * a.k, a.rounds, a.shots, offset=r * a.shots,
                               p_1q=a.p1q, p_2q=a.p2q, p_meas=a.pmeas)
            hit = {'failures': int(fails), 'shots': int(a.shots), 'pL': fails / a.shots}
            runcache.save('07_memory_error_rate', params, hit)
        F += hit['failures']
        N += hit['shots']
    e = F / N
    return e, float(np.sqrt(max(e * (1 - e), 1e-12) / N))


def run_selene(a, d):
    from qec_timing.circuit_level import effective_noise
    eff = effective_noise.measure(a.code, d, with_p2q=bool(a.p2q))
    path = effective_noise.save(eff, runcache.DATA / "effective_noise_model")
    n1 = float(np.mean(eff['ops_per_qubit']))
    lam = a.k / (n1 + (0.8 * a.p2q / a.p1q * float(np.mean(eff['n2_per_qubit']))
                       if a.p2q else 0.0))
    print(f"  n_1 = {n1:.3f} ops/data qubit/round, ancilla {eff['ops_anc']:.3f}"
          + (f", n_2 = {np.mean(eff['n2_per_qubit']):.3f}" if a.p2q else "")
          + f"  ->  lambda = {lam:.4f}   [{path.name}]")
    G = sm.Geometry(a.code, d, constants=eff)
    rows = []
    for tau in a.taus:
        t0 = time.time()
        pL, s = selene_logical_error(G, a, tau)
        rows.append(score(d, tau, a.rounds, pL, s))
        report(rows[-1], t0)
    return pd.DataFrame(rows).assign(lam=lam)


# ---- stim --------------------------------------------------------------------------------
def run_stim(a, d):
    from qec_timing.circuit_level.stim_memory import StimMemory, idle_channel
    M = StimMemory(a.code, d, basis=a.basis, p1=a.p1q, p2=a.p2q, p_meas=a.pmeas,
                   p_init=a.pinit, idle_depol=a.idle_depol,
                   dephase_linear=a.dephase_linear, dephase_quad=a.dephase_quad)
    circ_dir = runcache.DATA / "stim_circuits"
    circ_dir.mkdir(parents=True, exist_ok=True)
    cpath = circ_dir / f"{tag(a)}_d{d}_tau{a.taus[0]}.stim"
    M.circuit(a.taus[0], a.rounds).to_file(cpath)
    print(f"  circuit-level Stim, basis {a.basis}   [{cpath.name}]")
    rows = []
    for tau in a.taus:
        t0 = time.time()
        F = N = 0
        for r in range(a.replicas):
            params = dict(code=a.code, distance=d, basis=a.basis, tau=tau, rounds=a.rounds,
                          shots=a.shots, seed=r + 1, p1=a.p1q, p2=a.p2q, p_meas=a.pmeas,
                          p_init=a.pinit, idle_depol=a.idle_depol,
                          dephase_linear=a.dephase_linear, dephase_quad=a.dephase_quad,
                          decoder='pymatching-dem')
            hit = runcache.load('stim_memory', params)
            if hit is None:
                _, _, fails = M.logical_error(tau, a.rounds, a.shots, seed=r + 1)
                hit = {'failures': fails, 'shots': a.shots, 'pL': fails / a.shots}
                runcache.save('stim_memory', params, hit)
            F += hit['failures']
            N += hit['shots']
        pL = F / N
        s = float(np.sqrt(max(pL * (1 - pL), 1e-12) / N))
        px, py, pz = idle_channel(tau, a.idle_depol, a.dephase_linear, a.dephase_quad)
        rows.append({**score(d, tau, a.rounds, pL, s),
                     'idle_px': px, 'idle_py': py, 'idle_pz': pz})
        report(rows[-1], t0)
    return pd.DataFrame(rows)


# ---- shared ------------------------------------------------------------------------------
def report(row, t0):
    print(f"    tau={row['tau']:3d}  T={row['T']:4d}  pL={row['pL']:.4f} +/- {row['sigma']:.4f}"
          f"  pL/T={row['pL_per_T']:.3e}{'  SATURATED' if row['saturated'] else ''}"
          f"  [{time.time()-t0:.0f}s]")


def tag(a):
    t = f"{a.backend}_{a.code}_R{a.rounds}_p1q{a.p1q:g}_pmeas{a.pmeas:g}"
    if a.p2q:
        t += f"_p2q{a.p2q:g}"
    if a.backend == "selene":
        return t + (f"_k{a.k}" if a.k != 1 else "")
    t += f"_{a.basis}_idep{a.idle_depol:g}"
    for name, v in (("lin", a.dephase_linear), ("quad", a.dephase_quad), ("pinit", a.pinit)):
        if v:
            t += f"_{name}{v:g}"
    return t


def describe(a):
    s = (f"{a.code} surface code, {a.backend}, R={a.rounds} fixed, p_1q={a.p1q:g}, "
         f"p_2q={a.p2q:g}, p_meas={a.pmeas:g}")
    if a.backend == "selene":
        return s + (f", k={a.k}" if a.k != 1 else "")
    s += f", basis {a.basis}\nidle per tau unit: depol {a.idle_depol:g}"
    if a.dephase_linear:
        s += f", linear dephasing {a.dephase_linear:g}"
    if a.dephase_quad:
        s += f", quadratic dephasing f={a.dephase_quad:g}"
    return s + (f", p_init={a.pinit:g}" if a.pinit else "")


def plot(df, a, path):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), constrained_layout=True)
    for (d, g), col in zip(df.groupby('d'), PALETTE):
        ok, bad = g[~g['saturated']], g[g['saturated']]
        for ax, y, ye in ((axes[0], 'pL_per_T', 'pL_per_T_sigma'),
                          (axes[1], 'Gamma', 'Gamma_sigma')):
            ax.errorbar(ok['tau'], ok[y], yerr=ok[ye], marker='o', ms=5, lw=1.8,
                        capsize=2, color=col, label=f'd={d}')
            if ax is axes[0] and len(bad):
                ax.plot(bad['tau'], bad[y], 'o', mfc='none', mec=col, ms=6)
            if len(ok):
                k = ok[y].idxmin()
                ax.scatter([ok.loc[k, 'tau']], [ok.loc[k, y]], s=170, facecolor='none',
                           edgecolor=col, linewidth=2, zorder=4)
    axes[0].set_ylabel('pL / T   (T = R * tau)')
    axes[0].set_title('the paper\'s figure of merit, pL / T\n'
                      'ringed: minimum; hollow: saturated (pL within 3 sigma of 0.5)')
    axes[1].set_ylabel('Gamma = -ln(1 - 2 pL) / 2T')
    axes[1].set_title('decay rate: equals pL/T while pL << 0.5,\n'
                      'and does not reward a failed code (see 07)')
    for ax in axes:
        ax.set_xscale('log', base=2)
        ax.set_yscale('log')
        ax.set_xlabel('syndrome interval tau')
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(plt.NullFormatter())
        ax.grid(alpha=.5, which='both')
        ax.set_axisbelow(True)
        ax.legend(frameon=False)
    fig.suptitle(f"{describe(a)}, {a.shots * a.replicas} shots/point",
                 fontsize=11, color=INK)
    PLOTS.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def main():
    a = parse()
    print(f"cache: {runcache.DATA}\n{describe(a)}")
    frames = []
    for d in a.distances:
        print(f"\n== {a.code} d={d} ==")
        frames.append(run_selene(a, d) if a.backend == "selene" else run_stim(a, d))
    df = pd.concat(frames, ignore_index=True)

    out = runcache.DATA / "interval_analysis"
    out.mkdir(parents=True, exist_ok=True)
    csv = out / f"{tag(a)}.csv"
    df.to_csv(csv, index=False)
    png = PLOTS / f"pL_per_T_{tag(a)}.png"
    plot(df, a, png)

    print("\nminima (unsaturated points only; an edge minimum is not an interior optimum)")
    for d, g in df.groupby('d'):
        ok = g[~g['saturated']]
        for col in ('pL_per_T', 'Gamma'):
            if len(ok) < 3:
                print(f"  d={d}  {col}: too few unsaturated points")
                continue
            i = ok[col].idxmin()
            edge = ok.loc[i, 'tau'] in (ok['tau'].min(), ok['tau'].max())
            print(f"  d={d}  {col:9s} min at tau={ok.loc[i, 'tau']:3d}"
                  f"  (pL={ok.loc[i, 'pL']:.4f}){'  EDGE' if edge else ''}")
    print(f"\nwrote {csv}\nwrote {png}")


if __name__ == "__main__":
    main()
