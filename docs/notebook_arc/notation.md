# Notation

One symbol, one meaning, across the paper, the docs and notebooks 04–07. Where a notebook used
something else before, the old form is listed so older text can still be read.

Source column: **paper** = Haug, Bharti & Aolita (2026), equations in `paper_equations.md`;
**repo** = a quantity this repo defines.

## Time

| symbol | meaning | units | source | old forms |
| --- | --- | --- | --- | --- |
| $\Delta t$ | syndrome interval: time between the starts of consecutive rounds | see "time unit" below | paper | `τ`, `tau`, `dt` |
| $\Delta t^*$ | the optimal interval (minimises $R$) | same | paper (Eq. 6) | `τ*`, `tau*`, `dt*` |
| $T$ | total storage time, fixed in every sweep | same | paper | "total idle" (06) |
| $N = T/\Delta t$ | number of syndrome rounds | — | repo | `R`, `rounds` |
| $t_R$ | duration of one extraction round | — | repo (07) | |

**Time unit, by notebook.** The paper's time is dimensionless.
- 05: the paper's dimensionless time ($T = 96$).
- 06a/b/c: one idle operation per data qubit. $\Delta t$ = number of idle ops before each
  round, $T$ = 96 idle ops in total; the extraction round itself adds no time.
- 07: $t_R$. $\Delta t$ includes the round, so $\Delta t \ge 1$.

## Error rates

| symbol | meaning | source | old forms |
| --- | --- | --- | --- |
| $p$ | physical error scale; in 04–06c, the extraction cost of one round | paper | |
| $\lambda$ | idling rate | paper (Eq. 1) | |
| $p_\mathrm{idle} = 1 - e^{-p\lambda\Delta t}$ | idle bit-flip per interval | paper (Eq. 1) | |
| $p_\mathrm{stab} = p$ | fault on each data qubit per round, caused by measuring | paper (Eq. 2) | |
| $p_\mathrm{data} = 1-(1-p_\mathrm{stab})(1-p_\mathrm{idle})$ | total data bit-flip per round | paper (Eq. 2) | |
| $p_\mathrm{read} = b_\mathrm{read}\,p$ | flip of a recorded syndrome bit | paper (Eq. 3) | |
| $b_\mathrm{read}$ | readout coefficient; a free knob (paper uses 0.5, 1, 2) | paper | |
| $p_L$ | logical error probability over the whole run (time $T$) | paper | |
| $R$ | **logical error rate per unit time**, $R = -\ln(1-2p_L)/(2T) \approx p_L/T$ | paper (Eqs. 4, A3) | `r`, "rate" |

## Circuit-level noise (Selene, Stim)

| symbol | meaning | source | old forms |
| --- | --- | --- | --- |
| $p_{1q}$ | depolarizing error per single-qubit gate | repo | `p1` (07) |
| $p_{2q}$ | depolarizing error per two-qubit gate | repo | `p2` (07) |
| $p_\mathrm{meas}$ | flip of a recorded measurement outcome (every measurement, final data readout included) | repo | |
| $p_\mathrm{init}$ | preparation error | repo | |
| $f_1 = 2p_{1q}/3$ | bit-flip probability per charged single-qubit op (X and Y flip $\lvert 0\rangle$) | repo (03, 04) | `f1`, `F1` |
| $n_\mathrm{data}$, $n_\mathrm{anc}$ | charged ops per round on a data qubit / a Z-ancilla | repo (04) | `ops_mean`, `ops_anc` |
| $p_\mathrm{anc}$ | Z-ancilla bit-flip per round from its own gates | repo (04, 05, 06c) | |
| $a \oplus b = a + b - 2ab$ | two independent flips (they cancel when both occur) | — | |

Effective parameters (04 → 05, 06c): $p = $ accumulate$(n_\mathrm{data}, f_1)$,
$\lambda = 1/n_\mathrm{data}$, $p_\mathrm{read} = p_\mathrm{anc} \oplus p_\mathrm{meas}$ (04d §3;
05 and 06c used $1-(1-p_\mathrm{anc})(1-p_\mathrm{meas})$, which differs by $p_\mathrm{anc}p_\mathrm{meas}$).

## Dephasing model (07)

| symbol | meaning | source |
| --- | --- | --- |
| $\gamma$ | linear dephasing rate: Z error $\gamma\,\Delta t$ per interval | repo (07) |
| $c$ | quadratic dephasing: Z error $\sin^2(\sqrt{c}\,\Delta t) \approx c\,\Delta t^2$ | repo (07) |
| $q(\Delta t) = \gamma\Delta t + \sin^2(\sqrt{c}\,\Delta t)$ | idle Z error per interval | repo (07) |
| $p_e$ | effective extraction error in the 07 ansatz (the role $p_\mathrm{stab}$ plays in the paper) | repo (07) |

$\gamma$ and $c$ are quoted relative to $p_{2q}$, e.g. $c/p_{2q} = 0.01$.

## Ansatz

| symbol | meaning | source |
| --- | --- | --- |
| $d$ | code distance | paper |
| $\alpha = g(d+1)/2$ | exponent of the idle factor | paper (Eq. 5) |
| $g$, $A$, $\beta$, $p_\mathrm{th}$ | fitted constants of the rate ansatz | paper (Eq. 5) |

## Code identifiers

Code keeps ASCII names; printed tables show them. `tau`, `dt` = $\Delta t$;
`R`, `rounds` = $N$ (rounds); `r` (07) and the fitted rate in 05/06c = $R$;
`p1`, `p2` (07) = $p_{1q}$, $p_{2q}$; `F1`, `f1` = $f_1$.
