# Source paper

The PDF is **not redistributed in this repository**. Download it from arXiv:

> **https://arxiv.org/abs/2608.06242**

## Citation

> T. Haug, K. Bharti and L. Aolita,
> *Exponential logical-error reduction in quantum memories via optimal
> syndrome-measurement timing*,
> arXiv:2608.06242 [quant-ph] (2026).
> DOI: [10.48550/arXiv.2608.06242](https://doi.org/10.48550/arXiv.2608.06242)

Submitted 6 August 2026. No journal reference at the time of writing.

```bibtex
@misc{haug2026timing,
  title        = {Exponential logical-error reduction in quantum memories via
                  optimal syndrome-measurement timing},
  author       = {Haug, Tobias and Bharti, Kishor and Aolita, Leandro},
  year         = {2026},
  eprint       = {2608.06242},
  archivePrefix= {arXiv},
  primaryClass = {quant-ph},
  doi          = {10.48550/arXiv.2608.06242}
}
```

## Why the PDF is not here

The arXiv abstract page lists the licence as **"arXiv.org perpetual
non-exclusive license 1.0"**. That licence grants arXiv the right to distribute
the work; it is **not** a Creative Commons licence and does not grant third
parties permission to redistribute the PDF. Copyright remains with the authors.

So `docs/paper/` is git-ignored apart from this file. Download the PDF from the
link above and drop it in this directory if you want it locally — it will not be
committed.

## What you actually need

For day-to-day work you should not need the PDF. The equations used by this
repository are transcribed in **[`../paper_equations.md`](../paper_equations.md)**,
which is the canonical reference here (CLAUDE.md: prefer that file over the PDF
text for any formula). Equation numbers throughout the code and docs refer to it.
