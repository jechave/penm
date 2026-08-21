# docs

Notes and theory for `penm`. Migrated from `penmscan` on 2026-08-20
and reorganised the same day; nothing was rewritten in the move.

Dates below are the **true** last-modified dates, recovered from penmscan's git
history — filesystem dates all read 2026-08-20 because that is when the files
were copied.

## Where to start

- New to the model: `theory/lfenm.md`, then `theory/sc_lfenm.md`.
- The full development: `theory/dresden2017/` — read the **annotated** PDF.
- Open issues in the mutation models: `theory/NOTES.Rmd`.

---

## theory/ — the penm model

The elastic network, its perturbation by mutation, wild-type vs mutant.

| file | date | what it is |
|---|---|---|
| `lfenm.md` | 2020-05-14 | LFENM potential. States the key result: the mutant's surface is shifted, not rotated or deformed — same modes, same eigenvalues. |
| `sc_lfenm.md` | 2022-06-30 | SC-LFENM: why $K_{wt} \neq K(r_{mut})$ is inconsistent, and the fix. Ends with the open `#todo` that the model is **not reversible**. |
| `wt_mut_energies.md` | 2020-05-20 | Wild-type and mutant energies, thermodynamic functions. Gives **two ways to compute the mutant's minimum** — they should agree, which makes it a checkable invariant. |
| `NOTES.Rmd` | 2020-05-14 | Titled *"Notes"*; its section is *"Many ways of updating ENM..."*. Compares LFENM and SC-LFENM and raises five open issues; see below. |
| `sc_lfenm.pdf` | 2022-06-30 | rendered from `sc_lfenm.md` |
| `echave-fernandez-2010-lfenm-appendix.pdf` | 2020-05-25 | *(was `references/2010_proteins_appendix.pdf`)* **The published LFENM derivation** (Echave & Fernández 2010, appendix). The authoritative source for §4: expansion of $d_{ij}$ in displacements, $K_{ij} = -G_{ij}$ with $G_{ij}=k_{ij}\hat{d}_{ij}\hat{d}_{ij}^T$, the mutation as $\bar d_{ij} \to \bar d_{ij} + \delta\bar d_{ij}$, and $\delta\bar r = K^{-1}f$ with $K^{-1}$ **explicitly defined as the pseudo-inverse over the $3N-6$ non-zero modes**. |

### The five open issues in `NOTES.Rmd`

The file is short. It sets LFENM against SC-LFENM — the first changes no $K$, the
second rebuilds it — and raises five issues. Some concern SC-LFENM specifically,
others the mutation models generally:

1. **Which $K$ for $\delta r$** — $K_{wt}$ or $K_{mut}$? Raised nowhere else.
2. **Reversibility** — restated independently in `sc_lfenm.md` two years later.
3. **Path-independence** — does the order of two successive mutations matter?
   Raised *only* here.
4. **Connectivity changes** — topology change may cause large energy jumps.
5. **$V_0$ changes** — whether unrelaxable mutational stress enters the Hamiltonian.

It also distinguishes two models both called SC-LFENM: update $K$ and use it for
subsequent deformations (full, slow), or update it only to compute entropies and
dynamics (fast).

### theory/dresden2017/

The 2017 sabbatical research plan, *Biophysical origin of long-range functional
constraints in enzyme evolution*. Historical — kept as written, not maintained.

| file | what it is |
|---|---|
| `dresden2017-penm-theory-annotated.tex` / `.pdf` | **read this one.** The 2017 text, unmodified, plus commentary marked **[T]** (theory) and **[C]** (correspondence with the package). Strip it all by redefining `\claude` and `\claudecode` to `{}`. |
| `dresden2017-penm-theory.tex` / `.pdf` | the source as received (2026-08-20), unmodified, and its build |
| `dresden2017-build-2017-02-15.pdf` | *(was `references/msa_model_theory.pdf`)* **the original 2017 build**, 16 pp. Same document, earlier state — see below. |
| `dresden2017-offcut.tex` | a 13-line fragment cut from the main document: opens mid-sentence, has a `\subsection{?}`, and its `\ref`s resolve only against the main file |
| `dresden2017-offcut-standalone.tex` / `.pdf` | the fragment wrapped in a preamble so it compiles; its two `\ref`s render as `??` |

**The two builds differ, and the difference is the sign question.** The 2017-02-15
build (16 pp) prints the Hessian with $(1 - l_{ij}/d_{ij})$ and nothing else. The
current source (17 pp) adds a `\comment{}` block giving a second derivation with
$(l_{ij}/d_{ij} - 1)$, prefaced *"But beware of the sign"*, plus the trace
argument. So the doubt was introduced **after** February 2017. The annotated
build resolves it: the second form is correct as written.

Aside from that, the two builds differ only in how `\neq` renders (`6=` vs `≠`),
a font artifact. Verified by text extraction.

---

## scanning/ — response and mutational scanning

**This is penmscan's theory, not penm's.** Kept here because penm was extracted
from penmscan and the notes came with it. If penmscan is ever slimmed to
`Imports: penm`, this directory should follow the scanning code.

| file | date | what it is |
|---|---|---|
| `prs_theory.md` | 2020-08-05 | Structural response theory. The $\delta e / \delta f / \delta r$ scaling laws ($df^2 \propto 1/\sigma^2$, $de^2 \propto 1$, $dr^2 \propto \sigma^2$), the sum-over-edges → sum-over-sites identity, and a note that entropy separates rigorously by mode but not by site. Several results marked `#check`. |
| `superfast_response_theory.md` | 2022-06-30 | Analytical mean response without simulating mutations. The full derivation is `superfast-response-derivation.pdf` (below). |
| `superfast-response-derivation.pdf` | 2020-05-18 | *(was `references/superfast_response_calculation_theory.pdf`)* Derivation of the response formulae, ending in the "magic formula" $R_i^l = \sum_{k\sim l}\hat{d}_{kl}^T(A_{il}-A_{ik})^T(A_{il}-A_{ik})\hat{d}_{kl}$. Carries a **Warning** the `.md` omits: $C^{1/2} = U\Lambda^{-1/2}U^T$, the inverse square root. |
| `superfast_energy_response.md` | 2020-05-25 | Energy response. Contains the open `#todo` on **negative deformation energies** under frustration. |
| `deij2_constant.md` | 2020-05-15 | Conjecture that $\langle \delta e_i^2 \rangle_j \approx$ constant. Two unchecked boxes. |
| `influence_theory.md` | 2022-06-30 | Influence profiles. Three `#todo`s, no content yet. |
| `prs_convergence.md` | 2020-05-19 | **A measurement.** Simulation converges to the analytical formulae; inflection ~30 mutations/site; matrices converge poorly (~0.5 rel. err.), profiles well (<0.1 with ~20). |
| `double_mutation_response_scanning.md` / `.tex` / `.html` | 2022-06-30 | Double-mutation scanning; constrained maximisation of quadratic forms. |

---

## archive/ — superseded, kept for the record

| file | date | why archived |
|---|---|---|
| `LOGFILE.md` / `.Rmd` / `.pdf` | 2022-08-11 / 2023-04-12 | Development log, 2020–2023. Records two unresolved warnings: *"v_min changes between update_enm = F and update_enm = T"*, and *"mut_graph assumes frustrated = T? For now use frustrated = F everywhere, until I check it thoroughly"* — no later entry records the check. |
| `FUNCTIONS.Rmd` / `.pdf` | 2022-08-10 | **Stale.** A function-signature list; still shows `seed`, renamed to `ensemble` on 2026-08-19. Read the code, not this. |
| `daily_notes.md` | 2020-08-05 | Two entries on optimising response calculation. Useful line: the bottleneck turned out to be generating the forces, not the matrix multiply. |
| `NOTES-2020.pdf` | 2020-05-14 | rendered `NOTES.Rmd`; regenerate rather than trust |

---

## Deleted in the 2026-08-20 tidy

Build artifacts, duplicates and one derived file:

- `sc_lfenm.aux`, `double_-mutation_resonse_scanning.log`, two `.DS_Store`
- `double_-mutation_resonse_scanning.html` — regenerated into `scanning/` under a
  corrected filename; body content verified byte-identical before deletion
- `LFENM___Active_Site_Selection_Model.pdf` — a 2026-08-20 build of the same
  source as `dresden2017-penm-theory.pdf`; verified identical by text extraction
  apart from `\neq` glyph rendering
- `LFENM___Active_Site_Selection_Model.zip` — its two `.tex` files are extracted
  and in place

The `references/` directory was dissolved: none of its three files was an
external paper. All three are the author's own, and each now sits beside the
note it underpins.

## Conventions

- Notes are dated in this index, not in the files. Filesystem dates are wrong.
- PDFs with a `.md`/`.Rmd`/`.tex` sibling are derived; regenerate rather than edit.
- `#todo`, `#check` and `#warning` markers in the notes are the author's own open
  questions. They are real and mostly still open.
