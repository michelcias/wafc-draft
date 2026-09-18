# Statistica Sinica: instructions for authors

Transcription of the journal pages "Journal Information / Aims and Scope"
and "Instructions for Authors"
(<https://www3.stat.sinica.edu.tw/statistica/>), captured by the author as
PDF on 2026-09-18 and transcribed here. The text is kept in English, as on
the source page; quotation marks mark verbatim wording. The editorial
policy carries the date **August 1, 2020**; no revision date is given for
the manuscript-preparation section.

The templates are versioned in [`../manuscript/ss-template/`](../manuscript/ss-template/).
What the pages do **not** say is listed in §7.

---

## 1. Journal and scope

- Established 1991; co-sponsored by the Institute of Statistical Science,
  Academia Sinica (ISSAS), Taiwan, and the International Chinese Statistical
  Association; fully funded by ISSAS.
- Quarterly: January, April, July, October.
- Indexed in Scopus, Science Citation Index Expanded, Current Index to
  Statistics; included in JSTOR and EBSCO.
- Print ISSN 1017-0405; online ISSN 1996-8507.
- Scope: "innovative work of high quality in all areas of statistics and
  data science, including theory, methodology, and applications".

## 2. Editorial policy

- Papers must be "of significant interest and novelty, in terms of
  applications, methodology or theory"; "minor improvements of existing
  methods or routine applications are discouraged". Review papers of
  emerging areas are encouraged.
- No prior publication and no concurrent submission elsewhere. Conference
  proceedings are an exception, provided the submission goes "beyond the
  published work in much greater depth" and highlights the differences.
- A revised version of a manuscript previously rejected by the journal is
  declined without review; authors may appeal in writing, and a manuscript
  under appeal cannot be submitted elsewhere.
- **Screening.** Co-editors and screening associate editors screen every
  paper; about 40% survive and go to an associate editor, who pre-screens
  within two weeks and sends about 60% of those to referees.
- Every claim must be "supported by numerical evidence or theoretical
  analyses"; reported numerical evidence "must be reproducible", and the
  editors may ask for the software and the datasets that produced it.
- The paper must make clear how it advances the state of knowledge and must
  acknowledge predecessors.
- Target turnaround: 2-3 months for the initial review, 2 months for a
  revision.
- **Length: the main document, including references and any appendix, must
  not exceed 40 double-spaced pages in the journal template (12pt.)**
  Supplementary material is allowed but published online only. The editorial
  office checks template compliance on receipt and returns non-compliant
  papers immediately.
- Full descriptions of methods, theoretical results and data analysis belong
  in the main paper; "technical lemmas and proofs as well as details of
  simulations and tables will mostly be placed in the online Supplement".
  Simulations may be summarized briefly in the main paper; they "are not
  always required and do not need to be comprehensive".
- "Quality of writing is an important criterion for acceptance and poor
  quality of writing may be ground for rejection at any stage."

## 3. Layout

1. Double-spaced throughout, including references; at least 1.0 inch (25mm)
   margins on all four sides; type no smaller than 12 characters per inch
   (25mm); pages numbered consecutively; 40-page cap as in §2.
2. Short **Abstract** first, then **Key words and phrases in alphabetical
   order**. Avoid formulas in the abstract.
3. A short **running title of less than 45 characters**. Institutions and
   e-mail addresses of every author go **on the last page**.
4. Numbered sections with short titles. Subsections allowed and numbered;
   **sub-subsections are not**.
5. No footnotes except in tables.
6. Citations and references: see §4.
7. Acknowledgments at the end, as brief as possible. "Essential derivations
   should appear in a separate appendix, which follows the Acknowledgments
   but precedes the References." Each appendix gets its own title.

## 4. Citations and references

In-text (`natbib`, author-year):

- Up to three authors: list all, e.g. Brown, Ivanoff and Weber (1986).
- Four or more: first author plus et al., e.g. Hainzl et al. (2013).
- Books and chapters: cite the latest edition, normally with page, section
  or chapter, e.g. Box, Jenkins and Reinsel (1994, Chap. 9).

Reference list: it must match the in-text citations. For seven or more
authors, list the first six followed by "et al." Book and chapter entries
include the book title, the editors, the first and last page of the
chapter, the publisher and the place of publication. Journal names are
given in full, the volume in bold, and the year in parentheses after the
authors. Examples as printed on the page:

> Hall, P., Kerkyacharian, G. and Picard, D. (1999). On the minimax
> optimality of block thresholded wavelet estimators. Statistica Sinica 9,
> 33-49.
>
> Hastie, T., Tibshirani, R. and Friedman, J. (2009). The Elements of
> Statistical Learning: Data Mining, Inference, and Prediction. 2nd Edition.
> Springer, New York.
>
> Beadle, G. W. (1957). The role of the nucleus in heredity. In The Chemical
> Basis of Heredity (Edited by W. D. McElroy and B. Glass), 3-22. Johns
> Hopkins Press, Baltimore.
>
> Yaida, S. (2020). Non-Gaussian processes and neural networks at finite
> widths. In Proceedings of The First Mathematical and Scientific Machine
> Learning Conference (eds. J. Lu and R. Ward) PMLR 107, 165-192.

The manuscript template offers both routes: a hand-written
`thebibliography` (`SS-template.tex`) and BibTeX with
`\bibliographystyle{chicago}` (`SS-template-bib.tex`).

## 5. Formulas, equations, tables and figures

**Formulas.**

1. Make efficient use of space; avoid unnecessary displays. Short formulas
   stay in the text but must not exceed one line in height: write
   `\sum x_i` (when the limits are obvious) or `x_1 + ... + x_n`; a binomial
   coefficient must not be left in the text, use `a!/[b!(a-b)!]` instead.
   Avoid repeating lengthy expressions by introducing notation. Number only
   the equations referred to in the text, and place the number on the right.
2. Avoid multiple under- and overbars (such as a bar over a hat over `x`).
   Align subscripts and superscripts horizontally or mark them clearly;
   avoid higher-order sub- and superscripts.
3. Mark script letters clearly; distinguish easily confused characters, such
   as `w`/`omega`, `u`/`v`, `k`/`kappa`/`K`, `o`/`O`/`0`, `l`/`1`.
4. Use `E(X)`, not `EX` and not a script `E(X)`; `Var(X)`, not `Var X`;
   `x_1, ..., x_n`, not `x_1, x_2, ..., x_n`; `a/(bc)` or `a(bc)^{-1}`, not
   a stacked fraction in the text. Avoid a centered dot for a product.
5. Except where an operator requires otherwise, use the bracket hierarchy
   `[{()}]`, including for functions.
6. Do not start a sentence with a symbol. Avoid special symbols such as the
   universal quantifier and "=" topped by "def". Follow the usual
   conventions for brackets, exp and the square root sign, as in a recent
   issue; do not use a square root sign over a lengthy expression.

**Tables and figures.**

1. Numbered separately, in order of appearance, referred to by number in the
   text, and arranged to use the journal page efficiently. Avoid landscape.
2. Each one gets a short title or description, as self-explanatory as
   possible.
3. On acceptance, original artwork must be a printable PDF or PostScript
   file; photocopies are not accepted.
4. Artwork is typically black and white. Color in print is allowed only when
   it "provides a scientific benefit", and the corresponding author must
   alert the editors on submission; then (i) the same figure must work in
   color and in black and white, and (ii) captions and text must not refer
   to any color.

## 6. Supplementary material and submission

- Suitable content: detailed proofs, step-by-step computational procedures,
  extended data examples, useful software code. Online-only, linked to the
  online manuscript.
- **Add a section titled "Supplementary Material" as the last section of the
  paper, before the Acknowledgments**, with a brief description of the
  online material. (The manuscript template prints it as
  `\section*{Supplementary Materials}`, plural.)
- The supplement must be referred to properly in the main manuscript.
- It is submitted at the same time as the manuscript and is subject to
  review.
- **One single PDF file, at most 10 MB**, with the same title and author
  list as the main manuscript.
- The supplement is not copy-edited; it may be rejected or returned for
  revision.
- Submission through ScholarOne:
  <http://mc04.manuscriptcentral.com/statisticasinica>. Authors without
  access may mail three hard copies with a cover letter to Editors,
  Statistica Sinica, Institute of Statistical Science, Academia Sinica, 128
  Academia Road, Sec. 2, Nankang, Taipei 11529, Taiwan; submitted materials
  are not returned.
- "The Latex file must be provided after the manuscript is accepted." The
  page offers LaTeX templates for Windows and for Mac, plus a template for
  the supplement.
- For notation and abbreviations the page refers to an article by Prof.
  Jordan M. Stoyanov.

## 7. What the pages do not say

Gaps found while transcribing; each needs the journal site, the ScholarOne
system or a recent issue to settle.

1. **Review model.** Single or double blind is not stated anywhere, and
   there is no anonymization requirement. `alvo-revista.md` §3 still carries
   this row as "to confirm".
2. **Fees, copyright and licence.** Nothing on page charges, open access,
   copyright transfer or reuse.
3. **Data and code availability.** Beyond "software producing the evidence
   should be available for examination as well as pertinent datasets" in the
   editorial policy, there is no formal policy, no repository requirement
   and no availability statement.
4. **Cover letter, suggested editors, ORCID, funding and author
   contributions.** Not mentioned; these are probably fields in ScholarOne.
5. **Abstract length and number of key words.** No limit given.
6. **Which template.** The page distinguishes a Windows and a Mac LaTeX
   template; the zip archived here (`SS-template.zip`, downloaded
   2026-09-18) does not say which one it is. Both `.tex` files compile with
   the local TeX Live 2023.
7. **One lost symbol.** In the formulas section, item 2 reads "must be
   avoided, as must [image]": the second forbidden construct failed to load
   in the captured page. Presumably another stacked-accent example.
8. **Links not captured.** The URLs behind "LaTeX Template (for Manuscript)"
   and "(for Supplement)", and the reference to the Stoyanov notation
   article, are not in the PDF text.

## 8. What this changes for the project

- The cap is **40 double-spaced pages including references**, not the 30
  pages recorded from the web search in `alvo-revista.md` §3. The target
  structure in §5 of that document (30 pages) fits with room to spare, but
  the arithmetic there should be redone against the template.
- Proofs go to the supplement by the journal's own preference, which is what
  the plan already assumes (`supp_{k}.tex` from the start).
- Reproducibility is an editorial requirement, not a courtesy: the numbers
  in the paper must come from `wafc-studies`, and the code must be
  presentable on request even while `wafc/` stays private (D4).
- A "Supplementary Material" section before the Acknowledgments is part of
  the skeleton of `ms_1.tex` in E5a.
