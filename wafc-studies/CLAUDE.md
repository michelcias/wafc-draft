# wafc-studies

Reproducibility compendium of the simulation study and the data
applications of the WAFC article.
Read `README.md` (what it reproduces, how to run it) and `INSTRUCTIONS.md`
(how it is kept reproducible) before changing anything.

## Rules

1. **Commits never carry co-authorship.** No `Co-Authored-By:`, no
   "Generated with", no trailer naming the assistant. Commit and push only
   when the author asks; never rewrite history.
2. **Answers in the chat are as short as possible:** the verdict, the
   reason in a line or two, and the code when it helps.
3. **Run from the root of the compendium**, so that `renv` is active.
4. **No long run without the author's go-ahead.** A run of more than a few
   minutes is prepared, its cost estimated, and started only when the
   author says so.
5. **The configuration is data.** A grid, a cell, a seed or a method
   changes only when the author decides it; cached units that would mix
   two configurations are refitted, never kept beside each other.
6. **The WAFC code is not edited from here.** The compendium loads it
   (`code` in `config/study.yaml`); a change it needs is reported, not
   made.
7. **Documentation and comments are written as those of the final
   study,** in English: what the code does and why, not how it came to be.
8. Outputs and caches are not versioned (`outputs/`, `renv/library/`).
