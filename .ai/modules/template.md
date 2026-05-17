# `.ai/modules/template.md` — contest skeletons

Files in [`template/`](../../template/):

| File | Purpose |
|---|---|
| [`temp.cpp`](../../template/temp.cpp) | Default contest starter. |
| [`fhc.cpp`](../../template/fhc.cpp) | Facebook Hacker Cup variant (multi-test, output to `output.txt`). |
| [`leet.cpp`](../../template/leet.cpp) | LeetCode local harness — class skeleton + `main` driver. |
| [`Makefile`](../../template/Makefile) | One-shot `make foo` to build `foo.cpp` with the project's stock flags. |

## Module-specific conventions

- These files are **deliberately monolithic** — they include
  `<bits/stdc++.h>` and pull in the helpers commonly used by the author
  in contests. Do not refactor them into multiple translation units.
- The Makefile is the canonical way to compile a single contest source
  with the right flags. Mirror its CXXFLAGS if you write any external
  build glue.
- When you adapt a template into a new contest problem, **don't change
  the header section** (debug helpers, IO, y_combinator import, etc.) —
  add your problem-specific code below the `int main()` skeleton.
