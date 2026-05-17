---
trigger: always_on
description: cp_book — competitive-programming C++ library style and rules
---

# Windsurf rules — cp_book

> **Source of truth: `AGENTS.md` at the repo root and the matching
> `.ai/modules/<folder>.md` for whatever folder you're editing.**
> Read those *before* writing any code.

## Project shape

- C++20, GCC ≥ 11, Linux x86_64.
- Header-only competitive-programming library. Hand-tuned for speed.
- Build via CMake. Tests via Catch2 v3.

## Style rules — match exactly, do not "improve"

- **Hard tabs** (never spaces). K&R braces.
- **`snake_case`** everywhere; **`struct`** never `class`.
- **CAPS** for size/dimension template params and constants: `N`, `M`,
  `K`, `MOD`, `INF`. One-capital-letter for type params (`T`, `U`, `P`).
- Pass all inputs by **`const T&`** — *including primitives*.
- Prefer **inheritance from STL containers** for vector-shaped
  structures.
- Doc-block opens `/**` and closes `**/`. Title is `STRUCT_NAME !!!`.
- Keep the deliberate spaces inside `++ i`, `-- i`, `~ v`, `- 1`, and
  around `operator <`.
- Asserts wrapped in `#ifdef _GLIBCXX_DEBUG`.
- Always write `std::`. No `using namespace std;` in headers.
- Use `typename std::enable_if<…>::type` for SFINAE. **Do not** use
  C++20 `concepts` / `requires`.
- Headers don't include `<bits/stdc++.h>`. Tests include it once.

## Don't

- Don't reformat with `clang-format`.
- Don't add docstrings, comments, or type annotations to code you
  didn't change.
- Don't introduce `shared_ptr`, `<format>`, coroutines, ranges, or
  concepts on hot paths.
- Don't bypass `[!benchmark]`.
- Don't add `AGENTS.md` / `CLAUDE.md` / README files inside source
  folders. All AI metadata lives under `.ai/` or at the repo root.

## Do

- Read the matching `.ai/modules/<folder>.md` first.
- Add a regression test alongside any library patch (see
  `library_fixes.md`).
- Mirror source layout when adding tests
  (`ds/foo.hpp` ⇒ `tests/ds/foo.cpp`).
