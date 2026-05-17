# GitHub Copilot — cp_book

> **Source of truth: [`AGENTS.md`](../AGENTS.md) at the repo root and the
> matching [`.ai/modules/<folder>.md`](../.ai/modules/) for whatever
> folder you're editing.** Read those *before* writing any code.

This file is a short recap so the most important rules are always in
Copilot's context window.

## Project shape

- C++20, GCC ≥ 11, Linux x86_64.
- Header-only competitive-programming library. Hand-tuned for speed.
- Build via CMake (`cmake -S . -B build && cmake --build build -j`).
- Tests: Catch2 v3 (fetched). `ctest --test-dir build --output-on-failure -j`.

## Style rules — match exactly, do not "improve"

- **Hard tabs** (never spaces). K&R braces.
- **`snake_case`** everywhere; **`struct`** never `class`.
- **CAPS** for size/dimension template params and constants: `N`, `M`,
  `K`, `MOD`, `INF`. One-capital-letter for type params (`T`, `U`,
  `P`).
- Pass all inputs by **`const T&`** — *including primitives*.
  `graph(const int& N)` is intentional.
- Prefer **inheritance from STL containers** for vector-shaped
  structures (`struct foo : public std::vector<T> { using std::vector<T>::vector; ... };`).
- Doc-block at the top of every header opens with `/**` and closes
  with `**/` (two stars). Title is `STRUCT_NAME !!!`.
- Spaces inside unary `++`/`--`/`~`, around `operator <`, and in
  literals like `- 1`. Don't tighten them.
- Asserts wrapped in `#ifdef _GLIBCXX_DEBUG`.
- Always write `std::` — never `using namespace std;` in headers.
- Use `typename std::enable_if<…>::type` for SFINAE. **Do not** use
  C++20 `concepts` / `requires` (deliberate).
- Headers do **not** include `<bits/stdc++.h>`. Tests include it once.

## Don't

- Don't reformat with `clang-format`.
- Don't add docstrings, comments, or type annotations to code you
  didn't change.
- Don't introduce `shared_ptr`, `<format>`, coroutines, ranges, or
  concepts on hot paths.
- Don't bypass `[!benchmark]` (benchmarks must stay opt-in).
- Don't add `AGENTS.md` / `CLAUDE.md` / README files inside source
  folders (`dp/`, `ds/`, …). All AI metadata lives under `.ai/` or at
  the repo root.

## Do

- Read the matching [`.ai/modules/<folder>.md`](../.ai/modules/) first.
- Add a regression test alongside any library patch (see
  [`library_fixes.md`](../library_fixes.md)).
- Mirror source layout when adding tests
  (`ds/foo.hpp` ⇒ `tests/ds/foo.cpp`).
