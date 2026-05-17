# CLAUDE.md — cp_book

> **Source of truth: [`AGENTS.md`](AGENTS.md) at the repo root.** Read
> it in full before editing anything. Then read the matching
> [`.ai/modules/<folder>.md`](.ai/modules/) for whatever folder you're
> about to touch.
>
> Per-folder docs are **not** placed inside the source folders
> themselves — the source tree (`dp/`, `ds/`, `graph/`, …) is kept
> clean of AI metadata. All per-module notes live under
> [`.ai/modules/`](.ai/modules/).

## Quick links

- Master onboarding: [AGENTS.md](AGENTS.md)
- Per-module docs: [.ai/modules/](.ai/modules/)
  - [dp.md](.ai/modules/dp.md)
  - [ds.md](.ai/modules/ds.md) · [ds-bst.md](.ai/modules/ds-bst.md) ·
    [ds-heap.md](.ai/modules/ds-heap.md)
  - [fft.md](.ai/modules/fft.md)
  - [flow.md](.ai/modules/flow.md)
  - [geo.md](.ai/modules/geo.md)
  - [graph.md](.ai/modules/graph.md)
  - [math.md](.ai/modules/math.md)
  - [misc.md](.ai/modules/misc.md)
  - [mod.md](.ai/modules/mod.md)
  - [string.md](.ai/modules/string.md)
  - [template.md](.ai/modules/template.md)
- Library patches: [library_fixes.md](library_fixes.md)
- Benchmarks: [benchmarks.md](benchmarks.md)

## TL;DR style rules (full version in AGENTS.md §3)

- C++20, GCC ≥ 11, Linux. Header-only. Hand-tuned for speed.
- **Hard tabs**, K&R braces, `snake_case`, `struct` (never `class`).
- CAPS for dimension template parameters and constants (`N`, `M`, `K`,
  `MOD`, `INF`).
- All inputs by `const T&` — *including primitives*.
- Vector-shaped types inherit from `std::vector`.
- Doc-block opens `/**` and closes `**/` (two stars).
- Keep deliberate spaces inside `++ i`, `-- i`, `~ v`, `- 1`,
  `operator <`.
- Asserts wrapped in `#ifdef _GLIBCXX_DEBUG`.
- Always write `std::`. No `using namespace std;` in headers. No
  C++20 `concepts` / `requires`.

## Build / test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG"
cmake --build build -j
ctest --test-dir build --output-on-failure -j

# Benchmarks (opt-in, hidden by default)
./build/tests/bench/ds "[!benchmark]"
```

## Do / Don't (Claude-specific reminders)

- **Do** read the matching `.ai/modules/<folder>.md` before any edit.
- **Do** add a regression test next to any library patch.
- **Don't** add per-folder `CLAUDE.md` / `AGENTS.md` files inside
  source folders — keep `dp/`, `ds/`, `graph/`, … clean.
- **Don't** reformat with `clang-format` or "modernise" code that
  matches the existing style.
- **Don't** run destructive git operations (`push --force`,
  `reset --hard`, `rebase -i` on shared history) without explicit
  permission.
