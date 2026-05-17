# AGENTS.md — cp_book

> Master onboarding for AI coding agents (Claude, Cursor, Copilot,
> Windsurf, Codex…). **Read this entire file before editing anything.**
> Then read the `AGENTS.md` inside the specific folder you're about to
> touch.

---

## 1. Project at a glance

`cp_book` is a personal **competitive-programming C++ algorithm
library**. It is header-only, hand-tuned for speed, and follows a very
specific in-house style. The library author considers the style
**canonical** — your generated code MUST match it exactly. Do not
"modernize", "refactor", or "clean up" the existing code unless the
user explicitly asks.

- **Language**: C++20, GCC ≥ 11, Linux x86_64.
- **Build**: CMake ≥ 3.20 (see [`CMakeLists.txt`](CMakeLists.txt)).
- **Tests**: Catch2 v3 (fetched via `FetchContent`).
- **pb_ds** (`<ext/pb_ds/...>`) is used freely (Linux-only).

---

## 2. Repository map

| Source folder | Purpose | Per-folder doc |
|---|---|---|
| `dp/` | DP optimisations (CHT, monotonic structures) | [.ai/modules/dp.md](.ai/modules/dp.md) |
| `ds/` | Data structures (BIT, seg tree, sparse table, sets, …) | [.ai/modules/ds.md](.ai/modules/ds.md) |
| `ds/bst/` | Balanced BSTs (treap, splay, link-cut) | [.ai/modules/ds-bst.md](.ai/modules/ds-bst.md) |
| `ds/heap/` | Priority-queue variants (radix, skew) | [.ai/modules/ds-heap.md](.ai/modules/ds-heap.md) |
| `fft/` | FFT / NTT / Berlekamp–Massey / Kitamasa | [.ai/modules/fft.md](.ai/modules/fft.md) |
| `flow/` | Max flow + min-cost flow (`flow_graph` based) | [.ai/modules/flow.md](.ai/modules/flow.md) |
| `geo/` | Computational geometry | [.ai/modules/geo.md](.ai/modules/geo.md) |
| `graph/` | Graph algorithms (SCC, matching, MST, 2-SAT, …) | [.ai/modules/graph.md](.ai/modules/graph.md) |
| `math/` | Number theory & matrix | [.ai/modules/math.md](.ai/modules/math.md) |
| `misc/` | Utilities (compressor, y_combinator, tensor, …) | [.ai/modules/misc.md](.ai/modules/misc.md) |
| `mod/` | `modnum` + multiplier policies (Barrett, naive) | [.ai/modules/mod.md](.ai/modules/mod.md) |
| `string/` | Suffix array, Z, KMP, manacher, hashing | [.ai/modules/string.md](.ai/modules/string.md) |
| `template/` | Contest skeleton `.cpp` files | [.ai/modules/template.md](.ai/modules/template.md) |
| `tests/` | Catch2 unit + stress + benchmark sources | — (mirrors source layout) |
| `debug.hpp` | Pretty-printers used inside tests | — |

When you need to edit a file under any folder above, **first open the
matching doc in [.ai/modules/](.ai/modules/)**. It lists every file in
the folder, the top public name(s), and any quirks specific to that
module. The `.ai/` directory is the *only* place AI context lives —
the source folders themselves are kept clean.

---

## 3. Code style — non-negotiable

The library style is consistent across every file. Match it byte for
byte. New code that violates these rules will be rejected.

### 3.1 Formatting

- **Hard tabs** for indentation (never spaces). One tab = one level.
- **K&R braces**: opening brace on the same line.
- Single-statement bodies frequently kept on one line:
  `if (x < md) return T(x);`
- **No trailing whitespace, no final newline noise.**
- Headers do **not** include `<bits/stdc++.h>`; tests do (exactly once
  at the top).

### 3.2 Naming

- **`snake_case`** for absolutely everything: types, functions,
  variables, members, namespaces.
- **No `class`** — every aggregate is `struct` (all members public).
- **CAPS for size/dimension template parameters and constants**:
  `N`, `M`, `V`, `K`, `MOD`, `BARRETT_M`, `INF`. Lowercase letters for
  the type parameter itself (`T`, `num`, `Vector`, `Fun`, …).
- **Template parameters**: one-capital-letter is the norm
  (`T`, `U`, `P`, `Q`, `F`, `K`), with descriptive lowercase camel only
  when meaningful (`Vector`, `String`, `Order`, `Network`, `Fun`).
- **Namespaces are lowercase**, ending comment is required:
  `} // namespace foo`.

### 3.3 The "spaces inside operators" rule (deliberate)

You will see — and must reproduce — these spacings throughout:

```cpp
for (int i = N - 2; i >= 0; -- i) { ... }   // space inside -- / ++
sa[first_endpoint = -- tmp[c1]] = i + 1;
return n - 1;                                // not n-1
res[1][res[0].back()] = N + 1;
~ v                                          // unary ~ with a space
- 1                                          // bare -1 literal often spaced
friend inline bool operator < (const A&, const B&);   // spaced operator
```

This is **on purpose** — do not "tighten" it. Both `++i` and `++ i`
appear in older code; for *new* code prefer the spaced form to match
the dominant style.

### 3.4 Parameters: `const T&` everywhere (yes, even ints)

```cpp
graph(const int& N) : V(N), adj(N) {}
inline point_t c(bool z) const { return point_t((a << 1) | z); }
indexed_set(const int& M) { build(M); }
inline void insert(int i) { ... }     // by value when mutated locally
```

Pass by `const T&` for inputs (including primitives), by value when the
parameter is mutated inside the function. Do **not** insert
`const int` (no reference) for primitives — the codebase uses
`const int&` to keep one uniform style.

### 3.5 Inheritance from STL containers (CRTP-style reuse)

A very common pattern — reuse it for any new "vector-like" structure:

```cpp
template <typename T> struct binary_indexed_tree : public std::vector<T> {
	using std::vector<T>::vector;     // inherit ctors
	// ...
};

struct suffix_array : public std::vector<int32_t> { ... };
template <typename num> struct combinatorics : std::vector<std::array<num, 2>> { ... };
template <typename Line, typename Comp = std::less<>>
struct line_multiset : public std::multiset<Line, Comp> { ... };
```

Pros: free `.size()`, iterators, range-for, ctor forwarding. Use this
unless there's a real reason not to.

### 3.6 Doc-block at the top of every header

```cpp
/**
 * STRUCT_NAME !!!
 *
 * Description: short paragraph.
 *
 * Usage:
 *   * step 1: ...
 *   * step 2: ...
 *
 * NOTE: any sharp edges.
**/
```

The opening is `/**` and the **closing is `**/`** (two stars, not one).
Title is the struct/function name in ALL CAPS with `!!!` (or `!!!!`).
Reference URLs go inline as comments. Time complexity goes as
`// Time: O(...)` either above the struct or above each non-trivial
method.

### 3.7 Operator overloads

`friend` definitions inside the struct are preferred when the operator
needs access to private state of multiple types:

```cpp
friend inline bool operator < (const line& v, const T& o) { return v.x < o; }
friend inline bool operator < (const T& o, const line& v) { return o < v.x; }
```

Always provide both directions for asymmetric comparisons.

### 3.8 Asserts and debug

Wrap asserts behind `_GLIBCXX_DEBUG` so Release builds stay free of
runtime checks:

```cpp
#ifdef _GLIBCXX_DEBUG
	assert(0 <= int(N) && int(N) < int(this->size()));
#endif
```

Plain unconditional `assert(...)` is reserved for invariants that are
truly free at runtime.

### 3.9 `std::` is always written out

Headers never `using namespace std;`. Tests may, but the existing tests
do not — match local convention.

### 3.10 SFINAE / template constraints

Use the classic `typename std::enable_if<..., R>::type` return-type
trick. C++20 `concepts` / `requires` are intentionally **not** used in
this codebase yet — do not introduce them.

```cpp
template <typename U, typename P, typename Q, typename K = T>
constexpr static typename std::enable_if<sizeof(K) == 4, T>::type
multiply(const P& a, const Q& b, U md) { ... }
```

---

## 4. Build & test commands

```sh
# Release (for benchmarks)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG"
cmake --build build -j

# Run every unit + stress test
ctest --test-dir build --output-on-failure -j

# Run a single test binary
./build/tests/ds/swag

# Benchmarks (opt-in via tag, not run by default ctest)
./build/tests/bench/ds "[!benchmark]"
```

Catch2 tags used in this repo:

- `[<module>]` — primary module tag (e.g. `[bit]`, `[sa]`, `[scc]`)
- `[stress]` — randomised cross-check vs brute force
- `[!benchmark]` — Catch2 `BENCHMARK` blocks, **hidden by default**
- `[.skip]` — intentionally skipped (documents a known library bug)

---

## 5. Pre-existing fixes and known bugs

[`library_fixes.md`](library_fixes.md) lists every library patch that
has been applied and why. [`benchmarks.md`](benchmarks.md) lists every
benchmark and its measured time. Read both before claiming "this code
is wrong" — the answer is almost always already documented.

---

## 6. Operational guardrails for the agent

- **Don't** add docstrings, comments, or type annotations to code you
  didn't change.
- **Don't** reformat with `clang-format` — this code has its own dense
  style.
- **Don't** introduce `std::shared_ptr`, `<format>`, coroutines,
  ranges, or `concepts` into hot paths.
- **Don't** replace `<bits/stdc++.h>` in tests with explicit includes.
- **Don't** bypass the `[!benchmark]` tag — benchmarks must stay opt-in.
- **Don't** run `git push`, `git reset --hard`, `git rebase -i` on
  shared history, or any destructive git op without an explicit ask.
- **Do** add a regression test next to any library patch.
- **Do** mirror the source folder layout when adding tests
  (`ds/foo.hpp` ⇒ `tests/ds/foo.cpp`).
- **Do** read the matching [.ai/modules/](.ai/modules/) doc before
  editing any file in a source folder.
- **Do not** add `AGENTS.md`, `CLAUDE.md`, `README`, or any other meta
  file inside the source folders themselves (`dp/`, `ds/`, `graph/`,
  …) — keep them clean. All AI/agent metadata lives under `.ai/` or
  at the repo root.

---

## 7. Tool routing

| Tool | Config file | Status |
|---|---|---|
| GitHub Copilot | `.github/copilot-instructions.md` | mirrors this file |
| Cursor | `.cursor/rules/cp_book.mdc` | `alwaysApply: true` |
| Windsurf | `.windsurf/rules/cp_book.md` + `.windsurfrules` | always-on |
| Claude Code | `CLAUDE.md` (root) + `.ai/modules/*.md` | reads automatically |
| Codex / generic | `AGENTS.md` (this file) + `.ai/modules/*.md` | reads automatically |

**No `AGENTS.md` lives inside the source folders themselves** — every
per-module doc is under [.ai/modules/](.ai/modules/). If you are a tool
not listed above, default to reading **this file** plus the matching
`.ai/modules/<folder>.md` for whatever folder you are editing.
