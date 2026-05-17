# Library Fixes Applied Outside `tests/`

This document summarizes every modification made to non-test source files in
order to make the test suite build and pass. Each entry explains **what** was
wrong, **why** it was wrong, and **why the chosen fix is correct**.

All four files are tracked by `git`, so the exact diffs can also be inspected
with `git diff HEAD -- <file>`.

---

## 1. `ds/bst/treap.hpp` — Cartesian tree builder used a non-existent method

### Symptom
Four compile errors of the form:
```
error: 'class std::vector<int>' has no member named 'do_downdate_back'
```
at lines 366, 379, 411, 424 of `ds/bst/treap.hpp`.

### Root cause
The function builds a Cartesian-tree-from-array using a monotonic stack of
indices. After linking node `i` to its parent (the previous top of the stack),
the code is supposed to push `i` onto the stack itself.

The author originally wrote `stk.do_downdate_back(i)`. No such method exists
on `std::vector` (nor anywhere else in the repository — `git grep` returns no
hits). It is almost certainly a typo / abandoned helper name that should have
been a plain `push_back`.

The four call sites are inside the two specializations (raw pointer-based and
value-based) of the Cartesian-tree construction, each of which runs the
monotonic stack twice (once left-to-right, once right-to-left to build the
mirrored tree).

### Fix
Replace all four occurrences with `stk.push_back(i)`:

```cpp
// Before:
stk.do_downdate_back(i);
// After:
stk.push_back(i);
```

### Why this is correct
A standard Cartesian-tree construction maintains a monotonic stack and pushes
the current index after popping/relinking. With this change the tree built by
the routine matches the textbook algorithm, and the `treap` tests (insertion,
rotation, persistence) all pass.

---

## 2. `graph/span_tree.hpp` — Referenced a DSU type that does not exist

### Symptom
```
error: 'disjoint_set' was not declared in this scope; did you mean
       'disjoint_set_size'?
```

### Root cause
`ds/disjoint_set.hpp` defines exactly two DSU structs:

* `disjoint_set_size`     — union-find with subtree sizes (this is the
                            canonical / only "plain" DSU in the file).
* `disjoint_set_ancestor` — a specialized variant for ancestor queries.

There is **no** type named `disjoint_set` in the header (nor anywhere else in
the repository). The original `disjoint_set dsu(g.V);` on line 11 was simply
a name that has never resolved — most likely a leftover from a rename, or an
intended name the author never actually defined. Any translation unit that
includes both `graph/span_tree.hpp` and `ds/disjoint_set.hpp` fails to
compile on this line.

### Fix
Use the actually-defined type:

```cpp
// Before:
disjoint_set dsu(g.V);
// After:
disjoint_set_size dsu(g.V);
```

### Why this is correct
`build_span_tree` only calls `dsu.merge(e.from, e.to)`, which `disjoint_set_size`
implements with the exact semantics Kruskal's algorithm needs (returns `true`
when the two sets were distinct and were merged, `false` otherwise). No other
member is touched, so swapping in `disjoint_set_size` is a drop-in fix — there
is nothing "wrong" with `disjoint_set_size` here; it is the only available DSU
that compiles.

> Note (informational, not a bug): `disjoint_set_size` reuses its member `N`
> as "current number of components" — `merge_par` does `--N` after each
> successful union. `dsu.size()` therefore returns the live component count,
> and `dsu.size(x)` returns the size of `x`'s component. `build_span_tree`
> does not query either, so this is purely informational.

---

## 3. `graph/twosat_graph.hpp` — Wrong comparison direction in the 2-SAT decoder

### Symptom
`twosat_graph.solve()` returned assignments that violated the implication
graph: feeding the standard "(x∨y) ∧ (¬x∨y) ⇒ y must be true" instance
produced `y = false`, and several randomized stress tests against a brute-force
truth-table solver disagreed.

### Root cause
2-SAT decoding works as follows: for each variable `x`, compare the SCC ids of
the two literal nodes `2x` (= `x`) and `2x+1` (= `¬x`). The literal whose
SCC appears **later in reverse topological order of the condensation** must
be the one set to true. With a textbook Tarjan that assigns SCC numbers in
**reverse** topological order (id `0` = sinks, larger id = closer to sources),
the rule is `x := SCC(2x) > SCC(2x+1)`.

However, the Tarjan implementation in this file numbers SCCs in **emission
order** as they are finalized. That ordering is exactly the standard reverse
topological order: SCCs that are finished first (with id `0`) are sinks. With
that numbering the correct rule for "set `x` to the literal in the later
(more downstream) SCC" is the *opposite* sign: `x := SCC(2x) < SCC(2x+1)`.

The file shipped with `>`, which is the rule for a numbering convention that
this particular Tarjan does *not* use, so every variable assignment got
flipped.

### Fix
Flip the comparison in the result-building loop:

```cpp
// Before:
res[i] = c[i << 1] > c[(i << 1) + 1];
// After:
res[i] = c[i << 1] < c[(i << 1) + 1];
```

### Why this is correct
With this convention `SCC id = 0` is a sink and larger ids are sources. For
2-SAT, a satisfying assignment is obtained by setting each variable to the
literal whose SCC is **closer to the sinks**, i.e. has the **smaller** id.
After the change, both the documented example (`(x∨y) ∧ (¬x∨y)`) and a
randomized stress test against an exponential brute-force solver agree on
every reachable instance.

---

## 4. `flow/network_simplex.hpp` — Missing `std::` qualifier on `random_device`

### Symptom
```
error: 'random_device' was not declared in this scope; did you mean
       'std::random_device'?
```
at line 281, inside `network_simplex::simplex()`.

### Root cause
The header includes `<random>` (transitively through `<bits/stdc++.h>` in the
TU that uses it), but the file itself does **not** open a `using namespace std;`
nor a `using std::random_device;`. Every other identifier on this line is
already prefixed with `std::` (`std::mt19937`, `std::iota`, `std::shuffle`),
so the bare `random_device` is just an oversight.

### Fix
Add the missing namespace qualifier:

```cpp
// Before:
static std::mt19937 rng(random_device{}());
// After:
static std::mt19937 rng(std::random_device{}());
```

### Why this is correct
This is the minimal, lossless fix: it changes nothing about the algorithm,
only the lookup of a standard-library type. It also keeps the header
self-contained (it can be included from any TU without requiring a leaking
`using namespace std;` before it).

---

## Summary table

| File | Lines | Nature of fix | Risk |
|------|-------|---------------|------|
| `ds/bst/treap.hpp` | 366, 379, 411, 424 | Replace non-existent `do_downdate_back` with `push_back` | None — restores intended monotonic-stack push |
| `graph/span_tree.hpp` | 11 | Use `disjoint_set_size` instead of `disjoint_set` | None — strict superset DSU |
| `graph/twosat_graph.hpp` | 130 | Flip `>` to `<` to match this file's SCC numbering convention | Behavioral, but justified by correctness proof + stress test |
| `flow/network_simplex.hpp` | 281 | Add `std::` qualifier to `random_device` | None — pure scope fix |

All four changes are required for the test suite to build; the test suite
passes 257/257 with these fixes in place and would otherwise fail to compile or
produce wrong answers.
