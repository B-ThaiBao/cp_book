# `.ai/modules/ds-bst.md` — balanced BSTs

Files in [`ds/bst/`](../../ds/bst/):

| File | Top-level name | One-liner |
|---|---|---|
| [`treap.hpp`](../../ds/bst/treap.hpp) | `treap` family | Implicit treap (sequence) with `push_back`/`merge`/`split` and lazy stack. Inspired by Tourist / Ecnerwala / Benq submissions referenced in the header. **Patch history**: `stk.do_downdate_back(i)` → `stk.push_back(i)` at four call sites (see [library_fixes.md](../../library_fixes.md)). |
| [`splay_tree.hpp`](../../ds/bst/splay_tree.hpp) | `splay_tree` | Splay tree used as a building block for `link_cut_tree`. |
| [`link_cut_tree.hpp`](../../ds/bst/link_cut_tree.hpp) | `link_cut_tree` | LCT on top of splay tree. Path aggregate + subtree query (see Codeforces blog #80383, #67637, #80145 linked in the doc-block). |

## Module-specific conventions

- BSTs are intrusive — users derive a node struct from a base type
  exported by the header. Don't try to add a "generic value_type" — the
  whole point is letting the user store custom payload + lazy info per
  node.
- Operations are written as **free / static functions** that take
  raw node pointers, not member functions. Pattern follows Ecnerwala /
  Tourist competitive-programming style.
- Tests in `tests/ds/bst/` define their own node struct and exercise
  the public free functions. Use them as the example when adding a new
  BST.
