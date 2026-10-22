# Crossword generator using a genetic algorithm

`Crossword.py` searches for word placements on a square grid. It requires Python
3 and uses only the standard library (`random`, `sys`, and `os`); no packages need
installing. The repository does not include input datasets.

## Setup and batch execution

1. Install Python 3 and open a terminal in the repository directory.
2. Create `inputs/input1.txt` through `inputs/input100.txt`.
3. Ensure `outputs` does **not** already exist: the script calls
   `os.makedirs(output_dir)` without `exist_ok=True`.
4. Run `python Crossword.py`.

Paths are relative to the current working directory. The script prints generation
progress and an `Output for file N generated.` message for each completed file.
It writes `outputs/output1.txt` through `outputs/output100.txt` sequentially.
An exception interrupts the batch, leaving any completed files behind. Move or
remove that output directory before rerunning.

There are no command-line options. Change these settings in `Crossword.py`:

| Setting | Default | Meaning |
| --- | --- | --- |
| `population_size` | `2000` | Initial candidate count |
| `number_of_generations` | `200` | Reproduction rounds per file |
| `mutation_rate` | `0.15` | Mutation probability per child |
| `grid_size` | `20` | Both grid dimensions |
| `input_dir` | `'inputs'` | Input directory |
| `output_dir` | `'outputs'` | New output directory |
| `number_of_files` | `100` | Consecutively numbered input count |

The final `main(...)` call is unconditional: importing the module starts the
batch too. `main(input_dir, output_dir, num_files)` ignores `num_files` and uses
the global `number_of_files` instead.

## Input and output

Each input contains one word per line. `read_files` strips surrounding whitespace
from every line, retaining blank lines as empty strings. It does not change
case, deduplicate words, or validate characters. Files use Python's default text
encoding.

The original input description specifies 5–20 words per file and at most 20
letters per word; these limits are not validated. Use nonempty words, at least
four words for crossover, and word lengths no greater than `grid_size`. Avoid
`-`, which represents empty grid cells.

Example `inputs/input1.txt`:

```text
CAT
CAR
ART
TAR
RAT
```

Each output line contains `x y direction`, in input order:

- `x`: zero-based **row**, increasing downward.
- `y`: zero-based **column**, increasing to the right.
- `0`: **horizontal**; `1`: **vertical**.

The output contains neither words nor a rendered grid. Match each row with the
corresponding input line. The reduced example below produced:

```text
3 0 0
1 3 0
3 4 1
2 1 1
6 1 0
```

Here `CAT` starts at row 3, column 0 and extends horizontally. These placements
render on an 8 by 8 grid as follows. `build_grid` and `print_grid` provide this
view, but normal batch output does not include it.

```text
--------
---CAR--
-T------
CAT-A---
-R--R---
----T---
-RAT----
--------
```

## Chromosome representation and initialization

A chromosome is a list of `(x, y, direction, word)` tuples, one per input line,
in input order. Internal directions are `'horizontal'` and `'vertical'`, rather
than the output file's numbers.

Initialization samples random rows, columns, and orientations, accepting a word
if its endpoint fits in the grid. It tries up to 100 times per word. If any word
fails, it discards and restarts the entire candidate, continuing until the
population is filled. It does not check intersections or adjacency.

`build_grid` starts with `-` cells and writes words in chromosome order. Later
words overwrite earlier letters at overlaps. It skips words extending beyond
the positive grid edge; it is not a general validator for arbitrary tuples.

## Fitness calculation

Higher scores are preferred. `calculate_fitness(chromosome, grid_size, word_list)`
does not use `word_list`. For each word, it builds a grid of all the other words
and examines each letter:

- A matching occupied cell adds **1**.
- A mismatching occupied cell counts as an invalid intersection, penalized **3**.
- Only when the cell is empty in that other-word grid, each occupied neighbor
  perpendicular to the word counts as a neighbor, penalized **4**.

Using the full rendered grid, an occupied cell immediately before or after each
word along its direction incurs **8** penalty points, if that cell is in bounds.
The score gains **10** if there are no invalid intersections and **10** if there
are no counted perpendicular neighbors.

A recursive depth-first search counts connected components of occupied cells
using four-way adjacency. Exactly one component adds **20**; otherwise the score
loses **2 raised to the number of components**. This counts touching cells,
including conflicts and adjacency, rather than valid word crossings.

With matching-cell count `M`, invalid count `I`, neighbor count `N`, occupied
word-end count `E`, and component count `C`:

```text
fitness = M - 3*I - 4*N - 8*E
          + (10 if I == 0 else 0)
          + (10 if N == 0 else 0)
          + (20 if C == 1 else -(2**C))
```

Counts come from per-word scans, not unique pairs: a simple two-word matching
crossing is normally counted twice. When several words occupy one cell, grid
overwriting affects which letter is compared.

## Selection, crossover, and mutation

Each generation sorts candidates by descending fitness and retains the first
`elitism` candidates unchanged. Batch defaults retain 1,000 of 2,000 candidates.
Each of 1,000 new children uses two independently, uniformly chosen parents from
these survivors; both choices can be the same parent. Selection within the
surviving group is not weighted by fitness.

Crossover samples two distinct indices from `range(1, len(parent_1) - 1)`, sorts
them, and copies the slice between them from parent 2 into parent 1. It produces
one child, preserving input order. The first and last genes always come from
parent 1. Sampling requires at least four genes.

With probability `mutation_rate`, a child has one randomly chosen gene relocated.
Mutation samples a new row, column, and orientation, making up to 100 attempts
to find an in-bounds endpoint. If none succeeds, the gene stays unchanged. The
word itself is preserved; mutation may choose its existing placement. Neither
crossover nor mutation repairs conflicts or adjacency.

For direct calls with a different population size, pass `elitism` explicitly.
Its default is evaluated as `population_size // 2` when the function is defined
(1,000 with original settings), not from the call's `population_size` argument.
Each round retains up to `elitism` survivors and creates exactly `elitism`
children; custom values need not preserve the initial population size.

## Stopping criteria and small example

The search performs exactly `number_of_generations` rounds, then sorts the final
population and returns its highest-scoring chromosome. There is no early exit
for validity, target score, convergence, or elapsed time. Normal batch execution
does not seed the random generator, so results vary.

For a small run without editing the script or starting its batch, save this as
a separate Python file in the repository and run it. It loads the definitions
while omitting only the final batch call:

```python
import ast
from pathlib import Path
import random

tree = ast.parse(Path('Crossword.py').read_text())
tree.body.pop()  # Omit the final, unconditional main(...) call.
api = {}
exec(compile(tree, 'Crossword.py', 'exec'), api)

random.seed(0)
words = ['CAT', 'CAR', 'ART', 'TAR', 'RAT']
best = api['genetic_algorithm'](
    words, population_size=20, number_of_generations=5,
    mutation_rate=0.15, grid_size=8, elitism=10,
)
print('Fitness:', api['calculate_fitness'](best, 8, words))
for x, y, direction, word in best:
    print(x, y, 0 if direction == 'horizontal' else 1)
api['print_grid'](api['build_grid'](best, 8))
```

This example was run successfully, producing the placements and grid above with
fitness **6.0**. An isolated run of the original batch entry point, with the same
seed and reduced settings and `number_of_files = 1`, wrote the same five rows.
This verifies a functional example; it is not a performance benchmark or evidence
of a valid crossword. Exact seeded results can depend on Python's random
implementation/version.

## Limitations

- Fitness is a heuristic. Returned layouts can contain disconnected words,
  mismatched overlaps, touching words, or same-direction overlaps. There is no
  final validity check or guarantee of finding a valid crossword.
- Impossible inputs, such as a word longer than the grid, can cause initialization
  to restart indefinitely. There is no overall retry limit or input validation.
  Fewer than four words cause crossover to fail when reproduction runs.
- Case is significant. Duplicate lines remain separate genes. Blanks and `-`
  characters have unintended interactions with placement and empty-cell handling.
- Fitness repeatedly rebuilds grids and recursively traverses occupied cells.
  Larger configurations may be costly or reach Python's recursion limit. No
  timing measurements or performance guarantees are provided.
- Batch paths, counts, and parameters require source edits. There is no import
  guard, checkpointing, exception recovery, or automated test suite in the repo.
