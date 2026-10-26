# Crossword Generator with a Genetic Algorithm

A Python implementation of a genetic algorithm that arranges a list of words into a crossword-style layout on a 20×20 grid. Each candidate layout is scored on valid intersections, conflicts, touching words, and connectivity, and the population evolves through elitist selection, crossover, and mutation.

Pure Python 3, standard library only — no dependencies to install.

## Quick start

```bash
git clone https://github.com/AyoubElleuch/Crossword-Genetic-Algorithm.git
cd Crossword-Genetic-Algorithm
python Crossword.py
```

The repository ships with 100 sample word lists in [`inputs/`](inputs) (`input1.txt` … `input100.txt`, 5–20 words each). The script processes them in order and writes one result per file to `outputs/output1.txt` … `outputs/output100.txt`, printing generation progress as it goes.

> **Runtime:** with the default settings (population 2,000, 200 generations) a single word list takes roughly 1 minute (5 words) to 5 minutes (20 words) on a typical laptop, so the full 100-file batch is a long run. For a quick demo, lower `population_size` and `number_of_generations` at the top of `Crossword.py`, or set `number_of_files` near the bottom.

`outputs/` must not exist before a run (the script creates it), so delete it before running again. It is listed in `.gitignore`.

## Input and output

**Input** — one word per line:

```text
cat
car
art
tar
rat
```

**Output** — one line per input word, in the same order, formatted as `row column direction`:

- `row`, `column`: zero-based position of the first letter (row grows downward, column to the right)
- `direction`: `0` = horizontal, `1` = vertical

Example output for the list above and the grid it describes (one possible result — runs are random; shown on an 8×8 grid for readability):

```text
3 0 0        --------
1 3 0        ---car--
3 4 1        -t------
2 1 1        cat-a---
6 1 0        -r--r---
             ----t---
             -rat----
             --------
```

## How it works

**Chromosome.** A candidate solution is a list of `(row, column, direction, word)` genes, one per input word. The initial population places every word at a random in-bounds position and orientation.

**Fitness.** Higher is better. With `M` matching shared cells, `I` conflicting cells, `N` touching (parallel-adjacent) cells, `E` occupied cells directly before or after a word, and `C` connected components on the grid:

```text
fitness = M - 3·I - 4·N - 8·E
          + 10  if I == 0
          + 10  if N == 0
          + 20  if C == 1, otherwise - 2^C
```

**Evolution.** Each generation keeps the top half of the population by fitness (elitism), then refills it with children produced by two-point crossover between randomly chosen survivors. Each child has a 15% chance of mutation, which relocates one word to a new random position and orientation. After the final generation the fittest layout is written to the output file.

## Configuration

All settings are module-level variables in [`Crossword.py`](Crossword.py):

| Setting | Default | Meaning |
| --- | --- | --- |
| `population_size` | `2000` | Candidate layouts per generation |
| `number_of_generations` | `200` | Evolution rounds per word list |
| `mutation_rate` | `0.15` | Mutation probability per child |
| `grid_size` | `20` | Width and height of the grid |
| `input_dir` | `'inputs'` | Folder containing `input<N>.txt` |
| `output_dir` | `'outputs'` | Folder created for `output<N>.txt` |
| `number_of_files` | `100` | Number of consecutively numbered inputs to process |

Word lists should contain at least four words (two-point crossover needs them), and no word should be longer than `grid_size`.

## License

[MIT](LICENSE)
