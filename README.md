# Nonogram Solver

[Nonograms](https://en.wikipedia.org/wiki/Nonogram) are puzzles that follow a set of simple rules and require intermediate pattern recognition skills to solve. This solver eschews those typical patterns by explicitly checking all feasible solutions and reducing the domain of solutions until the correct one is found. The solver has minimal complexity, i.e., no lookahead. It only employs a basic search heuristic based on the confidence of the value of a particular cell (calculated using the entropy).

The code right now is a mess but works well enough. I find that puzzles up to around 40x40 with a "dense" solution (half or more of the cells are black) take less than a minute to solve on my old laptop.

I was hoping to employ probabilistic programming in the future to augment this solver considerably, but that depends on my level of patience with working with Julia.

## Input format

The input file's format is as follows:

R C

r1

r2

r3

...

c1

c2

c3

...

where R and C are the numbers of rows and columns respectively, and r#/c# are the specifications of "cell runs" for each row and column (read left-to-right and top-to-bottom). I have included a test file with moderate difficulty, [`shark.txt`](https://github.com/lukasmericle/nonogram-solver/blob/master/shark.txt). 

## Running the code

To run the Julia code, I usually just JIT compile by typing `julia nonogram.jl shark.txt` into the terminal. The result will print after every update and so takes longer than it probably could.
