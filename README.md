# code-breaking-games-optimal-solver
# serkan gur 2025

An ultra fast C++ Mastermind game optimal solution tree generator that uses a highly optimized recursive DFS algorithm to find the strategy with the lowest possible average number of moves (4.34028). Capable of solving any other combination and starting a game from any level!

Some features:

- intelligently prune non-promising branches of the search tree with Maxparts method
- pre-computed score table
- super fast (thanks to Clang compiler) hybrid polynomial dot product & hash table equivalence pruning
- custom tuned vectorization of loops

How to compile and run:

Single .cpp file

Platform toolset: LLVM (clang-cl)

All options:

/GS /W3 /Gy /Zi /O2 /Ob2 /fp:fast /WX- /arch:AVX2 /Gd /Oy /Oi /MD /std:c++20 /Fa






Example end of output:

...

112330  124520  114340  :3

112330  124521  112540  :3

112340  :1

5625 / 1296 = 4.34028 average moves to finish!

16.0356 seconds

Press Enter to exit...
