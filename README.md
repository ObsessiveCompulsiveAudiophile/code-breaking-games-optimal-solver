# code-breaking-games-optimal-solver
# serkan gur 2025

An even faster C++ classical Mastermind game optimal solver that uses a highly optimized recursive DFS algorithm to find the strategy with the lowest possible average number of moves (4.34028). Code can also be easily customized for different combinations of code breaking games other than 4 digits of 6 colors [4,6]

Some features:

- intelligently prune non-promising branches of the search tree with Maxparts method
- pre-computed score table
- super fast (thanks to Clang) polynomial dot product hashing for equivalence pruning
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

19.9356 seconds

Press Enter to exit...
