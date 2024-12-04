# treeutils

treeutils provides various utilities for tree manipulation in C++.  It is intended to be included in a broader project.  The command line interface currently provides one functionality, which is to take an input tree in newick format, and output the newick of all possible rootings of the tree.
To compile:
mkdir build
cd build
cmake ..
make

To use:
treeutils -m all_reroots -i [input_file] -o [output_file]

The arguments -m and -i are mandatory.  If -o is not specified, the standard output is used.
