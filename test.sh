#!/bin/sh
script_dir=$(dirname $0)
bin_dir="$script_dir/../../bin"

cd "$script_dir"

make clean
make

("$bin_dir/meshgen_hex" 2 2 2 1 1 1 -d "$script_dir")

# ("$script_dir/block_partitioner" 2 2 2 -ig "node.dat" -e 1e-6)