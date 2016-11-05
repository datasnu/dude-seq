#!/bin/bash

fasta_file=$1
stub=${fasta_file%.fasta}

k_1=5 # substitution errors


# correctiong substitution errors
echo "correcting substitution errors..."
cmd_str="bin/DUDE-Seq-1 -k $k_1 -i $fasta_file -o ${stub}_S${k_1}.fasta"
echo "$cmd_str"
$cmd_str
