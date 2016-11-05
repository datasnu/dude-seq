#!/bin/bash

dat_file=$1
stub="${dat_file%.dat}"

k_2=2 # homopolymer errors
k_1=5 # substitution errors


# formatting
cmd_str="Scripts/DAT_INT.pl < $dat_file > ${stub}.dat.int"
echo "$cmd_str"
$cmd_str

# correcting homopolymer errors
echo "correcting homopolymer errors..."
cmd_str="bin/DUDE-Seq-2 -k $k_2 -i ${stub}.dat.int -o ${stub}_H${k_2}.dat.int"
echo "$cmd_str"
$cmd_str

# base-calling
fasta_file="${stub}_H${k_2}.fasta"
cmd_str="Scripts/DAT_FASTA.pl ${stub}.dat.id ${stub}_H${k_2}.dat.int"
echo "$cmd_str"
$cmd_str > $fasta_file

# correctiong substitution errors
echo "correcting substitution errors..."
cmd_str="bin/DUDE-Seq-1 -k $k_1 -i $fasta_file -o ${fasta_file%.fasta}_S${k_1}.fasta"
echo "$cmd_str"
$cmd_str
