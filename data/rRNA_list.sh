#!/bin/bash

rna_types=("5S" "5_8S" "18S" "28S")

col3_values=("0.8" "0.85" "0.9")

for col2_val in $(seq 55 65); do
  for rna in "${rna_types[@]}"; do
    for col3_val in "${col3_values[@]}"; do
      echo -e "${rna}\t${col2_val}\t${col3_val}\t${col2_val}"
    done
  done
done

