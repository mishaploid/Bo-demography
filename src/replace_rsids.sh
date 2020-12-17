#!/bin/bash

# quick awk command to get bim file formatted with fake rsids
# helpful post here https://www.biostars.org/p/251598/

awk 'BEGIN{FS=OFS="\t"}{$2=$1":"$4":"$5":"$6;print}' $1
