#!/bin/bash

## linearlizes fasta to a tab delimited output, with headers in first column,
## and sequence in second. Removes empty lines as well. Output to stdout by 
## default.
cat $1 | sed 's/ /,/g' | \
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' \
| awk 'BEGIN{RS=">"}{print $1"\t"$2;}' | sed -r '/^\s*$/d'

