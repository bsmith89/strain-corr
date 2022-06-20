#!/usr/bin/env bash

library=$1
cat \
    <(bzcat "$2" | sed '1,1d') \
    <(bzcat "$3" | sed '1,1d') \
    | cut -f1,2 \
    | sort \
    | uniq \
    | cut -f1 \
    | uniq -c \
    | awk -v OFS='\t' -v library=$library 'NR>1{print library,$2,$1}'
