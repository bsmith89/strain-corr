#!/usr/bin/env bash
#
bzcat $3 \
    | awk -v OFS='\t' -v species=$1 -v library=$2 '$1==species {print library,$0}'
