#!/bin/bash
awk -F'\t' 'NR>1{print "data/marker_pat/"$1".pat"}' meta/samples.tsv
ls data/mix/*__d*.pat
