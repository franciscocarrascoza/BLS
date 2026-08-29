#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 metrics.csv [output.png]" >&2
  exit 1
fi

INPUT_CSV="$1"
OUTPUT_PNG="${2:-metrics.png}"

if [[ ! -f "$INPUT_CSV" ]]; then
  echo "Error: metrics file '$INPUT_CSV' not found" >&2
  exit 1
fi

gnuplot <<GNUPLOT
set datafile separator ','
set terminal pngcairo size 1200,800 enhanced
set output "${OUTPUT_PNG}"
set multiplot layout 2,1 title "BLS-SkipDFS Frame Metrics"
set xlabel "Frame"
set ylabel "Clusters / Seeds"
set key left top
plot "${INPUT_CSV}" using 1:12 with lines lw 2 lc rgb "#1b9e77" title "Number of clusters", \
     "${INPUT_CSV}" using 1:13 with lines lw 2 lc rgb "#d95f02" title "Max cluster size", \
     "${INPUT_CSV}" using 1:11 with lines lw 2 lc rgb "#7570b3" title "Seed hits", \
     "${INPUT_CSV}" using 1:10 with lines lw 1 lc rgb "#e7298a" title "Seeds tested"

set xlabel "Frame"
set ylabel "Voxels / Time (ms)"
plot "${INPUT_CSV}" using 1:14 with lines lw 2 lc rgb "#66a61e" title "Refined voxels", \
     "${INPUT_CSV}" using 1:15 with lines lw 2 lc rgb "#e6ab02" title "Elapsed ms"
unset multiplot
GNUPLOT

echo "Plot written to ${OUTPUT_PNG}"

