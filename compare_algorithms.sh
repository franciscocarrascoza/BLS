#!/bin/bash
# Compare all clustering algorithms on validation dataset

TRAJ="validation/g1000_20Ih.pdb"
CONF="validation/opt.in"
OUT_DIR="/tmp/algorithm_comparison"

mkdir -p $OUT_DIR

echo "=== Algorithm Comparison on 20-cluster Ice Ih Dataset ==="
echo ""

algorithms=("bls" "traditional_dfs" "skip_dfs" "dbscan" "hierarchical" "kmeans" "spectral" "gcbd" "cc3d" "hdbscan")

for algo in "${algorithms[@]}"; do
    echo "Testing: $algo"
    output="$OUT_DIR/${algo}_results.csv"

    case $algo in
        "hdbscan")
            ./build/bls_analyze --traj $TRAJ --conf $CONF --algo $algo \
                --algo-minclustersize 10 --algo-minsamples 5 --out $output --quiet
            ;;
        "cc3d")
            ./build/bls_analyze --traj $TRAJ --conf $CONF --algo $algo \
                --algo-connectivity 6 --out $output --quiet
            ;;
        *)
            ./build/bls_analyze --traj $TRAJ --conf $CONF --algo $algo \
                --out $output --quiet
            ;;
    esac

    if [ $? -eq 0 ]; then
        result=$(tail -1 $output | awk -F, '{printf "  Clusters: %-3s | Max: %-4s | Voxels: %-5s | Time: %s ms\n", $12, $13, $14, $15}')
        echo "$result"
    else
        echo "  FAILED"
    fi
    echo ""
done

echo "Results saved to: $OUT_DIR"
