#!/bin/bash
# Test script for single-chromosome testing with whole-genome statistics
#
# Usage:
#   Step 1: Run stats-only mode on full genome to get WG stats
#   Step 2: Run single-chrom with the extracted WG stats
#
# This enables fast iterative testing on a single chromosome while using
# statistics derived from the whole genome.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HMCNC="${SCRIPT_DIR}/../src/hmcnc"
TEST_DATA="${SCRIPT_DIR}/data"

# Check that hmcnc is built
if [ ! -x "$HMCNC" ]; then
    echo "ERROR: hmcnc not found at $HMCNC"
    echo "Please build first: cd src && make"
    exit 1
fi

# Example 1: Extract stats from coverage bed (stats-only mode)
echo "=== Example 1: Stats-only mode ==="
echo "Extracting whole-genome statistics from chr22 test data..."
echo ""

$HMCNC "${TEST_DATA}/chr22.fa" \
    -b "${TEST_DATA}/chr22.bed" \
    -s "${TEST_DATA}/chr22.snv" \
    --stats-only \
    -P results/test_stats.param \
    2>&1

echo ""
echo "Parameter file written to results/test_stats.param"
echo ""

# Example 2: Run HMM with pre-computed stats
echo "=== Example 2: Run with pre-computed WG stats ==="
echo "Running HMM on chr22 with explicitly provided stats..."
echo ""

# These would be extracted from a whole-genome run
WG_MEAN=30.0
WG_VAR=50.0
WG_CLIP_MEAN=0.5
WG_CLIP_VAR=2.0

$HMCNC "${TEST_DATA}/chr22.fa" \
    -b "${TEST_DATA}/chr22.bed" \
    -s "${TEST_DATA}/chr22.snv" \
    --wg-mean $WG_MEAN \
    --wg-var $WG_VAR \
    --wg-clip-mean $WG_CLIP_MEAN \
    --wg-clip-var $WG_CLIP_VAR \
    -o results/test_output.vcf \
    -P results/test_trained.param \
    2>&1 | head -50

echo ""
echo "VCF output written to results/test_output.vcf"
echo "Trained parameters written to results/test_trained.param"
echo ""

# Example 3: Workflow - full genome stats, then single chrom run
echo "=== Example 3: Two-step workflow ==="
echo "Step 1: Extract stats from full genome coverage bed"
echo "  hmcnc ref.fa -b full_genome.cov.bed --stats-only"
echo ""
echo "Step 2: Run HMM on single chromosome with those stats"
echo "  hmcnc chr22.fa -b chr22.cov.bed \\"
echo "    --wg-mean <mean> --wg-var <var> \\"
echo "    --wg-clip-mean <clip_mean> --wg-clip-var <clip_var> \\"
echo "    -o chr22.vcf"
echo ""
echo "This allows iterative development on a single chromosome"
echo "while maintaining statistics from the whole genome."

echo ""
echo "=== All tests completed ==="
