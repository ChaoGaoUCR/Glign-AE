#!/bin/bash

GRAPH_PATH=/home/cgao037/graph/lj/lj_snap.adj.weight
QUERY_FILE=/home/cgao037/Glign-AE/query_input/LJ_queries.txt
# 只保留尚未执行的组合（从你的清单中提取）
COMBINATIONS=(
  "64 16"
  "32 16"
  "16 16"
  "256 8"
  "128 8"
  "32 8"
  "16 8"
  "256 4"
  "128 4"
  "32 4"
  "16 4"
)

for combo in "${COMBINATIONS[@]}"; do
  read -r batch version <<< "$combo"
  echo "Running: batch=$batch, parallelVersion=$version"
  ./fusion_sssp \
    -option glign \
    -batch $batch \
    -max_combination 256 \
    -mode 2 \
    -delay \
    -qf $QUERY_FILE \
    -batchNum 16 \
    -batchRatio 0.01 \
    -parallelVersion $version \
    $GRAPH_PATH
  echo "Finished: batch=$batch, version=$version"
  echo "----------------------------------------"
done
