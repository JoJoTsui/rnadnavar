#!/usr/bin/env bash



SEQ2NEO=/t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seq2neo
CFG=$SEQ2NEO/config/runner.yaml



# Run set2 including extra samples (first pair used automatically)
python3 $SEQ2NEO/scripts/run_batch_from_json.py --config $CFG --set 2