#!/usr/bin/env bash

# mkdirs
mkdir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/ 
mkdir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/seq2neo
mkdir /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801

# nf-core configs
rsync -avP /t9k/mnt/hdd/work/Vax/pipeline/configs/ /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/configs/
# rnadnavar related databases
rsync -avP ~/joey/bio_db/ ~/WorkSpace/data/ngs/xuzhenyu/bio_db/
# deepvariant related binaries & models
rsync -avP /t9k/mnt/hdd/work/Vax/deepvariant/DeepVariant-1.9.0 /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/
rsync -avP /t9k/mnt/hdd/work/Vax/deepvariant/models /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/dv/
# rnadnavar test dataset
rsync -avP /t9k/mnt/hdd/work/Vax/sequencing/aim_exp/rdv_test/C008801/input/ /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/work/rnadnavar_test/C008801/input/
# rnadnavar pipeline
rsync -avP /t9k/mnt/hdd/work/Vax/pipeline/rnadnavar/ /t9k/mnt/WorkSpace/data/ngs/xuzhenyu/pipeline/rnadnavar/