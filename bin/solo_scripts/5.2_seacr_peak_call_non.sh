#!/bin/bash
#5.2_seacr_peak_call_non.sh

data=$1
results=$2
seacr=$3

bedgraph_output=$results/bedgraph

# -------------------------------------------------------------------------------------------------
# 4: Run the peak caller SEACR
# -------------------------------------------------------------------------------------------------
seacr_output=$results/peak_calls_noctrl
[ ! -d $seacr_output ] && mkdir -p $seacr_output

#Bedgraph conversion loop.
while IFS= read -r sample_name; do
	bash $seacr/SEACR_1.3.sh $bedgraph_output/${sample_name}_fragments.bedgraph 0.01 non stringent $seacr_output/${sample_name}
done < $data/all_samplelist.txt