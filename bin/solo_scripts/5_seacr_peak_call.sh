#!/bin/bash
#5_seacr_peak_call.sh

data=$1
results=$2
control_bedgraph=$3

# -------------------------------------------------------------------------------------------------
# 4: Run the peak caller SEACR
# -------------------------------------------------------------------------------------------------
seacr_output=$results/peak_calls
[ ! -d $seacr_output ] && mkdir -p $seacr_output

#Bedgraph conversion loop.
while IFS= read -r sample_name; do
	bash $seacr $bedgraph_output/${sample_name}_fragments.bedgraph $control_bedgraph 0.01 norm stringent $seacr_output/${sample_name}
done < $data/all_samplelist.txt