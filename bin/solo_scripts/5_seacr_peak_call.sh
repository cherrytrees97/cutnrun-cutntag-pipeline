#!/bin/bash
#5_seacr_peak_call.sh

data=$1
results=$2
control_bedgraph=$3
seacr=$4

bedgraph_output=$results/bedgraph

# -------------------------------------------------------------------------------------------------
# 4: Run the peak caller SEACR
# -------------------------------------------------------------------------------------------------
seacr_output=$results/peak_calls
[ ! -d $seacr_output ] && mkdir -p $seacr_output

#Bedgraph conversion loop.
while IFS= read -r sample_name; do
	bash $seacr/SEACR_1.3.sh $bedgraph_output/${sample_name}_fragments.bedgraph $control_bedgraph norm stringent $seacr_output/${sample_name}
done < $data/all_samplelist.txt