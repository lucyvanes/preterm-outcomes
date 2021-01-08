#!/bin/bash


cd /data/project/BIPP/Analyses/NNMF15_AAL/rank15_AAL_times100000
maps=$(ls NNMF_network* | sed 's/.nii.gz$//')


for m in ${maps}
do
echo ${m}


echo "applying weights"
fslmaths ../allSubsJac2mm.nii.gz -mul ${m}.nii.gz allSubsJac2mm_${m}.nii.gz

fslstats -t allSubsJac2mm_${m}.nii.gz -M > mean_weighted_${m}.txt


done 