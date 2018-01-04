#!/bin/bash
R=35
nItr=4000
nGene=1000
nMethy=2000
signal=weak
nCausalGene=10

mkdir  /home1/fangz/BI/BICV_05_0
mkdir  /home1/fangz/BI/BICV_0_05
mkdir  /home1/fangz/BI/BICV_025_025

for r in $(seq 1 $R); do
    Rscript --vanilla BIrunCV.R  /home1/fangz/BI/BICV_05_0 1 15217 $nItr 0.5 0 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait

for r in $(seq 1 $R); do
    Rscript --vanilla BIrunCV.R  /home1/fangz/BI/BICV_0_05 1 15217 $nItr 0 0.5 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait

for r in $(seq 1 $R); do
    Rscript --vanilla BIrunCV.R  /home1/fangz/BI/BICV_025_025 1 15217 $nItr 0.25 0.25 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait

Rscript --vanilla BIsummaryCV.R  /home1/fangz/BI/BICV_0_05 1 15217 $nItr 0 0.5 full $nGene $nCausalGene $nMethy $R $signal &
Rscript --vanilla BIsummaryCV.R  /home1/fangz/BI/BICV_05_0 1 15217 $nItr 0.5 0 full $nGene $nCausalGene $nMethy $R $signal &
Rscript --vanilla BIsummaryCV.R  /home1/fangz/BI/BICV_025_025 1 15217 $nItr 0.25 0.25 full $nGene $nCausalGene $nMethy $R $signal &
wait

mkdir /home1/fangz/BI/BICV
Rscript --vanilla BI3by3CV.R  /home1/fangz/BI/ 1 15217 $nItr 0.5 0 full $nGene $nCausalGene $nMethy $R $signal

for theta in "05_0" "0_05" "025_025"
do
    cp /home1/fangz/BI/BICV_$theta/decision1 /home1/fangz/BI/BICV/decision1_$theta
    cp /home1/fangz/BI/BICV_$theta/decision2 /home1/fangz/BI/BICV/decision2_$theta
    cp /home1/fangz/BI/BICV_$theta/decision3 /home1/fangz/BI/BICV/decision3_$theta
    cp /home1/fangz/BI/BICV_$theta/decision4 /home1/fangz/BI/BICV/decision4_$theta

    cp /home1/fangz/BI/BICV_$theta/decisionAdjust1 /home1/fangz/BI/BICV/decisionAdjust1_$theta
    cp /home1/fangz/BI/BICV_$theta/decisionAdjust2 /home1/fangz/BI/BICV/decisionAdjust2_$theta
    cp /home1/fangz/BI/BICV_$theta/decisionAdjust3 /home1/fangz/BI/BICV/decisionAdjust3_$theta
    cp /home1/fangz/BI/BICV_$theta/decisionAdjust4 /home1/fangz/BI/BICV/decisionAdjust4_$theta

    cp /home1/fangz/BI/BICV_$theta/err /home1/fangz/BI/BICV/err_$theta
    cp /home1/fangz/BI/BICV_$theta/adjusterr /home1/fangz/BI/BICV/adjusterr_$theta

    cp /home1/fangz/BI/BICV_$theta/mse.RData /home1/fangz/BI/BICV/mse_$theta.rData
done
