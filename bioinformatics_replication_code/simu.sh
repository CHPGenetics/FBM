#!/bin/bash
R=35
nItr=200
nGene=50
nMethy=80
signal=weak
nCausalGene=5

mkdir  /home1/fangz/BI/BI_05_0
for r in $(seq 1 $R); do
    Rscript --vanilla BIrunStable.R  /home1/fangz/BI/BI_05_0 1 15217 $nItr 0.5 0 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait
Rscript --vanilla BIsummaryStable.R  /home1/fangz/BI/BI_05_0 1 15217 $nItr 0.5 0 full $nGene $nCausalGene $nMethy $R $signal

mkdir  /home1/fangz/BI/BI_0_05
mkdir  /home1/fangz/BI/BI_02_02
mkdir  /home1/fangz/BI/BI_01_01
mkdir  /home1/fangz/BI/BI_0_01
mkdir  /home1/fangz/BI/BI_025_0
mkdir  /home1/fangz/BI/BI_0_025
mkdir  /home1/fangz/BI/BI_025_025
mkdir  /home1/fangz/BI/BI_03_03
mkdir  /home1/fangz/BI/BI_01_0

cp -r /home1/fangz/BI/BI_05_0/*Full  /home1/fangz/BI/BI_0_05/
cp -r /home1/fangz/BI/BI_05_0/*Full  /home1/fangz/BI/BI_02_02/
cp -r /home1/fangz/BI/BI_05_0/*Full  /home1/fangz/BI/BI_01_01/
cp -r /home1/fangz/BI/BI_05_0/*Full  /home1/fangz/BI/BI_0_01/
cp -r /home1/fangz/BI/BI_05_0/*Full  /home1/fangz/BI/BI_025_0/
cp -r /home1/fangz/BI/BI_05_0/*Full  /home1/fangz/BI/BI_0_025/
cp -r /home1/fangz/BI/BI_05_0/*Full  /home1/fangz/BI/BI_025_025/
cp -r /home1/fangz/BI/BI_05_0/*Full  /home1/fangz/BI/BI_03_03/
cp -r /home1/fangz/BI/BI_05_0/*Full  /home1/fangz/BI/BI_01_0/

for r in $(seq 1 $R); do
    Rscript --vanilla BIrunStable.R  /home1/fangz/BI/BI_0_05 1 15217 $nItr 0 0.5 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait
Rscript --vanilla BIsummaryStable.R  /home1/fangz/BI/BI_0_05 1 15217 $nItr 0 0.5 full $nGene $nCausalGene $nMethy $R $signal

for r in $(seq 1 $R); do
    Rscript --vanilla BIrunStable.R  /home1/fangz/BI/BI_02_02 1 15217 $nItr 0.2 0.2 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait
Rscript --vanilla BIsummaryStable.R  /home1/fangz/BI/BI_02_02 1 15217 $nItr 0.2 0.2 full $nGene $nCausalGene $nMethy $R $signal

for r in $(seq 1 $R); do
    Rscript --vanilla BIrunStable.R  /home1/fangz/BI/BI_01_01 1 15217 $nItr 0.1 0.1 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait
Rscript --vanilla BIsummaryStable.R  /home1/fangz/BI/BI_01_01 1 15217 $nItr 0.1 0.1 full $nGene $nCausalGene $nMethy $R $signal

for r in $(seq 1 $R); do
    Rscript --vanilla BIrunStable.R  /home1/fangz/BI/BI_0_01 1 15217 $nItr 0 0.1 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait
Rscript --vanilla BIsummaryStable.R  /home1/fangz/BI/BI_0_01 1 15217 $nItr 0 0.1 full $nGene $nCausalGene $nMethy $R $signal

for r in $(seq 1 $R); do
    Rscript --vanilla BIrunStable.R  /home1/fangz/BI/BI_025_0 1 15217 $nItr 0.25 0 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait
Rscript --vanilla BIsummaryStable.R  /home1/fangz/BI/BI_025_0 1 15217 $nItr 0.25 0 full $nGene $nCausalGene $nMethy $R $signal

for r in $(seq 1 $R); do
    Rscript --vanilla BIrunStable.R  /home1/fangz/BI/BI_0_025 1 15217 $nItr 0 0.25 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait
Rscript --vanilla BIsummaryStable.R  /home1/fangz/BI/BI_0_025 1 15217 $nItr 0 0.25 full $nGene $nCausalGene $nMethy $R $signal


for r in $(seq 1 $R); do
    Rscript --vanilla BIrunStable.R  /home1/fangz/BI/BI_025_025 1 15217 $nItr 0.25 0.25 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait
Rscript --vanilla BIsummaryStable.R  /home1/fangz/BI/BI_025_025 1 15217 $nItr 0.25 0.25 full $nGene $nCausalGene $nMethy $R $signal

for r in $(seq 1 $R); do
    Rscript --vanilla BIrunStable.R  /home1/fangz/BI/BI_03_03 1 15217 $nItr 0.3 0.3 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait
Rscript --vanilla BIsummaryStable.R  /home1/fangz/BI/BI_03_03 1 15217 $nItr 0.3 0.3 full $nGene $nCausalGene $nMethy $R $signal

for r in $(seq 1 $R); do
    Rscript --vanilla BIrunStable.R  /home1/fangz/BI/BI_01_0 1 15217 $nItr 0.1 0 full $nGene $nCausalGene $nMethy $R $signal $r &
done
wait
Rscript --vanilla BIsummaryStable.R  /home1/fangz/BI/BI_01_0 1 15217 $nItr 0.1 0 full $nGene $nCausalGene $nMethy $R $signal


############ wait for all
mkdir /home1/fangz/BI/BI
for theta in "05_0" "0_05" "025_0" "0_025" "01_0" "0_01" "01_01" "02_02" "025_025" "03_03"
do
    cp /home1/fangz/BI/BI_$theta/CI.pdf /home1/fangz/BI/BI/CI_$theta.pdf
    cp /home1/fangz/BI/BI_$theta/MSEyHat.pdf /home1/fangz/BI/BI/MSEyHat_$theta.pdf
    cp /home1/fangz/BI/BI_$theta/FS.pdf /home1/fangz/BI/BI/FS_$theta.pdf
   cp /home1/fangz/BI/BI_$theta/ROC.pdf /home1/fangz/BI/BI/ROC_$theta.pdf
done
