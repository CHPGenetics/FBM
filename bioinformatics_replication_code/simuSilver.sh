#!/bin/bash
R=10
nItr=3000
nGene=1000
nMethy=2000
signal=weak
nCausalM=10
nCausalMbar=10
fdrcut=0.2

# mkdir  /home1/fangz/BI/Silver05_0
# mkdir  /home1/fangz/BI/Silver0_05
# mkdir  /home1/fangz/BI/Silver02_02

# cp /home1/fangz/BI/S0/dataSilver.rdata /home1/fangz/BI/Silver05_0/
# cp /home1/fangz/BI/S0/dataSilver.rdata /home1/fangz/BI/Silver0_05/
# cp /home1/fangz/BI/S0/dataSilver.rdata /home1/fangz/BI/Silver02_02/

# # # # Rscript --vanilla BIrunSilver.R  /home1/fangz/BI/Silver05_0 1 2334 $nItr 0.5 0 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal 0 $fdrcut &
# # # # Rscript --vanilla BIrunSilver.R  /home1/fangz/BI/Silver0_05 1 2334 $nItr 0 0.5 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal 0 $fdrcut &
# # # # Rscript --vanilla BIrunSilver.R  /home1/fangz/BI/Silver02_02 1 2334 $nItr 0.2 0.2 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal 0 $fdrcut &

# wait

# for r in $(seq 1 $R); do
#     Rscript --vanilla BIrunSilver.R  /home1/fangz/BI/Silver05_0 1 2334 $nItr 0.5 0 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $r $fdrcut &
# done
# for r in $(seq 1 $R); do
#     Rscript --vanilla BIrunSilver.R  /home1/fangz/BI/Silver0_05 1 2334 $nItr 0 0.5 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $r $fdrcut &
# done
# for r in $(seq 1 $R); do
#     Rscript --vanilla BIrunSilver.R  /home1/fangz/BI/Silver02_02 1 2334 $nItr 0.2 0.2 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $r $fdrcut &
# done

# wait

# Rscript --vanilla BIsummarySilver.R  /home1/fangz/BI/Silver05_0 1 2334 $nItr 0.5 0 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $fdrcut &
# Rscript --vanilla BIsummarySilver.R  /home1/fangz/BI/Silver0_05 1 2334 $nItr 0 05 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $fdrcut &
# Rscript --vanilla BIsummarySilver.R  /home1/fangz/BI/Silver02_02 1 2334 $nItr 0.2 0.2 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $fdrcut &

# wait
 

# mkdir ~/BI/SilverAll
# for theta in "05_0" "0_05" "02_02"
# do
#     cp /home1/fangz/BI/Silver$theta/TP.pdf /home1/fangz/BI/SilverAll/TP_$theta.pdf
#     cp /home1/fangz/BI/Silver$theta/MSEyHat.pdf /home1/fangz/BI/SilverAll/MSEyHat_$theta.pdf
# done




# mkdir  /home1/fangz/BI/S05_0
# mkdir  /home1/fangz/BI/S0_05
# mkdir  /home1/fangz/BI/S02_02

# cp /home1/fangz/BI/S0/dataSilver.rdata /home1/fangz/BI/S05_0/
# cp /home1/fangz/BI/S0/dataSilver.rdata /home1/fangz/BI/S0_05/
# cp /home1/fangz/BI/S0/dataSilver.rdata /home1/fangz/BI/S02_02/

# wait

# for r in $(seq 1 $R); do
#     Rscript --vanilla BIrunS.R  /home1/fangz/BI/S05_0 1 2334 $nItr 0.5 0 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $r $fdrcut &
# done
# for r in $(seq 1 $R); do
#     Rscript --vanilla BIrunS.R  /home1/fangz/BI/S0_05 1 2334 $nItr 0 0.5 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $r $fdrcut &
# done
# for r in $(seq 1 $R); do
#     Rscript --vanilla BIrunS.R  /home1/fangz/BI/S02_02 1 2334 $nItr 0.2 0.2 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $r $fdrcut &
# done

# wait

Rscript --vanilla BIsummarySilver.R  /home1/fangz/BI/S05_0 1 2334 $nItr 0.5 0 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $fdrcut &
Rscript --vanilla BIsummarySilver.R  /home1/fangz/BI/S0_05 1 2334 $nItr 0 05 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $fdrcut &
Rscript --vanilla BIsummarySilver.R  /home1/fangz/BI/S02_02 1 2334 $nItr 0.2 0.2 full $nGene $nCausalM $nCausalMbar $nMethy $R $signal $fdrcut &

wait
 

# mkdir ~/BI/SAll
# for theta in "05_0" "0_05" "02_02"
# do
#     cp /home1/fangz/BI/S$theta/TP.pdf /home1/fangz/BI/SAll/TP_$theta.pdf
# done


 Rscript --vanilla BIRMSEplot.R
Rscript --vanilla BIPPplot.R

for theta in "05_0" "0_05" "02_02"
do
    for r in {1..10}
    do
	cp /home1/fangz/BI/S$theta/n1_r$r/est /home1/fangz/BI/SAll/$theta\_r$r
	cp /home1/fangz/BI/S$theta/n1_r$r\_CC/estCC /home1/fangz/BI/SAll/$theta\_r$r\_CC
    done
done

theta="02_02"
for r in {1..10}
do
    cp /home1/fangz/BI/S$theta/n1_r$r/est /home1/fangz/BI/SAll/$theta\_r$r
    cp /home1/fangz/BI/S$theta/n1_r$r\_CC/estCC /home1/fangz/BI/SAll/$theta\_r$r\_CC
    cp /home1/fangz/BI/S$theta/n1_r$r\_Gene/estGene /home1/fangz/BI/SAll/$theta\_r$r\_Gene
done
