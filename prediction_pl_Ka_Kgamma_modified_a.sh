#!/bin/csh
# options
#$ -N prediction
#$ -m be
#$ -M mtrachse@nd.edu
#$ -pe smp 23
#$ -q *@@crc_gpu
#$ -l gpu=1
setenv OMP_NUM_THREADS 23

./pred_stepps_veg_no_time sample num_samples=2000 num_warmup=500 data file=prediction_15_taxa_6796_cells_260_knots_cal_pl_Ka_Kgamma_EPs_modified_a_88_sites_only_abies.dump output file=prediction_pl_Ka_Kgamma_modified_a_88_sites_nb.csv


