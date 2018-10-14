#!/bin/bash

# source sims_maf
source /global/common/software/lsst/cori-haswell-gcc/stack/setup_current_sims.sh
# set up my own version rn
setup sims_maf -r /global/homes/a/awan/LSST/lsstRepos/sims_maf
setup sims_maf_contrib -r /global/homes/a/awan/LSST/lsstRepos/sims_maf_contrib

##############################################################################
# run the script for all the cadences
#for yr_cut in 1 3 6 10
#do
#    python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 --yr_cut=$yr_cut \
#                                                            --outDir='/global/homes/a/awan/LSST/output/coadd_output_allwps_perNight/' &
#done

# run the script for the new cadences
#for yr_cut in 1 3 6 10
#do
#    for db in kraken_2042.db kraken_2044.db mothra_2049.db nexus_2097.db
#    do
#        python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 --yr_cut=$yr_cut \
#                                                            --specific_db=$db \
#                                                            --outDir='/global/homes/a/awan/LSST/output/coadd_output_allwps_perNight/' &
#    done
#done

##############################################################################
outdir=/global/cscratch1/sd/awan/lsst_output/coadd_output_allwps_slair_altsched/
# run the script for the slair and altsced cadences. different db_paths and outdirs.
mkdir -p ${outdir}sl_scripts    # create a directory for the slurm script if it doesnt exist already.
cd ${outdir}sl_scripts

for yr_cut in 1 3 6 10
do
    # ------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------
    # run analysis on slair outputs
    for db in cadence_roll_75_mix_rolling_mix_10yrs \
                roll_mix_100_rolling_mix_10yrs \
                roll_mix_rolling_mix_10yrs \
                rolling_10yrs  tms_roll_10yrs
    do
        # ------------------------------------------------------------------------------------
        # write a slurm script for this coadd run.
        cat > ${db}_${yr_cut}yrs.sl << EOF
#!/bin/bash -l
#SBATCH --nodes=1               # Use 1 node; has multiple cores
#SBATCH --qos=regular
#SBATCH -t 00:30:00             # Set 20 min time limit
#SBATCH --constraint=haswell    # Use Haswell nodes
#SBATCH --output=${outdir}/sbatch_output/calc_coadd_${db}_${yr_cut}yrs.%j.out
#SBATCH --job-name=calc_coadd
srun python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 --yr_cut=$yr_cut \
                                                --specific_db=${db}.db --no_dither --slair \
                                                --dbs_path='/global/cscratch1/sd/awan/dbs_wp_unzipped/slair_altsched' \
                                                --outDir=${outdir}
EOF
        # ------------------------------------------------------------------------------------
        # run the slurm script
        sbatch ${db}_${yr_cut}yrs.sl
        echo Job submitted for ${db}_${yr_cut}yrs
    done

    # ------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------
    # need to separate alt_sched outputs since slair set up is markedly different.
    for db in alt_sched  alt_sched_rolling
    do
        # ------------------------------------------------------------------------------------
        # write a slurm script for this coadd run.
        cat > ${db}_${yr_cut}yrs.sl << EOF
#!/bin/bash -l
#SBATCH --nodes=1               # Use 1 node; has multiple cores
#SBATCH --qos=regular
#SBATCH -t 00:30:00             # Set 20 min time limit
#SBATCH --constraint=haswell    # Use Haswell nodes
#SBATCH --output=${outdir}/sbatch_output/calc_coadd_${db}_${yr_cut}yrs.%j.out
#SBATCH --job-name=calc_coadd
srun python /global/homes/a/awan/LSST/lsstRepos/ObsStrat/code/modifyWFD/calc_coadd_depth.py --nside=256 --yr_cut=$yr_cut \
                                                --specific_db=${db}.db --no_dither \
                                                --dbs_path='/global/cscratch1/sd/awan/dbs_wp_unzipped/slair_altsched' \
                                                --outDir=${outdir}
EOF
        # run the slurm script
        sbatch ${db}_${yr_cut}yrs.sl
        echo Job submitted for ${db}_${yr_cut}yrs
    done
done