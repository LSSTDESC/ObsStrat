The goal here to construct a modified WFD footprint that is optimized for extragalactic science (by avoiding high-extinction areas) and for LSS (wider area coverage). We start with using the wider-coverage OpSim output, `pontus_2002`.

#### Code Details
`bash_calc_coadd.sh` is the script that loads the right sims stacks and runs `calc_coadd_depth.py` that calculates the 5$\sigma$ coadded depth for 1yr and 10yr surveys (using the code from [mafContrib](https://github.com/LSST-nonproject/sims_maf_contrib ); the latest code still needs to be merged with the master branch but the latest versions are in [my fork](https://github.com/humnaawan/sims_maf_contrib )).

`WFD_footprint_forOpsimTeam.ipynb` runs `implement_ebv_cut_only.py` to implement an EBV>0.2 cut on the database. Please see the notebook for more details.