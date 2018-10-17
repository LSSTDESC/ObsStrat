The goal here is to implement depth and extinctions cuts for each OpSim cadence to define a footprint that is usable for extragalactic science. 

#### Code Details
`bash_calc_coadd_all.sh` calculates the 5$\sigma$ coadded depth for each cadence at different times into the full survey, i.e. 1yr, 3yr, 6yr,  and finally 10yr. `bash_implement_cuts_all_cadences.sh` calls `implement_depth_ebv_cuts.py` which reads in the coadded depth data, implements various cuts, and saves the outputs.

Summary markdown tables are saved the `mdfiles` folder while the rest of the outputs are in `/global/homes/a/awan/desc/depth_data_outputs/` on NERSC and should be readable by anyone with `lsst` group affiliation. Also, `summarize_outputs.ipynb` makes some plots that summarize the trends in the some of statistics on the final footprints (saved in `/global/homes/a/awan/desc/depth_data_outputs/comparisons/`).