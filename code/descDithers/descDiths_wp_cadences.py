from save_csv_dithers import save_csv_dithers

dbs_path= '/global/cscratch1/sd/awan/dbs_wp_unzipped'
outDir = '/global/homes/a/awan/desc/wp_descDithers_csvs'

save_csv_dithers(dbs_path, outDir, rot_rand_seed=42, trans_rand_seed=42,
                 print_progress=False, show_diagnostic_plots=False, save_plots=True)
