# dictionary mapping the overall categories of the new fbs outputs and the db filenames to look for for each group.
folder_map = {'mw_heavy': ['mw_heavy_v1.6_10yrs.db', 'mw_heavy_nexp2_v1.6_10yrs.db'],
 'dm_heavy': ['dm_heavy_v1.6_10yrs.db', 'dm_heavy_nexp2_v1.6_10yrs.db'],
 'ddf_heavy': ['ddf_heavy_v1.6_10yrs.db', 'ddf_heavy_nexp2_v1.6_10yrs.db'],
 'rolling_fpo': ['rolling_fpo_2nslice0.8_v1.6_10yrs.db',
  'rolling_fpo_2nslice0.9_v1.6_10yrs.db',
  'rolling_fpo_2nslice1.0_v1.6_10yrs.db',
  'rolling_fpo_3nslice0.8_v1.6_10yrs.db',
  'rolling_fpo_3nslice0.9_v1.6_10yrs.db',
  'rolling_fpo_3nslice1.0_v1.6_10yrs.db',
  'rolling_fpo_6nslice0.8_v1.6_10yrs.db',
  'rolling_fpo_6nslice0.9_v1.6_10yrs.db',
  'rolling_fpo_6nslice1.0_v1.6_10yrs.db'],
 'rolling_exgal': ['rolling_exgal_mod2_dust_sdf_0.80_v1.6_10yrs.db',
  'rolling_exgal_mod2_dust_sdf_0.80_nexp2_v1.6_10yrs.db'],
 'ss_heavy': ['ss_heavy_v1.6_10yrs.db', 'ss_heavy_nexp2_v1.6_10yrs.db'],
 'combo_dust': ['combo_dust_v1.6_10yrs.db', 'combo_dust_nexp2_v1.6_10yrs.db'],
 'baseline': ['baseline_nexp1_v1.6_10yrs.db',
  'baseline_nexp2_v1.6_10yrs.db',
  'baseline_nexp2_scaleddown_v1.6_10yrs.db'],
 'even_filters': ['even_filtersv1.6_10yrs.db',
  'even_filters_altv1.6_10yrs.db',
  'even_filters_alt_g_v1.6_10yrs.db',
  'even_filters_g_v1.6_10yrs.db'],
 'barebones': ['barebones_v1.6_10yrs.db', 'barebones_nexp2_v1.6_10yrs.db']}

for group in folder_map:
    for i, key in enumerate( folder_map[group] ):
        folder_map[group][i] = key.split('.db')[0]