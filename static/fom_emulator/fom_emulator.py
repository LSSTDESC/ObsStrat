"""
The goal of this script is to enable FoM emulation for WL+CL+LSS given a set of FoMs on a grid of
(area, depth).

Still to do:
- Get files from Tim.
- Check file read-in for consistency with area, depth.
- Inspect plots.
- Set up actual emulation.
- Look at outputs.
- Do the depth optimization.
"""
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
cm = plt.cm.get_cmap('RdYlBu')

# Set some defaults.
year_vals = ['Y1', 'Y3', 'Y6', 'Y10']

areas = np.zeros((3,4)) # 3 grid values, 4 years
areas[:,0] = np.array([7.5, 13., 16.]) * 1000.0 # deg^2 for Y1
areas[:,1] = np.array([10., 15., 20.]) * 1000.0 # for Y3, Y6, Y10
areas[:,2] = np.array([10., 15., 20.]) * 1000.0 # for Y3, Y6, Y10
areas[:,3] = np.array([10., 15., 20.]) * 1000.0 # for Y3, Y6, Y10

depths = np.zeros((3,4)) # 3 grid values, 4 years
depths[:,0] = np.array([24.9, 25.2, 25.5])
depths[:,1] = np.array([25.5, 25.8, 26.1])
depths[:,2] = np.array([25.9, 26.1, 26.3])
depths[:,3] = np.array([26.3, 26.5, 26.7])

area_mid = areas[1,1] # arbitrary area for rescaling

def load_foms(dir='FoM', prior=False):
    """Script to load some FoM values on a 3x3 grid in (area, depth)."""
    ### TO CHECK: ordering of FoM values in the file - is it area then depth or ...?

    # Get list of files in directory.
    file_list = glob.glob(os.path.join(dir, '*'))

    # Initialize FoM arrays.  The first dimension is area, second is depth, 3rd is year.
    fom_arr = np.zeros((3,3,4))

    # Find the ones with Y1, Y3, Y6, Y10.  For each, get the right types of FoMs in order.
    for year_str in year_vals:
        # Find the one that has that year string and no other, otherwise Y1 and Y10 are confused.
        my_file = None
        for file_name in file_list:
            # Find the one that has that year string.
            if year_str in file_name:
                # Then make sure none of the others is there.
                other_found = False
                for tmp_year_str in year_vals:
                    if tmp_year_str != year_str and not tmp_year_str in year_str:
                        if tmp_year_str in file_name:
                            other_found = True
                if not other_found:
                    if my_file is None:
                        my_file = file_name
                    else:
                        if my_file != file_name:
                            raise ValueError(
                                "Year %s found more than once in dir %s!"%\
                                (year_str, dir))
        if my_file is None:
            raise ValueError("Year %s not found in dir %s!"%(year_str, dir))
        print 'Found file %s for year %s in dir %s'%(my_file, year_str, dir)

        # Now read in the relevant parts of the file.
        lines = [line.rstrip('\n') for line in open(my_file)]

        # Find the right lines based on whether they include/exclude the Stage III prior.
        if prior:
            grep_str = 'incl'
        else:
            grep_str = 'excl'


        n_found = 0
        for line in lines:
            if grep_str in line:
                tmp_fom = float(line.split('=')[1])
                fom_arr[n_found % 3, n_found / 3, year_vals.index(year_str)] = tmp_fom
                n_found += 1
        print n_found,' relevant lines found'

    return fom_arr

def load_strategy_table(year_str = 'Y1'):
    infile = './strategy_table_%s.txt'%year_str
    if not os.path.exists(infile):
        raise ValueError("Cannot find file %s for year %s!"%(infile, year_str))

    # Read it in
    lines = [line.rstrip('\n') for line in open(infile)]

    strat_name = []
    area = []
    median_depth = []
    for line in lines:
        tmp_vals = line.split("|")
        strat_name.append(tmp_vals[1].strip())
        area.append(float(tmp_vals[3].strip()))
        median_depth.append(float(tmp_vals[4].strip()))
    print 'Strategy table loaded: %d lines in %s for year %s!'%\
        (len(strat_name), infile, year_str)
    return strat_name, area, median_depth

def emulate_fom(area_vals, depth_vals, grid_area_vals, grid_depth_vals, grid_fom_vals):
    """Fake this up for now."""
    fom_vals = area_vals/np.max(area_vals) + depth_vals/np.max(depth_vals)
    return fom_vals

# Get FoM values.  
foms_prior = load_foms(prior=True)
print foms_prior[:,:,0]
print foms_prior[:,:,3]
foms_noprior = load_foms(prior=False)
print foms_noprior[:,:,0]
print foms_noprior[:,:,3]
# We will pretend these are (area, depth, year). Need to confirm.

# Make some basic plots:
# In each year, make a 2D color plot of the FoM vs. depth and area without and with prior, without
# and with removal of area scaling.
for year_ind in range(len(year_vals)):
    plot_areas = np.repeat(areas[:,year_ind],3).reshape(areas[:,year_ind].shape[0],3)
    plot_depths = np.repeat(depths[:,year_ind],3).reshape(depths[:,year_ind].shape[0],3).transpose()
    fom_max = max(np.max(foms_prior[:,:,year_ind]), np.max(foms_noprior[:,:,year_ind]))
    fom_max_rescale = max(np.max(foms_prior[:,:,year_ind]*area_mid/plot_areas),
                             np.max(foms_noprior[:,:,year_ind]*area_mid/plot_areas))
    max_val = max(fom_max, fom_max_rescale)

    # No prior, no area rescaling.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sc = ax.scatter(plot_areas, plot_depths, c=foms_prior[:,:,year_ind]/max_val,
                    cmap=cm, s=80, edgecolors='none')
    plt.colorbar(sc)
    plt.title('FoM/%f (with Stage III prior)'%max_val)
    plt.xlabel('Area [sq. deg.]')
    plt.ylabel('Median i-band depth')
    plt.savefig('figs/fom_emulator_%s_prior.pdf'%year_vals[year_ind])

    # Prior, no area rescaling.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sc = ax.scatter(plot_areas, plot_depths, c=foms_noprior[:,:,year_ind]/max_val,
                    cmap=cm, s=80, edgecolors='none')
    plt.colorbar(sc)
    plt.title('FoM/%f (no prior)'%max_val)
    plt.xlabel('Area [sq. deg.]')
    plt.ylabel('Median i-band depth')
    plt.savefig('figs/fom_emulator_%s_noprior.pdf'%year_vals[year_ind])

    # No prior, area rescaling.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sc = ax.scatter(plot_areas, plot_depths,
                    c=foms_prior[:,:,year_ind]*area_mid/plot_areas/max_val,
                    cmap=cm, s=80, edgecolors='none')
    plt.colorbar(sc)
    plt.title('(Rescaled FoM with Stage III prior)/%f'%max_val)
    plt.xlabel('Area [sq. deg.]')
    plt.ylabel('Median i-band depth')
    plt.savefig('figs/fom_emulator_%s_prior_rescaled.pdf'%year_vals[year_ind])

    # Prior, area rescaling.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sc = ax.scatter(plot_areas, plot_depths,
                    c=foms_noprior[:,:,year_ind]*area_mid/plot_areas/max_val,
                    cmap=cm, s=80, edgecolors='none')
    plt.colorbar(sc)
    plt.title('(Rescaled FoM without prior)/%f'%max_val)
    plt.xlabel('Area [sq. deg.]')
    plt.ylabel('Median i-band depth')
    plt.savefig('figs/fom_emulator_%s_noprior_rescaled.pdf'%year_vals[year_ind])

    # Evaluate various strategies using emulator as needed.

    # Load strategy tables:
    tmp_strat, tmp_area, tmp_depth = load_strategy_table(year_vals[year_ind])
    emulated_fom_prior = emulate_fom(tmp_area, tmp_depth, plot_areas, plot_depths,
                                     foms_prior[:,:,year_ind])
    emulated_fom_noprior = emulate_fom(tmp_area, tmp_depth, plot_areas, plot_depths,
                                       foms_noprior[:,:,year_ind])
    inds = emulated_fom_noprior.argsort()[::-1]
    print ''
    print 'Emulated from best to worst in year %s'%year_vals[year_ind]
    print 'Strategy, FoM without prior, FoM with prior'
    for ind in inds:
        print tmp_strat[ind], emulated_fom_noprior[ind], emulated_fom_prior[ind]
    print ''

# Depth optimization.
