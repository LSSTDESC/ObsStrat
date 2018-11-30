########################################################################################################################
# Consolidation of some helper functions in various notebooks.
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################################################################
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib.lines import Line2D

import matplotlib as mpl
fontsize = 16
mpl.rcParams['figure.figsize'] = (10, 6)
mpl.rcParams['axes.labelsize'] = fontsize
mpl.rcParams['xtick.labelsize'] = fontsize-2
mpl.rcParams['ytick.labelsize'] = fontsize-2
mpl.rcParams['legend.fontsize'] = fontsize-2
mpl.rcParams['axes.titlesize'] = fontsize
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['axes.grid'] = True
mpl.rcParams['figure.titlesize'] = fontsize

__all__ = ['get_area', 'pix_to_radec', 'plot_skymap', 'plot_skymap_somepix', 'plot_matplot']

########################################################################################################################
def get_area(pix_arr, nside):
    return len(pix_arr)*hp.nside2pixarea(nside=nside, degrees=True)

########################################################################################################################
def pix_to_radec(pixels, nside):
    lat, ra = hp.pix2ang(ipix=pixels, nside=nside)
    ra = np.remainder(ra+np.pi*2, np.pi*2)
    dec = np.pi/2.0 - lat
    return np.rad2deg(ra), np.rad2deg(dec)

########################################################################################################################
def plot_skymap(arr, title='', cmin=None, cmax=None):
    plt.clf()
    hp.mollview(arr, flip='astro', rot=(0,0,0), title=title, min=cmin, max=cmax)
    hp.graticule(dpar=20, dmer=20, verbose=False)
    plt.gcf().set_size_inches(6, 6)
    plt.show()

    ########################################################################################################################
def plot_skymap_somepix(pix_list, nside, title_append=None, pix_list_2=None):
    # set up the skymaps
    npix = hp.nside2npix(nside)

    val = np.zeros(npix)+1.
    footprint = val.view(np.ma.MaskedArray)
    
    footprint.mask = [True]*npix
    footprint.data[:] = 0
    footprint.fill_value = np.nan
    footprint.mask[pix_list] = False
    footprint.data[pix_list] = 500

    if pix_list_2 is not None:
        footprint.mask[pix_list_2] = False
        footprint.data[pix_list_2] += 250
    if title_append is not None:
        title = '%s\n%.2f deg2; nside %s'%(title_append,
                                           get_area(pix_list, nside),
                                           nside)
    else:
        title = '%.2f deg2; nside %s'%(get_area(pix_list, nside),
                                       nside)
    # plot the skymap
    plt.clf()
    hp.mollview(footprint, flip='astro', rot=(0,0,0), title=title, cbar=False)
    hp.graticule(dpar=20, dmer=20, verbose=False)
    plt.gcf().set_size_inches(6, 6)
    plt.show()

########################################################################################################################
def plot_matplot(nside, pixels_dict, colors, alphas, titles, syms,
                 saveplot=False, fname=None, add_ec_mw=False, xcoord=None, ycoord=None):
    # plots the skymap using matplotlib mollweide projection since it labels RA, Dec lines
    plt.figure()
    plt.subplot(111, projection="mollweide")

    # set up some plotting params
    size = 4.0
    
    if syms is None: syms = ['o']*len(colors)
    # ras, decs = {}, {}
    custom_lines, labels = [], []
    # loop over the pixels
    for i, pix_type in enumerate(pixels_dict):
        # convert healpix pixels to ra, dec values
        lat, ra = hp.pix2ang(ipix=pixels_dict[pix_type], nside=nside)
        ra = np.remainder(ra+np.pi*2, np.pi*2)
        dec = np.pi/2.0 - lat
        c = SkyCoord(ra=ra * u.radian, dec=dec * u.radian, frame='icrs')
        ra_rad = c.ra.wrap_at(180. * u.deg).radian
        dec_rad = c.dec.radian

        print('%s: %.2f <= dec <= %.2f'%(pix_type, min(np.degrees(dec_rad)),
                                         max(np.degrees(dec_rad))))
        # plot the ra, dec
        plt.scatter(-ra_rad, dec_rad, s=size, marker=syms[i], color=colors[i], alpha=alphas[i])
        
        custom_lines.append(Line2D([0], [0], color=colors[i], lw=10))
        labels.append('%s (%.f deg$^2$)'%(titles[i], get_area(pixels_dict[pix_type], nside)))

    if add_ec_mw:
        # code adapted from lsst sims_maf
        # https://github.com/lsst/sims_maf/blob/ff7bec6daa7d4291a1b87f63aca6930ba8bfedcc/python/lsst/sims/maf/plots/spatialPlotters.py#L429-L460
        # ------------------------------------------------------------------------
        # plot the eccliptic spur
        ra_center = 0
        peak_width = np.radians(10.)
        taper_length = np.radians(80.)
        ecinc = 23.439291 * (np.pi / 180.0)
        ra_ec = np.arange(0, np.pi * 2., (np.pi * 2. / 360.))
        dec_ec = np.sin(ra_ec) * ecinc
        lon = -(ra_ec - ra_center - np.pi) % (np.pi * 2) - np.pi
        # plot
        plt.plot(lon, dec_ec, 'r.', markersize=2, alpha=0.6)
        # add legend
        custom_lines.append(Line2D([0], [0], color='r', lw=2))
        labels.append('Ecliptic Spur')
        # ------------------------------------------------------------------------
        # Calculate galactic coordinates for mw location.
        step = 0.02
        gal_l = np.arange(-np.pi, np.pi + step / 2., step)
        val = peak_width * np.cos(gal_l / taper_length * np.pi / 2.)
        gal_b1 = np.where(np.abs(gal_l) <= taper_length, val, 0)
        gal_b2 = np.where(np.abs(gal_l) <= taper_length, -val, 0)
        # Convert to ra/dec.
        # Convert to lon/lat and plot.
        for gal_b in [gal_b1, gal_b2]:
            c = SkyCoord(l=gal_l*u.radian, b=gal_b*u.radian, frame='galactic')
            c = c.transform_to(frame='icrs')
            ra, dec = c.ra.radian, c.dec.radian
            lon = -(ra - ra_center - np.pi) % (np.pi * 2) - np.pi
            plt.plot(lon, dec, 'b.', markersize=2, alpha=0.6)
        # add legend
        custom_lines.append(Line2D([0], [0], color='b', lw=2))
        labels.append('Galactic Plane')
    # ---------------------------------------------
    # add legend/title
    if i==0:
        plt.title('%s (%.f deg$^2$)'%(titles[i],
                                      get_area(pixels_dict[pix_type], nside)),
                  y=1.05)
    else:
        if add_ec_mw:
            del_x, del_y = 0.1, 0.2
            labelspacing = 0.5
        else:
            del_x, del_y = 0., 0.
            labelspacing = 0.001
        if xcoord is None or ycoord is None:
            if i==1:
                xcoord, ycoord = 0.95+del_x, 1.25+del_y
            elif i==2:
                xcoord, ycoord = 0.95+del_x, 1.35+del_y
            elif i==3:
                xcoord, ycoord = 0.9+del_x, 1.4+del_y
            else:
                xcoord, ycoord = 1.+del_x, 1.25+del_y
        ax = plt.gca()
        ax.legend(custom_lines, labels, bbox_to_anchor=(xcoord, ycoord),
                  frameon=False, labelspacing=labelspacing)
    # ---------------------------------------------
    # relabel stuff: want RA increasing left-right
    locs, label = plt.xticks()
    labels = np.rad2deg(locs)
    print('\nOld labels: %s'%labels)
    labels[labels!=0] = -labels[labels!=0]
    print('News labels: %s'%labels)
    labels = [r'%.f$\degree$'%(f) for f in labels]
    plt.gca().set_xticklabels(labels)
    
    if saveplot:
        if fname is None:
            raise ValueError('Need fname to save plot.')
        else:
            plt.savefig(filename=fname, format='png', bbox_inches='tight')
            print('\n## Saved %s'%fname.split('/')[-1])
    
    plt.show()