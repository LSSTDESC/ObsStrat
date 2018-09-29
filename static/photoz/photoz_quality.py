### Let's make a new metric that computes the coadded depth at a point on the sky, in every filter.
### Includes option to correct for MW dust extinction, which has been commented out because we don't want that.
import lsst.sims.maf.metrics as metrics
import numpy as np

class Photo_z_quality(metrics.BaseMetric):
    """Calculate extragalacitc coadded depth in each filter.
    """
    def __init__(self, m5Col='fiveSigmaDepth', units='XXXsome kind of unit', maps=['DustMap'],
                 wavelen_min=None , wavelen_max=None , wavelen_step=1., filterCol='filter',
                 implement_depth_ebv_cut=False, i_lim_mag=26.0, **kwargs):
        self.filternames = ['u', 'g', 'r', 'i', 'z', 'y']
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.implement_depth_ebv_cut = implement_depth_ebv_cut
        self.i_lim_mag = i_lim_mag
        ### Without the corrention for MW extinction
        self.coaddSimple = metrics.Coaddm5Metric()
        ### With a correction for MW extinction
#         self.coadd_with_dust = {}
#         for filtername in self.filternames:
#            self.coadd_with_dust[filtername] = metrics.ExgalM5(lsstFilter=filtername, 
#                                                          m5Col=m5Col, units=units, 
#                                                          maps=maps, wavelen_min=wavelen_min,
#                                                          wavelen_max=wavelen_max, 
#                                                          wavelen_step=1.,**kwargs)
        # set up for i-band extincted coadded depth
        if self.implement_depth_ebv_cut:
            self.coadd_iband_with_dust = metrics.ExgalM5(lsstFilter='i',
                                                         m5Col=self.m5Col, units=units,
                                                         maps=maps, **kwargs)

        super(Photo_z_quality, self).__init__(col=[self.m5Col, self.filterCol],
                                              maps=maps, units=units, **kwargs)
        self.metricDtype = 'object'
        
    def run(self, dataslice, slicePoint=None):
        coadd_depths = []
        nfilters = 0
        for filtername in self.filternames:
            in_filt = np.where(dataslice[self.filterCol] == filtername)[0]
            if len(in_filt) > 0:
                ### Without the correction for MW extinction
                coadd_depths.append(self.coaddSimple.run(dataslice[in_filt]))               
                ### With the correction for MW extinction
                ### NOTE that this APPENDS and so you may print out BOTH,
                ###  but then this will change the format of the data files
#                 coadd_depths.append(self.coadd_with_dust[filtername].run(dataslice[in_filt],
#                                                                        slicePoint=slicePoint))
                nfilters += 1
            else:
                coadd_depths.append(self.badval)

        if self.implement_depth_ebv_cut:
            # find the i-band extincted coadded depth
            ext_iband_coadd = self.coadd_iband_with_dust.run(dataSlice=dataslice,
                                                             slicePoint=slicePoint)

            # get ebv
            ebv = slicePoint['ebv']

        # figure out the conditions
        discard_condition = (nfilters !=6) # want coverage in all 6 filters
        if self.implement_depth_ebv_cut:
            discard_condition = discard_condition or (ebv>0.2) or (ext_iband_coadd<self.i_lim_mag)

        if discard_condition: # want coverage in all 6 filters
            # return a single badval which will mask the datapoint in the bundle.metricValues.
            coadd_depths = self.badval
            
        return coadd_depths