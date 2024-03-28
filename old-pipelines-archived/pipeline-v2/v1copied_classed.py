from astroquery.mast import Observations
from astropy import units as u
from astropy.io import fits
from astropy.timeseries import LombScargle

import matplotlib.pyplot as plt
import numpy as np


class LC:
    # Initialize with TIC ID, sector, and data source
    def __init__(self, tid, sector, source):
        self.tid = tid
        self.sector = sector
        if source not in ["spoc-lc", "ffi-lc", "ffi-only"]:
            raise AttributeError(f"source must be one of ['spoc-lc', 'ffi-lc', 'ffi-only'], input was \"{source}\"")
        self.source = source
        # Set a lot of variables
        self.get_fits_path()
        self.parse_quality()
        self.make_pgrm()

    def get_fits_path(self):
        ''' Using the TIC ID, sector, and source, download the file
        '''
        if self.source == "spoc-lc":
            ob = Observations.query_criteria(target_name=self.tid
                                             ,sequence_number=self.sector
                                             ,provenance_name="SPOC")
        elif self.source == "ffi-lc":
            ob = Observations.query_criteria(target_name=self.tid
                                             ,sequence_number=self.sector
                                             ,provenance_name="TESS-SPOC")
        # TO-DO: how will we support FFIs?
        prod = Observations.get_product_list(ob)
        lc = Observations.download_products(prod, productSubGroupDescription="LC")
        fits_path = lc['Local Path'][0]
        self.fits = fits_path


    def parse_quality(self):
        ''' Read in a fits path, output the times / fluxes
        (remove the "bad quality" values from both)
        Also remove NaNs, or else we can't generate periodograms
        '''
        # Open file, read in the times, SAP flux, quality flags, and cadence
        with fits.open(self.fits, mode="readonly") as hdulist:
            tmz = hdulist[1].data['TIME']*u.day
            pdc_fluxes = hdulist[1].data['PDCSAP_FLUX']
            quality = hdulist[1].data['QUALITY']
            pdc_error = hdulist[1].data['PDCSAP_FLUX_ERR']
            timedel = hdulist[1].header['TIMEDEL']*u.day
            elapse = hdulist[1].header['TELAPSE']*u.day
        # Convert flags to binary
        bad_bits = np.array([1,2,3,4,5,6,8,10,12,13,15])
        value = 0
        for v in bad_bits:
            value = value + 2**(v-1)
        
        # "bad data" will match the quality flags
        bad_quality = np.bitwise_and(quality, value) >= 1
        
        # Include any NaNs in bad data, 
        nan_data = np.isnan(pdc_fluxes)
        
        # Concat bad quality and NaN indexes
        bad_data = np.logical_or(bad_quality, nan_data)

        # Set all of the class variables
        self.times = tmz[~bad_data]
        self.flux = pdc_fluxes[~bad_data]
        self.flux_err = pdc_error[~bad_data]
        self.timedel = timedel
        self.elapsed = elapse

    
    def make_pgrm(self):
        t = self.times
        # Normalize these values
        f_norm = self.flux/np.median(self.flux)
        f_err = self.flux_err/np.median(self.flux_err)
        
        # Nyquist limits
        max_f = (1/(2*self.timedel)).to(u.microhertz)
        min_f = (2/self.elapsed).to(u.microhertz)
        samples = int(5*self.elapsed*max_f)
        freqs = np.linspace(min_f, max_f, samples)
        
        # do psd lombscargle
        LS = LombScargle(t, f_norm, normalization='psd')
        pwrs = LS.power(freqs)
        max_p = max(pwrs)
        
        # get false alarm probability for maximum peak
        LS_err = LombScargle(t, f_norm, f_err, normalization='psd')
        max_pwr_err = LS_err.power(freqs).max()
        sigma3 = float(LS_err.false_alarm_level(0.997))*max_p/max_pwr_err
        
        # special correction to convert to amplitude
        pwrs = np.sqrt(pwrs) * np.sqrt(4.0 / len(t))
        sigma3 = np.sqrt(sigma3) * np.sqrt(4.0 / len(t))

        self.freqs=freqs
        self.amps=pwrs 
        self.sigma3=sigma3 
 
    
    def make_plots(self):
        f_norm = self.flux/np.nanmedian(self.flux)
        fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,5), gridspec_kw={'width_ratios':[3,1]})
        # plot lc
        ax1.scatter(self.times ,(f_norm-1)*100 ,c='black',s=2)
        ax1.set_xlabel('TESS day')
        ax1.set_ylabel('amplitude [%]')
        
        # plot prgm
        ax2.plot(self.freqs, self.amps*100, c='black')
        ax2.set_xlabel('Frequency [$\mu$Hz]')
        ax2.set_ylabel('amplitude [%]')
        ax2.hlines(self.sigma3, self.freqs[0].value, self.freqs[-1].value, linestyles='dashed', colors='tab:orange')
        
        # Set up period at top
        def invert(x):
            # 1/x with special treatment of x == 0
            x = np.array(x).astype(float)
            near_zero = np.isclose(x, 0)
            x[near_zero] = np.inf
            x[~near_zero] = 277.77778 / x[~near_zero]
            return np.around(x, decimals=1)
        per_ax2 = ax2.secondary_xaxis('top', functions=(invert, invert))
        per_ax2.set_xticks([1, 0.25])
        per_ax2.set_xlabel('Period [h]')
    
    def get_stats(self):
        # No variability detected
        if self.sigma3 > max(self.amps*100):
            # Maximum amplitude under one day
            m1d = max(self.amps[np.argwhere(lc.freqs>11.574*u.microhertz)])
            # Envelope containing 99.7% of the data
            v99 = np.nanpercentile(self.flux/np.nanmedian(self.flux), 99.865) - np.nanpercentile(self.flux/np.nanmedian(self.flux), 0.135)
            print(f"No variability detected\nm1d: {m1d}\nv99: {v99}")
        else:
            max_amp = max(self.amps)
            max_freq = self.freqs[np.argwhere(self.amps==max_amp)]
            max_per = (1/max_freq).to(u.day)
            print(f"Max freq: {max_freq} = {max_per}\n Amp @ max period: {max_amp*100}")