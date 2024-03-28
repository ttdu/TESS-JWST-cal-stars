from astroquery.mast import Observations
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import LombScargle, TimeSeries

import matplotlib.pyplot as plt
import numpy as np

# Suppress astroquery download messages
from contextlib import redirect_stdout
import io

# Suppress sqrt warning
import warnings

plt.rcParams['figure.dpi'] = 300



class LC:
    # Initialize with fitspath
    def __init__(self, fitspath):
        self.fits = fitspath
        # Sets a lot of variables
        if self.fits:
            self.parse_quality()
            self.make_pgrm()
            self.get_stats()       

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
            self.timedel = hdulist[1].header['TIMEDEL']*u.day
            self.elapsed = hdulist[1].header['TELAPSE']*u.day
            self.crowd = hdulist[1].header['CROWDSAP']
            self.tic = hdulist[0].header['OBJECT']
            self.sector = hdulist[0].header['SECTOR']
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
        self.fluxnorm = self.flux/np.nanmedian(self.flux)
        self.flux_err = pdc_error[~bad_data]

    
    def make_pgrm(self):
        t = self.times
        # Normalize these values
        f_err = self.flux_err/np.median(self.flux_err)
        
        # Nyquist limits
        max_f = (1/(2*self.timedel)).to(u.microhertz)
        min_f = (2/self.elapsed).to(u.microhertz)
        samples = int(5*self.elapsed*max_f)
        freqs = np.linspace(min_f, max_f, samples)
        
        # do psd lombscargle
        LS = LombScargle(t, self.fluxnorm, normalization='psd')
        pwrs = LS.power(freqs)
        max_p = max(pwrs)
        
        # get false alarm probability for maximum peak
        LS_err = LombScargle(t, self.fluxnorm, f_err, normalization='psd')
        max_pwr_err = LS_err.power(freqs).max()
        sigma3 = float(LS_err.false_alarm_level(0.997))*max_p/max_pwr_err
        
        # special correction to convert to amplitude
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pwrs = np.sqrt(pwrs) * np.sqrt(4.0 / len(t))
        sigma3 = np.sqrt(sigma3) * np.sqrt(4.0 / len(t))

        self.freqs=freqs
        self.amps=pwrs 
        self.sigma3=sigma3 
 
    
    def make_plots(self):
        fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,5), gridspec_kw={'width_ratios':[3,1]})
        fig.suptitle(f"{self.tic} : sector {self.sector}")
        # plot lc
        if self.isvar:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ts = TimeSeries({"time":Time(self.times, format='mjd'), "flux":(self.fluxnorm-1)*100}).fold(period=self.max_per)
            ax1.scatter(ts['time'].value, ts['flux'].value, c='black', s=2, alpha=0.2)
            ax1.set_xlabel('TESS day')
            ax1.set_ylabel('amplitude [%]')
        else:
            ax1.scatter(self.times ,(self.fluxnorm-1)*100 ,c='black',s=2)
            ax1.set_xlabel('TESS day')
            ax1.set_ylabel('amplitude [%]')
        
        # plot prgm
        ax2.plot(self.freqs, self.amps*100, c='black')
        ax2.set_xlabel('Frequency [$\mu$Hz]')
        ax2.set_ylabel('amplitude [%]')
        ax2.hlines(self.sigma3, self.freqs[0].value, self.freqs[-1].value, linestyles='dashed', colors='tab:orange')
        if self.isvar:
            df = 150 * u.microhertz
            ax2.set_xlim((self.max_freq-df)/u.microhertz, (self.max_freq+df)/u.microhertz)
        
        # Set up period at top
        # def invert(x):
        #     # 1/x with special treatment of x == 0
        #     x = np.array(x).astype(float)
        #     near_zero = np.isclose(x, 0)
        #     x[near_zero] = np.inf
        #     x[~near_zero] = 277.77778 / x[~near_zero]
        #     return np.around(x, decimals=1)
        # per_ax2 = ax2.secondary_xaxis('top', functions=(invert, invert))
        # per_ax2.set_xticks([1, 0.25])
        # per_ax2.set_xlabel('Period [h]')
        if self.isvar:
            plt.savefig(f'./var/{self.tic}-{self.sector}.png')
        else:
            plt.savefig(f'./nonvar/{self.tic}-{self.sector}.png')
        plt.close()
    
    def get_stats(self):
        # No variability detected
        if self.sigma3 > max(self.amps*100):
            self.isvar = False
            # Maximum amplitude under one day
            self.m1d = max(self.amps[np.argwhere(self.freqs>11.574*u.microhertz)])*100
            # Envelope containing 99.7% of the data
            self.v99 = (np.nanpercentile(self.fluxnorm, 99.865) - np.nanpercentile(self.fluxnorm, 0.135))*100
        # Variability detected
        else:
            self.isvar = True
            self.v95 = (np.nanpercentile(self.fluxnorm, 97.5) - np.nanpercentile(self.fluxnorm, 2.5))*100
            self.max_amp = max(self.amps)*100
            self.max_freq = self.freqs[np.argwhere(self.amps==max(self.amps))][0][0]
            self.max_per = (1/self.max_freq).to(u.day)
        