from astroquery.mast import Observations
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import LombScargle, TimeSeries
from numba_fft.main import false_alarms
from find_spoc import get_spoc_obs
from vo_magic import get_name_type

import matplotlib.pyplot as plt
import numpy as np

# Suppress astroquery download messages
from contextlib import redirect_stdout
import io

# Suppress sqrt warning
import warnings

class LC:
    # Initialize with fitspath
    def __init__(self, ticid=None, sector=None, ob=None):
        self.ob = ob
        if not self.ob:
            self.ob = get_spoc_obs(ticid, sector)
        # Set a lot of variables
        self.get_fits_path()
        # Get the "real" name from simbad
        self.name, self.otype, self.stype = get_name_type(ticid)
        if self.fits:
            self.parse_quality()
            self.make_pgrm()
            self.get_stats()

    def get_fits_path(self):
        ''' Use the ob to download the file, get the fits path
        '''
        prod = Observations.get_product_list(self.ob)
        prod = Observations.filter_products(prod, productSubGroupDescription=["LC","FAST-LC"])
        if len(prod) == 0:
            self.fits = None
        else:
            with redirect_stdout(io.StringIO()) as f:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    lc  = Observations.download_products(prod)
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


    def make_pgrm(self, sg=False, lse=False):
        # Nyquist limits
        max_f = (1/(2*self.timedel)).to(u.microhertz)
        min_f = (2/self.elapsed).to(u.microhertz)
        samples = int(5*self.elapsed*max_f)
        self.freqs = np.linspace(min_f, max_f, samples)

        # do PSD lombscargle, get powers
        LS = LombScargle(self.times, self.fluxnorm, normalization="psd")
        self.pwrs = LS.power(self.freqs)

        # Get the false alarm levels
        #sigmas = np.array([.682689, .954499, .997300])
        self.alarmz = false_alarms(np.array(self.times).astype(np.float64), np.array(self.fluxnorm).astype(np.float64))
        self.alarmz = 100*self.alarmz

        # special correction to convert from power to amplitude
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            factor = np.sqrt(4.0 / len(self.times))
            self.amps = np.sqrt(self.pwrs) * factor
            self.alarms = np.sqrt(self.alarmz) * factor

    def make_plots(self, savefig=False):
        fig, ax = plt.subplot_mosaic(mosaic=[["A","A"],["C", "B"]], figsize=(12,7), width_ratios=[1,3], dpi=300)
        plt.subplots_adjust(hspace=0.3)
        try:
            fig.suptitle(f"{self.tic}, sector {self.sector}, crowd factor: {self.crowd:.3f}\n{self.name} | {self.otype} | {self.stype}")
        except:
            print(f"{self.ob} has no target name?")
        # plot lc
        #if self.isvar:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ts = TimeSeries({"time":Time(self.times, format='mjd'), "flux":(self.fluxnorm-1)*100}).fold(period=self.max_per)

        # Full light curve
        ax['A'].scatter(self.times, (self.fluxnorm-1)*100, c='black', s=2, alpha=0.2)
        ax['A'].set_xlabel("TESS day (BTJD)")
        ax['A'].set_ylabel("amplitude [%]")

        # Folded light curve
        ax['B'].scatter(ts['time'].value, ts['flux'].value, c='black', s=2, alpha=0.2)
        ax['B'].set_xlabel('Days')
        ax['B'].set_ylabel('amplitude [%]')
        ax['B'].text(0.05,0.1, f"Period: {self.max_per:.2f}",bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5}, transform=ax['B'].transAxes)

        # Periodogram, plotted in % amplitudes
        ax['C'].plot(self.freqs, self.amps*100, c='black')
        ax['C'].set_xlabel('Frequency [$\mu$Hz]')
        ax['C'].set_ylabel('amplitude [%]')
        for alarm in self.alarmz:
            ax['C'].axhline(alarm, 0.05, 0.95, c='g')

        #ax['C'].legend()
        plt.plot()
        # if self.isvar:
        #     # Only plot a "focused" area around the maximum peak
        #     df = 150 * u.microhertz
        #     xmin = (self.max_freq-df)/u.microhertz
        #     xmax = (self.max_freq+df)/u.microhertz
        #     # Lowest/highest limits shouldn't go beyond bounds of data
        #     if xmin < 0:
        #         xmin = 0
        #     if xmax > max(self.freqs)/u.microhertz:
        #         xmax = max(self.freqs)/u.microhertz
        #     ax['C'].set_xlim(xmin, xmax)

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
        if savefig:
            if self.isvar:
                plt.savefig(f'./var/{self.tic}-{self.sector}.png')
            else:
                plt.savefig(f'./nonvar/{self.tic}-{self.sector}.png')
            plt.close()

    def get_stats(self):
        self.max_freq = self.freqs[np.argwhere(self.amps==max(self.amps))][0][0]
        self.max_per = (1/self.max_freq).to(u.day)
        self.max_amp = max(self.amps)*100
        if np.any(self.amps>self.alarms[2]):
            # Variability detected
            self.isvar = True
            self.v95 = (np.nanpercentile(self.fluxnorm, 97.5) - np.nanpercentile(self.fluxnorm, 2.5))*100

        else:
            # No variability detected
            self.isvar = False
            # Envelope containing 99.7% of the data
            self.v99 = (np.nanpercentile(self.fluxnorm, 99.865) - np.nanpercentile(self.fluxnorm, 0.135))*100

    def plot_pgrm(self, xmin=0, xmax=None):
        fig, ax = plt.subplots()
        ax.plot(self.freqs, self.amps*100, c='black')
        ax.set_xlabel('Frequency [$\mu$Hz]')
        ax.set_ylabel('amplitude [%]')
        if not xmax:
            xmax = max(self.freqs)/u.microhertz
        ax.set_xlim(xmin, xmax)
