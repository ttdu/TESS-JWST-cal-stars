from astroquery.mast import Observations
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from astropy.timeseries import LombScargle, TimeSeries, aggregate_downsample
from numba_fft.main import false_alarms
from vo_magic import get_name_type

import matplotlib.pyplot as plt
import numpy as np

# Suppress astroquery download messages
from contextlib import redirect_stdout
import io

# Suppress sqrt warning
import warnings

class LC:
    """Lightcurve object; may hold a single sector or multiple sectors worth
    of data.
    """
    def __init__(self, obs, ticid):
        self.obs = obs
        self.ticid = ticid
        
        # Create some numpy arrays to hold these values
        self.times = np.array([])
        self.fluxnorm = np.array([])

        # Create some empty lists to hold these values
        self.sectors = []
        self.timedel = []
        self.elapsed = []
        
        # download the data, parse it, make periodograms, calculate statistics
        self.get_fits_path()
        if any(self.fits):
            self.parse_quality()
            self.make_pgrm()
            self.get_stats()

            # Get the "real" name from simbad
            self.name, self.otype, self.stype = get_name_type(self.ticid)
            self.make_plots(savefig=True)


    def get_fits_path(self):
        ''' Use the ob to download the file, get the fits path
        '''
        # get list of product IDs to download
        prod = Observations.get_product_list(self.obs)

        # logic in get_spoc_obs prevents us from getting duplicates here
        prod = Observations.filter_products(prod, productSubGroupDescription=["LC","FAST-LC"])
        if len(prod) == 0:
            # Make this a list so we can check this later
            self.fits = [None]
        else:
            with redirect_stdout(io.StringIO()) as f:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    lc  = Observations.download_products(prod)
            self.fits = lc['Local Path'].value


    def parse_quality(self):
        ''' Read in a fits path, output the times / fluxes
        (remove the "bad quality" values from both)
        Also remove NaNs, or else we can't generate periodograms
        '''
        for path in self.fits:
            # Open file, read in the times, SAP flux, quality flags, and cadence
            with fits.open(path, mode="readonly") as hdulist:
                tmz = hdulist[1].data['TIME']
                pdc_fluxes = hdulist[1].data['PDCSAP_FLUX']
                quality = hdulist[1].data['QUALITY']
                self.crowd = hdulist[1].header['CROWDSAP']
                self.tic = hdulist[0].header['OBJECT']
                # There are a few special values that we need to preserve
                self.sectors.append(hdulist[0].header['SECTOR'])
                self.timedel.append(hdulist[1].header['TIMEDEL']*u.day)
                self.elapsed.append(hdulist[1].header['TELAPSE']*u.day)
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

            # Get time, flux, normalized flux
            self.times = np.append(self.times, tmz[~bad_data])
            
            # keep only the "good" values, convert to normalized flux
            flux = pdc_fluxes[~bad_data]
            self.fluxnorm = np.append(self.fluxnorm, flux/np.nanmedian(flux))
        
        # for the nyquist frequency, we need the highest timedel
        self.timedel = max(self.timedel)
        # total elapsed time is, happily, just a sum
        self.elapsed = sum(self.elapsed)

    def make_pgrm(self, get_alarms=True):
        # Nyquist limits
        max_f = (1/(2*self.timedel)).to(u.microhertz)
        min_f = (2/self.elapsed).to(u.microhertz)
        samples = int(5*self.elapsed*max_f)
        self.freqs = np.linspace(min_f, max_f, samples)

        # do PSD lombscargle, get powers
        LS = LombScargle(self.times*u.day, self.fluxnorm, normalization="psd")
        self.pwrs = LS.power(self.freqs)
        # ignore weird values in the square root
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.amps = 100*np.sqrt(4*self.pwrs/len(self.times))

        # Get the false alarm levels
        if get_alarms:
            self.alarms = false_alarms(datax = self.times, datay = self.fluxnorm)
            self.alarms = 100*self.alarms

    def make_plots(self, savefig=False, bin=True, n_bins=50):
        fig, ax = plt.subplot_mosaic(mosaic=[["A","A"],["C", "B"]], figsize=(12,7), width_ratios=[1,3], dpi=300)
        plt.subplots_adjust(hspace=0.3)
        try:
            fig.suptitle(f"{self.tic}, sector {self.sectors}, crowd factor: {self.crowd:.3f}\n{self.name} | {self.otype} | {self.stype}")
        except:
            print(f"{self.tic} has no target name?")
        # plot lc
        #if self.isvar:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ts = TimeSeries({"time":Time(self.times, format='mjd'), "flux":(self.fluxnorm-1)*100}).fold(period=self.max_per)

        # Plot A: full, flux-normalized lightcurve, plotted in % amplitudes
        ax['A'].scatter(self.times, (self.fluxnorm-1)*100, c='black', s=2, alpha=0.2)
        ax['A'].set_xlabel("TESS day (BTJD)")
        ax['A'].set_ylabel("amplitude [%]")

        # Plot B: folded light curve, folded at period with maximum power
        #
        # report sensible units for the period
        if self.max_per < 1*u.day:
            if self.max_per < 1*u.hour:
                per = self.max_per.to(u.min)
            else:
                per = self.max_per.to(u.hour)
            self.max_per = per
        
        # bin down to 50 points by default
        if bin:
            bs = aggregate_downsample(ts, n_bins=n_bins)
            time, flux = bs['time_bin_start'].value, bs['flux'].value
            sz = None
            alp = 1
        else:
            time, flux = ts['time'].value, ts['flux'].value
            sz = 2
            alp = 0.2
        # Plot B: actually make the plot
        ax['B'].scatter(time, flux, c='black', s=sz, alpha=alp)
        ax['B'].set_xlabel('Days')
        ax['B'].set_ylabel('amplitude [%]')
        ax['B'].text(0.05,0.1, f"Period: {self.max_per:.2f}",bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 5}, transform=ax['B'].transAxes)

        # Plot C: periodogram, plotted in % amplitudes
        ax['C'].plot(self.freqs, self.amps, c='black')
        ax['C'].set_xlabel('Frequency [$\mu$Hz]')
        ax['C'].set_ylabel('amplitude [%]')
        for alarm in self.alarms:
            ax['C'].axhline(alarm, 0.05, 0.95, c='g')

        if self.isvar:
            # Only plot a "focused" area around the maximum peak
            df = 300 * u.microhertz
            xmin = (self.max_freq-df)/u.microhertz
            xmax = (self.max_freq+df)/u.microhertz
            # Lowest/highest limits shouldn't go beyond bounds of data
            if xmin < 0:
                xmin = 0
            if xmax > max(self.freqs)/u.microhertz:
                xmax = max(self.freqs)/u.microhertz
            ax['C'].set_xlim(xmin, xmax)

        plt.plot()

        # Set up period at top -- commented out because the inversion
        # produces a broad range of values (e.g. 0.1, 15, 240) that are not super useful
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
            plt.savefig(f"./diagnostic-plots/{self.tic}-{self.sectors}.png")
            plt.close()

    def get_stats(self):
        self.max_amp = max(self.amps)
        self.max_freq = self.freqs[np.argwhere(self.amps==self.max_amp)][0][0]
        self.max_per = (1/self.max_freq).to(u.day)
        
        # Only mark as variable if we pass the 3 sigma mark
        if np.any(self.amps>self.alarms[2]):
            self.isvar = True
            
        # No variability detected
        else:
            self.isvar = False
