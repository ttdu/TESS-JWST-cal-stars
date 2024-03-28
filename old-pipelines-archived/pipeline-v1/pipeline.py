from astropy import units as u
from astropy.io import fits
from astropy.timeseries import TimeSeries, LombScargle, aggregate_downsample
from astroquery.mast import Observations
from astropy.time import Time

import matplotlib.pyplot as plt
import numpy as np

def main(tid, sector):
    # Get the path to the fits file
    fits_path = get_fits_path(tid, sector)
    
    # Open the fits file, parse quality flags, and return time, flux, and error
    t, flux, err = parse_quality(fits_path)
    
    # Bin the light curve
    tb, fluxb, errb = bin_lc(t, flux, err)
    
    # Run a periodogram, return the frequencies, amplitude, and 3-sigma noise threshold
    freqs, amps, sigma3 = make_pgrm(tb, fluxb, errb)
    
    make_plots(tb, fluxb, freqs, amps, sigma3)


def get_fits_path(tid, sector):
    ''' Given a TIC ID and sector,
    download the LC file and return the local path
    '''
    obs = Observations.query_criteria(target_name=tid
                                      , provenance_name='SPOC'
                                      , sequence_number=sector)
    prod = Observations.get_product_list(obs)
    lc = Observations.download_products(prod[0])

    fits_path = lc['Local Path'][0]
    return fits_path

def parse_quality(fits_path):
    ''' Read in a fits path, output the times / fluxes
    (remove the "bad quality" values from both)
    Also remove NaNs, or else we can't generate periodograms
    '''
    # Open file, read in the times, SAP flux, and quality flags
    with fits.open(fits_path, mode="readonly") as hdulist:
        tmz = hdulist[1].data['TIME']
        pdc_fluxes = hdulist[1].data['PDCSAP_FLUX']
        quality = hdulist[1].data['QUALITY']
        pdc_error = hdulist[1].data['PDCSAP_FLUX_ERR']
    
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
    
    return tmz[~bad_data], pdc_fluxes[~bad_data], pdc_error[~bad_data]

def bin_lc(t, f, f_err, bs=6*u.min):
    t = Time(t, format='jd')
    ts = TimeSeries({'time':t, 'flux':f, 'flux_err':f_err})
    ts_bin = aggregate_downsample(ts, time_bin_size=bs)
    t_bin = ts_bin['time_bin_start']
    f_bin = ts_bin['flux']
    ferr_bin = ts_bin['flux_err']
    return t_bin[~f_bin.data.mask].value, f_bin.compressed(), ferr_bin.compressed()
    

def make_pgrm(t, f, f_err):
    f_norm = f/np.median(f)
    f_err = f_err/np.median(f)
    
    # make input freqs
    t = t*u.day
    freqs = np.linspace(0.5, 60, 50000) * u.microhertz
    
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
    return freqs, pwrs, sigma3
    

def make_plots(t, f, freqs, pwrs, sigma3):
    f_norm = f/np.nanmedian(f)
    fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,5), gridspec_kw={'width_ratios':[3,1]})
    # plot lc
    ax1.scatter(t ,(f_norm-1)*100 ,c='black',s=2)
    ax1.set_xlabel('TESS day')
    ax1.set_ylabel('amplitude [%]')
    
    # plot prgm
    ax2.plot(freqs, pwrs*100, c='black')
    ax2.set_xlabel('Frequency [$\mu$Hz]')
    ax2.set_ylabel('amplitude [%]')
    ax2.hlines(sigma3*100, freqs[0].value, freqs[-1].value, linestyles='dashed', colors='tab:orange')
    
    # Set up period at top
    def invert(x):
        # 1/x with special treatment of x == 0
        x = np.array(x).astype(float)
        near_zero = np.isclose(x, 0)
        x[near_zero] = np.inf
        x[~near_zero] = 277.77778 / x[~near_zero]
        return np.around(x, decimals=1)
    per_ax2 = ax2.secondary_xaxis('top', functions=(invert, invert))
    per_ax2.set_xticks(invert(freqs)[2000::11999])
    per_ax2.set_xlabel('Period [h]')
