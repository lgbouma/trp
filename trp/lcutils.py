import numpy as np

def p2p_rms(flux):
    """
    e.g., section 2.6 of Nardiello+2020:
        The point-to-point RMS is not sensitive to intrinsic stellar
        variability.  It's obtained by calculating the 84th-50th percentile of
        the distribution of the sorted residuals from the median value of
        δF_j = F_{j} - F_{j+1}, where j is an epoch index.
    """

    dflux = np.diff(flux)

    med_dflux = np.nanmedian(dflux)

    up_p2p = (
        np.nanpercentile( np.sort(dflux-med_dflux), 84 )
        -
        np.nanpercentile( np.sort(dflux-med_dflux), 50 )
    )
    lo_p2p = (
        np.nanpercentile( np.sort(dflux-med_dflux), 50 )
        -
        np.nanpercentile( np.sort(dflux-med_dflux), 16 )
    )

    p2p = np.mean([up_p2p, lo_p2p])

    return p2p / np.sqrt(2)



