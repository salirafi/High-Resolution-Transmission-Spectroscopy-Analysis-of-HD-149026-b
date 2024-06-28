import numpy as np
from scipy.interpolate import interp1d
import astropy.constants as const

# function for high-pass filter
def high_pass_filter(spectrum_input, binsize, maxval=False):

    binned, x = [], []
    for i in range(0,spectrum_input.size - binsize, binsize):
        
        x.append(i + binsize / 2.) # define the bin center
        if maxval: binned.append(np.max(spectrum_input[i:i+binsize])) # taking maximum value in each bin
        else: binned.append(np.median(spectrum_input[i:i+binsize])) # taking median value in each bin

    filtered = interp1d(np.asarray(x),np.asarray(binned),kind='linear',bounds_error=False,
                        fill_value='extrapolate')(np.arange(spectrum_input.size))

    return filtered

# function for extracting the +-1 sigma and the median value of a normal distribution
def med_n_lim(sample, dist):

    # calculating the cumulative distribution function ('dist' should be a PDF-like distribution)
    cdf = np.cumsum(dist)

    # normalizing CDF by its summed value or, equivalently, its last value
    cdf_norm = cdf / cdf[-1]

    # extracting the median-1sigma, median, median+1sigma values
    median, one_sigma = 0.5, 0.6827 # median and one-sigma values of a CDF of a Gaussian-like distribution
    low1sig, x_med, high1sig = interp1d(cdf_norm,sample,\
                                        fill_value='extrapolate')\
                                        (np.array([median-one_sigma*median,median,median+one_sigma*median]))
    
    return low1sig-x_med, x_med, high1sig-x_med

# function for removing masked pixels
def mask_pix(wave, flux, error):

    """ data must have dimensions of (time x wavelength) """

    maskpix = error[0] != np.inf
    fl = np.zeros((flux.shape[0], flux[0][maskpix].shape[0]))
    er = fl.copy()
    wv = fl.copy()
    for frame in range(flux.shape[0]):
        fl[frame] = flux[frame][maskpix]
        er[frame] = error[frame][maskpix]
        wv[frame] = wave[frame][maskpix]

    return wv, fl, er, maskpix

# function for constructing a Gaussian kernel centered in 0 km/s
def gaussian_kernel(FWHM, sampling_step, kernel_width=5):

    """
    Gaussian kernel evenly sampled in velocity where
    the input should also be evenly sampled in velocity or log wavenumber

    Args:
        FWHM            : Full-Width-at-Half-Maximum of the Gaussian kernel (km/s)
        sampling_step   : sampling step of the input (km/s)
        kernel_width    : kernel width (sigma)

    Return:
        velocity of the kernel, gaussian kernel
    """

    # calculate the standard deviation of the Gaussian kernel (km/s)
    sigma = FWHM / (np.sqrt(8. * np.log(2.)))

    # define the width of the kernel in Gaussian sigma
    kernel_width = kernel_width

    # kernel x sample
    x_sample = np.arange(-int(kernel_width * sigma / sampling_step) * sampling_step,
                         int(kernel_width * sigma / sampling_step) * sampling_step + sampling_step,
                         sampling_step) # this is to ensure that the kernel is centered at 0 km/s
    
    # Gaussian kernel
    gauss = np.exp(-0.5 * (x_sample)**2 / (sigma**2))

    # normalising the kernel
    gauss /= np.nansum(gauss) # normalise by the sum of the kernel because the sum of area inside the kernel has to be 1 (since the kernel is a PDF) --> preserve the line shape and depth

    return x_sample, gauss

# function for convolving data with a Gaussian kernel
def inst_gaussian(wvnumber, flux, kernel_FWHM):

    """
    Convolve a Gaussian kernel to an evenly sampled in velocity

    Args:
        wvnumber    : input wavenumber evenly-spaced in log/in constant velocity (cm-1)
        flux        : input flux in any unit
        kernel_FWHM : Full-Width-at-Half-Maximum of the Gaussian kernel (km/s)

    Return:
        convolved flux
    """

    # check if the input is evenly spaced in log
    if np.max(np.diff(np.diff(np.log(wvnumber)))) > 1e-5:
        raise ("Sample the input in a evenly log spaced wavenumber or constant velocity")

    c = const.c.value / 1000.  # velocity of light in vacuum (km/s)
    sampling_step = np.mean((1 / wvnumber[:-1] - 1 / wvnumber[1:]) / (1 / wvnumber[1:]) * c) # sampling step of the input in km/s
    
    vel_sample, norm_gauss = gaussian_kernel(kernel_FWHM, sampling_step)

    # handling of the edges
    nf = len(flux)
    flux = np.concatenate((np.ones(nf) * flux[0], flux, np.ones(nf) * flux[-1]))

    return np.convolve(flux, norm_gauss, mode='same')[nf:-nf]

# function for constructing rotational broadening kernel
def rotation_kernel(vsini, u1, u2, sampling_step):

    """
    Solid body rotational kernel evenly sampled in velocity based on equation (55) in Kawahara et al. 2022,
       Gray 2005, and Exojax's exojax.response.rotkernel. 
       The input should be evenly sampled in velocity or log wavenumber.
       This function assumes quadratic limb darkening model with coefficients u1 and u2.
    Args:
        vsini        : projected rotational velocity (km/s)
        u1           : limb-darkening coefficient 1
        u2           : limb-darkening coefficient 2      
        sampling_step: sampling step of the input (km/s)
    Return:
        velocity of the kernel, rotational kernel
    """ 

    # kernel width in km/s
    kernel_width = vsini * 2

    # kernel x sample
    x_sample = np.arange(-int(kernel_width / sampling_step) * sampling_step,
                       int(kernel_width / sampling_step) * sampling_step + sampling_step,
                       sampling_step) # this is to ensure that the kernel is centered at 0 km/s
    
    # rotation kernel
    eta = x_sample**2 / vsini**2
    kernel = np.where(eta <= 1.0, np.pi / 2.0 * u1 * (1.0 - eta) - 2.0 / 3.0 * np.sqrt(1.0 - eta) * (3.0*u1 + u2 + 2.0*u2*eta**2 - 3.0), 0.0) # accounting NaN values caused by the square root
    
    # normalising the kernel
    kernel /= np.nansum(kernel) # normalise by the sum of the kernel because the sum of area inside the kernel has to be 1 (since the kernel is a PDF) --> preserve the line shape and depth

    return x_sample, kernel

# def rotation_kernel_transm(R_s, R_p, P_rot, u1, u2, sampling_step):

#     """
#     Solid body rotational kernel evenly sampled in velocity based on equation (55) in Kawahara et al. 2022,
#        Gray 2005, and Exojax's exojax.response.rotkernel. 
#        The input should be evenly sampled in velocity or log wavenumber.
#        This function assumes quadratic limb darkening model with coefficients u1 and u2.
#     Args:
#         vsini        : projected rotational velocity (km/s)
#         u1           : limb-darkening coefficient 1
#         u2           : limb-darkening coefficient 2      
#         sampling_step: sampling step of the input (km/s)
#     Return:
#         velocity of the kernel, rotational kernel
#     """ 

#     # kernel width in km/s
#     dx = R_p/(50*R_s)
#     vsini = 2*np.pi*np.arange(-1,1,dx)*R_s/P_rot

#     kernel = np.exp(-0.5*()/(3.73/(2*np.sqrt(2*np.log(2))))) # accounting NaN values caused by the square root
    
#     # normalising the kernel
#     kernel /= np.nansum(kernel) # normalise by the sum of the kernel because the sum of area inside the kernel has to be 1 (since the kernel is a PDF) --> preserve the line shape and depth

#     return x_sample, kernel

def rotation_flux(wvnumber, flux, vsini, u1, u2):
    """Convolve a rigid body rotational kernel to a evenly sampled in velocity

    Args:
        wvnumber: input wavenumber evenly-spaced in log/in constant velocity (cm-1)
        flux    : input flux in any unit
        vsini   : projected rotational velocity (km/s)
        u1      : limb-darkening coefficient 1
        u2      : limb-darkening coefficient 2   

    Return:
        convolved flux
    """ 
    
    #check if the input is evenly spaced in log
    if np.max(np.diff(np.diff(np.log(wvnumber)))) >1e-5:
        raise ("Sample the input in a evenly log spaced wavenumber or constant velocity")
        
    c = const.c.value / 1000.  # velocity of light in vacuum (km/s)
    
    sampling_step = np.mean(
        (1 / wvnumber[:-1] - 1 / wvnumber[1:]) / (1 / wvnumber[1:]) * c) # sampling step of the input in km/s
    
    vel_sample, norm_rot = rotation_kernel(vsini,u1,u2,sampling_step)

    # handling of the edges
    nf = len(flux)
    flux = np.concatenate((np.ones(nf) * flux[0], flux, np.ones(nf) * flux[-1]))

    return np.convolve(flux, norm_rot, mode="same")[nf:-nf]

def solve_E(M, e, tol=1e-12, max_iter=1000):
    """
    Solves for E in the equation M = E - e*sin(E) using the Newton-Raphson method.
    """
    # Initial guess for E
    E = M
    
    # Iterate using the Newton-Raphson method until convergence or max_iter is reached
    for i in range(max_iter):
        f = E - e*np.sin(E) - M
        df = 1 - e*np.cos(E)
        dE = -f / df
        E += dE
        
        # Check for convergence
        if np.all(np.abs(dE) < tol):
            break
    
    return E

def find_nearest(array, value):
    return (np.abs(array - value)).argmin()

def split_matrix_for_dot(a,b,split):
    b_split = np.empty((int(b.shape[1]/split),b.shape[0],split))
    b_split[0] = b[:,1*split-split:1*split]
    a_final = np.dot(a,b_split[0])
    for i in range(1+1,int(b.shape[1]/split)+1,1):
        b_split[i-1] = b[:,i*split-split:i*split]
        a_dot = np.dot(a,b_split[i-1])
        a_final = np.concatenate((a_final,a_dot),axis=1)
    return a_final