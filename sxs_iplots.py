import re
import sxs
import numpy as np
import bilby
import scipy.interpolate
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import numpy as np
import ipywidgets
import mpl_interactions.ipyplot as iplt

G = 6.67430e-11
c = 3e8
M = 6.563e+31
r = 3.086e24

def load_strain(h):
    simulations = sxs.load("simulations")
    sxs_bbh_2139 = sxs.load("SXS:BBH:"+h)
    metadata = sxs_bbh_2139.metadata
    h = sxs_bbh_2139.h

    reference_index = h.index_closest_to(metadata.reference_time)
    h=h[reference_index:]
    print(f"Mass ratio: {metadata.reference_mass_ratio}")
    return metadata, h

def initialize(h, G, c, M, r):
    h.time = h.time * G * (M/(c**3))
    h = h * (M/r) * (G/(c**2))
    t = h.t
    h22 = h.data[:, h.index(2,2)]
    A22 = np.abs(h22)
    phi22 = np.unwrap(np.angle(h22))
    phi22_t = scipy.interpolate.CubicSpline(t, phi22)
    f22 = phi22_t.derivative()(t) / (2*np.pi)
    df22dt = phi22_t.derivative(2)(t) / (2*np.pi)
    dt = 1e-4 #np.min(np.diff(h.t))
    i1 = h.index_closest_to(0.5)
    i2 = h.max_norm_index()
    return h, t, f22, df22dt, dt, i1, i2
    
def convertlm(m, f22, df22dt):
    flm = f22 * (m / 2)
    dflmdt = df22dt * (m / 2)  
    return flm, dflmdt

def SPA_fft_calc(l,m, h, t, f22, df22dt, dt, i1, i2, cut):
    """
    Calculate the SPA and FFT of user selected mode of strain h
    """

    #calculate for SPA
    hlm = h.data[:, h.index(l,m)]
    Alm = np.abs(hlm)
    philm = np.unwrap(np.angle(hlm))
    philm_t = scipy.interpolate.CubicSpline(t, philm)
    flm, dflmdt = convertlm(m, f22, df22dt)

    #calculate for FFT
    h_lm = h[:, h.index(l,m)]
    h_lm_interpolated = h_lm.interpolate(np.arange(h_lm.t[0], h_lm.t[-1], dt))
    hlm_tapered = h_lm_interpolated.taper(0, h.t[0]+1000*(G*(M/(c**3))))
    hlm_transitioned = hlm_tapered.transition_to_constant(h.max_norm_time()+100*(G*(M/(c**3))), h.max_norm_time()+200*(G*(M/(c**3))))
    hlm_padded = hlm_transitioned.pad(6000.0*(G*(M/(c**3))))
    hlm_line_subtracted = hlm_padded.line_subtraction()

    #print(hlm_line_subtracted.ndarray)
    htilde_lm = np.fft.rfft(hlm_line_subtracted)*dt
    frequencies_lm = np.fft.rfftfreq(len(hlm_line_subtracted), dt)
    htilde_lm_scaled = 2*(np.abs(htilde_lm))*(np.sqrt(frequencies_lm))

    #amplitude and frequency scaling
    f= -flm[i1:i2]
    amp = 2*np.abs(((1/2)*np.abs(hlm[i1:i2])) / np.sqrt(np.abs(dflmdt[i1:i2])))*np.sqrt(np.abs(flm[i1:i2]))

    if cut == True:
        #cutting off lower freq for frequencies_lm, htilde_lm_scaled
        frequencies_lm = frequencies_lm[60:]
        htilde_lm_scaled = htilde_lm_scaled[60:]
        relmax, relmin = find_extrema_indices(htilde_lm_scaled, 2.0e-25)
        frequencies_lm[:relmax[0]] = np.nan
        htilde_lm_scaled[:relmax[0]] = np.nan

    return f, amp, frequencies_lm, htilde_lm_scaled

def find_extrema_indices(data, threshold):
    """Finds indices of relative extrema greater than a threshold."""

    # Find indices of local maxima
    max_indices = argrelextrema(data, np.greater)[0]

    # Filter for maxima greater than the threshold
    filtered_max_indices = max_indices[data[max_indices] > threshold]

    # Find indices of local minima
    min_indices = argrelextrema(data, np.less)[0]

    # Filter for minima greater than the threshold
    filtered_min_indices = min_indices[data[min_indices] > threshold]

    return filtered_max_indices, filtered_min_indices

def create_functions(h, t, f22, df22dt, dt, i1, i2, cut=True):

    hlm = []
    for ell in range(2, h.ell_max + 1):
        for m in range(-ell, ell + 1):
            if ell < 7:
                if m > 0 and m >= ell - 2:
                    hlm.append([ell, m])
            else:
                if m == ell:
                    hlm.append([ell, m])
    
    f=[]; amp=[]; frequencies_lm=[]; htilde_lm_scaled=[];
    for i in hlm:
        f_i, amp_i, frequencies_lm_i, htilde_lm_scaled_i = SPA_fft_calc(i[0], i[1], h, t, f22, df22dt, dt, i1, i2, cut);
        f.append(np.array(f_i))
        amp.append(np.array(amp_i))
        frequencies_lm.append(np.array(frequencies_lm_i))
        htilde_lm_scaled.append(np.array(htilde_lm_scaled_i))
    f=np.array(f); amp=np.array(amp); frequencies_lm=np.array(frequencies_lm, dtype=object); htilde_lm_scaled=np.array(htilde_lm_scaled, dtype=object)

    t1 = h.max_norm_time()
    i1 = h.index_closest_to(h.max_norm_time())
    print(i1)
    
    def f_x_SPA(Mass, Distance):
        mscale = (Mass*1.989e+30)/M
        dscale = (Distance*3.086e+22)/r
        return np.array((f/mscale)).T[::12]

    def f_y_SPA(x, Mass, Distance):
        mscale = (Mass*1.989e+30)/M
        dscale = (Distance*3.086e+22)/r
        return np.array((mscale**(3/2) * amp / dscale)).T[::12]

    def f_x_strain(Mass, Distance):
        mscale = (Mass*1.989e+30)/M
        dscale = (Distance*3.086e+22)/r
        return np.array((frequencies_lm/mscale), dtype=object).T[::12]
        
    def f_y_strain(x, Mass, Distance):
        mscale = (Mass*1.989e+30)/M
        dscale = (Distance*3.086e+22)/r
        return np.array((mscale**(3/2) * htilde_lm_scaled / dscale), dtype=object).T[::12]
        
    return f_x_SPA, f_y_SPA, f_x_strain, f_y_strain

def iplt_lm(h, t, ratio, f22, df22dt, dt, i1, i2, cut=True):

    f_x_SPA_test, f_y_SPA_test, f_x_strain_test, f_y_strain_test = create_functions(h, t, f22, df22dt, dt, i1, i2, cut)
    
    Mass = np.linspace(33, 165)
    Distance = np.linspace(100, 500)

    fig, ax = plt.subplots()
    plt.xscale("log")
    plt.yscale("log")
    controls = iplt.plot(
        f_x_SPA_test, 
        f_y_SPA_test,
        Mass=Mass,
        Distance=Distance,
        xlim="fixed",
        ylim="fixed",
        color="lightblue"
    )
    with controls:
        iplt.plot(f_x_strain_test, f_y_strain_test, label="strain", xlim="fixed", ylim="fixed", ls="--", color="pink")
    plt.xlim(3, 5e3)
    plt.ylim(1e-25, 5e-21)
    plt.title("Mass Ratio: "+str(ratio))