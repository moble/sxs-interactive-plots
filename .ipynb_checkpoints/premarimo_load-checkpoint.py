import sxs
import math
import numpy as np
import pandas as pd
import scipy.interpolate
import csv

G = 6.67430e-11
c = 3e8
M = 6.563e+31
r = 3.086e24
dt = 1e-4 #np.min(np.diff(h.t))


df = sxs.load("dataframe", tag="3.0.0")

def sort(binary_id):
    for i in range(len(df[:-20])):
        if df.iloc[i].name[-8:] == binary_id:
            return i

def load_strain(h_id):
    simulations = sxs.load("simulations")
    sxs_bbh_n = sxs.load(h_id)
    metadata = sxs_bbh_n.metadata
    h = sxs_bbh_n.h
    reference_index = h.index_closest_to(metadata.reference_time)
    h=h[reference_index:]
    return metadata, h

def dimensionalize(h):
    h.time = h.time * G * (M/(c**3))
    #print(h)
    h = h * (M/r) * (G/(c**2))
    t = h.t
    return h, t

def SPA_fft_calc(l,m, h, t, metadata):
    """
    Calculate the SPA and FFT of user selected mode of strain h
    """

    #calculate for FFT
    h_lm = h[:, h.index(l,m)]
    h_lm_interpolated = h_lm.interpolate(np.arange(h_lm.t[0], h_lm.t[-1], dt))
    hlm_tapered = h_lm_interpolated.taper(0, h.t[0]+1000*(G*(M/(c**3))))
    #hlm_transitioned = hlm_tapered.transition_to_constant(h.t[h.max_norm_index()+100])#, h.max_norm_time()+200*(G*(M/(c**3))))
    hlm_transitioned = hlm_tapered.transition_to_constant(h.max_norm_time()+100*(G*(M/(c**3))))#, h.max_norm_time()+200*(G*(M/(c**3))))
    if type(metadata.reference_eccentricity) == float and ((metadata.reference_eccentricity) > 0.3):
        hlm_padded = hlm_transitioned.pad(100000*(G*(M/(c**3))))
    else:
        hlm_padded = hlm_transitioned.pad(25000*(G*(M/(c**3))))
    hlm_line_subtracted = hlm_padded.line_subtraction().real
    #print(type(hlm_line_subtracted.ndarray))
    htilde_lm = np.abs(np.fft.rfft(hlm_line_subtracted.ndarray.astype(float))*dt)
    frequencies_lm = np.abs(np.fft.rfftfreq(len(hlm_line_subtracted.ndarray.astype(float)), dt))
    htilde_lm_scaled = 2*(np.abs(htilde_lm))*np.abs(np.sqrt(frequencies_lm))

    return frequencies_lm, htilde_lm_scaled #f, amp

def find_index(data, value):
    array = np.asarray(data)
    idx = (np.abs(data - value)).argmin()
    return idx


def create_functions(h, t, hlm, metadata):
    
    frequencies_lm=[]; htilde_lm_scaled=[];
    fin_freq = (np.linalg.norm(h.angular_velocity[h.max_norm_index()]) / (2*np.pi))
    for i in hlm:
        frequencies_lm_i, htilde_lm_scaled_i = SPA_fft_calc(i[0], i[1], h, t, metadata);
        ini_freq_m = ((metadata.initial_orbital_frequency / (2*np.pi)) * i[-1]) * c**3/(G*M) * 1.15
        fin_freq_m = fin_freq * i[-1]
        ini_index_strain = find_index(frequencies_lm_i, ini_freq_m)
        fin_index_strain = find_index(frequencies_lm_i, fin_freq_m)
        cutoff_amp = htilde_lm_scaled_i[fin_index_strain] * 1e-2
        cutoff_index = find_index(htilde_lm_scaled_i[fin_index_strain:], cutoff_amp)
        #test_index = find_index(htilde_lm_scaled_i, amp_i[0])
        frequencies_lm.append(np.array(frequencies_lm_i)[ini_index_strain:fin_index_strain+cutoff_index][::12])
        htilde_lm_scaled.append(np.array(htilde_lm_scaled_i)[ini_index_strain:fin_index_strain+cutoff_index][::12])
    #print(len(f))
    frequencies_lm=np.array(frequencies_lm, dtype=object); htilde_lm_scaled=np.array(htilde_lm_scaled, dtype=object)
        
    return frequencies_lm, htilde_lm_scaled

def create_files():
    hlm = []
    for ell in range(2, 9):
        for m in range(-ell, ell + 1):
            if ell < 7:
                if m > 0 and m >= ell - 2:
                    hlm.append([ell, m])
            else:
                if m == ell:
                    hlm.append([ell, m])
                    
    binaries = ["SXS:BBH:2378", "SXS:BBH:2516", "SXS:BBH:3551", "SXS:BBH:2524", "SXS:BBH:2139", "SXS:BBH:1154", "SXS:BBH:1441", "SXS:BBH:1107", "SXS:BBH:2527", "SXS:BBH:3946", "SXS:BBH:2550", "SXS:BBH:2557", "SXS:BBH:2442", "SXS:BBH:2443", "SXS:BBH:0832"]

    strain_data = [np.array(hlm)]
    for h_id in binaries:
        metadata, h = load_strain(h_id)
        h, t = dimensionalize(h)
        frequencies, htilde = create_functions(h, t, hlm, metadata)
        metadata_info = [["Number of orbits", "Mass Ratio", "eccentricity", "chi1", "chi2", "chi1_perp", "chi2_perp"],[metadata.number_of_orbits, metadata.reference_mass_ratio, metadata.reference_eccentricity, [{", ".join(f"{c:.3g}" for c in metadata.reference_dimensionless_spin1)}], [{", ".join(f"{c:.3g}" for c in metadata.reference_dimensionless_spin2)}], metadata.reference_chi1_perp, metadata.reference_chi2_perp]]
        strain_info = [h_id, np.array(frequencies), np.array(htilde), np.array(metadata_info, dtype=object)]
        strain_data.append(np.array(strain_info, dtype=object))

    #maybe add catalog to data file too?

    return np.array(strain_data, dtype=object)
    
data = create_files()

np.savez_compressed("marimo_data", data)
    