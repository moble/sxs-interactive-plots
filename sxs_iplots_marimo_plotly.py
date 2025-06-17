import re
import sxs
import numpy as np
import pandas as pd
import bilby
import scipy.interpolate
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from matplotlib.widgets import CheckButtons
import ipywidgets as widgets
import marimo as mo
from IPython.display import display
import mpl_interactions.ipyplot as iplt
import sxs_iplots as isxs
import plotly.io as pio
pio.renderers.default = 'iframe'

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
    sxs_bbh_n = sxs.load("SXS:"+h_id)
    metadata = sxs_bbh_n.metadata
    h = sxs_bbh_n.h

    reference_index = h.index_closest_to(metadata.reference_time)
    h=h[reference_index:]
    print(f"Mass ratio: {metadata.reference_mass_ratio} \nReference eccentricity: {metadata.reference_eccentricity} \nReference Chi1_Perp (Precession): {metadata.reference_chi1_perp} \nReference Chi2_Perp (Precession): {metadata.reference_chi2_perp}")
    #print(int(h_id))
    #print(type(int(h_id)))
    n = sort(h_id)
    display(df[n:n+1])
    return metadata, h

def dimensionalize(h, G, c, M, r):
    h.time = h.time * G * (M/(c**3))
    h = h * (M/r) * (G/(c**2))
    t = h.t
    return h, t
    
def convertlm(m, f22, df22dt):
    flm = f22 * (m / 2)
    dflmdt = df22dt * (m / 2)  
    return flm, dflmdt

def SPA_fft_calc(l,m, h, t, metadata):
    """
    Calculate the SPA and FFT of user selected mode of strain h
    """
    
    #calculate for SPA
    h22 = h.data[:, h.index(2,2)]
    phi22 = np.unwrap(np.angle(h22))
    phi22_t = scipy.interpolate.CubicSpline(t, phi22)
    f22 = phi22_t.derivative()(t) / (2*np.pi)
    df22dt = phi22_t.derivative(2)(t) / (2*np.pi)
    hlm = h.data[:, h.index(l,m)]
    i1 = h.index_closest_to(0.5)
    i2 = h.max_norm_index()
    dt = 1e-4 #np.min(np.diff(h.t))
    Alm = np.abs(hlm)
    philm = np.unwrap(np.angle(hlm))
    philm_t = scipy.interpolate.CubicSpline(t, philm)
    flm, dflmdt = convertlm(m, f22, df22dt)

    #calculate for FFT
    h_lm = h[:, h.index(l,m)]
    h_lm_interpolated = h_lm.interpolate(np.arange(h_lm.t[0], h_lm.t[-1], dt))
    hlm_tapered = h_lm_interpolated.taper(0, h.t[0]+1000*(G*(M/(c**3))))
    hlm_transitioned = hlm_tapered.transition_to_constant(h.max_norm_time()+100*(G*(M/(c**3))), h.max_norm_time()+200*(G*(M/(c**3))))
    if type(metadata.reference_eccentricity) == float and ((metadata.reference_eccentricity) > 0.3):
        hlm_padded = hlm_transitioned.pad(100000*(G*(M/(c**3))))
    else:
        hlm_padded = hlm_transitioned.pad(25000*(G*(M/(c**3))))
    hlm_line_subtracted = hlm_padded.line_subtraction()
    #print(type(hlm_line_subtracted.ndarray))
    htilde_lm = np.fft.rfft(hlm_line_subtracted.ndarray.astype(float))*dt
    frequencies_lm = np.fft.rfftfreq(len(hlm_line_subtracted.ndarray.astype(float)), dt)
    htilde_lm_scaled = 2*(np.abs(htilde_lm))*(np.sqrt(frequencies_lm))

    #amplitude and frequency scaling
    f= -flm[i1:i2]
    amp = 2*np.abs(((1/2)*np.abs(hlm[i1:i2])) / np.sqrt(np.abs(dflmdt[i1:i2])))*np.sqrt(np.abs(flm[i1:i2]))

    return f, amp, frequencies_lm, htilde_lm_scaled

def find_index(data, value):
    array = np.asarray(data)
    idx = (np.abs(data - value)).argmin()
    return idx

def create_functions(h, t, metadata):

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
        f_i, amp_i, frequencies_lm_i, htilde_lm_scaled_i = SPA_fft_calc(i[0], i[1], h, t, metadata);
        ini_freq_m = (metadata.initial_orbital_frequency / (2*np.pi)) * i[-1]
        print(ini_freq_m)
        ini_index_SPA = find_index(f_i, ini_freq_m)
        ini_index_strain = find_index(frequencies_lm_i, ini_freq_m)
        f.append(np.array(f_i)[ini_index_SPA:][::12])
        amp.append(np.array(amp_i)[ini_index_strain:][::12])
        frequencies_lm.append(np.array(frequencies_lm_i)[::12])
        htilde_lm_scaled.append(np.array(htilde_lm_scaled_i)[::12])
    #print(len(f))
    f=np.array(f); amp=np.array(amp); frequencies_lm=np.array(frequencies_lm, dtype=object); htilde_lm_scaled=np.array(htilde_lm_scaled, dtype=object)
        
    return f, amp, frequencies_lm, htilde_lm_scaled, hlm
"""
ideally function looks like this to reduce bloat
def iplt_lm(strain, **kwargs)
    ...
    
where **kwargs includes ratio, and maybe SPA_enable if managing to make an interactive plot enable/disable is impossible

strain would be the class that includes all of the parameters needed to calculate the strain and SPA

"""
def load_plots(h, t, metadata):
    f_x_SPA_test, f_y_SPA_test, f_x_strain_test, f_y_strain_test, hlm = create_functions(h, t, metadata)
    """
    print(f_x_SPA_test(Mass, Distance))
    print(f_x_SPA_test(Mass, Distance)[0])
    print(f_x_SPA_test(Mass, Distance)[1])
    print(len(f_x_SPA_test(Mass, Distance)))
    print(len(f_y_SPA_test(Mass, Distance)))
    print(len(f_x_strain_test(Mass, Distance)))
    print(len(f_y_strain_test(Mass, Distance)))
    """
    return f_x_SPA_test, f_y_SPA_test, f_x_strain_test, f_y_strain_test, hlm
    
def iplt_lm(xSPA, ySPA, xStrain, yStrain, hlm, Mass, Distance):
    mscale = (Mass*1.989e+30)/M
    dscale = (Distance*3.086e+22)/r
    f_x_SPA = xSPA/mscale
    f_y_SPA = mscale**(3/2) * ySPA / dscale
    f_x_strain = xStrain/mscale
    f_y_strain = mscale**(3/2) * yStrain / dscale
    
    fig = go.Figure()
    for i in range(len(f_x_SPA)):
        fig.add_trace(go.Scatter(x=f_x_SPA[i], y=f_y_SPA[i],
                         line=dict(color='firebrick', width=1),
                         name=f"{hlm[i]} SPA"))
        fig.add_trace(go.Scatter(x=f_x_strain[i], y=f_y_strain[i],
                         line=dict(color='royalblue', width=1),
                         name=f"{hlm[i]} Strain"))
        #print(i)

    # Edit the layout
    fig.update_layout(
        title=dict(
            text='CE vs aLIGO comparison'
        ),
        xaxis=dict(
            range=[np.log10(3), np.log10(5e3)],
            title=dict(
                 text='frequency'
            )
        ),
        xaxis_type="log",
        yaxis=dict(
            range=[np.log10(1e-25), np.log10(5e-21)],
            title=dict(
                text='amplitude'
             ),
        ),
        yaxis_type="log"
    )
    return fig
    
def chooseplot():
    dropdown = widgets.Dropdown(options=['BBH', 'BHNS', 'NSNS'], layout={'width':'72px'})
    text_input = widgets.Text(description=":", style={'description_width': '5px'}, layout={'width':'72px'})
    button = widgets.Button(description="Plot!", layout={'width':'72px'})
    output = widgets.Output() # For displaying results

    def binary_id(b):
        with output:
            output.clear_output()
            binary_val = dropdown.value+ ":" + text_input.value
            metadatan, hn = isxs.load_strain(binary_val)
            hn, tn = isxs.dimensionalize(hn, G, c, M, r)
            isxs.iplt_lm(h=hn, t=tn, metadata=metadatan)
            #display(filtered_df)

    # Connect the button click to the function
    button.on_click(binary_id)

    display(widgets.HBox([dropdown, text_input, button]))
    display(output)

#def change_padding(h, t, metadata):
    