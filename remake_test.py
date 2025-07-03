import marimo

__generated_with = "0.14.0"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    import re
    import sxs
    import math
    import numpy as np
    import pandas as pd
    import bilby
    import scipy.interpolate  
    import matplotlib.pyplot as plt
    from IPython.display import display
    import isxs_marimo as isxs
    import plotly.graph_objects as go
    import plotly.io as pio
    pio.renderers.default = 'iframe'
    return bilby, go, isxs, mo, re


@app.cell(hide_code=True)
def _(bilby, re):
    ce = bilby.gw.detector.psd.PowerSpectralDensity(psd_file = 'CE_psd.txt', asd_file = 'CE_wb_asd.txt')
    ligo_o4 = bilby.gw.detector.psd.PowerSpectralDensity(asd_file = 'aLIGO_O4_high_asd.txt')

    ce_asd_file = open(ce.asd_file, "r")
    ce_asd = ce_asd_file.readlines()
    ce_asd_file.close()
    ce_asd_amplitude = []; ce_asd_frequency = []
    for i in ce_asd:
        split_line = re.split(" |\n", i)
        ce_asd_amplitude.append(float(split_line[0]))
        ce_asd_frequency.append(float(split_line[1]))
    #print(ce_asd_frequency)

    ligo_o4_asd_file = open(ligo_o4.asd_file, "r")
    ligo_o4_asd = ligo_o4_asd_file.readlines()
    ligo_o4_asd_file.close()
    ligo_o4_asd_amplitude = []; ligo_o4_asd_frequency = []
    for i in ligo_o4_asd:
        split_line = re.split(" |\n", i)
        ligo_o4_asd_amplitude.append(float(split_line[0]))
        ligo_o4_asd_frequency.append(float(split_line[1]))
    #print(ligo_o4_asd_frequency)
    return (
        ce_asd_amplitude,
        ce_asd_frequency,
        ligo_o4_asd_amplitude,
        ligo_o4_asd_frequency,
    )


@app.cell
def _(isxs):
    hlm, h_id_list, strain_data, metadata_list = isxs.load_data()
    return h_id_list, hlm, strain_data


@app.cell
def _(h_id_list, isxs):
    idx = isxs.load_index(h_id_list, "SXS:BBH:2378")
    return (idx,)


@app.cell
def _(idx, isxs, strain_data):
    frequencies, htildes = isxs.load_plots(strain_data, idx)
    return frequencies, htildes


@app.cell
def _(mo):
    Distance = mo.ui.slider(100.0,10000.0,10.0, label="Distance (Mpc)", include_input=True)
    Mass = mo.ui.slider(5.0,10000.0,1.0, label="Mass (Solar Mass M☉)", value=33 ,include_input=True)
    return Distance, Mass


@app.cell
def _(
    Distance,
    Mass,
    ce_asd_amplitude,
    ce_asd_frequency,
    frequencies,
    go,
    hlm,
    htildes,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    mo,
):
    fig = isxs.iplt_lm(frequencies, htildes, hlm, Mass.value, Distance.value)
    fig.add_trace(go.Scatter(x=ce_asd_amplitude, y=ce_asd_frequency,
                             line=dict(color='orange', width=2),
                             name="CE Noise Curve"))
    fig.add_trace(go.Scatter(x=ligo_o4_asd_amplitude, y=ligo_o4_asd_frequency,
                             line=dict(color='orchid', width=2),
                             name="aLIGO Noise Curve"))
    mo.vstack([Distance, Mass, fig])
    return


@app.cell
def _(frequencies):
    print(frequencies[0][0])
    return


@app.cell
def _(mo):
    dropdown_MR = mo.ui.dropdown(
        options=["SXS:BBH:1154 (MR:1)", "SXS:BBH:2139 (MR:3)", "SXS:BBH:1441 (MR:8)", "SXS:BBH:1107 (MR:10)"],
        value="SXS:BBH:1154 (MR:1)",
        label="Choose a Mass Ratio",
        searchable=True,
    )
    dropdown_MR
    return (dropdown_MR,)


@app.cell
def _(dropdown_MR, h_id_list, isxs):
    idx_MR = isxs.load_index(h_id_list, dropdown_MR.value[:12])
    return (idx_MR,)


@app.cell
def _(idx_MR, isxs, strain_data):
    frequencies_MR, htildes_MR = isxs.load_plots(strain_data, idx_MR)
    return frequencies_MR, htildes_MR


@app.cell
def _(mo):
    Distance_MR = mo.ui.slider(100.0,10000.0,10.0, label="Distance (Mpc)", include_input=True)
    Mass_MR = mo.ui.slider(5.0,10000.0,1.0, label="Mass (Solar Mass M☉)", value=33 ,include_input=True)
    return Distance_MR, Mass_MR


@app.cell
def _(
    Distance_MR,
    Mass_MR,
    ce_asd_amplitude,
    ce_asd_frequency,
    frequencies_MR,
    go,
    hlm,
    htildes_MR,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    mo,
):
    fig_MR = isxs.iplt_lm(frequencies_MR, htildes_MR, hlm, Mass_MR.value, Distance_MR.value)
    fig_MR.add_trace(go.Scatter(x=ce_asd_amplitude, y=ce_asd_frequency,
                             line=dict(color='orange', width=2),
                             name="CE Noise Curve"))
    fig_MR.add_trace(go.Scatter(x=ligo_o4_asd_amplitude, y=ligo_o4_asd_frequency,
                             line=dict(color='orchid', width=2),
                             name="aLIGO Noise Curve"))
    mo.vstack([Distance_MR, Mass_MR, fig_MR])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
