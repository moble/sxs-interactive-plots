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
    return bilby, go, isxs, mo, np, re


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
    return (
        ce_asd_amplitude,
        ce_asd_frequency,
        ligo_o4_asd_amplitude,
        ligo_o4_asd_frequency,
    )


@app.cell
def _(isxs):
    hlm, h_id_list, strain_data, metadata_list = isxs.load_data()
    return h_id_list, hlm, metadata_list, strain_data


@app.cell
def _(h_id_list, isxs):
    idx = isxs.load_index(h_id_list, "SXS:BBH:2442")
    return (idx,)


@app.cell
def _(idx, isxs, strain_data):
    frequencies, htildes = isxs.load_plots(strain_data, idx)
    return frequencies, htildes


@app.cell
def _(hlm):
    hlm
    return


@app.cell
def _(frequencies):
    len(frequencies[0])
    return


@app.cell
def _(frequencies, np):
    np.logspace(0,np.log10(len(frequencies[0])))
    return


@app.cell
def _(frequencies):
    frequencies[0]
    return


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
    #dropdown_MR
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
    dropdown_MR,
    frequencies_MR,
    go,
    hlm,
    htildes_MR,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    markdown,
    mo,
):
    fig_MR = isxs.iplt_lm(frequencies_MR, htildes_MR, hlm, Mass_MR.value, Distance_MR.value)
    fig_MR.add_trace(go.Scatter(x=ce_asd_amplitude, y=ce_asd_frequency,
                             line=dict(color='orange', width=2),
                             name="CE Noise Curve"))
    fig_MR.add_trace(go.Scatter(x=ligo_o4_asd_amplitude, y=ligo_o4_asd_frequency,
                             line=dict(color='orchid', width=2),
                             name="aLIGO Noise Curve"))
    mo.vstack([dropdown_MR, mo.md("-----------------------------"), Distance_MR, Mass_MR, fig_MR, markdown])
    return


@app.cell
def _(frequencies_MR):
    len(frequencies_MR[1])
    return


@app.cell
def _(metadata_list):
    metadata_list[0]
    return


@app.cell
def _(dropdown_MR, idx_MR, metadata_list, mo):
    line1 = mo.md(f"""<h1 style="font-size: 24px;">{dropdown_MR.value[:12]} Metadata Info:</h1>""")
    line3 = mo.md(
        f"""
        n orbits: {metadata_list[idx_MR][1][0]:.3g}  
        mass ratio: {metadata_list[idx_MR][1][1]:.3g}  
        eccentricity: {metadata_list[idx_MR][1][2]:.3g}  
        chi1: {metadata_list[idx_MR][1][3]}  
        chi2: {metadata_list[idx_MR][1][4]}  
        chi1_perp: {metadata_list[idx_MR][1][5]:.3g}  
        chi2_perp: {metadata_list[idx_MR][1][6]:.3g}
        """
    )
    markdown = mo.vstack([line1, mo.md("-----------"), line3])
    return (markdown,)


@app.cell
def _(mo):
    dropdown_ecc = mo.ui.dropdown(
        options=["SXS:BBH:2527 (MR:1)", "SXS:BBH:3946 (MR:2)", "SXS:BBH:2550 (MR:4)", "SXS:BBH:2557 (MR:6)"],
        value="SXS:BBH:2527 (MR:1)",
        label="Choose a system:",
        searchable=True,
    )
    return (dropdown_ecc,)


@app.cell
def _(dropdown_ecc, h_id_list, isxs):
    idx_ecc = isxs.load_index(h_id_list, dropdown_ecc.value[:12])
    return (idx_ecc,)


@app.cell
def _(idx_ecc, isxs, strain_data):
    frequencies_ecc, htildes_ecc = isxs.load_plots(strain_data, idx_ecc)
    return frequencies_ecc, htildes_ecc


@app.cell
def _(frequencies_ecc):
    len(frequencies_ecc[1])
    return


@app.cell
def _(mo):
    Distance_ecc = mo.ui.slider(100.0,10000.0,10.0, label="Distance (Mpc)", include_input=True)
    Mass_ecc = mo.ui.slider(5.0,10000.0,1.0, label="Mass (Solar Mass M☉)", value=33 ,include_input=True)
    return Distance_ecc, Mass_ecc


@app.cell
def _(
    Distance_ecc,
    Mass_ecc,
    ce_asd_amplitude,
    ce_asd_frequency,
    dropdown_ecc,
    frequencies_ecc,
    go,
    hlm,
    htildes_ecc,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    mo,
):
    fig_ecc = isxs.iplt_lm(frequencies_ecc, htildes_ecc, hlm, Mass_ecc.value, Distance_ecc.value)
    fig_ecc.add_trace(go.Scatter(x=ce_asd_amplitude, y=ce_asd_frequency,
                             line=dict(color='orange', width=2),
                             name="CE Noise Curve"))
    fig_ecc.add_trace(go.Scatter(x=ligo_o4_asd_amplitude, y=ligo_o4_asd_frequency,
                             line=dict(color='orchid', width=2),
                             name="aLIGO Noise Curve"))
    #fig_ecc.add_vline(x=xvalue)
    mo.vstack([dropdown_ecc, mo.md("-----------------------------"), Distance_ecc, Mass_ecc, fig_ecc])
    return


@app.cell
def _():
    #arr = np.array([1, 1.5, 2, 2.4, 2.5, 2.8, 4, 5, 6, 10, 12, 23, 30, 31, 40])
    return


@app.cell
def _():
    """
    i = 0
    while arr[i] != arr[-1]:
        if arr[i+1] < arr[i] * 2:
            arr = np.delete(arr, i+1)
            htilde = np.delete(htilde, i+1)
        else:
            i+=1
    """
    return


@app.cell
def _(hlm):
    hlm
    return


@app.cell
def _(mo):
    dropdown_prec = mo.ui.dropdown(
        options=["SXS:BBH:2442 (MR:1)", "SXS:BBH:2443 (MR:1)", "SXS:BBH:0832 (MR:2)"],
        value="SXS:BBH:2442 (MR:1)",
        label="Choose a system:",
        searchable=True,
    )
    return (dropdown_prec,)


@app.cell
def _(mo):
    Distance_prec = mo.ui.slider(100.0,10000.0,10.0, label="Distance (Mpc)", include_input=True)
    Mass_prec = mo.ui.slider(5.0,10000.0,1.0, label="Mass (Solar Mass M☉)", value=33 ,include_input=True)
    return Distance_prec, Mass_prec


@app.cell
def _(
    Distance_prec,
    Mass_prec,
    dropdown_prec,
    h_id_list,
    hlm,
    isxs,
    metadata_list,
    strain_data,
):
    isxs.run(dropdown_prec.value[:12], h_id_list, strain_data, metadata_list, hlm, Mass_prec, Distance_prec, dropdown_prec)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
