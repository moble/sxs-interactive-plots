import marimo

__generated_with = "0.14.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import sxs
    return mo, np, plt, sxs


@app.cell
def _(mo):
    dropdown_ecc = mo.ui.dropdown(
        options=["SXS:BBH:2553 (MR:1)", "SXS:BBH:3946 (MR:2)", "SXS:BBH:2550 (MR:4)", "SXS:BBH:2557 (MR:6)"],
        value="SXS:BBH:2553 (MR:1)",
        label="Choose a system:",
        searchable=True,
    )
    dropdown_ecc
    return (dropdown_ecc,)


@app.cell
def _(dropdown_ecc, sxs):
    sxs_bbh_1124 = sxs.load(dropdown_ecc.value[:12], extrapolation="N4")
    w = sxs_bbh_1124.h
    return (w,)


@app.cell
def _(np, w):
    dt = np.min(np.diff(w.t))
    w_interpolated = w.interpolate(np.arange(w.t[0], w.t[-1], dt))
    return dt, w_interpolated


@app.cell
def _(np, w_interpolated):
    θ, ϕ = np.pi/2, 0.0
    h = w_interpolated.evaluate(θ, ϕ)
    return (h,)


@app.cell
def _(h, plt):
    h_tapered = h.taper(0, h.t[0]+1000)

    plt.figure()
    plt.plot(h.t, h.real, lw=3, label="original")
    plt.plot(h_tapered.t, h_tapered.real, ls="--", label="tapered")
    plt.xlabel("Time")
    plt.ylabel("Strain $h_+$")
    plt.legend();
    plt.show()
    return (h_tapered,)


@app.cell
def _(h_tapered, plt, w):
    h_transitioned = h_tapered.transition_to_constant(w.max_norm_time()+100)#, w.max_norm_time()+200)

    plt.figure()
    y_f = h_transitioned[-1].real
    plt.plot(h_tapered.t, h_tapered.real, lw=3, label="before transition")
    plt.plot(h_transitioned.t, h_transitioned.real, ls="--", label="after transition")
    plt.xlim(w.max_norm_time(), h_transitioned.t[-1])
    plt.ylim(y_f - 4e-3, y_f + 4e-3)
    plt.xlabel("Time")
    plt.ylabel("Strain amplitude")
    plt.legend();
    plt.show()
    return (h_transitioned,)


@app.cell
def _(h_transitioned, plt):
    h_padded = h_transitioned.pad(25000)

    plt.figure()
    plt.plot(h_transitioned.t, h_transitioned.real, lw=3, label="windowed and tapered")
    plt.plot(h_padded.t, h_padded.real, ls="--", label="padded")
    plt.xlabel("Time")
    plt.ylabel("Strain $h_+$")
    plt.legend();
    plt.show()
    return (h_padded,)


@app.cell
def _(h_padded, plt):
    h_line_subtracted = h_padded.line_subtraction()

    plt.figure()
    plt.plot(h_padded.t, h_padded.real, lw=3, label="padded")
    plt.plot(h_line_subtracted.t, h_line_subtracted.real, ls="--", label="line-subtracted")
    plt.xlabel("Time")
    plt.ylabel("Strain $h_+$")
    plt.legend();
    plt.show()
    return


@app.cell
def _(mo):
    dropdown_lm = mo.ui.dropdown(
        options=["2,1", "2,2", "3,1", "3,2", "3,3", "4,2", "4,3", "4,4", "5,3", "5,4", "5,5", "6,4", "6,5", "6,6", "7,7", "8,8"],
        value="2,2",
        label="Choose an l,m mode:",
        searchable=True,
    )
    dropdown_lm
    return (dropdown_lm,)


@app.cell
def _(dropdown_lm, dt, h, np, w):
    l = int(dropdown_lm.value[0])
    m = int(dropdown_lm.value[-1])
    h_lm = w[:, w.index(l,m)]
    h_lm_interpolated = h_lm.interpolate(np.arange(h_lm.t[0], h_lm.t[-1], dt))
    hlm_tapered = h_lm_interpolated.taper(0, h.t[0]+1000)
    hlm_transitioned = hlm_tapered.transition_to_constant(w.max_norm_time()+100)
    hlm_padded = hlm_transitioned.pad(250000)
    hlm_line_subtracted = hlm_padded.line_subtraction().real
    return hlm_line_subtracted, hlm_padded


@app.cell
def _(hlm_line_subtracted, hlm_padded, plt):
    plt.figure()
    plt.plot(hlm_padded.t, hlm_padded.real, lw=3, label="padded")
    plt.plot(hlm_line_subtracted.t, hlm_line_subtracted.real, ls="--", label="line-subtracted")
    plt.xlabel("Time")
    plt.ylabel("Strain $h_+$")
    plt.legend();
    plt.show()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
