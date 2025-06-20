import marimo

__generated_with = "0.14.0"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Motivation!""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""The development of third generation detectors like Cosmic Explorer (CE) and Einstein Telescope (ET) marks an exciting new chapter in gravitational wave detection. For one, Cosmic Explorer has a lower sensitivity curve than VIRGO or Advanced LIGO and is capable of detecting gravitational waves at a slightly lower frequency and in a much wider frequency band (here are some numbers, if you're curious: Cosmic Explorer's frequency band starts at $5$ $Hz$ compared to aLIGO/VIRGO/KAGRA's lower limit at $10$ $Hz$. CE's strain sensitivity is better than ~$10^{-23} Hz^{-½}$ starting at ~$5.7$ $Hz$ compared to Advanced LIGO achieving ~$10^{-23} Hz^{-½}$ at ~$100$ $Hz$). These lower frequencies made accesible by CE is |important for detecting intermediate sized BH binaries ranging from ~100M☉ to ~1000M☉ (also include a piece regarding neutron star binaries and maybe black hole neutron star binaries if able)""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""So the main questions this notebook is trying to address are with regards to both the range/limits of detectable signals by Cosmic Explorer and more generally, the limits of current SXS simulations. For example, for a variety of different spins, eccentricities, etc characterized in various BHBH system, to what limits/ranges of mass and distance are we able to see them with CE? Are current models sufficient enough for matched filtering with possible signals, or are we missing information? SXS currently decomposes gravitational waves strains up to their (8,8) modes - would we need more modes to model a more massive BHBH system with sufficient enough accuracy (and if so, how many more?)? Cosmic Explorer's frequency band starts at a lower frequency - do our SXS simulations need more orbits or are our modelled inspirals long enough? Since the early inspiral is better modeled using Post-Newtonian, more orbits would provide more motivation for joining Post-Newtonian with Numerical Relativity for our simulations.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""This notebook is designed to let you see some examples and play around with the masses and distances of various BHBH systems, to see what these BHBH simulations with various mass ratios, eccentricities, and precession look like with respect to the CE noise curve. So let's delve right in!""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Importing and Initializing""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Let's start by just importing the necessary packages used in this notebook to run these examples...""")
    return


@app.cell
def _():
    # '%matplotlib ipympl' command supported automatically in marimo

    import re
    import sxs
    import numpy as np
    import bilby
    import scipy.interpolate  
    from scipy.signal import argrelextrema
    import matplotlib.pyplot as plt
    from matplotlib.widgets import CheckButtons
    import ipywidgets as widgets
    from IPython.display import display
    import mpl_interactions.ipyplot as iplt
    import sxs_iplots_marimo_plotly as isxs
    import plotly.graph_objects as go
    import plotly.io as pio
    pio.renderers.default = 'iframe'
    return bilby, go, isxs, np, plt, re, sxs


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""We should also set some basic constants we need here for all of our calculations later on (G is the gravitational constant, c is the speed of light, M is the mass of the whole binary in kg and r is the distance from the detector to the the binary in meters). Though since kg and m are just a tad inappropriate for the context we are in, we scale M and r to be in the more astronomical units of solar masses (M☉) and megaparsecs (Mpc) respectively (though the 'units' are kept in kg and meters for calculation purposes only). We initialize M to be 33 M☉ and r to be 100 Mpc""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    <font color='orange'>note to myself: what if we just change the initialized M to be something like 1 for just simplification purposes? Since binary total mass likely can't be as small as 1 solar mass, we could also maybe set it to 5? Unless you want to include examples for NSNS systems in the future, but this is just food for thought. 33 solar masses as the initial mass does feel a little too arbitrary and restrictive, since there's definitely BHBH binaries with a smaller total mass </font>
    <font color='magenta'>-j 5/29
    """
    )
    return


@app.cell
def _():
    # constants
    G = 6.67430e-11
    c = 3e8
    M = 6.563e+31
    r = 3.086e24
    dt = 1e-4 #np.min(np.diff(h.t))
    return G, M, c, r


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Import and create noise curves for LIGO O4 and Cosmic Explorer""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Here we download and extract some important data for the Cosmic Explorer noise curve in the frequency domain. We also download the noise curve of the most recent observation run of Advanced Ligo just to compare the two""")
    return


@app.cell
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


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Varying <span style="color:red">Mass Ratios</span> Examples #""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Okay now that that's out of the way let's delve into our examples! This first one will include a brief walk through of what's going on with the calculations behind the code and their (*brief*) explanations. The subsequent examples follow the same process, just using a different simulation for the purposes of seeing what changes with different mass ratios/eccentricities/precessions. At the very bottom there's an included 'choose your own waveform' for you to test out, if you're curious about any particular binary mass ratio, eccentricity, precession, etc (or any of their combinations)! <font color='orange'> note to myself: please don't forget to implement this! <font color='magenta'> -j 5/29""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Example 1(a): BBH:2139 (Mass Ratio: 3)""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""So first to start off, let's just first load the strain in and cut off the first bit of junk radiation before the surrogate signal even starts. Let's also get some basic information like the mass ratio, precession and eccentricity about the binary system. <font color='orange'> note to myself: maybe it would be better and more clear/useful to include all of the information of the simulated system in a table? <font color='magenta'> -j 5/29""")
    return


@app.cell
def _(isxs):
    metadata2139, h2139 = isxs.load_strain("BBH:2139")
    return h2139, metadata2139


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Just for the sake of the first example, it may be useful to have a visualization of what's going on so far. We've just loaded in a simulation of interest, so let's see what this simulation looks like in the time domain!""")
    return


@app.cell
def _(h2139, plt):
    plt.figure()
    plt.plot(h2139.t, h2139.data.view(float))
    plt.title(f"Extrapolated waveform")
    plt.xlabel(r"$(t_{\mathrm{corr}} - r_\ast)/M$")
    plt.ylabel(r"$r\, h^{(\ell,m)}/M$");
    plt.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Great! Now notice that the units for the x and y axes aren't quite the same as the units for time nor strain amplitude. You can see in the SXS metadata that there is data for the strain's amplitude and 'time', but they're labeled as something like rhOverM. This is because these 'units' are *actually* unitless. Since it's much more useful to generalize these for any particular strain, time and amplitude are made and saved as dimensionless in the metadata. But for our purposes, we definitely do care about the units so let's revert these back and make them dimensionful.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    So for a very very brief explanation on how to do this, let's see what's going on first: When measuring gravitational waves, it's useful to note that the gravitational wave strain $h$ is inversely proportional to the distance $r$ just by $$h \propto \frac{1}{r}$$
    But since the universe is quite large, $r$ can be assumed to $$r \rightarrow \infty$$
    So, we tend to measure and use the value $rh$ instead (which has units of distance).
    Now for simplicity and accuracy, we work in geometrized units, where distance $d$ and time $t$ are written as in terms of constants like $G$, $c$, and $M☉$, so $$d \propto \frac{GM☉}{c^2}$$ and $$t \propto \frac{GM☉}{c^3}$$
    Thus, while what LIGO physically sees and measures are indeed $rh$ and $t$, the amplitude strain $H$ and time $T$ that's in the metadata saved in SXS are unitless, where $$H \propto \frac{rh}{\frac{GM☉}{c^2}} \qquad \text{or,} \qquad h \propto \frac{\frac{GM☉}{c^2}}{r} * H$$
    and similarly, $$T \propto \frac{t}{\frac{GM☉}{c^3}} \qquad \text{or,} \qquad t \propto \frac{GM☉}{c^3} * T$$
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Okay great! Now that we know this, let's dimensionalize our strain's amplitude and time!""")
    return


@app.cell
def _(G, M, c, h2139, isxs, r):
    h2139_1, t2139 = isxs.dimensionalize(h2139, G, c, M, r)
    return h2139_1, t2139


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Now we can *really* get started with some calculations. Our goal is to compare our surrogate models with aLIGO O4 and Cosmic Explorer's noise curves, so we care more about what our signal looks like in the frequency domain rather than the time domain. So the easiest way to achieve this is to take the Fast Fourier Transform of our signal in the time domain. But before we do that, let's decompose our strain into its (l,m) modes. Because gravitational waves propagate in three dimensions, we model them using spherical harmonics. The lowest mode that SXS goes down to is (8,8) so we decompose our strain into all modes from (8,8) to (2,2) (except the m=0 modes).""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Now we *can* (and will!) take the FFT of all this raw data of every (l,m) mode, but we let's take another additional step for simplification and clarity purposes. The raw data of all these modes together can look quite chaotic and it's nontrivial to slice off all the noise from each mode, so it can be useful to have a way of looking at them without all of the messiness. To get a feel for what the amplitude of each (l,m) mode's Fourier transform looks like, specifically for their amplitudes during inspiral, we can also apply the method of stationary phase to each mode. You can apply the stationary phase approximation (SPA) method whenever you have some arbitrary function $B(t)$ of the form $B(t) = A(t)*\cos(\phi(t))$ where $\dfrac{d\ln A}{dt} \ll \dfrac{d\phi}{dt}$ (or in other words, $A$ varies much slower than $\phi$) and $\omega = 2\pi f = \dfrac{d\phi}{dt}$ is strictly monotonic. Then, you can approximate the Fourier transform of $B(t)$ as $$B(t) \approx A(t(f)) * (\dfrac{df}{dt})^{-\frac{1}{2}}$$
    (Or more intuitively, you can think of it as if you have a sinusoid with a rapidly changing phase. If this sinusoid has a region where the phase change is 'slower' - i.e. where the 'stationary phase' is - then this 'slower' region will dominate the Fourier transform integral while the faster, more rapidly changing regions will cancel out and can effectively be ignored.)
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Note that this approximation is good for modeling the early stages of a binary's inspiral, but starts to breaks down as the binary approach merger, since the frequency starts to increase more and more dramatically as the two black holes get closer and closer to coalescing. The early phases of the inspiral has a slower orbit, and thus is the location of their stationary phase.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Okay! So to recap: We're interested in looking at and comparing the (l,m) modes of our strain in the frequency domain to the noise curves of aLIGO and CE, so we break down the strain into the (2,2) to (8,8) modes and take the FFT of each. Since the *characteristic* strain is more useful when comparing to detector noise curves, this is what we will plot (we calculate it by multiplying the absolute value of our strain FFT by $\sqrt{f}$). We also find the SPA of all of these modes for a good approximation of what each mode looks like early on in the inspiral. Let's see what all of this finally looks like combined with the noise curves. The code cell below will complete all of the calculations behind the scenes, and run a plot of every mode together with its SPA and the aLIGO/CE noise curves.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    <font color='red'>(AS OF 5/29) Next steps: Finish explaining decomposing the strain into (2,2) to (8,8) modes. Mention how we calculate (2,2) mode first, and scale the frequency and amplitude of all other modes from that one. Also should convince yourself and explain why we don't include the m=0 modes. You can include a plot of all the (l,m) modes in the time domain for visualization. After that, mention how for working in the frequency domain, it's useful to model and measure the amplitude of our modes using the SPA of their respective modes. In particular, it's useful to approximate the earlier inspiral using SPA since the orbits are less frequent. Here you can give a brief explanation of how SPA works. After all this, we take the FFTs of both the actual strain and the SPA of each (l,m) mode and plot them together with also the two noise curves. This is how we see how our current models hold up to the next-gen noise curves (how much can be detected, if not all, and if we need more information like more modes or more orbits).  Then, talk about mass and distance scaling. Then in the next sets of examples you can include respective explanations about precession and eccentricity </font>
    <font color='magenta'>-j 5/29
    """
    )
    return


@app.cell
def _(mo):
    Distance = mo.ui.slider(100.0,10000.0,10.0, label="Distance")
    Distance
    return (Distance,)


@app.cell
def _(mo):
    Mass = mo.ui.slider(33.0,10000.0,10.0, label="Mass")
    Mass
    return (Mass,)


@app.cell(hide_code=True)
def _(
    Distance,
    Mass,
    ce_asd_amplitude,
    ce_asd_frequency,
    go,
    h2139_1,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    metadata2139,
    t2139,
):
    xSPA, ySPA, xStrain, yStrain, hlm = isxs.load_plots(h=h2139_1, t=t2139, metadata=metadata2139)
    fig = isxs.iplt_lm(xSPA, ySPA, xStrain, yStrain, hlm, Mass.value, Distance.value)
    fig.add_trace(go.Scatter(x=ce_asd_amplitude, y=ce_asd_frequency,
                             line=dict(color='orange', width=2),
                             name="CE Noise Curve"))
    fig.add_trace(go.Scatter(x=ligo_o4_asd_amplitude, y=ligo_o4_asd_frequency,
                             line=dict(color='orchid', width=2),
                             name="aLIGO Noise Curve"))
    fig
    return (xStrain,)


@app.cell
def _(xStrain):
    xStrain
    return


@app.cell
def _(h2139):
    h2139.data
    return


@app.cell
def _(h2139):
    h2139.t
    return


@app.cell
def _(metadata2139):
    omegai = metadata2139.initial_orbital_frequency
    omegai
    return (omegai,)


@app.cell
def _(np, omegai):
    fi = omegai/(2*np.pi)
    fi
    return (fi,)


@app.cell
def _(fi):
    1/fi
    return


@app.cell
def _(metadata2139):
    metadata2139.reference_orbital_frequency
    return


@app.cell
def _(metadata2139, np):
    np.linalg.norm(metadata2139.reference_orbital_frequency)
    return


@app.cell
def _(h2139, np):
    (np.linalg.norm(h2139.angular_velocity[h2139.max_norm_index()]) / (2*np.pi)) 
    return


@app.cell
def _(h2139):
    h2139.max_norm_index()+100
    return


@app.cell
def _(h2139):
    h2139.angular_velocity[14967]
    return


@app.cell
def _(np):
    np.linalg.norm([-5.01125441e-05, -3.49534150e-04,  9.38034656e+02])
    return


@app.cell
def _(G, M, c):
    c**3 / (G*M)
    return


@app.cell
def _(h2139):
    h2139.angular_velocity
    return


@app.cell
def _(sxs):
    h4001 = sxs.load("SXS:BBH:4001").h
    return (h4001,)


@app.cell
def _(h4001, np):
    (np.linalg.norm(h4001.angular_velocity[h4001.max_norm_index()]) / (2*np.pi)) 
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""<font color=orange> I remember this being pretty annoying to try to implement, but I still think it will be useful to include captions or some sort of label for each (l,m) modes so you can actually see which is which <font color="magenta"> -j 6/3""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""This plot is interactive, so you can change the mass and distance of this example binary via their respective sliders above the plot! If you're curious about looking at specifically the SPA's or the strains, you can do so via toggling the buttons on the top right. For a very very brief explanation on how the mass and distance scaling works, we let $a$ be the scaling factor for mass $M$ and $b$ be the scaling factor for distance $R$ such that $$M' \rightarrow a*M$$ and $$R' \rightarrow b*R$$""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    How the frequency of the strain/SPA changes by these scaling factors is relatively simple to see - from our geometrized units time is measured as $$t \propto \frac{GM☉}{c^3}$$
    so as $M' \rightarrow a*M$, $$t' = a*t$$
    But frequency is the inverse of time, so $$\boxed{f' = \frac{1}{a} *f}$$
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    To look at the amplitude, recall that $$h = \frac{\frac{GM☉}{c^2}}{r}*H \qquad \text{and} \qquad t = \frac{GM☉}{c^3}*T$$
    And that the Fourier transform of strain h is given by $$\tilde{h} = \int h(t)*e^{2\pi ift}dt$$
    A factor of $\frac{a^2}{b}$ comes from the $h*dt$ terms. But, since we plot the characteristic strain, we plot $|{\tilde{h}}|\sqrt{f}$, so the final scaling factor ends up being $$\boxed{A' = \frac{a^{3/2}}{b}*A}$$
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Please feel free to play around with these plots and the sliders to test out various masses and distances! There are also plenty of examples down below (including ones with high precession or high eccentricity), so if you're interested, please check them out! Enjoy!""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Example 1(b): BBH1441 (Mass Ratio: 8)""")
    return


@app.cell
def _(isxs):
    metadata1441, h1441 = isxs.load_strain("BBH:1441")
    return h1441, metadata1441


@app.cell
def _(G, M, c, h1441, isxs, r):
    h1441_1, t1441 = isxs.dimensionalize(h1441, G, c, M, r)
    return h1441_1, t1441


@app.cell
def _(
    ce_asd_amplitude,
    ce_asd_frequency,
    h1441_1,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    metadata1441,
    plt,
    t1441,
):
    isxs.iplt_lm(h=h1441_1, t=t1441, metadata=metadata1441)
    plt.loglog(ce_asd_amplitude, ce_asd_frequency, label='CE', color='orange')
    plt.loglog(ligo_o4_asd_amplitude, ligo_o4_asd_frequency, label='LIGO O4', color='red')
    plt.draw()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Example 1(c): BBH1107 (Mass Ratio:10)""")
    return


@app.cell
def _(isxs):
    metadata1107, h1107 = isxs.load_strain("BBH:1107")
    return h1107, metadata1107


@app.cell
def _(metadata1107):
    metadata1107.reference_eccentricity
    return


@app.cell
def _(G, M, c, h1107, isxs, r):
    h1107_1, t1107 = isxs.dimensionalize(h1107, G, c, M, r)
    return h1107_1, t1107


@app.cell
def _(
    ce_asd_amplitude,
    ce_asd_frequency,
    h1107_1,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    metadata1107,
    plt,
    t1107,
):
    isxs.iplt_lm(h=h1107_1, t=t1107, metadata=metadata1107)
    plt.loglog(ce_asd_amplitude, ce_asd_frequency, label='CE', color='orange')
    plt.loglog(ligo_o4_asd_amplitude, ligo_o4_asd_frequency, label='LIGO O4', color='red')
    plt.draw()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Example 1(d): BBH1154 (Mass Ratio: 1)""")
    return


@app.cell
def _(isxs):
    metadata1154, h1154 = isxs.load_strain("BBH:1154")
    return h1154, metadata1154


@app.cell
def _(metadata1154):
    metadata1154.reference_eccentricity
    return


@app.cell
def _(G, M, c, h1154, isxs, r):
    h1154_1, t1154 = isxs.dimensionalize(h1154, G, c, M, r)
    return h1154_1, t1154


@app.cell
def _(
    ce_asd_amplitude,
    ce_asd_frequency,
    h1154_1,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    metadata1154,
    plt,
    t1154,
):
    isxs.iplt_lm(h=h1154_1, t=t1154, metadata=metadata1154, cut=False)
    plt.loglog(ce_asd_amplitude, ce_asd_frequency, label='CE', color='orange')
    plt.loglog(ligo_o4_asd_amplitude, ligo_o4_asd_frequency, label='LIGO O4', color='red')
    plt.draw()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# High <span style="color:green">Eccentricity</span> Examples. #""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Here are some examples of binaries with high eccentricity. You'll notice that the SPA here isn't particularly useful and can be quite chaotic (or may not even show up). This is because the SPA can only be used if the rate of change of the phase $\phi$ is strictly monotonic (if $f$ is monotonic, we interchange the phase $\phi$ and time $t$ in these cases, so in other words if $t$ is strictly monotonic). But here, frequency will no longer be strictly monotonic, as eccentric orbits imply that the frequency of these orbits at periapsis will be relatively higher than the frequency of apoapsis (so frequencies will vary between these points and there will be more than one time $t$ for a given frequency $f$.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""For some of these plots, you'll notice that the SPA starts to match the strain at higher frequencies close to merger. This is because as the binary approaches merger, the orbits will continue to lose energy to the point where the orbits become less eccentric and more circular.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Example 2(a): BBH1361 (Mass Ratio: 6; High Eccentricity)""")
    return


@app.cell
def _(isxs):
    metadata2527, h2527 = isxs.load_strain("BBH:2527")
    return h2527, metadata2527


@app.cell
def _(h2527, mo, plt):
    #plt.figure()
    plt.plot(h2527.t, h2527.data.view(float))
    plt.title(f"Extrapolated waveform")
    plt.xlabel(r"$(t_{\mathrm{corr}} - r_\ast)/M$")
    plt.ylabel(r"$r\, h^{(\ell,m)}/M$");
    mo.mpl.interactive(plt.gcf())
    #plt.show()
    return


@app.cell
def _(G, M, c, h2527, isxs, r):
    h2527_1, t2527 = isxs.dimensionalize(h2527, G, c, M, r)
    return h2527_1, t2527


@app.cell
def _(Mass):
    print(type(Mass.value))
    return


@app.cell
def _():
    return


@app.cell
def _(
    Distance,
    Mass,
    ce_asd_amplitude,
    ce_asd_frequency,
    h2527_1,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    metadata2527,
    plt,
    t2527,
):
    isxs.iplt_lm(h=h2527_1, t=t2527, metadata=metadata2527, Mass=Mass.value, Distance=Distance.value, cut=False)
    plt.loglog(ce_asd_amplitude, ce_asd_frequency, label='CE', color='orange')
    plt.loglog(ligo_o4_asd_amplitude, ligo_o4_asd_frequency, label='LIGO O4', color='red')
    plt.draw()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Example 2(b): BBH1360 (Mass Ratio: 1; High Eccentricity)""")
    return


@app.cell
def _(isxs):
    metadata1360, h1360 = isxs.load_strain("BBH:1360")
    return h1360, metadata1360


@app.cell
def _(G, M, c, h1360, isxs, r):
    h1360_1, t1360 = isxs.dimensionalize(h1360, G, c, M, r)
    return h1360_1, t1360


@app.cell
def _(
    ce_asd_amplitude,
    ce_asd_frequency,
    h1360_1,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    metadata1360,
    plt,
    t1360,
):
    isxs.iplt_lm(h=h1360_1, t=t1360, metadata=metadata1360, cut=False)
    plt.loglog(ce_asd_amplitude, ce_asd_frequency, label='CE', color='orange')
    plt.loglog(ligo_o4_asd_amplitude, ligo_o4_asd_frequency, label='LIGO O4', color='red')
    plt.draw()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""where did the SPA go..""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Example 2(c): BBH2612 (Mass Ratio: 1 ; High Eccentricity)""")
    return


@app.cell
def _(isxs):
    metadata2612, h2612 = isxs.load_strain("BBH:2612")
    return h2612, metadata2612


@app.cell
def _(G, M, c, h2612, isxs, r):
    h2612_1, t2612 = isxs.dimensionalize(h2612, G, c, M, r)
    return h2612_1, t2612


@app.cell
def _(
    ce_asd_amplitude,
    ce_asd_frequency,
    h2612_1,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    metadata2612,
    plt,
    t2612,
):
    isxs.iplt_lm(h=h2612_1, t=t2612, metadata=metadata2612, cut=False)
    plt.loglog(ce_asd_amplitude, ce_asd_frequency, label='CE', color='orange')
    plt.loglog(ligo_o4_asd_amplitude, ligo_o4_asd_frequency, label='LIGO O4', color='red')
    plt.draw()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# High <span style="color:lightblue">Precession</span> Examples #""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Similarly to the high eccentricity plots, the SPA's here also won't be too helpful.""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Example 3(a): BBH2442 (Mass Ratio: 1; High Precession)""")
    return


@app.cell
def _(isxs):
    metadata2442, h2442 = isxs.load_strain("BBH:2442")
    return h2442, metadata2442


@app.cell
def _(G, M, c, h2442, isxs, r):
    h2442_1, t2442 = isxs.dimensionalize(h2442, G, c, M, r)
    return h2442_1, t2442


@app.cell
def _(
    ce_asd_amplitude,
    ce_asd_frequency,
    h2442_1,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    metadata2442,
    plt,
    t2442,
):
    isxs.iplt_lm(h=h2442_1, t=t2442, metadata=metadata2442, cut=False)
    plt.loglog(ce_asd_amplitude, ce_asd_frequency, label='CE', color='orange')
    plt.loglog(ligo_o4_asd_amplitude, ligo_o4_asd_frequency, label='LIGO O4', color='red')
    plt.draw()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    SPA goes crazy! why? 
    (lm modes are crazy too!)
    """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Example 3(b): BBH0832 (Mass Ratio: 2; High Precession)""")
    return


@app.cell
def _(isxs):
    metadata0832, h0832 = isxs.load_strain("BBH:0832")
    return h0832, metadata0832


@app.cell
def _(G, M, c, h0832, isxs, r):
    h0832_1, t0832 = isxs.dimensionalize(h0832, G, c, M, r)
    return h0832_1, t0832


@app.cell
def _(
    ce_asd_amplitude,
    ce_asd_frequency,
    h0832_1,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    metadata0832,
    plt,
    t0832,
):
    isxs.iplt_lm(h=h0832_1, t=t0832, metadata=metadata0832, cut=False)
    plt.loglog(ce_asd_amplitude, ce_asd_frequency, label='CE', color='orange')
    plt.loglog(ligo_o4_asd_amplitude, ligo_o4_asd_frequency, label='LIGO O4', color='red')
    plt.draw()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Example 3(c): BBH2443 (Mass Ratio: 1; High Precession)""")
    return


@app.cell
def _(isxs):
    metadata2443, h2443 = isxs.load_strain("BBH:2443")
    return h2443, metadata2443


@app.cell
def _(G, M, c, h2443, isxs, r):
    h2443_1, t2443 = isxs.dimensionalize(h2443, G, c, M, r)
    return h2443_1, t2443


@app.cell
def _(
    ce_asd_amplitude,
    ce_asd_frequency,
    h2443_1,
    isxs,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    metadata2443,
    plt,
    t2443,
):
    isxs.iplt_lm(h=h2443_1, t=t2443, metadata=metadata2443, cut=False)
    plt.loglog(ce_asd_amplitude, ce_asd_frequency, label='CE', color='orange')
    plt.loglog(ligo_o4_asd_amplitude, ligo_o4_asd_frequency, label='LIGO O4', color='red')
    plt.draw()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Can use the metadata to make an if case for the eccentric case for extra padding""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Choose your own plot!""")
    return


@app.cell
def _(sxs):
    df = sxs.load("dataframe", tag="3.0.0")
    df
    return


@app.cell
def _(isxs):
    isxs.chooseplot()
    return


@app.cell
def _(
    ce_asd_amplitude,
    ce_asd_frequency,
    ligo_o4_asd_amplitude,
    ligo_o4_asd_frequency,
    plt,
):
    # Add the noise curves to the interactive plot
    plt.loglog(ce_asd_amplitude, ce_asd_frequency, label="CE", color = "orange")
    plt.loglog(ligo_o4_asd_amplitude, ligo_o4_asd_frequency, label="LIGO O4", color = "red")
    plt.draw()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""<font color='orange'> as of 6/5: First draft of the explanations I would say is more or less complete. I'm going to push everything onto github for review, and we'll see what to change/fix/add. In terms of other changes I'm thinking of adding: captions/labels of all the modes will be helpful, though I'm sure a pain to implement. Also make sure to just re-add all the other missing modes, not just the ones that'll appear above the noise curves initially. You can also change the initial mass/add more masses. I want to also figure out how to deal with all the deprecated simulations - be able to plot them or no? I also want to make the table interactive, but maybe we can worry about that when/if we push everything on marimo?""")
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
