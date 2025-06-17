# test code cell
"""
xSPA, ySPA, xStrain, yStrain, hlm = isxs.load_plots(h=h2139_1, t=t2139, metadata=metadata2139)
mscale = (Mass.value*1.989e+30)/M
dscale = (Distance.value*3.086e+22)/r
f_x_SPA = xSPA/mscale
f_y_SPA = mscale**(3/2) * ySPA / dscale
f_x_strain = xStrain/mscale
f_y_strain = mscale**(3/2) * yStrain / dscale

fig = go.Figure()
for j in range(len(f_x_SPA)):
    fig.add_trace(go.Scatter(x=f_x_SPA[j], y=f_y_SPA[j],
                        line=dict(color='firebrick', width=1),
                        name=f"{hlm[j]} SPA"))
    fig.add_trace(go.Scatter(x=f_x_strain[j], y=f_y_strain[j],
                        line=dict(color='royalblue', width=1),
                        name=f"{hlm[j]} Strain"))

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
fig.add_trace(go.Scatter(x=ce_asd_amplitude, y=ce_asd_frequency,
                         line=dict(color='orange', width=2),
                         name="CE Noise Curve"))
fig.add_trace(go.Scatter(x=ligo_o4_asd_amplitude, y=ligo_o4_asd_frequency,
                         line=dict(color='orchid', width=2),
                         name="aLIGO Noise Curve"))
fig
#plt.loglog(ligo_o4_asd_amplitude, ligo_o4_asd_frequency, label='LIGO O4', color='red')
#plt.draw()
#plt.show()
"""