import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from matplotlib.widgets import Slider
from matplotlib.lines import Line2D


dfEXC = pd.read_csv("EXC_res_run-20250804T111207.csv")
dfGRO = pd.read_csv("GRO_res_run-20250804T111334.csv")


# Complete Binding Energy surface
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.plot_trisurf(
    dfEXC["Vc"],
    dfEXC["Vso"],
    dfEXC["energia"],
    cmap="viridis",
    edgecolor="none",
)
ax.plot_trisurf(
    dfGRO["Vc"],
    dfGRO["Vso"],
    dfGRO["energia"],
    cmap="copper",
    edgecolor="none",
)

custom_lines = [
    # approximated color to cmap viridis
    Line2D([0], [0], color="mediumseagreen", lw=4, label="Excited"),
    # approximated color to cmap copper
    Line2D([0], [0], color="peru", lw=4, label="Ground"),
]

ax.legend(handles=custom_lines)

ax.set_xlabel("Vc")
ax.set_ylabel("Vso")
ax.set_zlabel("Eb")
plt.title("Energy of bound state")
plt.show()


# Fit regression plane
eX = dfEXC[["Vc", "Vso"]]  # features
ez = dfEXC["energia"]  # target
model_e = LinearRegression()
model_e.fit(eX, ez)
print(model_e.score(eX, ez))
print("E coef:", model_e.coef_)
we = np.array(model_e.coef_)

gX = dfGRO[["Vc", "Vso"]]  # features
gz = dfGRO["energia"]  # target
model_g = LinearRegression()
model_g.fit(gX, gz)
print(model_g.score(gX, gz))
print("G coef:", model_g.coef_)
wg = np.array(model_g.coef_)

Hessian = (np.outer(we, we) + np.outer(wg, wg)) / 2
print(np.linalg.eig(Hessian))
"""
eigvals = [0.07855039, 0.00281386]
eigvect = 
    [[ 0.99786511, -0.06530865],
    [ 0.06530865,  0.99786511]]

upper limit = 2/lmb_max = 25.46136308171099
1/lmb_min = 355.3837077893001
"""

fig = plt.figure()
lr_ = np.linspace(0, 500)
plt.plot(lr_, 1 - lr_ * 0.07855039, label="max eigval")
plt.plot(lr_, 1 - lr_ * 0.00281386, label="min eigval")
plt.plot(lr_, [0] * len(lr_), "--")
plt.ylim((-1, 1))
plt.show()
print("L2 optim:", (0.07855039 + 0.00281386) / (0.07855039**2 + 0.00281386**2))

# Band from sampled surface
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

b = 1e-1
dfEXC_band = dfEXC[abs(dfEXC["energia"] + 1.990390) < b]
dfGRO_band = dfGRO[abs(dfGRO["energia"] + 2.467000) < b]

print(dfEXC_band.head())
print(dfGRO_band.head())

ax.plot(
    dfEXC_band["Vc"],
    dfEXC_band["Vso"],
    dfEXC_band["energia"],
    # cmap="viridis",
    # edgecolor="none",
)
ax.plot(
    dfGRO_band["Vc"],
    dfGRO_band["Vso"],
    dfGRO_band["energia"],
    # cmap="copper",
    # edgecolor="none",
)
ax.set_xlabel("Vc")
ax.set_ylabel("Vso")
ax.set_zlabel("Eb")
plt.show()


# Loss functions
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

c = 1
LossE = np.log10(c + (dfEXC["energia"] + 1.990390) ** 2)
LossG = np.log10(c + (dfGRO["energia"] + 2.467000) ** 2)
# (dfEXC["energia"] + 1.990390) ** 2 + (dfGRO["energia"] + 2.467000) ** 2

ax.plot_trisurf(dfGRO["Vc"], dfGRO["Vso"], LossE, cmap="viridis", edgecolor="none")
ax.plot_trisurf(dfGRO["Vc"], dfGRO["Vso"], LossG, cmap="copper", edgecolor="none")

ax.set_xlabel("Vc")
ax.set_ylabel("Vso")
ax.set_zlabel("Loss")
plt.title("Loss functions")
plt.show()


# Loss function combined in homotopy
fig = plt.figure()
plt.subplots_adjust(bottom=0.25)
ax = fig.add_subplot(111, projection="3d")

a = 1 / 2
Loss = a * LossE + (1 - a) * LossG

surf = ax.plot_trisurf(
    dfGRO["Vc"], dfGRO["Vso"], Loss, cmap="gist_stern", edgecolor="none"
)

ax_slider = plt.axes([0.2, 0.1, 0.65, 0.03])
slider = Slider(ax_slider, label="a", valmin=0.0, valmax=1.0, valinit=a)


def update(val):
    global surf
    a = slider.val

    surf.remove()

    Loss_new = a * LossE + (1 - a) * LossG
    surf = ax.plot_trisurf(dfGRO["Vc"], dfGRO["Vso"], Loss_new, cmap="gist_stern")
    fig.canvas.draw_idle()


ax.set_xlabel("Vc")
ax.set_ylabel("Vso")
ax.set_zlabel("Loss")
plt.title("Loss function")
slider.on_changed(update)
plt.show()
