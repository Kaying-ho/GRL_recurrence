#%% Theoretical Distribution of recurrence with varying φ/α (α denoted as Pk here in code)
'''
Figure 4 (b-c): Theoretical Distribution of recurrence with varying φ/α
This script plots theoretical distribution of recurrence with varying φ (at fixed α = 0.05) and varying α (at fixed φ = 0.8)
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import os
basepath = os.path.expanduser("~/Github")
savingpath = f"{basepath}/plots/Fig4"

x = np.arange(1, 100)

# First plot: varying phi, fixed α
phis = np.arange(0.1, 1.0, 0.1)
Pk_fixed = 0.05

# Second plot: fixed phi, varying α
phi_fixed = 0.8
Pks = np.arange(0.01, 0.1, 0.01)

# Create figure with 2 subplots (independent y-axes)
fig, axs = plt.subplots(1, 2, figsize=(18, 4.5))

# -- Left plot: Varying phi --
cmap_phi = get_cmap('Blues')
norm_phi = Normalize(vmin=min(phis), vmax=max(phis))

axs[0].plot([], [],  color='w', label=f"α = {Pk_fixed}")

for phi in phis:
    P = (1 - phi**x) * Pk_fixed
    y = P * ((1 - P) ** (x - 1))
    y /= y.sum()
    max_x = np.where(y == np.max(y))[0][0]
    mean_y = np.sum(x * y)
    sd_y = np.sqrt(np.sum(y * (x - mean_y) ** 2))
    color = cmap_phi(norm_phi(phi))
    Label_text = fr"$\phi={phi:.1f}$"
    axs[0].plot(x, y, linewidth=2, color=color, label=f"{Label_text}")

axs[0].set_xlabel("Return Period (day)", fontsize=12)
axs[0].set_yticks(np.arange(0, 0.051, 0.01))
axs[0].set_ylabel("Probability Density", fontsize=12)
axs[0].set_title(f"Return Period Distribution with Varying φ", fontsize=14)
axs[0].legend(fontsize=8)
axs[0].grid(True, linestyle='--', alpha=0.5)

# -- Right plot: Varying α --
cmap_pk = get_cmap('Reds')
norm_pk = Normalize(vmin=min(Pks), vmax=max(Pks))

axs[1].plot([], [],  color='w', label=f"φ = {phi_fixed}")

for Pk in Pks:
    P = (1 - phi_fixed**x) * Pk
    y = P * ((1 - P) ** (x - 1))
    y /= y.sum()
    max_x = np.where(y == np.max(y))[0][0]
    mean_y = np.sum(x * y)
    sd_y = np.sqrt(np.sum(y * (x - mean_y) ** 2))
    color = cmap_pk(norm_pk(Pk))
    Label_text = fr"$\alpha={Pk:.2f}$"
    axs[1].plot(x, y, linewidth=2, color=color, label=f"{Label_text}")

axs[1].set_xlabel("Return Period (day)", fontsize=12)
axs[1].set_yticks(np.arange(0, 0.061, 0.01))
axs[1].set_ylabel("Probability Density", fontsize=12)
axs[1].set_title(f"Return Period Distribution with Varying α", fontsize=14)
axs[1].legend(fontsize=8)
axs[1].grid(True, linestyle='--', alpha=0.5)

plt.tight_layout(pad=2)
fig.subplots_adjust(wspace=0.15) 
plt.savefig(f"{savingpath}/Distribution.png", dpi = 600)
