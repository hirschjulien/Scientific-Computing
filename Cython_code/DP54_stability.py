# %%
from DP54 import DP54_solver
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import time
from matplotlib.gridspec import GridSpec
import numba
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy.ma as ma

# %%
plt.style.use('seaborn-v0_8-whitegrid')  
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 13
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['legend.fontsize'] = 13

# %%
real_range = np.linspace(-5.0, 2.0, 1000)  # Re range
imag_range = np.linspace(-5.0, 5.0, 1000)  # Im range
real_grid, imag_grid = np.meshgrid(real_range, imag_range)

amplification_factors = np.zeros_like(real_grid)

t_range = np.array((0.0, 1e0 ))  # range of time
h_size = np.array((1e0, 1e0))  # the maximun and minimun step(initial, the algorithm will refresh)
x0 = np.array((1.0, 0.0, 0.0)) ### initial value, if the variable is less than three, just set the later as 0, make sure the array is have 3 variable.
reps = 1e-2
aeps = 1e-3
model = 7  # model 7 is stability test
MAX_STEPS = int(50000)  ### most of iteration are 1 step, so MAX_STEP should be 1+1 = 2, but there have a bug, only if MAX_STEPS is a large number can run successfuly

for i in range(len(real_range)):
    for j in range(len(imag_range)):
        argv = np.array((real_range[i], imag_range[j]))      
        # use DP54 solver
        t_history, x_history, t_adaptive_size, t_adaptive = DP54_solver(
            x0, argv, t_range, h_size, reps, aeps, model, MAX_STEPS) 
        # |y_new|/|y_initial|
        y_new_magnitude = np.sqrt(x_history[-1, 0]**2 + x_history[-1, 1]**2)
        y_initial_magnitude = np.sqrt(x0[0]**2 + x0[1]**2)  
        amplification_factors[j, i] = y_new_magnitude / y_initial_magnitude

# %%
fig, pic1 = plt.subplots(figsize=(8, 6))
plt.style.use('seaborn-v0_8-whitegrid') 
clipped_factors = np.minimum(amplification_factors, 2.0)
mask = amplification_factors > 1.0  
masked_data = ma.masked_array(amplification_factors, mask=mask)
mesh = plt.pcolormesh(real_grid, imag_grid, masked_data, cmap="bone_r", vmin=0, vmax=1)

pic1.axhline(y=0, color='k', linestyle='-', alpha=0.3)
pic1.axvline(x=0, color='k', linestyle='-', alpha=0.3)
pic1.set_title('Stability plot for DP54 Method', fontsize=20, fontweight='bold', color='#333333', y=0.98)
pic1.set_xlabel('Re(z)', fontsize=12)
pic1.set_ylabel('Im(z)', fontsize=12)
pic1.grid(True, alpha=0.3)
pic1.set_facecolor('#f9f9f9')
cbar = plt.colorbar(mesh, label='R(z)')
plt.tight_layout()
plt.savefig("DP54_stability.png")
plt.show()


