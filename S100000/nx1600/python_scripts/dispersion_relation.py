import numpy as np
from matplotlib import pyplot as plt
import pyvista as pv
import os
import matplotlib.colors as colors



# --- Parameters to set ---

gamma = 1.5 #1.66667 #1.66667
beta = 1000 #0.001, 1, 1000
S = 10000

nx = 1600 #6400 #3200 # 6400
ny = nx
nz = 1  # 2D simulation

#grid_cells_per_a = 13.3  # a_values * nx_values / L
L = 6.0
Lz = 1.0

# --- Derived quantities ---
#a0 = grid_cells_per_a * L / nx
a0 = 0.05
grid_cells_per_a = nx * a0 / L

p0 = 1 / gamma
rho0 = p0 ** (1 / gamma)
alpha = np.sqrt(2 / beta)
b0 = alpha / np.sqrt(gamma)
va = np.sqrt(b0 ** 2 / rho0)
cs = np.sqrt(gamma * p0/rho0)
cms = np.sqrt(va**2 + cs**2)
eta = L * va / S
epsp = 0.0 #0.3 * p0 #0.15 * p0


# set t_limit and dt using t_crossing
t_c = L / cms
#t_limit = 1000 * t_c
t_limit = 100 * t_c # increase to 100 or 1000
#t_limit = 10 * t_c
dt_vtk = 0.1 * t_c
dt_hst = 1 * t_c
dt_restart = 1 * t_c


raw_name = f"gamma{gamma:.2f}_beta{beta}_S{S}_gpa{grid_cells_per_a:.2f}_nx{nx}"
safe_name = raw_name.replace(".", "_")
jobname = safe_name



# params
time_steps = int(t_limit/dt_vtk)
N_grid = nx   
dpi = 150
cut_idx = 100 


# convert to physical units
t_array     = np.arange(0, time_steps) * dt_vtk


t_physical  = t_array / t_c
delta_t_physical = t_physical[1] - t_physical[0]
    #print('dt:', delta_t_physical)

#cs = np.sqrt(gamma*(1/gamma)**((gamma - 1)/gamma))

# grid in x (for 1D cut)
x = np.linspace(-L/2, L/2, N_grid)
x_physical = x/L
delta_x = x_physical[1] - x_physical[0]


# make 2D array for v(x,t)
v_xt = np.zeros((len(t_physical), len(x)))


# make sure output folder exists
os.makedirs("./plots_init", exist_ok=True)

for idx in range(time_steps):

    # filepaths
    #filename_bcc = f"../vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk"
    #filename_jz  = f"../vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
    filename_w   = f"../vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"


    # load data
    #data_jz  = pv.read(filename_jz)
    data_w   = pv.read(filename_w)
    #data_bcc = pv.read(filename_bcc)

    # extract vars
    vx        = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))

    vx_line = vx[cut_idx]
    v_xt[idx] = vx_line


# compute fft
v_tilde = np.fft.fftshift(np.fft.fft2(v_xt,norm="forward"))
v_tilde_sqr = np.abs(v_tilde)**2

# get frequency and wavenumber axes
omega = np.fft.fftshift(np.fft.fftfreq(len(t_physical), d = 1./delta_t_physical))  # normalized frequency
k = np.fft.fftshift(np.fft.fftfreq(len(x/L), d = 1./delta_x))              # normalized wavenumber

# select only positive quadrant
# get indices for positive quadrant
omega_idx = np.where(omega >= 0)[0]
k_idx = np.where(k >= 0)[0]

# slice arrays
v_tilde_sqr = v_tilde_sqr[np.ix_(omega_idx, k_idx)]
omega = omega[omega_idx]
k = k[k_idx]


# fastest ms wave
cms = np.sqrt(cs**2 + va**2)*t_c/L
omega_fast = cms*k
omega_a = va*t_c/L   *k
omega_s = cs*t_c/L   *k
print('v_a', va*t_c/L)
print('cms', cms)


# normalize using integral over the plane
# spacing in k and omega
delta_k = np.abs(k[1] - k[0])
delta_omega = np.abs(omega[1] - omega[0])

# normalize by total integral over k, omega plane
total_power = np.sum(v_tilde_sqr) * delta_k * delta_omega
v_tilde_sqr_norm = v_tilde_sqr / total_power


# plot power spectrum

fig, ax = plt.subplots(figsize=(8,6))
#vmin = np.min(v_tilde_sqr_norm[v_tilde_sqr_norm > 0])  # smallest positive entry
#vmax = np.max(v_tilde_sqr_norm)
#vmin = 10**(-7)
#vmin = 10**(-2) #10**(-7)
#vmax = 10**(2)    #10**(2)


# Plot full power spectrum
pcm = ax.pcolormesh(
    k, omega, v_tilde_sqr_norm,
    shading='auto',
    cmap='inferno',
    norm=colors.LogNorm()#vmin=vmin, vmax=vmax)
)

#pcm = ax.pcolormesh(
#    k[:20], omega, v_tilde_sqr[:,:20],
#    shading='auto',
#    cmap='inferno',
#)

# overlay line
#cutoff_omega_fast = np.where(omega_fast <= np.max(omega))
#cutoff_omega_a = np.where(omega_a <= np.max(omega))
#cutoff_omega_s = np.where(omega_s <= np.max(omega))
      
#ax.plot(k[cutoff_omega_a], omega_a[cutoff_omega_a], color='blue', lw=0.6, label="alfven wave", linestyle = '-')

#ax.plot(k[cutoff_omega_s], omega_s[cutoff_omega_s], color='white', lw=0.6, label="sound speed", linestyle = '-')

#ax.plot(k[cutoff_omega_fast], omega_fast[cutoff_omega_fast], color='cyan', lw=0.6, label="fast magnetosonic wave", linestyle = '-')

# extend curves to full extent of plotted k range


k_plot = np.linspace(0, np.max(k), 100)

omega_a_plot = va * t_c / L * k_plot
omega_s_plot = cs * t_c / L * k_plot
omega_fast_plot = cms * k_plot

ax.plot(k_plot, omega_a_plot, color='blue', lw=0.6, label="Alfvén wave", linestyle='-')
ax.plot(k_plot, omega_s_plot, color='white', lw=0.6, label="Sound speed", linestyle='-')
ax.plot(k_plot, omega_fast_plot, color='cyan', lw=0.6, label="Fast magnetosonic wave", linestyle='-')

ax.set_ylim(min(omega), max(omega))


ax.set_xlabel("$k L$")
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylabel("$\omega t_c$")
ax.set_title(f"$\\beta = $ {beta}")

#ax.set_xlim(0, 20)         # only show k up to 20
ax.set_xlim(0, 0.00005)

ax.legend(loc='lower right')
fig.colorbar(pcm, ax=ax, label=r"$|\tilde{v}(k,\omega)|^2$")


outname = f"./dispersion_plots/{jobname}_dispersion.png"
plt.savefig(outname, dpi=dpi)
plt.close(fig)


