import numpy as np
from matplotlib import pyplot as plt
import pyvista as pv
import os
import math
from celluloid import Camera

########## convert this file to read bin files ###############

# -------------------------
# Simulation parameters
# -------------------------
# --- Parameters to set ---

gamma = 1.5 #1.66667 #1.66667
beta = 1 # 0.001, 1, 1000
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

time_steps = int(t_limit/dt_vtk)
t_physical = np.arange(time_steps) * dt_vtk / t_c
N_grid = nx


print(f"Derived params: a0={a0:.5e}, p0={p0:.5e}, rho0={rho0:.5e}, b0={b0:.5e}, va={va:.5e}")
print(f"t_c={t_c:.5e}, t_limit={t_limit:.5e}, dt_vtk={dt_vtk:.5e}")

# -------------------------
# Filenames / jobname
# -------------------------
raw_name = f"gamma{gamma:.2f}_beta{beta}_S{S}_gpa{grid_cells_per_a:.2f}_nx{nx}"
safe_name = raw_name.replace(".", "_")   # matches your template naming
print(safe_name)
#safe_name = 'CS'
vtk_dir = "../vtk"                       # change if needed
jobname = safe_name

# grid in x (for 1D cut)
x = np.linspace(-L/2, L/2, N_grid)


# initialize normalization
rho_max_t0, P_max_t0, jz_max_t0 = 0, 0, 0

# set up figure with 5 subplots
dpi = 150
fig, axes = plt.subplots(1, 5, figsize=(20, 4), dpi=dpi)
ax_jz, ax_rho, ax_P, ax_vx, ax_va = axes



for idx in np.arange(0,time_steps,100):
#for idx in np.arange(0,49,1):

    # filepaths
    filename_bcc = f"../vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk"
    filename_jz  = f"../vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
    filename_w   = f"../vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"

    if not (os.path.exists(filename_bcc) and
        os.path.exists(filename_jz) and
        os.path.exists(filename_w)):
        print(f"  Skipping timestep {idx} (files not found)")
        continue

    # load data
    data_jz  = pv.read(filename_jz)
    data_w   = pv.read(filename_w)
    data_bcc = pv.read(filename_bcc)

    # extract vars
    density   = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid))
    internalE = np.array(data_w.cell_data["eint"]).reshape((N_grid, N_grid))
    jz        = np.array(data_jz.cell_data["jz"]).reshape((N_grid, N_grid))
    vx        = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))
    B         = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))

    va_t        = B/np.sqrt(density)

    # derived pressure
    pressure = (gamma - 1.0) * internalE

    # take 1D cut
    cut_idx = N_grid // 2
    rho_line = density[cut_idx]
    P_line   = pressure[cut_idx]
    jz_line  = jz[cut_idx]
    vx_line = vx[cut_idx]
    va_t_line = va_t[cut_idx]

    # normalization
    if idx == 0:
        rho_max_t0 = np.max(rho_line)
        P_max_t0   = np.max(P_line)
        jz_max_t0  = np.max(jz_line)


    rho_line /= 1 #rho_max_t0
    P_line   /= 1 #P_max_t0 : p is 0 when gamma = 1
    jz_line  /= 1 #jz_max_t0

    # restrict range
    x_lim = 0.5 #0.05
                        
    mask = (x/L >= -x_lim) & (x/L <= x_lim)
    idx_range = np.where(mask)[0]

    # plot lines
    ax_jz.plot((x/L)[idx_range], jz_line[idx_range],
                                       label=f"$t/t_c$={t_physical[idx]:.1f}")
    ax_rho.plot((x/L)[idx_range], rho_line[idx_range],
                                        label=f"$t/t_c$={t_physical[idx]:.1f}")
    ax_P.plot((x/L)[idx_range], P_line[idx_range],
                                      label=f"$t/t_c$={t_physical[idx]:.1f}")
    #ax_vx.plot((x/L)[idx_range], vx_line[idx_range],
    #          label=f"$t/t_c$={t_physical[idx]:.1f}")
    ax_va.plot((x/L)[idx_range], va_t_line[idx_range],
                                      label=f"$t/t_c$={t_physical[idx]:.1f}")
    ax_vx.plot((x/L), vx_line,
                                      label=f"$t/t_c$={t_physical[idx]:.1f}")
    #ax_va.plot((x/L), va_line,
                            #          label=f"$t/t_c$={t_physical[idx]:.1f}")

# add titles & labels
ax_jz.set_title("$j_z$")
ax_rho.set_title("$\\rho$")
ax_P.set_title("$P$")
ax_vx.set_title("$v_x$")
ax_va.set_title("$v_a = B/\sqrt{\\rho}$")

for ax in axes:
    ax.set_xlabel("x/L")
    ax.legend(fontsize=6, ncol=2)


fig.suptitle(
    f"Equilibrium (S = {S:.0f}, Nx = {N_grid}, $\\gamma$ = {gamma}, $\\beta$ = {beta})",
    fontsize=14
)

#fig.tight_layout()
outname = f"./plots_init/{jobname}_jz_rho_P.png"
plt.savefig(outname, dpi=dpi)
plt.close(fig)






#############
# animation
#############

fig_anim, axes_anim = plt.subplots(1, 5, figsize=(20, 4), dpi=dpi)
ax_jz, ax_rho, ax_P, ax_vx, ax_va = axes_anim
camera = Camera(fig_anim)

plt.subplots_adjust(top=0.75)

#for idx in np.arange(0,time_steps,100):
for idx in range(time_steps):

    filename_bcc = f"../vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk"
    filename_jz  = f"../vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
    filename_w   = f"../vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"

    if not (os.path.exists(filename_bcc) and
            os.path.exists(filename_jz) and
            os.path.exists(filename_w)):
        print(f"Skipping timestep {idx}")
        continue

    # --- Load data ---
    data_jz  = pv.read(filename_jz)
    data_w   = pv.read(filename_w)
    data_bcc = pv.read(filename_bcc)

    # --- Extract variables ---
    density   = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid))
    internalE = np.array(data_w.cell_data["eint"]).reshape((N_grid, N_grid))
    jz        = np.array(data_jz.cell_data["jz"]).reshape((N_grid, N_grid))
    vx        = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))
    B         = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))
    va_t      = B / np.sqrt(density)
    pressure  = (gamma - 1.0) * internalE

    # --- 1D cut ---
    cut_idx = N_grid // 2
    rho_line = density[cut_idx]
    P_line   = pressure[cut_idx]
    jz_line  = jz[cut_idx]
    vx_line  = vx[cut_idx]
    va_t_line = va_t[cut_idx]

    # --- Range restriction ---
    x_lim = 0.5 #0.05
    mask = (x/L >= -x_lim) & (x/L <= x_lim)
    idx_range = np.where(mask)[0]

    # --- Plot for animation ---
    ax_jz.plot((x/L)[idx_range], jz_line[idx_range])
    ax_rho.plot((x/L)[idx_range], rho_line[idx_range])
    ax_P.plot((x/L)[idx_range], P_line[idx_range])
    ax_vx.plot((x/L)[idx_range], vx_line[idx_range])
    ax_va.plot((x/L)[idx_range], va_t_line[idx_range])

    # add titles & labels
    ax_jz.set_title("$j_z$")
    ax_rho.set_title("$\\rho$")
    ax_P.set_title("$P$")
    ax_vx.set_title("$v_x$")
    ax_va.set_title("$v_a = B/\sqrt{\\rho}$")


    # --- Static figure title (does not change) ---
    fig_anim.suptitle(
        f"$S$={S},  $\\gamma$={gamma},  $\\beta$={beta},  $n_x$={N_grid}",
        fontsize=14,
        y=0.90
    )


    # --- Add time label to one subplot (e.g. ax_jz) ---
    t_label = ax_jz.text(
        0.5, 1.10,
        f"$t/t_c =$ {idx*dt_vtk/t_c:.2f}",
        ha="center", va="bottom",
        transform=ax_jz.transAxes,
        fontsize=12
    )

    # --- Capture frame ---
    camera.snap()
    #fig_time.remove()   



# save animation 
animation = camera.animate()#(interval=200, blit=True)
os.makedirs("./plots_init", exist_ok=True)
animation.save(f"./plots_init/{jobname}_1Danimation.mp4", writer='ffmpeg')
plt.close(fig_anim)