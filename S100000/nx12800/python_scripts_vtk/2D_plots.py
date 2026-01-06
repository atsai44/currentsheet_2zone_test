import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import cmasher as cmr
from matplotlib.colors import TwoSlopeNorm
from celluloid import Camera

# --- Parameters to set ---

gamma = 1.66667 #1.66667
beta = 1 # 1, 1000
S = 100000

nx = 6400 #12800 #6400 #3200 # 6400
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
epsp = 0.2 * p0 #0.15 * p0


# set t_limit and dt using t_crossing
t_c = L / cms
#t_limit = 1000 * t_c
t_limit = 100 * t_c # increase to 100 or 1000
#t_limit = 10 * t_c
dt_vtk = 0.1 * t_c
dt_hst = 1 * t_c
dt_restart = 1 * t_c



# --- File naming ---
raw_name = f"gamma{gamma:.2f}_beta{beta}_S{S}_gpa{grid_cells_per_a:.2f}_nx{nx}_epsp{epsp:.4f}"
#raw_name = f"rst1_gamma{gamma:.2f}_beta{beta}_S{S}_gpa{grid_cells_per_a}_nx{nx}"
safe_name = raw_name.replace(".", "_")
jobname = safe_name
#jobname = "gamma1_67_beta1_S10000_gpa13_3_nx1600"
#jobname = "nx_1600_S10000_a0p05_g1p6667_b1"

nsteps = int(t_limit/dt_vtk)  #100   # how many timesteps
N_grid = nx

# loop over first 5 timesteps
#for idx in range(nsteps):
#for idx in np.arange(0,100,5):
#for idx in np.arange(0,10,5):
#    filename_jz = f"../vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
#    filename_w = f"../vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"
#    filename_bcc = f"../vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk"

#    data_w = pv.read(filename_w)
#    data_jz = pv.read(filename_jz)
#    data_bcc = pv.read(filename_bcc)

#    # extract density
#    rho = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid)) # test density cuts

#    jz = np.array(data_jz.cell_data["jz"]).reshape((N_grid, N_grid))
#    By  = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))
#    print('max rho', np.max(rho))
#    print('max jz',np.max(jz))
#    print('max By',np.max(By))

#    # specify what is plotted
#    quantity = jz
#    quantity_name = 'jz'

    # --- take slices ---
#    q_x1 = quantity   # z=0 slice


#    # --- global normalization (shared color scale) ---
#    vmin = -100 #quantity.min()
#    vmax = 100 #quantity.max()
#    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

#    # --- plot three slices side by side ---
#    fig, axes = plt.subplots(1, 1, figsize=(6, 6), dpi=150)

#    im0 = axes.imshow(q_x1.T, origin="lower", extent=[-3,3,-3,3],
#                         cmap=cmr.redshift, norm=norm)
#    axes.set_title("x1")
#    axes.set_xlabel("x"); axes.set_ylabel("y")


#    # --- single shared colorbar ---
#    cbar = fig.colorbar(im0, ax=axes, orientation='vertical',
#                        fraction=0.046, pad=0.04, label=quantity_name)

#    plt.suptitle(f"Timestep {idx}", fontsize=14)
#    outname = f"./plots_init/2D_{jobname}_{idx:05d}.png"
#    plt.savefig(outname, dpi=150)
#    plt.close()




# --- setup ---
# animation

fig, ax = plt.subplots(figsize=(6,6), dpi=150)
camera = Camera(fig)


vmin, vmax = -45, 45         # for jz
#vmin, vmax = -1.5, 1.5
norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

nsteps = int(t_limit/dt_vtk)
#for idx in range(nsteps):
for idx in np.arange(0,78,1):
    filename_jz = f"../vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
    filename_w = f"../vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"
    filename_bcc = f"../vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk"

    data_w = pv.read(filename_w)
    data_jz = pv.read(filename_jz)
    data_bcc = pv.read(filename_bcc)

    # extract density
    rho = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid)) # test density cuts

    jz = np.array(data_jz.cell_data["jz"]).reshape((N_grid, N_grid))
    By  = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))

    internalE = np.array(data_w.cell_data["eint"]).reshape((N_grid, N_grid))

    pressure = (gamma - 1.0) * internalE


    print('max rho', np.max(rho))
    print('max jz',np.max(jz))
    print('max By',np.max(By))
    print('min By',np.min(By))
    print('max pressure',np.max(pressure))

    # specify what is plotted
    quantity = jz
    quantity_name = 'jz'
    #quantity = By
    #quantity_name = 'By'
    #quantity = rho
    #quantity_name = 'rho'
    #quantity = pressure
    #quantity_name = 'pressure'

    #im = ax.imshow(jz.T, origin="lower", extent=[-3,3,-3,3],
    #               cmap=cmr.redshift, norm=norm)
    #im = ax.imshow(quantity.T, origin="lower", extent=[-3,3,-3,3],
    #               cmap=plt.cm.plasma, norm=norm)
    im = ax.imshow(quantity.T, origin="lower", extent=[-3,3,-1,1],
                   cmap=plt.cm.plasma, norm=norm)
    title = ax.text(
        0.5, 1.02,                      # position (x,y) in axes coords
        f"{quantity_name}, $t/t_c =$ {idx*dt_vtk/t_c:.2f}, $S$={S},  $\\gamma$={gamma},  $\\beta$={beta},  $n_x$={N_grid}",
        transform=ax.transAxes,
        ha="center", va="bottom",
        fontsize=10
    )
    ax.set_xlabel("x"); ax.set_ylabel("y")
    camera.snap()

animation = camera.animate()#(interval=200, blit=True)
#animation.save(f"./plots_init/{jobname}_{quantity_name}_animation.gif", writer='pillow')
animation.save(f"./plots_init/{jobname}_{quantity_name}_animation.mp4", writer='ffmpeg')
plt.close()