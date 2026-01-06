import numpy as np
import matplotlib.pyplot as plt
import cmasher as cmr
from matplotlib.colors import TwoSlopeNorm
from celluloid import Camera
from read_binaries import*


# --- Parameters to set ---

gamma = 1.66667
beta = 1 # 0.001, 1, 1000
S = 100000 #100000

nx = 12800 #12800 #6400 #6400 #3200 # 6400
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
epsp = 0.2 * p0 #0.2 * p0 #0.5 * p0 #0.15 * p0


# set t_limit and dt using t_crossing
t_c = L / cms
#t_limit = 1000 * t_c
t_limit = 100 * t_c # increase to 100 or 1000
#t_limit = 10 * t_c
dt_vtk = 0.1 * t_c
dt_hst = 1 * t_c
dt_restart = 0.5 * t_c



# --- File naming ---
raw_name = f"gamma{gamma:.2f}_beta{beta}_S{S}_gpa{grid_cells_per_a:.2f}_nx{nx}_epsp{epsp:.4f}"
#raw_name = f"rst1_gamma{gamma:.2f}_beta{beta}_S{S}_gpa{grid_cells_per_a}_nx{nx}"
safe_name = raw_name.replace(".", "_")
jobname = safe_name
#jobname = "gamma1_67_beta1_S10000_gpa13_3_nx1600"
#jobname = "nx_1600_S10000_a0p05_g1p6667_b1"

nsteps = int(t_limit/dt_vtk)  #100   # how many timesteps
N_grid = nx



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
    data = ReadAthenaBinData(file_number=idx, basename=safe_name, data_dir="../bin")
    
    data.read("mag")
    data.read("vel")
    data.read("dens")
    data.read("press")
    data.read("cur")

    # Extract variables
    rho = data.dens[0, :, :, 0]
    jz  = data.cur[0, :, :, 0]
    pressure = data.press[0, :, :, 0]
    By = data.mag[1, :, :, 0]
    #internalE = 

    #pressure = (gamma - 1.0) * internalE


    print('max rho', np.max(rho))
    print('max jz',np.max(jz))
    #print('max By',np.max(By))
    #print('min By',np.min(By))
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
    im = ax.imshow(quantity, origin="lower", extent=[-3,3,-1,1],
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