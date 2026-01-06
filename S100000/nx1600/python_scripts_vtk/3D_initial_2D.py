import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import cmasher as cmr
from matplotlib.colors import TwoSlopeNorm

# --- parameters ---
N_grid = 200
jobname = "CS3D"
nsteps = 5   # how many timesteps

# loop over first 5 timesteps
for idx in range(nsteps):
    filename_jz = f"../vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
    filename_w = f"../vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"
    filename_bcc = f"../vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk"

    data_w = pv.read(filename_w)
    data_jz = pv.read(filename_jz)
    data_bcc = pv.read(filename_bcc)

    # extract density
    rho = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid, N_grid)) # test density cuts

    jz = np.array(data_jz.cell_data["jz"]).reshape((N_grid, N_grid, N_grid))
    By  = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid, N_grid))
    print('max rho', np.max(rho))
    print('max jz',np.max(jz))
    print('max By',np.max(By))

    print('min rho', np.min(rho))
    print('min jz',np.min(jz))
    print('min By',np.min(By))

    # specify what is plotted
    quantity = jz
    quantity_name = 'jz'

    # --- take slices ---
    q_xy = quantity[:, :, N_grid//2]   # z=0 slice
    q_yz = quantity[N_grid//2, :, :]   # x=0 slice
    q_xz = quantity[:, N_grid//2, :]   # y=0 slice

    # --- global normalization (shared color scale) ---
    vmin = -1 #quantity.min()
    vmax = 1 #quantity.max()
    norm = TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)

    # --- plot three slices side by side ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=150)

    im0 = axes[0].imshow(q_xy.T, origin="lower", extent=[-3,3,-3,3],
                         cmap=cmr.redshift, norm=norm)
    axes[0].set_title("XY plane (z=0)")
    axes[0].set_xlabel("x"); axes[0].set_ylabel("y")

    im1 = axes[1].imshow(q_yz.T, origin="lower", extent=[-3,3,-3,3],
                         cmap=cmr.redshift, norm=norm)
    axes[1].set_title("YZ plane (x=0)")
    axes[1].set_xlabel("y"); axes[1].set_ylabel("z")

    im2 = axes[2].imshow(q_xz.T, origin="lower", extent=[-3,3,-3,3],
                         cmap=cmr.redshift, norm=norm)
    axes[2].set_title("XZ plane (y=0)")
    axes[2].set_xlabel("x"); axes[2].set_ylabel("z")

    # --- single shared colorbar ---
    cbar = fig.colorbar(im2, ax=axes, orientation='vertical',
                        fraction=0.046, pad=0.04, label=quantity_name)

    plt.suptitle(f"Timestep {idx}", fontsize=14)
    outname = f"./plots_init/box_slices_{idx:05d}.png"
    plt.savefig(outname, dpi=150)
    plt.close()

