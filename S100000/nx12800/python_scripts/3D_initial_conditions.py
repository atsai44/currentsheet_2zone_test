import numpy as np
from matplotlib import pyplot as plt
import pyvista as pv
import os

if __name__ == '__main__':

    # parameter lists
    res_list   = [200]        # grid resolutions in x/y
    S_list     = [10000]      # resistivities
    a0_list    = [0.05]       # sheet half-thickness
    gamma_list = [1.6667]
    beta_list  = [1/5]        # plasma betas    

    # fixed domain size
    L = 6.0

    # plotting params
    time_steps = 5
    dpi = 150

    # make sure output folder exists
    os.makedirs("./plots_init", exist_ok=True)

    # loop over parameter combinations
    for res in res_list:
        N_grid = res  # only used for x-axis plotting

        for S in S_list:
            for a0 in a0_list:
                for gamma_val in gamma_list:
                    for beta in beta_list:

                        # ---- derived quantities ----
                        alpha  = np.sqrt(2.0 / beta)
                        b0     = alpha / np.sqrt(gamma_val)
                        rho_0  = (1.0 / gamma_val) ** (1.0 / gamma_val)
                        v_a    = np.sqrt(b0**2 / rho_0)
                        t_char = a0 / v_a
                        t_array     = np.arange(0, time_steps) * 0.1
                        t_physical  = t_array / t_char

                        # grid in x (for 1D cut)
                        x = np.linspace(-L/2, L/2, N_grid)

                        jobname = "CS3D"
                        print(f"Processing {jobname}...")

                        # initialize normalization
                        rho_max_t0, P_max_t0, jz_max_t0 = 0, 0, 0

                        # set up figure with 5 subplots
                        fig, axes = plt.subplots(1, 5, figsize=(20, 4), dpi=dpi)
                        ax_jz, ax_rho, ax_P, ax_vx, ax_va = axes

                        for idx in range(time_steps):

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

                            # get 3D cell dimensions
                            nx, ny, nz = data_w.dimensions
                            nx, ny, nz = nx - 1, ny - 1, nz - 1 
                            # PyVista uses Fortran ordering for cell_data
                            order = 'F'

                            # choose mid-plane in z (index nz//2)
                            k_mid = nz // 2

                            def extract_midplane(var):
                                arr3d = np.array(var, order='F').reshape((nx, ny, nz), order='F')
                                return arr3d[:, :, k_mid]  # slice at constant z

                            density   = extract_midplane(data_w.cell_data["dens"])
                            internalE = extract_midplane(data_w.cell_data["eint"])
                            jz        = extract_midplane(data_jz.cell_data["jz"])
                            vx        = extract_midplane(data_w.cell_data["velx"])
                            B         = extract_midplane(data_bcc.cell_data["bcc2"])

                            va        = B/np.sqrt(density)

                            # derived pressure
                            pressure = (gamma_val - 1.0) * internalE

                            # take 1D cut along x through y mid-plane
                            j_mid = ny // 2
                            rho_line = density[:, j_mid]
                            P_line   = pressure[:, j_mid]
                            jz_line  = jz[:, j_mid]
                            vx_line  = vx[:, j_mid]
                            va_line  = va[:, j_mid]

                            # normalization at t=0
                            if idx == 0:
                                rho_max_t0 = np.max(rho_line)
                                P_max_t0   = np.max(P_line)
                                jz_max_t0  = np.max(jz_line)

                            # (Keep raw values or normalize if desired)
                            rho_line /= 1  # e.g. rho_max_t0
                            P_line   /= 1  # e.g. P_max_t0
                            jz_line  /= 1  # e.g. jz_max_t0

                            # restrict range for zoom
                            x_lim = 0.05
                            mask = (x/L >= -x_lim) & (x/L <= x_lim)
                            idx_range = np.where(mask)[0]

                            # plot lines
                            label = f"$t/t_c$={t_physical[idx]:.1f}"
                            ax_jz.plot((x/L)[idx_range], jz_line[idx_range], label=label)
                            ax_rho.plot((x/L)[idx_range], rho_line[idx_range], label=label)
                            ax_P.plot((x/L)[idx_range], P_line[idx_range], label=label)
                            ax_va.plot((x/L)[idx_range], va_line[idx_range], label=label)
                            ax_vx.plot((x/L), vx_line, label=label)

                        # add titles & labels
                        ax_jz.set_title("Current density $j_z$")
                        ax_rho.set_title("Density $\\rho$")
                        ax_P.set_title("Pressure $P$")
                        ax_vx.set_title("Velocity $v_x$")
                        ax_va.set_title("Alfvén speed $v_a$")

                        for ax in axes:
                            ax.set_xlabel("x/L")
                            ax.legend(fontsize=6, ncol=2)

                        fig.suptitle(
                            f"Equilibrium (S = {S:.0f}, Nx = {N_grid}, $\\gamma$ = {gamma_val}, $\\beta$ = {beta})",
                            fontsize=14
                        )
                        fig.tight_layout()
                        outname = f"./plots_init/{jobname}_jz_rho_P.png"
                        plt.savefig(outname, dpi=dpi)
                        plt.close(fig)