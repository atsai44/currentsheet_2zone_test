import numpy as np
from matplotlib import pyplot as plt
import pyvista as pv
import os

if __name__ == '__main__':

    # parameter lists
    res_list   = [1600]                    # grid resolutions
    S_list   = [10000]                  # resistivities
    a0_list    = [0.05]                    # sheet half-thickness
    gamma_list=[1, 1.5, 2]
    beta_list=[0.001, 1, 1000]     # plasma betas    

    # fixed domain size
    L = 6.0

    # plotting params
    time_steps = 5
    dpi = 150

    # make sure output folder exists
    os.makedirs("./plots_init", exist_ok=True)

    # loop over parameter combinations
    for res in res_list:
        N_grid = res

        for S in S_list:
            for a0 in a0_list:
                for gamma_val in gamma_list:
                    for beta in beta_list:

                        # ---- derived quantities ----
                        alpha = np.sqrt(2.0 / beta)
                        b0 = alpha / np.sqrt(gamma_val)
                        rho_0 = (1.0 / gamma_val) ** (1.0 / gamma_val)
                        v_a = np.sqrt(b0**2 / rho_0)
                        t_char = a0 / v_a
                        t_array = np.arange(0, time_steps) * 0.1
                        t_physical = t_array / t_char

                        # grid in x (for 1D cut)
                        x = np.linspace(-L/2, L/2, N_grid)

                        # labels for filenames
                        S_label   = str(S).replace(".", "p")
                        a0_label    = str(a0).replace(".", "p")


                        gamma_label = f"g{str(gamma_val).replace('.', 'p')}"


                        beta_label = f"b{str(beta).replace('.', 'p')}"

                        jobname = f"nx_{res}_S{S_label}_a{a0_label}_{gamma_label}_{beta_label}"
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

                            # extract vars
                            density   = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid))
                            internalE = np.array(data_w.cell_data["eint"]).reshape((N_grid, N_grid))
                            jz        = np.array(data_jz.cell_data["jz"]).reshape((N_grid, N_grid))
                            vx        = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))
                            B         = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))

                            va        = B/np.sqrt(density)

                            # derived pressure
                            pressure = (gamma_val - 1.0) * internalE

                            # take 1D cut
                            cut_idx = N_grid // 2
                            rho_line = density[cut_idx]
                            P_line   = pressure[cut_idx]
                            jz_line  = jz[cut_idx]
                            vx_line = vx[cut_idx]
                            va_line = va[cut_idx]

                            # normalization
                            if idx == 0:
                                rho_max_t0 = np.max(rho_line)
                                P_max_t0   = np.max(P_line)
                                jz_max_t0  = np.max(jz_line)


                            rho_line /= 1 #rho_max_t0
                            P_line   /= 1 #P_max_t0 : p is 0 when gamma = 1
                            jz_line  /= 1 #jz_max_t0

                            # restrict range
                            x_lim = 0.05
                        
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
                            ax_va.plot((x/L)[idx_range], va_line[idx_range],
                                      label=f"$t/t_c$={t_physical[idx]:.1f}")
                            ax_vx.plot((x/L), vx_line,
                                      label=f"$t/t_c$={t_physical[idx]:.1f}")
                            #ax_va.plot((x/L), va_line,
                            #          label=f"$t/t_c$={t_physical[idx]:.1f}")

                        # add titles & labels
                        ax_jz.set_title("Current density $j_z$")
                        ax_rho.set_title("Density $\\rho$")
                        ax_P.set_title("Pressure $P$")
                        ax_vx.set_title("velocity $v_x$")
                        ax_va.set_title("Alven velocity $v_a = B/\sqrt{\\rho}$")

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
