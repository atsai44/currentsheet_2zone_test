import numpy as np
from matplotlib import pyplot as plt
import pyvista as pv
import cmasher as cmr
from matplotlib.colors import TwoSlopeNorm
from PLASMAtools.aux_funcs.derived_var_funcs import DerivedVars as DV
from celluloid import Camera # for animating plots

if __name__ == '__main__':

    jobname = 'gamma1_67_beta1_S10000_gpa13_3_nx1600'

    dpi = 150 #300 #150
    N_grid = 1600 #1024
    
    dv = DV(bcs="00",
            num_of_dims=2,
            L=[6,6])


    # set up figure
    f, ax = plt.subplots(1, 1, figsize=(N_grid/dpi, N_grid/dpi), dpi=dpi)
    camera = Camera(f)

    # Initialize lists to store max values over time
    max_vals = {
    "density": [],
    "vx": [],
    "vy": [],
    "Bx": [],
    "By": [],
    "pressure": [],
    "temperature": [],
    "Eint": [],
    "jz": []
    }



    for idx in range(0, 200):


        filename_bcc = f"./vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk" # b field
        filename_jz = f"./vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
        filename_w = f"./vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"

        # Load the .vtk file
        data_bcc = pv.read(filename_bcc)
        data_jz = pv.read(filename_jz)
        data_w = pv.read(filename_w)

        print(data_bcc.cell_data.keys()) # print vtk labels
        print(data_jz.cell_data.keys()) 
        print(data_w.cell_data.keys()) 

        # Extract cell data or point data

        density = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid))
        vx = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))
        vy = np.array(data_w.cell_data["vely"]).reshape((N_grid, N_grid))
        Bx = np.array(data_bcc.cell_data["bcc1"]).reshape((N_grid, N_grid))
        By = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))
        jz = np.array(data_jz.cell_data["jz"]).reshape((N_grid, N_grid))
        internalE = np.array(data_w.cell_data["eint"]).reshape((N_grid, N_grid))

        # derived quantities
        gamma = 1.6667
        pressure = (gamma - 1.0) * internalE

        temp = pressure/density


        # store max values
        max_vals["density"].append(np.max(density))
        max_vals["vx"].append(np.max(vx))
        max_vals["vy"].append(np.max(vy))
        max_vals["Bx"].append(np.max(Bx))
        max_vals["By"].append(np.max(By))
        max_vals["pressure"].append(np.max(pressure))
        max_vals["temperature"].append(np.max(temp))
        max_vals["Eint"].append(np.max(internalE))
        max_vals["jz"].append(np.max(jz))

        # plot
        im = ax.imshow(density,
                norm = TwoSlopeNorm(vmin=0, vcenter=35, vmax=70),
                extent=[-3, 3, -3, 3],
                cmap=cmr.cm.redshift,
                origin="lower"
                )
        #plt.colorbar(im, ax=ax)
        ax.axis("off")
        ax.set_title('Density')
        camera.snap()



    animation = camera.animate()
    animation.save(f'./animations/{jobname}_density.mp4')



# velx
    # set up figure
    f, ax = plt.subplots(1, 1, figsize=(N_grid/dpi, N_grid/dpi), dpi=dpi)
    camera = Camera(f)

    for idx in range(0, 200):


        filename_bcc = f"./vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk" # b field
        filename_jz = f"./vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
        filename_w = f"./vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"

        # Load the .vtk file
        data_bcc = pv.read(filename_bcc)
        data_jz = pv.read(filename_jz)
        data_w = pv.read(filename_w)

        print(data_bcc.cell_data.keys()) # print vtk labels
        print(data_jz.cell_data.keys()) 
        print(data_w.cell_data.keys()) 

        # Extract cell data or point data

        density = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid))
        vx = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))
        vy = np.array(data_w.cell_data["vely"]).reshape((N_grid, N_grid))
        Bx = np.array(data_bcc.cell_data["bcc1"]).reshape((N_grid, N_grid))
        By = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))

        internalE = np.array(data_w.cell_data["eint"]).reshape((N_grid, N_grid))

        # derived quantities
        gamma = 1.6667
        pressure = (gamma - 1.0) * internalE

        temp = pressure/density


        # plot
        im = ax.imshow(vx,
                norm = TwoSlopeNorm(vmin=0, vcenter=35, vmax=70),
                extent=[-3, 3, -3, 3],
                cmap=cmr.cm.redshift,
                origin="lower"
                )
        #plt.colorbar(im, ax=ax)
        ax.axis("off")
        ax.set_title('Velocity x')
        camera.snap()



    animation = camera.animate()
    animation.save(f'./animations/{jobname}_velx.mp4')


# Vel y
    # set up figure
    f, ax = plt.subplots(1, 1, figsize=(N_grid/dpi, N_grid/dpi), dpi=dpi)
    camera = Camera(f)

    for idx in range(0, 200):


        filename_bcc = f"./vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk" # b field
        filename_jz = f"./vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
        filename_w = f"./vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"

        # Load the .vtk file
        data_bcc = pv.read(filename_bcc)
        data_jz = pv.read(filename_jz)
        data_w = pv.read(filename_w)

        print(data_bcc.cell_data.keys()) # print vtk labels
        print(data_jz.cell_data.keys()) 
        print(data_w.cell_data.keys()) 

        # Extract cell data or point data

        density = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid))
        vx = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))
        vy = np.array(data_w.cell_data["vely"]).reshape((N_grid, N_grid))
        Bx = np.array(data_bcc.cell_data["bcc1"]).reshape((N_grid, N_grid))
        By = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))

        internalE = np.array(data_w.cell_data["eint"]).reshape((N_grid, N_grid))

        # derived quantities
        gamma = 1.6667
        pressure = (gamma - 1.0) * internalE

        temp = pressure/density


        # plot
        im = ax.imshow(vy,
                norm = TwoSlopeNorm(vmin=0, vcenter=35, vmax=70),
                extent=[-3, 3, -3, 3],
                cmap=cmr.cm.redshift,
                origin="lower"
                )
        #plt.colorbar(im, ax=ax)
        ax.axis("off")
        ax.set_title('Velocity y')
        camera.snap()



    animation = camera.animate()
    animation.save(f'./animations/{jobname}_vely.mp4')



# magx
    # set up figure
    f, ax = plt.subplots(1, 1, figsize=(N_grid/dpi, N_grid/dpi), dpi=dpi)
    camera = Camera(f)

    for idx in range(0, 200):


        filename_bcc = f"./vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk" # b field
        filename_jz = f"./vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
        filename_w = f"./vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"

        # Load the .vtk file
        data_bcc = pv.read(filename_bcc)
        data_jz = pv.read(filename_jz)
        data_w = pv.read(filename_w)

        print(data_bcc.cell_data.keys()) # print vtk labels
        print(data_jz.cell_data.keys()) 
        print(data_w.cell_data.keys()) 

        # Extract cell data or point data

        density = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid))
        vx = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))
        vy = np.array(data_w.cell_data["vely"]).reshape((N_grid, N_grid))
        Bx = np.array(data_bcc.cell_data["bcc1"]).reshape((N_grid, N_grid))
        By = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))

        internalE = np.array(data_w.cell_data["eint"]).reshape((N_grid, N_grid))

        # derived quantities
        gamma = 1.6667
        pressure = (gamma - 1.0) * internalE

        temp = pressure/density


        # plot
        im = ax.imshow(Bx,
                norm = TwoSlopeNorm(vmin=0, vcenter=35, vmax=70),
                extent=[-3, 3, -3, 3],
                cmap=cmr.cm.redshift,
                origin="lower"
                )
        #plt.colorbar(im, ax=ax)
        ax.axis("off")
        ax.set_title('Bx')
        camera.snap()


    animation = camera.animate()
    animation.save(f'./animations/{jobname}_Bx.mp4')



# By
    # set up figure
    f, ax = plt.subplots(1, 1, figsize=(N_grid/dpi, N_grid/dpi), dpi=dpi)
    camera = Camera(f)

    for idx in range(0, 200):


        filename_bcc = f"./vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk" # b field
        filename_jz = f"./vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
        filename_w = f"./vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"

        # Load the .vtk file
        data_bcc = pv.read(filename_bcc)
        data_jz = pv.read(filename_jz)
        data_w = pv.read(filename_w)

        print(data_bcc.cell_data.keys()) # print vtk labels
        print(data_jz.cell_data.keys()) 
        print(data_w.cell_data.keys()) 

        # Extract cell data or point data

        density = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid))
        vx = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))
        vy = np.array(data_w.cell_data["vely"]).reshape((N_grid, N_grid))
        Bx = np.array(data_bcc.cell_data["bcc1"]).reshape((N_grid, N_grid))
        By = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))

        internalE = np.array(data_w.cell_data["eint"]).reshape((N_grid, N_grid))

        # derived quantities
        gamma = 1.6667
        pressure = (gamma - 1.0) * internalE

        temp = pressure/density


        # plot
        im = ax.imshow(By,
                norm = TwoSlopeNorm(vmin=0, vcenter=35, vmax=70),
                extent=[-3, 3, -3, 3],
                cmap=cmr.cm.redshift,
                origin="lower"
                )
        #plt.colorbar(im, ax=ax)
        ax.axis("off")
        ax.set_title('By')
        camera.snap()



    animation = camera.animate()
    animation.save(f'./animations/{jobname}_By.mp4')
#plt.savefig(f"./density_plots/CS_{str(idx).zfill(5)}.png", dpi=dpi)


# pressure
    # set up figure
    f, ax = plt.subplots(1, 1, figsize=(N_grid/dpi, N_grid/dpi), dpi=dpi)
    camera = Camera(f)

    for idx in range(0, 200):


        filename_bcc = f"./vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk" # b field
        filename_jz = f"./vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
        filename_w = f"./vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"

        # Load the .vtk file
        data_bcc = pv.read(filename_bcc)
        data_jz = pv.read(filename_jz)
        data_w = pv.read(filename_w)

        print(data_bcc.cell_data.keys()) # print vtk labels
        print(data_jz.cell_data.keys()) 
        print(data_w.cell_data.keys()) 

        # Extract cell data or point data

        density = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid))
        vx = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))
        vy = np.array(data_w.cell_data["vely"]).reshape((N_grid, N_grid))
        Bx = np.array(data_bcc.cell_data["bcc1"]).reshape((N_grid, N_grid))
        By = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))

        internalE = np.array(data_w.cell_data["eint"]).reshape((N_grid, N_grid))

        # derived quantities
        gamma = 1.6667
        pressure = (gamma - 1.0) * internalE

        temp = pressure/density


        # plot
        im = ax.imshow(pressure,
                norm = TwoSlopeNorm(vmin=0, vcenter=35, vmax=70),
                extent=[-3, 3, -3, 3],
                cmap=cmr.cm.redshift,
                origin="lower"
                )
        #plt.colorbar(im, ax=ax)
        ax.axis("off")
        ax.set_title('Pressure')
        camera.snap()



    animation = camera.animate()
    animation.save(f'./animations/{jobname}_pressure.mp4')



# temp
    # set up figure
    f, ax = plt.subplots(1, 1, figsize=(N_grid/dpi, N_grid/dpi), dpi=dpi)
    camera = Camera(f)

    for idx in range(0, 200):


        filename_bcc = f"./vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk" # b field
        filename_jz = f"./vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
        filename_w = f"./vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"

        # Load the .vtk file
        data_bcc = pv.read(filename_bcc)
        data_jz = pv.read(filename_jz)
        data_w = pv.read(filename_w)

        print(data_bcc.cell_data.keys()) # print vtk labels
        print(data_jz.cell_data.keys()) 
        print(data_w.cell_data.keys()) 

        # Extract cell data or point data

        density = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid))
        vx = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))
        vy = np.array(data_w.cell_data["vely"]).reshape((N_grid, N_grid))
        Bx = np.array(data_bcc.cell_data["bcc1"]).reshape((N_grid, N_grid))
        By = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))

        internalE = np.array(data_w.cell_data["eint"]).reshape((N_grid, N_grid))

        # derived quantities
        gamma = 1.6667
        pressure = (gamma - 1.0) * internalE

        temp = pressure/density


        # plot
        im = ax.imshow(temp,
                norm = TwoSlopeNorm(vmin=0, vcenter=35, vmax=70),
                extent=[-3, 3, -3, 3],
                cmap=cmr.cm.redshift,
                origin="lower"
                )
        #plt.colorbar(im, ax=ax)
        ax.axis("off")
        ax.set_title('Temperature')
        camera.snap()



    animation = camera.animate()
    animation.save(f'./animations/{jobname}_temp.mp4')






# jz
    # set up figure
    f, ax = plt.subplots(1, 1, figsize=(N_grid/dpi, N_grid/dpi), dpi=dpi)
    camera = Camera(f)

    for idx in range(0, 200):


        filename_bcc = f"./vtk/{jobname}.mhd_bcc.{str(idx).zfill(5)}.vtk" # b field
        filename_jz = f"./vtk/{jobname}.mhd_jz.{str(idx).zfill(5)}.vtk"
        filename_w = f"./vtk/{jobname}.mhd_w.{str(idx).zfill(5)}.vtk"

        # Load the .vtk file
        data_bcc = pv.read(filename_bcc)
        data_jz = pv.read(filename_jz)
        data_w = pv.read(filename_w)

        print(data_bcc.cell_data.keys()) # print vtk labels
        print(data_jz.cell_data.keys()) 
        print(data_w.cell_data.keys()) 

        # Extract cell data or point data

        density = np.array(data_w.cell_data["dens"]).reshape((N_grid, N_grid))
        vx = np.array(data_w.cell_data["velx"]).reshape((N_grid, N_grid))
        vy = np.array(data_w.cell_data["vely"]).reshape((N_grid, N_grid))
        Bx = np.array(data_bcc.cell_data["bcc1"]).reshape((N_grid, N_grid))
        By = np.array(data_bcc.cell_data["bcc2"]).reshape((N_grid, N_grid))
        jz = np.array(data_jz.cell_data["jz"]).reshape((N_grid, N_grid))


        # derived quantities
        gamma = 1.6667
        pressure = (gamma - 1.0) * internalE

        temp = pressure/density


        # plot
        im = ax.imshow(jz,
                norm = TwoSlopeNorm(vmin=0, vcenter=35, vmax=70),
                extent=[-3, 3, -3, 3],
                cmap=cmr.cm.redshift,
                origin="lower"
                )
        #plt.colorbar(im, ax=ax)
        ax.axis("off")
        ax.set_title('Jz')
        camera.snap()



    animation = camera.animate()
    animation.save(f'./animations/jz.mp4')




# Create time array
time_steps = list(range(0, len(max_vals["density"])))

# Plot max values
plt.figure(figsize=(10, 6))
for key in max_vals:
    plt.plot(time_steps, max_vals[key], label=key)

plt.xlabel("Timestep Index")
plt.ylabel("Maximum Value")
plt.title("Maximum Value of Quantities Over Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"./animations/{jobname}_max_values_over_time.png", dpi=150)
plt.show()