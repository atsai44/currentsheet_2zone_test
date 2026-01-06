import os
import numpy as np

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
epsp = 0.2*p0  #0.2 * p0 #0.5 * p0 #0.15 * p0


# set t_limit and dt using t_crossing
t_c = L / cms
#t_limit = 1000 * t_c
t_limit = 100 * t_c # increase to 100 or 1000
#t_limit = 10 * t_c
dt_vtk = 0.1 * t_c
dt_hst = 1 * t_c
dt_restart = 0.5 * t_c

print(f't limit: {t_limit:0.6f}')
print(f'dt vtk: {dt_vtk:0.6f}')
print(f'dt hst: {dt_hst:0.6f}')
print(f'dt restart: {dt_restart:0.6f}')

# --- Read template ---
with open("input_temp.athinput") as f:
    template = f.read()

# --- File naming ---
raw_name = f"gamma{gamma:.2f}_beta{beta}_S{S}_gpa{grid_cells_per_a:.2f}_nx{nx}_epsp{epsp:.4f}"
#raw_name = f"rst1_gamma{gamma:.2f}_beta{beta}_S{S}_gpa{grid_cells_per_a}_nx{nx}"
safe_name = raw_name.replace(".", "_")
filename = safe_name + ".athinput"

# --- Print summary ---
print("\n" + "=" * 60)
print(f"Creating: {filename}")
print(f" gamma = {gamma:.2f}, beta = {beta}, S = {S}, grid_cells_per_a = {grid_cells_per_a:.2f}, epsp{epsp:.4f}")
print(f" Grid: nx = {nx}, ny = {ny}, L = {L}, a0 = {a0:.2f}")
print(f" Derived: alpha = {alpha:.2f}, b0 = {b0:.2f}, eta = {eta:.2e}")
print("=" * 60)

# --- Replace placeholders ---
text = (template
    .replace("BASE_NAME", safe_name)
    .replace("N_X", str(nx))
    .replace("N_Y", str(ny))
    .replace("N_Z", str(nz))
    .replace("GAMMA", f"{gamma:.4f}")
    .replace("T_LIMIT", f"{t_limit:.4f}")
    .replace("DT_VTK", f"{dt_vtk:.6f}")
    .replace("DT_HST", f"{dt_hst:.4f}")
    .replace("DT_RESTART", f"{dt_restart:.4f}")
    .replace("EPSP", f"{epsp:.4f}")
    .replace("A0", f"{a0:.4f}")
    .replace("ETA_OHM", f"{eta:.4e}")
    .replace("B0", f"{b0:.4f}")
    .replace("L2X", f"{L/2:.4f}")
    .replace("L2Z", f"{Lz/2:.4f}")
    .replace("P0", f"{p0:.4f}")
)

# --- Write input file to current directory ---
with open(filename, "w") as out:
    out.write(text)

print(f"Wrote {filename} successfully\n")

