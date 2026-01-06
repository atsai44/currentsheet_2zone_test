import numpy as np
import matplotlib.pyplot as plt

# --- Load the data ---
#jobname = 'gamma1_50_beta0_001_S10000_gpa13_33_nx1600'
#jobname = 'gamma1_50_beta1_S10000_gpa13_33_nx1600'
jobname = 'gamma1_50_beta1000_S10000_gpa13_33_nx1600'
filename = f"../{jobname}.mhd.hst"
data = np.loadtxt(filename, comments='#')

# --- Extract columns ---
time = data[:, 0]
mass = data[:, 2]
mom1 = data[:, 3]
mom2 = data[:, 4]
mom3 = data[:, 5]
totE = data[:, 6]
KE1  = data[:, 7]
KE2  = data[:, 8]
KE3  = data[:, 9]
ME1  = data[:, 10]
ME2  = data[:, 11]
ME3  = data[:, 12]

# --- Plot everything ---
plt.figure(figsize=(10, 6))

plt.plot(time, mass, label='Mass')
plt.plot(time, totE, label='Total Energy')
plt.plot(time, KE1, '--', label='KE1')
plt.plot(time, KE2, '--', label='KE2')
plt.plot(time, KE3, '--', label='KE3')
plt.plot(time, ME1, ':', label='ME1')
plt.plot(time, ME2, ':', label='ME2')
plt.plot(time, ME3, ':', label='ME3')

plt.xlabel('Time')
plt.xscale('log')
plt.yscale('log')
plt.title(f'Quantities from HST file: {jobname}')
plt.legend()
plt.grid(True)
plt.tight_layout()


outname = f"./plots_init/{jobname}_Energy.png"
plt.savefig(outname)
plt.close()