import h5py
import numpy as np

fname = "../h5py/gamma1_67_beta1_S10000_gpa13_33_nx1600_epsp0_0000.mhd_w.00000.athdf"
with h5py.File(fname, "r") as f:
    data = f["uov"][:]
    print("Shape of uov:", data.shape)
    print("Min/max per channel:")
    for i in range(data.shape[2]):
        print(f"Channel {i}: min={np.min(data[:,:,i])}, max={np.max(data[:,:,i])}")


# Open file in read mode
with h5py.File("../h5py/gamma1_67_beta1_S10000_gpa13_33_nx1600_epsp0_0000.mhd_w.00000.athdf", "r") as f:
    # List top-level groups/datasets
    print(list(f.keys()))
    # Access a dataset
    uov = f["uov"][:]  # reads the full array into memory
    # Access an attribute
    time = f.attrs["Time"]