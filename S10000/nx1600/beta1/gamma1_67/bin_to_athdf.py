import os
import glob
from bin_convert import convert_file
from bin_convert import read_binary
from bin_convert import write_athdf
from bin_convert import write_xdmf_for

# --- Directories ---
bin_dir = "./bin"
athdf_dir = "./athdf"

# Create h5py folder if it doesn't exist
os.makedirs(athdf_dir, exist_ok=True)

# Get all .bin files in bin_dir
bin_files = sorted(glob.glob(os.path.join(bin_dir, "*.bin")))

if not bin_files:
    print(f"No .bin files found in {bin_dir}")
else:
    print(f"Found {len(bin_files)} .bin files. Converting...")

# Loop over files and convert
for bin_file in bin_files:
    # Read binary
    filedata = read_binary(bin_file)
    
    # Construct output filenames in h5py folder
    base_name = os.path.basename(bin_file).replace(".bin", "")
    athdf_fname = os.path.join(athdf_dir, base_name + ".athdf")
    xdmf_fname  = os.path.join(athdf_dir, base_name + ".athdf.xdmf")
    
    # Write HDF5 and XDMF directly to h5py folder
    write_athdf(athdf_fname, filedata)
    write_xdmf_for(xdmf_fname, os.path.basename(athdf_fname), filedata)

print("Conversion complete. .athdf files are in the athdf/ folder.")