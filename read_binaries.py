import bin_convert
import os
import numpy as np
import matplotlib.pyplot as plt

class ReadAthenaBinData:
    """
    Container for stitched Athena bin outputs with convenient component access.

    - vel: array shape (3, nx, ny, nz) with components [vx, vy, vz]
    - mag: array shape (3, nx, ny, nz) with components [Bx, By, Bz] (if provided)
    - cur: array shape (3, nx, ny, nz) with components [Jx, Jy, Jz] (if provided)
    - dens: array shape (1, nx, ny, nz)
    - press: array shape (1, nx, ny, nz)
    
    In general, scalar fields have the shape (1, nx, ny, nz) for consistency.
    3-component vector fields have shape (3, nx, ny, nz) with components along axis 0.
    
    Usage:
        data = AthenaBinData(file_number=40, basename="CS", data_dir=".")
        data.read("mag")
        data.read("vel")
        
    Author: 
        James R. Beattie
    
    """

    def __init__(self, 
                 file_number: int, 
                 basename: str = "CS", 
                 data_dir: str = ".",
                 w_binary_fname: str | None = None, 
                 b_binary_fname: str | None = None,
                 write_hdf: bool = False) -> None:
        
        self.file_number = file_number
        self.basename = basename
        self.data_dir = data_dir
        self.w_binary_fname = w_binary_fname or self._default_fname("mhd_w")
        self.b_binary_fname = b_binary_fname or self._default_fname("mhd_bcc")
        self.write_hdf = write_hdf
        self.meta_w = None
        self.meta_b = None
        self.dens = None
        self.vel = None
        self.press = None
        self.mag = None
        self.cur = None
        self._loaded_w = set()
        self._cache_w = {}

    def _default_fname(self, 
                       kind: str) -> str:
        return os.path.join(
            self.data_dir,
            f"{self.basename}.{kind}.{self.file_number:05d}.bin",
        )

    def _convert_and_stack(self, 
                           binary_fname: str, 
                           variables=None) -> dict:
        """
        Convert binary output to athdf/xdmf and return requested variables stitched as (x,y,z).
        Only the requested variables are read and stitched.
        """
        filedata = bin_convert.read_binary(binary_fname)
        athdf_fname = None
        xdmf_fname = None
        if self.write_hdf:
            athdf_fname = binary_fname.replace(".bin", ".athdf")
            xdmf_fname = athdf_fname + ".xdmf"
            print(f"[convert] Writing HDF5/XDMF for {binary_fname} -> {athdf_fname}")
            bin_convert.write_athdf(athdf_fname, filedata)
            bin_convert.write_xdmf_for(xdmf_fname, os.path.basename(athdf_fname), filedata)
        else:
            print(f"[convert] Skipping HDF5/XDMF write for {binary_fname}")

        levels = np.unique(filedata["mb_logical"][:, 3])
        if len(levels) != 1:
            raise ValueError(f"mixed refinement levels not supported (levels={levels})")

        nx1_mb = filedata["nx1_out_mb"]
        nx2_mb = filedata["nx2_out_mb"]
        nx3_mb = filedata["nx3_out_mb"]

        lx1 = filedata["mb_logical"][:, 0]
        lx2 = filedata["mb_logical"][:, 1]
        lx3 = filedata["mb_logical"][:, 2]

        lx1_min, lx1_max = lx1.min(), lx1.max()
        lx2_min, lx2_max = lx2.min(), lx2.max()
        lx3_min, lx3_max = lx3.min(), lx3.max()

        nx1_global = (lx1_max - lx1_min + 1) * nx1_mb
        nx2_global = (lx2_max - lx2_min + 1) * nx2_mb
        nx3_global = (lx3_max - lx3_min + 1) * nx3_mb

        vars_to_use = variables if variables is not None else filedata["var_names"]
        missing = [v for v in vars_to_use if v not in filedata["var_names"]]
        if missing:
            raise ValueError(f"Requested variables not in file: {missing}")

        # allocate in (x, y, z) order to avoid an extra full-array transpose later
        stacked = {}
        for name in vars_to_use:
            dtype = filedata["mb_data"][name][0].dtype
            stacked[name] = np.zeros((nx1_global, nx2_global, nx3_global), dtype=dtype)

        for mb_id in range(filedata["n_mbs"]):
            lx1_idx = filedata["mb_logical"][mb_id][0] - lx1_min
            lx2_idx = filedata["mb_logical"][mb_id][1] - lx2_min
            lx3_idx = filedata["mb_logical"][mb_id][2] - lx3_min

            i_start = lx1_idx * nx1_mb
            j_start = lx2_idx * nx2_mb
            k_start = lx3_idx * nx3_mb

            i_end = i_start + nx1_mb
            j_end = j_start + nx2_mb
            k_end = k_start + nx3_mb

            for name in vars_to_use:
                block = filedata["mb_data"][name][mb_id]  # (k, j, i)
                stacked[name][i_start:i_end, j_start:j_end, k_start:k_end] = np.transpose(block, (2, 1, 0))

        return {
            "data": stacked,
            "athdf": athdf_fname,
            "xdmf": xdmf_fname,
            "filedata": filedata,
        }

    def _convert_and_stack_b(self,
                             binary_fname: str):
        """Convenience wrapper for magnetic components; returns Bx, By, Bz."""
        b_names = ["bcc1", "bcc2", "bcc3"]
        res = self._convert_and_stack(binary_fname,
                                      variables=b_names)
        return {
            "Bx": res["data"][b_names[0]],
            "By": res["data"][b_names[1]],
            "Bz": res["data"][b_names[2]],
            "athdf": res["athdf"],
            "xdmf": res["xdmf"],
            "filedata": res["filedata"],
        }

    @staticmethod
    def _get_param_from_header(filedata,
                               block_name: str, 
                               key: str, 
                               cast=float) -> float:
        """Minimal parser for values in the stored header."""
        block = None
        for line in filedata["header"]:
            if line.startswith("<"):
                block = line
                continue
            if block == f"<{block_name}>":
                if line.strip().startswith(f"{key}"):
                    _, val = line.split("=")
                    return cast(val.strip())
        raise KeyError(f"{block_name}/{key} not found in header")

    def read(self, 
             what: str):
        """
        Lazily read specified field(s) into the object.

        what: 'vel', 'dens', 'press', 'mag', 'cur'
        """
        if what in ("vel", "dens", "press"):
            # Determine which hydro variables are missing
            needed = []
            if what in ("dens", "press") and "dens" not in self._loaded_w:
                needed.append("dens")
            if what == "press" and "eint" not in self._loaded_w:
                needed.append("eint")
            if what == "vel":
                for vname in ("velx", "vely", "velz"):
                    if vname not in self._loaded_w:
                        needed.append(vname)

            if needed:
                w_res = self._convert_and_stack(self.w_binary_fname, variables=needed)
                if self.meta_w is None:
                    self.meta_w = w_res["filedata"]
                self._loaded_w.update(needed)
                self._cache_w.update(w_res["data"])

            if what in ("dens", "press") and self.dens is None and "dens" in self._cache_w:
                self.dens = self._cache_w["dens"][None, ...]
            if what == "vel" and self.vel is None:
                if all(v in self._cache_w for v in ("velx", "vely", "velz")):
                    self.vel = np.stack(
                        [
                            self._cache_w["velx"],
                            self._cache_w["vely"],
                            self._cache_w["velz"],
                        ],
                        axis=0,
                    )
            if what == "press" and self.press is None:
                if self.dens is not None and "eint" in self._cache_w:
                    gamma = self._get_param_from_header(self.meta_w, "mhd", "gamma")
                    self.press = (gamma - 1.0) * self.dens * self._cache_w["eint"][None, ...]
        if what == "mag":
            if self.b_binary_fname is None:
                raise ValueError("No magnetic field binary file provided.")
            if self.meta_b is None or self.mag is None:
                b_res = self._convert_and_stack_b(self.b_binary_fname)
                self.meta_b = b_res["filedata"]
                self.mag = np.stack(
                    [
                        b_res["Bx"],
                        b_res["By"],
                        b_res["Bz"]
                    ],
                    axis=0,
                )
        if what == "cur":
            cur_fname = self._default_fname("mhd_jz")
            c_res = self._convert_and_stack(cur_fname)
            # Expect curx/cury/curz or a single-component current (e.g., jz)
            names = [n for n in ("curx", "cury", "curz") if n in c_res["data"]]
            if not names:
                names = list(c_res["data"].keys())
            self.cur = np.stack([c_res["data"][n] for n in names], axis=0)
            
if __name__ == "__main__":
    
    data = ReadAthenaBinData(file_number=40, basename="CS", data_dir=".")
    data.read("mag")
    data.read("vel")
    data.read("dens")
    data.read("press")
    data.read("cur")
    print("mag shape", data.mag.shape)
    print("vel shape", data.vel.shape)
    print("dens shape", data.dens.shape)
    print("press shape", data.press.shape)
    print("cur shape", data.cur.shape)
    