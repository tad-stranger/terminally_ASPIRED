import json
from pathlib import Path
import pandas as pd
import numpy as np
from astropy.io import fits
from terminally_ASPIRED import image_reduction, spectral_reduction
# from aspired import image_reduction, spectral_reduction
from astroscrappy import detect_cosmics
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.widgets import RectangleSelector
import warnings
import sys
import shutil
from importlib import resources

# This is a class of my 1.9M pipeline to be used to make terminally_ASPIRED

class Tee:
    """Redirects writes to both console and file"""
    def __init__(self, logfile, stream = sys.stdout):
        self.stream = stream
        self.log = open(logfile, 'w')

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
        self.log.write(data)
        self.log.flush()

    def flush(self):
        self.stream.flush()
        self.log.flush()


class SpectralReductionPipeline:
    def __init__(self, science_file, arc_file, std_file, std_arc_file, config_path="",
                 grating="7", bias_path=None, flat_path=None, interactive_trim=False, show_plots=False, smooth=1, verbose=False,
                 no_warnings=True, output_dir_name=None, sky=False, error=False):

        # Hack to get depreciated pkg_resources to work by tricking it into believing the aspired package is present.
        # an alternative is to replace every call to:
        #   pkg_resources.resource_filename('aspired', ...)
        # with the corrected:
        #   from importlib import resources
        #   filename = resources.files("terminally_aspired.vendor.aspired_forked").joinpath("data/some_file.dat")
        # Not sure how many calls to pkg_resources there are so went with hacky fix that works ¯\_(ツ)_/¯
        import sys
        import terminally_ASPIRED.vendor.aspired_fork as aspired
        sys.modules['aspired'] = aspired

        self.science_path = Path(science_file)
        self.arc_path = Path(arc_file)
        self.std_path = Path(std_file)
        self.arc_std_path = Path(std_arc_file)
        self.grating = grating
        self.show_sky = sky
        self.show_error = error
        self.smoothing_value = smooth
        self.verbose = verbose
        self.master_bias = None
        self.interactive_trim = interactive_trim
        # Deduce base observation directory (e.g., ".../0503")
        # self.obs_base = self.science_path.parents[2]
        if no_warnings:
            warnings.filterwarnings("ignore")
        if bias_path == "":
            self.bias_folder = Path("DO_NOT_USE_BIAS")
        else:
            self.bias_folder = Path(bias_path)

        if flat_path == "":
            self.flat_folder = Path("DO_NOT_USE_FLATS")
        else:
            self.flat_folder = Path(flat_path)
        if config_path == "":
            with resources.open_text("terminally_ASPIRED.config_files", "defaults.json") as f:
                self.config = json.load(f)
        else:
            self.config = self._load_config(config_path)

        self.show_plots = show_plots
        self.config["display_plots"] = show_plots

        # Placeholder attributes
        self.object_name = None
        self.output_dir = output_dir_name
        self.cleaned = {}

        # FITS data
        self.sci_data = None
        self.arc_data = None
        self.std_data = None
        self.arc_std_data = None
        self.hdr_sci = None
        self.hdr_arc = None
        self.hdr_std = None



    def _load_config(self, path):
        with open(path, 'r') as f:
            return json.load(f)

    def _trim_raw_data(self):
        bounds = self.config["trim_bounds"]

        def trim(data):
            return data[
                   bounds["y_min"]:bounds["y_max"],
                   bounds["x_min"]:bounds["x_max"]
                   ]

        self.sci_data = trim(self.sci_data)
        self.arc_data = trim(self.arc_data)
        self.std_data = trim(self.std_data)
        self.arc_std_data = trim(self.arc_std_data)

    def extract_data(self):
        def _extract(path):
            data, header = fits.getdata(path, header=True)
            return np.flip(data, axis=1), header  # flip for ASPIRED format

        sci_path = self.science_path
        arc_path = self.arc_path
        std_path = self.std_path
        arc_std_path = self.arc_std_path

        self.sci_data, self.hdr_sci = _extract(sci_path)
        self.arc_data, self.hdr_arc = _extract(arc_path)
        self.std_data, self.hdr_std = _extract(std_path)
        self.arc_std_data, self.hdr_arc_std = _extract(arc_std_path)

        # Check if raw files all have the same shape.
        shapes  = {"science": self.sci_data.shape,
                   "arc": self.arc_data.shape,
                   "std": self.std_data.shape,
                   "arc_std": self.arc_std_data.shape}

        if len(set(shapes.values())) != 1:
            raise ValueError(f"Shape mismatch in raw frames before trimming: {shapes}")


        # Trim Raw files
        if self.interactive_trim:
            self.interactive_trim_function(tag = "science")
            self._trim_raw_data()
        else:
            self._trim_raw_data()

        # Load and trim bias frames (if the bias folder is set and exists)
        if self.bias_folder != Path("DO_NOT_USE_BIAS"):
            bias_files = sorted(self.bias_folder.glob("*.fits"))
            self.bias_data_list = []
            for bf in bias_files:
                data, header = fits.getdata(bf, header=True)
                data = np.flip(data, axis=1)
                # Trim using the same bounds
                data = data[
                       self.config["trim_bounds"]["y_min"]:self.config["trim_bounds"]["y_max"],
                       self.config["trim_bounds"]["x_min"]:self.config["trim_bounds"]["x_max"]
                       ]
                self.bias_data_list.append(data)

        # Load and trim flat frames (if a flat folder is set and exists)
        if self.flat_folder != Path("DO_NOT_USE_FLATS"):
            flat_files = sorted(self.flat_folder.glob("*.fits"))
            self._original_flat_files = flat_files # Keep for later saving
            self.flat_data_list = []
            for ff in flat_files:
                data, header = fits.getdata(ff, header=True)
                data = np.flip(data, axis=1)
                # Trim using the same bounds
                data = data[
                       self.config["trim_bounds"]["y_min"]:self.config["trim_bounds"]["y_max"],
                       self.config["trim_bounds"]["x_min"]:self.config["trim_bounds"]["x_max"]
                       ]
                self.flat_data_list.append(data)

        self.object_name = self.hdr_sci.get("OBJECT", "Unknown").replace(" ", "")
        if self.output_dir is None:
            self.output_dir = Path(f"./ReducedSpectra/{self.object_name}")
        else:
            self.output_dir = Path(f"./ReducedSpectra/{self.output_dir}")
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def bias_subtract(self):
        if self.bias_folder == Path("DO_NOT_USE_BIAS") or not hasattr(self, "bias_data_list") or len(
                self.bias_data_list) == 0:
            if self.verbose:
                print("No Bias files supplied or loaded. Bias subtraction skipped.")
            return

        # Compute master bias from trimmed bias frames
        self.master_bias = np.median(self.bias_data_list, axis=0)
        self.master_bias = np.nan_to_num(self.master_bias, nan=0.0)

        def subtract_and_clean(image, label):
            image = np.nan_to_num(image, nan=0.0)
            subtracted_image = image - self.master_bias
            subtracted_image = np.nan_to_num(subtracted_image, nan=0.0)
            subtracted_image[subtracted_image < 0] = 0.0
            out_path = self.output_dir / f"{self.object_name}_{label}_bias_subtracted.fits"
            fits.writeto(out_path, subtracted_image.astype(np.float32), overwrite=True)
            return subtracted_image

        self.sci_data = subtract_and_clean(self.sci_data, 'sci')
        self.std_data = subtract_and_clean(self.std_data, 'std')
        self.arc_data = subtract_and_clean(self.arc_data, 'arc')
        self.arc_std_data = subtract_and_clean(self.arc_std_data, 'arc_std')

        # Also subtract bias from trimmed flats (if any)
        if self.flat_folder != Path("DO_NOT_USE_FLATS") and hasattr(self, "flat_data_list") and len(
                self.flat_data_list) > 0:
            bias_flat_dir = self.output_dir / "flat_bias_subtracted"
            bias_flat_dir.mkdir(parents=True, exist_ok=True)

            for i, flat_data in enumerate(self.flat_data_list):
                flat_data = np.nan_to_num(flat_data, nan=0.0)
                subtracted_flat = flat_data - self.master_bias
                subtracted_flat = np.nan_to_num(subtracted_flat, nan=0.0)
                subtracted_flat[subtracted_flat < 0] = 0.0

                # Save cleaned flat with the original filename
                flat_files = sorted(self.flat_folder.glob("*.fits"))
                out_path = bias_flat_dir / flat_files[i].name
                fits.writeto(out_path, subtracted_flat.astype(np.float32), overwrite=True)

            self.flat_folder = bias_flat_dir

    def _ensure_trimmed_flats_on_disk(self):
        """If no bias is used, write the in-memory trimmed flats to disk so the reducer
        gets correctly trimmed files with the right shape."""
        if self.flat_folder == Path("DO_NOT_USE_FLATS"):
            return
        if not hasattr(self, "flat_data_list") or len(self.flat_data_list) == 0:
            return

        # If bias_subtract() ran, it already wrote bias-subtracted (and trimmed) flats
        # and repointed self.flat_folder to that directory. In that case, do nothing.
        # We only need this when NO bias was used.
        if self.master_bias is not None:
            return

        # Write trimmed-only flats to disk, mirror original filenames
        out_dir = self.output_dir / "flat_trimmed"
        out_dir.mkdir(parents=True, exist_ok=True)

        # Use original filenames if available, else just enumerate
        if hasattr(self, "_original_flat_files"):
            names = [p.name for p in self._original_flat_files]
        else:
            names = [f"flat_{i:03d}.fits" for i in range(len(self.flat_data_list))]

        for arr, name in zip(self.flat_data_list, names):
            # Write with a minimal header; if you want to preserve headers, load them above and save here.
            fits.writeto(out_dir / name, np.asarray(arr, dtype=np.float32), overwrite=True)

        # Point flat_folder to the trimmed set on disk so filelists use the right files
        self.flat_folder = out_dir


    def write_filelists(self):

        # Ensure flats on disk have the same shape as trimmed science/arc
        self._ensure_trimmed_flats_on_disk()
        def write(filelist_path, science_file, arc_file):
            flat_files = sorted(self.flat_folder.glob("*.fits"))
            with open(filelist_path, "w") as f:
                for file in flat_files:
                    f.write(f"flat, {file.resolve()}\n")
                f.write(f"arc, {arc_file}\n")
                f.write(f"light, {science_file}\n")

        if self.bias_folder == Path("DO_NOT_USE_BIAS"):
            self.arc_std_file_resolved = self._save_fits(self.arc_std_data, self.hdr_arc_std,f"{self.object_name}_arc_std.fits")
            self.arc_file_resolved = self._save_fits(self.arc_data, self.hdr_arc,f"{self.object_name}_arc.fits")
            sci_resolved = self._save_fits(self.sci_data, self.hdr_sci, f"{self.object_name}_sci.fits")
            std_resolved = self._save_fits(self.std_data, self.hdr_std, f"{self.object_name}_std.fits")
        else:
            self.arc_std_file_resolved = self._save_fits(self.arc_std_data, self.hdr_arc_std,f"{self.object_name}_arc_std_bias_subtracted.fits")
            self.arc_file_resolved = self._save_fits(self.arc_data, self.hdr_arc,f"{self.object_name}_arc_bias_subtracted.fits")
            sci_resolved = self._save_fits(self.sci_data, self.hdr_sci, f"{self.object_name}_sci_bias_subtracted.fits")
            std_resolved = self._save_fits(self.std_data, self.hdr_std, f"{self.object_name}_std_bias_subtracted.fits")


        # Save FITS files first
        write(self.output_dir / f"{self.object_name}_science_file.list", sci_resolved, self.arc_file_resolved)
        write(self.output_dir / f"{self.object_name}_standard_file.list", std_resolved, self.arc_std_file_resolved)

    def _save_fits(self, data, header, filename):
        path = (self.output_dir / filename).resolve()
        fits.writeto(path, data, header, overwrite=True)
        return path


    def reduce_images(self):
        def reduce(list_file, tag):
            ir = image_reduction.ImageReduction(verbose=self.verbose)
            ir.add_filelist(str(list_file))
            ir.load_data()
            ir.reduce()
            ir.save_fits(str(self.output_dir / f"{self.object_name}_{tag}_image_reduced.fits"), overwrite=True)

        reduce(self.output_dir / f"{self.object_name}_science_file.list", tag="science")
        reduce(self.output_dir / f"{self.object_name}_standard_file.list", tag="standard")

    def cosmic_clean(self):
        bounds = self.config["trim_bounds"]

        def clean(data):
            _, cleaned = detect_cosmics(data, **self.config["cosmic_ray"], verbose=self.verbose)
            return cleaned

        self.cleaned["sci"] = np.nan_to_num(clean(self._load_image("science")), nan = 0.0)
        self.cleaned["arc"] = np.nan_to_num(clean(self.arc_data), nan = 0.0)
        self.cleaned["std"] = np.nan_to_num(clean(self._load_image("standard")), nan = 0.0)
        self.cleaned["arc_std"] = np.nan_to_num(clean(self.arc_std_data), nan = 0.0)

    def _load_image(self, tag):
        data, _ = fits.getdata(self.output_dir / f"{self.object_name}_{tag}_image_reduced.fits", header=True)
        return data

    def extract_2dspec(self):
        self.sci2d = self._extract_twodspec(self.cleaned["sci"], self.cleaned["arc"], self.hdr_sci, is_standard=False)
        self.std2d = self._extract_twodspec(self.cleaned["std"], self.cleaned["arc_std"], self.hdr_std, is_standard=True)

    def _extract_twodspec(self, data, arc, header, is_standard):
        twod = spectral_reduction.TwoDSpec(data, header=header, cosmicray=False, verbose=self.verbose)
        twod.set_properties(saxis=1, flip=False, cosmicray=False, verbose=self.verbose)

        trace_key = "standard" if is_standard else "science"
        trace_kwargs = self.config["trace_kwargs"][trace_key].copy()
        extract_kwargs = self.config["extract_kwargs"][trace_key].copy()

        # Inject display flag
        trace_kwargs["display"] = self.config.get("display_plots", False)
        extract_kwargs["display"] = self.config.get("display_plots", False)

        twod.ap_trace(**trace_kwargs)
        twod.ap_extract(**extract_kwargs)

        twod.add_arc(arc=arc)
        twod.extract_arc_spec()
        return twod

    def calibrate_wavelength(self):
        self.onedspec = spectral_reduction.OneDSpec(verbose=self.verbose)
        self.onedspec.from_twodspec(self.sci2d, stype="science")
        self.onedspec.from_twodspec(self.std2d, stype="standard")

        cal_cfg = self.config["wavelength_cal"]
        atlas = self.config["atlas_lines"][self.grating]
        element = ["CuAr"] * len(atlas)

        self.onedspec.find_arc_lines(
            prominence=cal_cfg["science"]["prominence"],
            refine=cal_cfg["science"]["refine"],
            display=self.config["display_plots"],
            stype='science'
        )

        self.onedspec.find_arc_lines(
            prominence=cal_cfg["standard"]["prominence"],
            refine=cal_cfg["standard"]["refine"],
            display=self.config["display_plots"],
            stype='standard'
        )

        self.onedspec.initialise_calibrator(stype ='science+standard')

        self.onedspec.add_user_atlas(
            elements=element,
            wavelengths=atlas,
            stype='science+standard'
        )

        self.onedspec.set_hough_properties(**self.config["hough"])

        self.onedspec.do_hough_transform(stype='science+standard')

        self.onedspec.verbose = True

        self.onedspec.set_ransac_properties(
            minimum_matches=cal_cfg["science"]["ransac_minimum_matches"],
            stype="science"
        )

        self.onedspec.set_ransac_properties(
            minimum_matches=cal_cfg["standard"]["ransac_minimum_matches"],
            stype="standard"
        )

        self.onedspec.fit(
            max_tries=cal_cfg["science"]["fit_max_tries"],
            fit_deg=cal_cfg["science"]["fit_deg"],
            fit_type=cal_cfg["science"]["fit_type"],
            display=self.config["display_plots"],
            stype='science'
        )

        self.onedspec.fit(
            max_tries=cal_cfg["standard"]["fit_max_tries"],
            fit_deg=cal_cfg["standard"]["fit_deg"],
            fit_type=cal_cfg["standard"]["fit_type"],
            display=self.config["display_plots"],
            stype='standard'
        )

        self.onedspec.apply_wavelength_calibration(stype='science')
        self.onedspec.apply_wavelength_calibration(stype='standard')

    def calibrate_flux(self):
        std_name = self.config.get("standard_name")
        self.onedspec.load_standard(target=std_name)
        if self.config["display_plots"]:
            self.onedspec.inspect_standard()
        self.onedspec.get_sensitivity()
        if self.config["display_plots"]:
            self.onedspec.inspect_sensitivity()
        self.onedspec.apply_flux_calibration()
        if self.config["display_plots"]:
            self.onedspec.inspect_reduced_spectrum()

    def save_final_spectrum(self):
        output_dir = self.output_dir
        object_name = self.object_name

        # Save just wavelength and flux
        self.onedspec.save_csv(
            stype='science',
            spec_id=0,
            output='wavelength+flux',
            filename=output_dir / object_name,
            overwrite=True
        )

        # Merge wavelength and flux
        wav_path = output_dir / f"{object_name}_science_wavelength.csv"
        flux_path = output_dir / f"{object_name}_science_flux.csv"

        wav = pd.read_csv(wav_path, skiprows=1, names=['wav'])
        flux = pd.read_csv(flux_path, skiprows=1, names=['flux', 'uflux', 'sky'])

        merged = pd.concat([wav, flux['flux'],flux['uflux'], flux['sky']], axis=1)
        merged.to_csv(output_dir / f"{object_name}_final.csv", header=True, index=False)

        # SNID-style space-separated output
        snid_path = output_dir / f"{object_name}_snid.csv"
        snid_file = pd.concat([wav, flux['flux']], axis=1)
        snid_file.to_csv(snid_path, sep=' ', header=True, index=False)


    def plot_final_spectrum(self):
        def box_smoothing(y_data, box_pts):
            box = np.ones(box_pts) / box_pts
            return np.convolve(y_data, box, mode='same')

        # Load y_data
        path = self.output_dir / f"{self.object_name}_final.csv"
        data = pd.read_csv(path, sep=',')
        wav, flux, uflux, sky = data.iloc[:, 0] /10, data.iloc[:, 1], data.iloc[:, 2], data.iloc[:, 3]

        # Define the line wavelengths and names
        lines = [656.279, 486.135, 434.0472, 397.0075, 388.9064, 383.5397, 420, 468.6]
        line_names = ['H-α', 'H-β', 'H-γ', 'H-δ', 'H-ζ', 'H-η', 'He I', 'He II']

        # Plot the spectrum
        plt.figure(figsize=(12, 6))
        plt.plot(wav, box_smoothing(flux, self.smoothing_value), color='black', linewidth=0.75, label='Final Spectrum')
        if self.show_sky:
            plt.plot(wav, sky, color='blue', label='Sky Flux')
        if self.show_error:
            plt.plot(wav, uflux, color='red', label='Error')
        plt.yscale('log')
        plt.xlabel("Wavelength (nm)", fontsize=12)
        plt.ylabel("Flux (arb)", fontsize=12)

        # Plot vertical lines and labels
        for line_wav, name in zip(lines, line_names):
            if wav.min() < line_wav < wav.max():
                # Find the nearest index for labelling height
                idx = (np.abs(wav - line_wav)).argmin()
                y = flux.iloc[idx]
                plt.axvline(line_wav, color='blue', linestyle='--', alpha=0.5)
                plt.text(line_wav + 2, y * 1.2, name, color='blue', fontsize=9, rotation=90, ha='left', va='bottom')

        plt.title(f"Final Spectrum: {self.object_name}", fontsize=14)
        plt.tight_layout()
        plt.legend()
        plt.savefig(self.output_dir / f"{self.object_name}_log_plot.png", dpi=300)
        plt.show()

    def interactive_trim_function(self, tag="science"):

        # Helper function to load raw, flipped data
        def _load_raw(path):
            data,_ = fits.getdata(path, header=True)
            return np.flip(data, axis = 1)

        # Load Raw images
        data = _load_raw(self.science_path)  # 2D reduced science image
        arc = _load_raw(self.arc_path)  # Raw or reduced arc frame (should already be flipped and loaded)

        fig, (ax2d, ax1d_sci, ax1d_arc) = plt.subplots(
            3, 1, figsize=(10, 10), gridspec_kw={'height_ratios': [3, 1, 1]}
        )

        ax2d.imshow(data, cmap='jet', origin='lower', aspect='auto', norm=LogNorm())
        ax2d.set_title(f"Select trimming region for {tag} image")
        ax1d_sci.set_title("1D Science Spectrum (sum over spatial axis)")
        ax1d_arc.set_title("1D Arc Spectrum (sum over spatial axis)")

        coords = {}

        def update_1d_plots(x1, x2, y1, y2):
            # Ensure trimming stays within bounds
            x1, x2 = np.clip([x1, x2], 0, data.shape[1] - 1)
            y1, y2 = np.clip([y1, y2], 0, data.shape[0] - 1)

            # Trim regions
            sci_sub = data[y1:y2, x1:x2]
            arc_sub = arc[y1:y2, x1:x2]

            # Clear previous plots
            ax1d_sci.cla()
            ax1d_arc.cla()

            if sci_sub.size > 0:
                sci_1d = np.sum(sci_sub, axis=0)
                ax1d_sci.plot(sci_1d, color='black')
                ax1d_sci.set_xlim(0, sci_1d.size)

            if arc_sub.size > 0:
                arc_1d = np.sum(arc_sub, axis=0)
                ax1d_arc.plot(arc_1d, color='black')
                ax1d_arc.set_xlim(0, arc_1d.size)

            fig.canvas.draw_idle()

        def onselect(eclick, erelease):
            x1, y1 = int(eclick.xdata), int(eclick.ydata)
            x2, y2 = int(erelease.xdata), int(erelease.ydata)

            x_min, x_max = sorted((x1, x2))
            y_min, y_max = sorted((y1, y2))

            coords.update({
                'x_min': x_min,
                'x_max': x_max,
                'y_min': y_min,
                'y_max': y_max
            })

            print(f"\nSelected region:\n{json.dumps(coords, indent=4)}")
            update_1d_plots(x_min, x_max, y_min, y_max)

        selector = RectangleSelector(
            ax2d, onselect,
            useblit=True,
            button=[1],
            minspanx=5, minspany=5,
            spancoords='pixels',
            interactive=True
        )

        plt.tight_layout()
        plt.show()

        if not coords:
            print("No region selected. Aborting trim update.")
            return

        save = input("Save this region as the new default trimming bounds? (y/n): ").strip().lower()
        if save == 'y':
            defaults_path = "config_files/defaults.json"
            try:
                with open(defaults_path, 'r') as f:
                    defaults = json.load(f)
            except FileNotFoundError:
                print(f"Could not find {defaults_path}.")

            defaults["trim_bounds"] = coords

            with open(defaults_path, 'w') as f:
                json.dump(defaults, f, indent=4)
            print(f"Trim bounds saved to {defaults_path}")
        else:
            print("New trimming bounds not saved. Using for current session only.")

        # update in-memory config
        self.config["trim_bounds"] = coords

    def clean_dir(self):
        """
        Move all intermediate files to a subdirectory 'temp',
        leaving only the final products: final.csv, snid.csv, and log file.
        """
        temp_dir = self.output_dir / "temp"
        temp_dir.mkdir(exist_ok=True)

        # Files we want to keep in the main directory
        final_files = [
            f"{self.object_name}_final.csv",
            f"{self.object_name}_snid.csv",
            f"{self.object_name}_log_file.txt"
        ]

        for file in self.output_dir.iterdir():
            # Skip final files and the temp directory itself
            if file.name in final_files or file == temp_dir:
                continue

            try:
                shutil.move(str(file), str(temp_dir / file.name))
            except Exception as e:
                print(f"Could not move {file.name}: {e}")

        print(f"Intermediate files moved to {temp_dir}, final products remain in {self.output_dir}")

    def run(self):
        self.extract_data()
        log_path = self.output_dir / f"{self.object_name}_log_file.txt"
        sys.stdout = Tee(log_path, sys.__stdout__)
        sys.stderr = Tee(log_path, sys.__stderr__)
        print(f"Logging everything to {log_path}")
        self.bias_subtract()
        self.write_filelists()
        self.reduce_images()
        self.cosmic_clean()
        self.extract_2dspec()
        self.calibrate_wavelength()
        self.calibrate_flux()
        self.save_final_spectrum()
        self.plot_final_spectrum()
        self.clean_dir()