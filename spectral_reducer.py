import json
from pathlib import Path
import pandas as pd
import numpy as np
from astropy.io import fits
from aspired import image_reduction, spectral_reduction
from astroscrappy import detect_cosmics
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.widgets import RectangleSelector

# This is a class of my 1.9M pipeline to be used to make terminally_ASPIRED
class SpectralReductionPipeline:
    def __init__(self, science_file, arc_file, std_file, config_path="config_files/defaults.json", use_bias_flats = False, bias_path = None, flat_path = None):
        self.science_path = Path(science_file)
        self.arc_path = Path(arc_file)
        self.std_path = Path(std_file)

        # Deduce base observation directory (e.g., ".../0503")
        self.obs_base = self.science_path.parents[2]

        if use_bias_flats and bias_path and flat_path:
            self.bias_folder = Path(bias_path)
            self.flat_folder = Path(flat_path)
        else:
            self.bias_folder = Path("DO_NOT_USE_BIAS")
            self.flat_folder = Path("DO_NOT_USE_FLATS")

        self.config = self._load_config(config_path)

        # Placeholder attributes
        self.object_name = None
        self.output_dir = None
        self.cleaned = {}

        # FITS data
        self.sci_data = None
        self.arc_data = None
        self.std_data = None
        self.hdr_sci = None
        self.hdr_arc = None
        self.hdr_std = None

    def _load_config(self, path):
        with open(path, 'r') as f:
            return json.load(f)

    def extract_data(self):
        def _extract(path):
            data, header = fits.getdata(path, header=True)
            return np.flip(data, axis=1), header  # flip for ASPIRED format

        sci_path = self.science_path
        arc_path = self.arc_path
        std_path = self.std_path

        self.sci_data, self.hdr_sci = _extract(sci_path)
        self.arc_data, self.hdr_arc = _extract(arc_path)
        self.std_data, self.hdr_std = _extract(std_path)

        self.object_name = self.hdr_sci.get("OBJECT", "Unknown").replace(" ", "")
        self.output_dir = Path(f"./ReducedSpectra/{self.object_name}")
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def write_filelists(self):
        def write(filelist_path, science_file):
            bias_files = sorted(self.bias_folder.glob("*.fits"))
            flat_files = sorted(self.flat_folder.glob("*.fits"))
            with open(filelist_path, "w") as f:
                for file in bias_files:
                    f.write(f"bias, {file.resolve()}\n")
                for file in flat_files:
                    f.write(f"flat, {file.resolve()}\n")
                f.write(f"arc, {self.arc_file_resolved}\n")
                f.write(f"light, {science_file}\n")

        # Save FITS files first
        self.arc_file_resolved = self._save_fits(self.arc_data, self.hdr_arc, f"{self.object_name}_arc.fits")
        sci_resolved = self._save_fits(self.sci_data, self.hdr_sci, f"{self.object_name}_sci.fits")
        std_resolved = self._save_fits(self.std_data, self.hdr_std, f"{self.object_name}_std.fits")

        write(self.output_dir / f"{self.object_name}_science_file.list", sci_resolved)
        write(self.output_dir / f"{self.object_name}_standard_file.list", std_resolved)

    def _save_fits(self, data, header, filename):
        path = (self.output_dir / filename).resolve()
        fits.writeto(path, data, header, overwrite=True)
        return path

    def reduce_images(self):
        def reduce(list_file, tag):
            ir = image_reduction.ImageReduction()
            ir.add_filelist(str(list_file))
            ir.load_data()
            ir.reduce()
            ir.save_fits(str(self.output_dir / f"{self.object_name}_{tag}_image_reduced.fits"), overwrite=True)

        reduce(self.output_dir / f"{self.object_name}_science_file.list", tag="science")
        reduce(self.output_dir / f"{self.object_name}_standard_file.list", tag="standard")

    def trim_and_clean(self):
        bounds = self.config["trim_bounds"]

        def trim(data):
            return data[
                   bounds["y_min"]:bounds["y_max"],
                   bounds["x_min"]:bounds["x_max"]
                   ]

        def clean(data):
            _, cleaned = detect_cosmics(data, **self.config["cosmic_ray"])
            return cleaned

        self.cleaned["sci"] = clean(trim(self._load_image("science")))
        self.cleaned["arc"] = clean(trim(self.arc_data))
        self.cleaned["std"] = clean(trim(self._load_image("standard")))

    def _load_image(self, tag):
        data, _ = fits.getdata(self.output_dir / f"{self.object_name}_{tag}_image_reduced.fits", header=True)
        return data

    def extract_2dspec(self):
        self.sci2d = self._extract_twodspec(self.cleaned["sci"], self.cleaned["arc"], self.hdr_sci, is_standard=False)
        self.std2d = self._extract_twodspec(self.cleaned["std"], self.cleaned["arc"], self.hdr_std, is_standard=True)

    def _extract_twodspec(self, data, arc, header, is_standard):
        twod = spectral_reduction.TwoDSpec(data, header=header, cosmicray=False)
        twod.set_properties(saxis=1, flip=False, cosmicray=False)
        twod.ap_trace(**self.config["trace_kwargs"]["standard" if is_standard else "science"])
        twod.ap_extract(**self.config["extract_kwargs"]["standard" if is_standard else "science"])
        twod.add_arc(arc=arc)
        twod.extract_arc_spec()
        return twod

    def calibrate_wavelength(self):
        self.onedspec = spectral_reduction.OneDSpec()
        self.onedspec.from_twodspec(self.sci2d, stype="science")
        self.onedspec.from_twodspec(self.std2d, stype="standard")

        atlas = self.config["atlas_lines"]
        element = ["CuAr"] * len(atlas)
        cal_cfg = self.config["wavelength_cal"]

        self.onedspec.find_arc_lines(
            prominence=cal_cfg["prominence"],
            refine=cal_cfg["refine"],
            display=True,
            stype='science+standard'
        )

        self.onedspec.initialise_calibrator(stype='science+standard')
        self.onedspec.add_user_atlas(
            elements=element,
            wavelengths=atlas,
            stype='science+standard'
        )

        self.onedspec.set_hough_properties(**self.config["hough"])

        self.onedspec.do_hough_transform(stype='science+standard')

        self.onedspec.set_ransac_properties(
            minimum_matches=cal_cfg["ransac_minimum_matches"]
        )

        self.onedspec.fit(
            max_tries=cal_cfg["fit_max_tries"],
            fit_deg=cal_cfg["fit_deg"],
            fit_type=cal_cfg["fit_type"],
            display=True,
            stype='science+standard'
        )

        self.onedspec.apply_wavelength_calibration(stype='science+standard')

    def calibrate_flux(self):
        std_name = self.config.get("standard_name", "hilt600")
        self.onedspec.load_standard(target=std_name)
        self.onedspec.inspect_standard()
        self.onedspec.get_sensitivity()
        self.onedspec.inspect_sensitivity()
        self.onedspec.apply_flux_calibration()

    def save_final_spectrum(self):
        output_dir = self.output_dir
        object_name = self.object_name

        # Save just wavelength + flux
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

        merged = pd.concat([wav, flux['flux']], axis=1)
        merged.to_csv(output_dir / f"{object_name}.csv", header=True, index=False)

        # SNID-style space-separated output
        snid_path = output_dir / f"{object_name}_final.csv"
        snid_file = pd.concat([wav, flux['flux']], axis=1)
        snid_file.to_csv(snid_path, sep=' ', header=True, index=False)


    def plot_final_spectrum(self):
        # Load data
        path = self.output_dir / f"{self.object_name}_final.csv"
        data = pd.read_csv(path, sep=' ')
        wav, flux = data.iloc[:, 0] /10, data.iloc[:, 1]

        # Define the line wavelengths and names
        lines = [656.279, 486.135, 434.0472, 397.0075, 388.9064, 383.5397, 420, 440]
        line_names = ['H-α', 'H-β', 'H-γ', 'H-δ', 'H-ζ', 'H-η', 'He I', 'He II']

        # Plot the spectrum
        plt.figure(figsize=(12, 6))
        plt.plot(wav, flux, color='black', linewidth=0.75, label='Final Spectrum')
        plt.yscale('log')
        plt.xlabel("Wavelength (nm)", fontsize=12)
        plt.ylabel("Flux (arb)", fontsize=12)

        # Plot vertical lines and labels
        for line_wav, name in zip(lines, line_names):
            if wav.min() < line_wav < wav.max():
                # Find nearest index for labeling height
                idx = (np.abs(wav - line_wav)).argmin()
                y = flux.iloc[idx]
                plt.axvline(line_wav, color='blue', linestyle='--', alpha=0.5)
                plt.text(line_wav + 2, y * 1.2, name, color='blue', fontsize=9, rotation=90, ha='left', va='bottom')

        plt.title(f"Final Spectrum: {self.object_name}", fontsize=14)
        plt.tight_layout()
        plt.legend()
        plt.show()

    def interactive_trim(self, tag="science"):
        data = self._load_image(tag)  # 2D reduced science image
        arc = self.arc_data  # Raw or reduced arc frame (should already be flipped and loaded)

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

        #plt.tight_layout()
        plt.show()

        if not coords:
            print("No region selected. Aborting trim update.")
            return

        save = input("Save this region as the new default trimming bounds? (y/n): ").strip().lower()
        if save == 'y':
            self.config["trim_bounds"] = coords
            save_path = "config_files/trim_bounds.json"
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            with open(save_path, 'w') as f:
                json.dump(coords, f, indent=4)
            print(f"New trimming bounds saved to {save_path}")
        else:
            print("New trimming bounds not saved. Using for current session only.")

        # Always update in-memory config
        self.config["trim_bounds"] = coords

    def run(self):
        self.extract_data()
        self.write_filelists()
        self.reduce_images()
        self.trim_and_clean()
        self.extract_2dspec()
        self.calibrate_wavelength()
        self.calibrate_flux()
        self.save_final_spectrum()
        self.plot_final_spectrum()

    def run_with_interactive_trim(self, tag = "science"):
        self.extract_data()
        self.write_filelists()
        self.reduce_images()
        self.interactive_trim(tag=tag)
        self.trim_and_clean()
        self.extract_2dspec()
        self.calibrate_wavelength()
        self.calibrate_flux()
        self.save_final_spectrum()
        self.plot_final_spectrum()