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
import lineid_plot

# This is a class of my 1.9M pipeline to be used to make terminally_ASPIRED
class SpectralReductionPipeline:
    def __init__(self, science_file, arc_file, std_file, config_path="config_files/defaults.json"):
        self.science_path = Path(science_file)
        self.arc_path = Path(arc_file)
        self.std_path = Path(std_file)

        # Deduce base observation directory (e.g., ".../0503")
        self.obs_base = self.science_path.parents[2]

        # BIASS and FLATS directories are at the same level
        self.bias_folder = self.obs_base / "BIASS"
        self.flats_folder = self.obs_base / "FLATS"

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
            flat_files = sorted(self.flats_folder.glob("*.fits"))
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

        self.onedspec.find_arc_lines(prominence=7, refine=True, display=True, stype='science+standard')
        self.onedspec.initialise_calibrator(stype='science+standard')
        self.onedspec.add_user_atlas(elements=element, wavelengths=atlas, stype='science+standard')
        self.onedspec.set_hough_properties(**self.config["hough"])
        self.onedspec.do_hough_transform(stype='science+standard')
        self.onedspec.set_ransac_properties(minimum_matches=10)
        self.onedspec.fit(max_tries=1000, fit_deg=3, fit_type='poly', display=True, stype='science+standard')
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
        path = self.output_dir / f"{self.object_name}_final.csv"
        data = pd.read_csv(path, sep = ' ')
        wav, flux = data.iloc[:, 0], data.iloc[:, 1]
        plt.plot(wav, flux, label="Final Spectrum")
        plt.xlabel("Wavelength (Å)")
        plt.ylabel("Flux (arb)")
        plt.legend()
        plt.show()

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



