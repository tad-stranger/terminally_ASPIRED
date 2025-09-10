<h1>
  <img src="assets/terminally_ASPIRED_logo_no_text.png" alt="Logo" width="60" style="vertical-align: middle; margin-right: 12px;">
  terminally_ASPIRED
</h1>


 A customizable, command-line-driven spectral reduction pipeline built around the [ASPIRED](https://github.com/cylammarco/ASPIRED) framework — tailored for data from the **SpUpNIC spectrograph** on the **SAAO 1.9m telescope**.
Designed for flexibility and reproducibility, `terminally_ASPIRED` combines robust automation with fine-grained user control via configuration files and interactive tools.
---

## ✨ Features

- Optional **bias and flat field correction**
- Fully configurable via a `defaults.json` file
- Interactive 2D trimming and live preview of spectrum extraction
- Custom **wavelength calibration** with user-defined atlas lines
- Built-in flux calibration and sensitivity inspection
- CLI wrapper for fast processing with reproducible configuration
- Automatic generation of final CSV output (science wavelength vs. flux) in [SNID](https://people.lam.fr/blondin.stephane/software/snid/) - ready format for fast transient classification.

---

## 📦 Installation

### Step 1: Clone and Download
Download or clone this repository, and ensure that you have the following key files:

- The pipeline class:`spectral_reducer.py`
- The command-line interface: `terminally_ASPIRED.py`
- The environment file: `environment.py`
- The `config_files` directory which should contain `defaults.json`. This directory can also be populated with your own presets (`my_preset.json`) for particular telescopes or gratings. 

### Step 2: Set up Conda Environment

We recommend using a virtual environment:
```bash
conda env create -f environment.yml
conda activate terminally_ASPIRED
```

---
## 🚀 CLI Usage
### Basic Usage 
```bash
python terminally_ASPIRED.py science.fits arc.fits standard.fits standard_arc.fits
```
### ⚙️ Required Arguments
| Argument            | Description                                                        |
| ------------------- | ----------------------------------------------------               |
| `science.fits`      | Path to the science target FITS file                               |
| `arc.fits`          | Path to the arc lamp FITS file for the science frame               |
| `standard.fits`     | Path to the standard star FITS file                                |
| `standard_arc.fits` | Path to the arc lamp FITS file for the standard star               |
| `-gr`,`--grating`   | Wavelength atlas for gratings of 1.9m. Options: `7`, `6`, `13` or `custom` |                                                     |


### 🛠️ Optional Arguments
| Flag                 | Description                                                          | Default                      |
| -------------------- | -------------------------------------------------------------------- | ---------------------------- |
| `--config`           | Path to JSON config file                                             | `config_files/defaults.json` |
| `-b`, `--bias`       | Path to directory with bias frames (skip if empty)                   | `""`                         |
| `-f`, `--flat-field` | Path to directory with flat fields (skip if empty)                   | `""`                         |
| `-t`,`--interactive-trim` | Enable interactive trimming of 2D spectra                            | Off                     |
| `--show-plots`       | Show intermediate plots during reduction                             | Off                          |
| `-s`, `--smooth`     | Smoothing box size for final 1D spectrum                             | `1` (no smoothing)           |
| `-v`, `--verbose`    | Enable verbose output                                                | Off                          |
| `--no-warnings`      | Suppress warning messages                                            | Off                          |
| `-O`, `--output-dir` | Custom output directory name (default: object name from FITS header) | `None`                       |
| `--show-sky`         | Show extracted sky spectrum in final plot                            | Off                          |
| `--show-error`       | Show extracted flux error in final plot                              | Off                          |

### ⚠️ Note on Image Reduction
Dark frame correction is not implemented. This is because the CCD used on the 1.9m telescope is cryogenically cooled to temperatures around 170 K, rendering the dark current negligible.

---

## ⚙️ Configuration: `defaults.json`
This pipeline uses a JSON configuration file to control all steps of spectral reduction and calibration. Below is a description of the parameters and their purpose.

### 1. Trim Bounds (`trim_bounds`)
Defines the region of the CCD image to keep after trimming the raw frame.

- `x_min`, `x_max`, `y_min`, `y_max`: pixel boundaries to crop the image.  
  Removes overscan regions and unused parts of the detector.


### 2. Display Plots (`display_plots`)
Boolean flag to toggle plotting of intermediate steps.  
Useful for debugging and visual verification.


### 3. Cosmic Ray Rejection (`cosmic_ray`)
Controls the removal of cosmic ray hits from science and calibration frames using a sigma-clipping algorithm.

- `sigclip`: sigma threshold above which pixels are flagged as cosmic rays.  
- `sigfrac`: fractional detection threshold used in identifying cosmic rays.  
- `objlim`: minimum contrast between object and cosmic ray to distinguish them.


### 4. Atlas Lines (`atlas_lines`)
Reference arc-lamp wavelengths used for wavelength calibration.

- Keys correspond to grating identifiers (`"7"`, `"6"`, `"13"`, `"custom"`).  
- Values are lists of known emission line wavelengths in Ångström for each grating.  
- `"custom"` can be used for user-supplied line lists.


### 5. Standard Star Name (`standard_name`)
Name of the spectrophotometric standard star used for flux calibration (e.g., `ltt6248`).


### 6. Tracing Parameters (`trace_kwargs`)
Controls how the spectrum is traced along the dispersion axis.  
Separate settings exist for **science** and **standard** frames.

- `nspec`: number of spectra to trace in the frame.  
- `nwindow`: half-window size (in pixels) used when locating the spectrum.  
- `trace_width`: width of the spectral trace in pixels.  
- `resample_factor`: factor by which the trace is oversampled before fitting.  
- `fit_deg`: degree of the polynomial used to fit the trace.  
- `display`: whether to show the trace fitting plot.  
- `width`, `height`: dimensions of the input image (used for normalization).


### 7. Extraction Parameters (`extract_kwargs`)
Controls how the 1D spectrum is extracted from the 2D frame.  
Separate settings exist for **science** and **standard** frames.

- `apwidth`: half-width of the extraction aperture (in pixels).  
- `optimal`: whether to use optimal extraction (weights by signal-to-noise).  
- `algorithm`: extraction algorithm, e.g., `"horne86"` (Horne 1986 optimal extraction).  
- `skysep`: minimum distance (in pixels) between spectrum and sky regions.  
- `skywidth`: width of sky background region (in pixels).  
- `skydeg`: polynomial degree for fitting the sky background.  
- `display`: whether to show the extraction step.


### 8. Hough Transform Parameters (`hough`)
Controls automated line identification using a Hough transform.

- `num_slopes`: number of slope values sampled in the transform.  
- `xbins`, `ybins`: binning for the transform parameter space.  
- `min_wavelength`, `max_wavelength`: expected wavelength range (Å).  
- `range_tolerance`: tolerance (Å) when matching atlas lines.


### 9. Wavelength Calibration (`wavelength_cal`)
Controls how arc lamp lines are matched and fitted to derive the wavelength solution.  
Separate settings exist for **science** and **standard** frames.

- `prominence`: required prominence of peaks for line detection.  
- `refine`: whether to iteratively refine the wavelength solution.  
- `ransac_minimum_matches`: minimum number of line matches required for RANSAC fitting.  
- `fit_max_tries`: maximum number of fitting attempts before failing.  
- `fit_deg`: polynomial degree of the wavelength solution.  
- `fit_type`: type of fit (e.g., `"poly"` for polynomial).

---
## 📤 Output
Several output files are created, but there are three main output files that are important:

- `object_name_final.csv` - final fully reduced csv file containing wavelength (Å), flux, flux error, sky flux
- `object_name_snid.csv` - SNID ready formatted csv file (space-seperated) for quick transient identification
- `object_name_log_file.txt` - log file containing all terminal output from ASPIRED. If run in verbose mode, you will get all of the useful statistics ASPIRED produces saved in this file. 

---

## 🧪 Developer Notes
This pipeline is built on the excellent ASPIRED library by Marco Lam. terminally_ASPIRED wraps its lower-level functionality to streamline workflows, reduce bugs, and standardize output for further scientific use.

---
## ☄️Acknowledgements
Developed by Francois Campher and Lloyd Landsberg as part of our Masters' Dissertations. We aim to provide a useful tool for quick spectral reduction for the transients and variable stars research teams within the BlackGEM and MeerLICHT consortia, The University of Cape Town (UCT) and the South African Astronomical Observatory (SAAO). We would like to also thank the exellent developers of the RASCAL and ASPIRED packages.

---



