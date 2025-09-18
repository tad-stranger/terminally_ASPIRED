<h1 align="left">
  <img src="assets/terminally_ASPIRED_logo_no_text.png" alt="Logo" width="80" style="vertical-align: middle; margin-right: 12px;">
  terminally_ASPIRED
</h1>

<p align="left">
  <a href="https://pypi.org/project/terminally-aspired/">
    <img src="https://img.shields.io/pypi/v/terminally-aspired?color=blue&label=PyPI" alt="PyPI version">
  </a>
  <a href="https://github.com/tad-stranger/terminally_ASPIRED/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/tad-stranger/terminally_ASPIRED" alt="License">
  </a>
  <img src="https://img.shields.io/badge/python-3.11%2B-blue" alt="Python">
  <a href="https://doi.org/10.5281/zenodo.17153499">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17153499.svg" alt="DOI">
  </a>
</p>

 A customizable, command-line-driven spectral reduction pipeline built around the [ASPIRED](https://github.com/cylammarco/ASPIRED) framework — tailored for data from the **SpUpNIC spectrograph** on the **SAAO 1.9m telescope**.
Designed for flexibility and reproducibility, `terminally_ASPIRED` combines robust automation with fine-grained user control via configuration files and interactive tools.
---

## ✨ Features

- Fully configurable via a `config.json` file
- Interactive 2D trimming and live preview of spectrum extraction
- Includes **wavelength atlas** for gratings 6 and 7 on the 1.9m
- Custom **wavelength calibration** with user-defined atlas lines
- Built-in flux calibration and sensitivity inspection
- CLI wrapper for fast processing with reproducible configuration
- Automatic generation of final CSV output (science wavelength vs. flux) in [SNID](https://people.lam.fr/blondin.stephane/software/snid/) - ready format for fast transient classification.

---

## 📦 Installation and Setup


### Recommended: Conda environment 🐍

#### Step 1: Create a Virtual Environment

```bash
conda create -n terminally_ASPIRED python=3.11
```

#### Step 2: Activate the Virtual Environment

```bash
conda activate terminally_ASPIRED
```

#### Step 3: Download the package from PyPI

```bash
pip install --upgrade setuptools jmespath -i https://pypi.org/simple && \
pip install terminally-aspired
```

#### Step 4: Verify the Installation

```bash
tA --help
```

#### Step 5: Copy `defaults.json` as template and save in working directory 

The main way you will interact with terminally_ASPIRED is through the config `.json` file, so you should copy the `defaults.json` file to use as a template and rename it to whatever is convenient. 
You can use a built-in function in terminally_ASPIRED to do this: 
```bash
make-config -n name_of_config_file.json -O directory_to_store_config
```
This will create a copy of the defaults.json file and rename it to what you choose (remeber to include the .json extension) as well as save it to a output directory of your choice. Then when calling terminally_ASPIRED, you can point to your own config file by making use of the `--config` argument. 


### Alternative: Python venv 🫙

#### Step 1: Create a Virtual Environment (Ensure that python >= 3.11)
If you wish to use a python venv, than make sure that the python version on your device is >= 3.11
```bash
python3 -m venv ~/.venv/terminally_ASPIRED
```

#### Step 2: Activate the Virtual Environment

```bash
source ~/.venv/terminally_ASPIRED/bin/activate
```

#### Step 3: Download the package from PyPI

```bash
pip install --upgrade setuptools jmespath -i https://pypi.org/simple && \
pip install terminally-aspired
```

#### Step 4: Verify the Installation

To make sure that the package installed correctly, run the --help argument. It sometimes takes a few minutes to run on a fresh install.
```bash
tA --help
```


#### Step 5: Copy `defaults.json` as template and save in working directory 

The main way you will interact with terminally_ASPIRED is thorugh the config `.json` file, so you should copy the `defaults.json` file to use as a template and rename it to whatever is convienient. 
You can use a built-in function in terminally_ASPIRED to do this: 
```bash
make-config -n name_of_config_file -O directory_to_store_config
```
This will create a copy of the defaults.json file and rename it to what you choose (remeber to include the .json extention) as well as save it to a output directory of your choice. Then when calling terminally_ASPIRED, you can point to your own config file by making use of the `--config` argument. 

### Alternative(Not Reccomended): Clone 

#### Step 1: Clone and Download
Download or clone this repository, and ensure that you have the following key files:

- The pipeline class:`spectral_reducer.py`
- The command-line interface: `terminally_ASPIRED.py`
- The environment file: `environment.py`
- The `config_files` directory which should contain `defaults.json`

#### Step 2: Set up Conda Environment

We recommend using a virtual environment:
```bash
conda env create -f environment.yml
conda activate terminally_ASPIRED
```

#### Step 3: Calls to terminally_ASPIRED.py
If you prefer to download the files like this, then you will need to call terminally_ASPIRED in the following way:
```bash
python terminally_ASPIRED.py arguments

```

---
## 🚀 CLI Usage
### Basic Usage 
```bash
tA science.fits arc.fits standard.fits standard_arc.fits -gr 7 -f flat_dir -b bias_dir --config path/to/config.json
```
This will call terminally_ASPIRED and reduce the science.fits file for data from grating 7. This will work well once you have gotten some settings in the config files just right for this particular set of data.

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
## 🪁Example - First time reduction with a new set of data

Often, when you have a fresh set of data you will need to carefully monitor how the pipeline deals with trace extraction & wavelength and flux calibration. Therefore for a first run with new data we reccomend to make use of the `--show-plots` argument which will show you the ASPIRED intermediate plots. 

### Step 1: Standard Star name in `config.json`

Remember to input your standard star name in lower case letters into your `config.json` file under the parameter `standard_name` (examples: hilt600, ltt6248 etc)

### Step 2: Run the pipeline with `--show-plots` enabled

```bash
tA science.fits arc.fits standard.fits standard_arc.fits -gr 7 -f flat_dir -b bias_dir --config path/to/config.json --show-plots
```

### Step 3: Wavelength Calibration 
A particular sticky point is typically wavelength calibration, due to the way in which the auto-detection of the peaks of the arc spectrum works. You can control the wavelength calibration in the `config.json` file, located under the header `hough` and  `wavelength_cal`. 

#### `hough`
Here you need to pay attention to `min_wavelength` and `max_wavelength`, as these will depend on the particular grating and grating angle you use. For:
- grating 6 : Use `min_wavelength` = 4500  , `max_wavelength` = 7000
- grating 7 : Use `min_wavelength` = 3500  , `max_wavelength` = 9000


#### `wavelength_cal`
Here you need to pay attention to the `prominence` parameter. This sets the level of relative flux above which the auto-peak detection will detect peaks in the arc spectrum. You must make sure that the `prominence` level is such that the peaks that are auto detected match up with the peaks within the wavelength atlas of the chosen grating. You can find the wavelength atlases for the different grating on the [SAAO SpUpNIC TOPS Wiki](https://topswiki.saao.ac.za/index.php/SPUPNIC). The pipeline by default contains atlases for gratings 6 and 7 (and soon 13), but you still need to make sure that they line up with the detected peaks correctly. You can also add your own lines to these atlases. Once this is done the pipeline should do the wavelength calibration quickly and well. 

### Step 4: Inspect output
Make sure that the spectrum that is returned looks reasonable. You can plot the error and sky flux by making use of the `--show-error` and `--show-sky` arguments for further confirmation. Alternatively you can use the interactive html plots produced by ASPIRED from the `--show-plots` argument.

### Step 5: Reduce the rest of the data.
Once you are happy with a particular spectrum from a night, you can process the rest of them by making use of a simple call to terminally_ASPIRED: 
```bash
tA science.fits science_arc.fits standard.fits standard_arc.fits -gr 7 --config config.json -f flat_dir -b bias_dir

```
This will simply process the images and plot the final spectrum for inspection. 

---

## ⚙️ Configuration: `config.json`
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
- `object_name_snid.csv` - SNID ready formatted csv file (space-separated) for quick transient identification
- `object_name_log_file.txt` - log file containing all terminal output from ASPIRED. If run in verbose mode, you will get all of the useful statistics ASPIRED produces saved in this file. 

---

## 📖 Citation  

If you use **`terminally_ASPIRED`** in your work, please make sure to cite both **ASPIRED** and **terminally_ASPIRED**.  

- **First, cite ASPIRED** (the underlying spectral reduction framework):  

```bibtex
@ARTICLE{2023AJ....166...13L,
       author = {{Lam}, Marco C. and {Smith}, Robert J. and {Arcavi}, Iair and {Steele}, Iain A. and {Veitch-Michaelis}, Josh and {Wyrzykowski}, Lukasz},
        title = "{Automated SpectroPhotometric Image REDuction (ASPIRED)}",
      journal = {AJ},
     year = 2023,
     volume = {166},
     number = {1},
     eid = {13},
     pages = {13},
     doi = {10.3847/1538-3881/acd75c}
}
```

- **Then, cite terminally_ASPIRED** (this pipeline):
  
```bibtex
@software{campher_landsberg2025terminallyASPIRED,
  author       = {Francois Campher and Lloyd Landsberg},
  title        = {{terminally_ASPIRED: A command-line spectral reduction pipeline built on ASPIRED}},
  year         = {2025},
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {10.5281/zenodo.17153500},
  url          = {https://doi.org/10.5281/zenodo.17153499}
}
```

## 🧪 Developer Notes
This pipeline is built on the excellent ASPIRED toolkit by Marco Lam and Robert Smith. terminally_ASPIRED wraps its lower-level functionality to streamline workflows, reduce bugs, and standardize output for further scientific use.

---
## ☄️Acknowledgements
Developed by Francois Campher and Lloyd Landsberg as part of our Masters' Dissertations. We aim to provide a useful tool for quick spectral reduction for the transients and variable stars research teams within the BlackGEM and MeerLICHT consortia, The University of Cape Town (UCT) and the South African Astronomical Observatory (SAAO). We would like to also thank the excellent developers of the RASCAL and ASPIRED packages.

---



