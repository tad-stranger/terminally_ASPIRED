# terminally_ASPIRED 🔭

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
- Automatic generation of final CSV output (science wavelength vs. flux)

---

## 📦 Installation

### Step 1: Clone and Download
Download or clone this repository, and ensure that you have the following key files:

- The pipeline class:`spectral_reducer.py`
- The command-line interface: `terminally_ASPIRED.py`
- The environment file: `environment.py`
- The `config_files` directory which should contain both `defaults.json` & `trim_bounds.json`

### Step 2: Set up Conda Environment

We recommend using a virtual environment:
```bash
conda env create -f environment.yml
conda activate terminally_ASPIRED
```
## ⚙️ Configuration: `defaults.json`
This file controls all pipeline behaviour.
Key sections include:
- **Trimming Bounds**(`trim_bounds`): pixel coordinates for cropping the 2D image
- **Cosmic ray cleaning**(`cosmic_ray`): parameters passed to AstroScrappy
- **Arc line atals**(`arc_lines`): wavelengths for calibration.
- **Extraction parameters**(`extract_kwargs`): per-object optimal extraction.
- **Tracing config**(`trace_kwargs`): for locating spectrum traces.
- **Wavelength solution fitting** (`wavelength_cal`, `hough`)
---

## 🚀 CLI Usage

---
## 📤 Output

---

## 🧪 Developer Notes
This pipeline is built on the excellent ASPIRED library by Marco Lam. terminally_ASPIRED wraps its lower-level functionality to streamline workflows, reduce bugs, and standardize output for further scientific use.

---
## ☄️Acknowledgements

Developed by Francois Campher and Lloyd Landsberg as part of our Masters' Dissertations. We aim to provide a useful tool for quick spectral reduction for the transients research team at The University of Cape Town (UCT)
and the South African Astronomical Observatory (SAAO). We would like to also thank the exellent developers of the RASCAL and ASPIRED packages.

---
