import os

from astropy.io import fits
from aspired import image_reduction, spectral_reduction
from astroscrappy import detect_cosmics
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import plotly.io as pio
from pathlib import Path
import pandas as pd
import glob
import lineid_plot

def extract_1D_spectrum(science_data):
    spectrum_1d = np.sum(science_data, axis = 0)
    return spectrum_1d

def extract_data(filename):
    hdul = fits.open(filename)
    header = hdul[0].header
    data = hdul[0].data
    hdul.close()
    return data, header

def plot_image(image_data, keyword, object_name):
    plt.imshow(image_data, cmap='jet', origin='lower', aspect='auto', norm=LogNorm())
    plt.colorbar()
    plt.title(f"{keyword} Spectrum of {object_name}")
    plt.xlabel("Pixel (Dispersion)")
    plt.ylabel("Spatial Axis")
    plt.show()

def plot_spectrum(spec_data, keyword, object_name):
    plt.plot(spec_data, label=f"{keyword} Spectrum")
    plt.xlabel("Pixel (Dispersion)")
    plt.ylabel("Counts (arb)")
    plt.title(f"{keyword} Spectrum of {object_name}")
    plt.legend()
    plt.show()

def write_filelist_for_aspired(bias_dir, flat_dir, science_file, arc_file, output_name="aspired_filelist.list"):
    '''This function loads the science and calibration frames in a filelist
    in a format which ASPIRED can read'''


    bias_files = sorted(glob.glob(f"{bias_dir}/*.fits"))
    flat_files = sorted(glob.glob(f"{flat_dir}/*.fits"))

    with open(output_name, "w") as f:
        for file in bias_files:
            f.write(f"bias, {Path(file).resolve()}\n")
        for file in flat_files:
            f.write(f"flat, {Path(file).resolve()}\n")
        f.write(f"arc, {Path(arc_file).resolve()}\n")
        f.write(f"light, {Path(science_file).resolve()}\n")

    print(f"Filelist written to {output_name}")

def plot_csv_final_spectrum(filename, scale):
    data = pd.read_csv(filename, sep = ' ')
    lines = [656.279, 486.135, 434.0472, 397.0075, 388.9064, 383.5397, 420, 440]
    line_names = ['H-α', 'H-β', 'H-γ', 'H-δ', 'H-ζ', 'H-η', 'HeI', 'HeII']

    wavelength = data['wav'] / 10
    flux = data['flux']


    lineid_plot.plot_line_ids(wavelength, flux, lines, line_names, color='blue', linewidth=0.75)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Flux [arb]')
    if scale == 'log':
        plt.yscale('log')
    plt.title(f'{object_name} Reduced Spectrum', pad=60)
    plt.savefig(f'{output_dir}/{object_name}_final_{scale}.png')
    plt.show()

def read_in_data(observations_path, science_file_name, arc_file_name, std_file_name):
    sci_data, hdr_sci = extract_data(f"{observations_path}SCIENCE/{science_file_name}")
    arc_data, hdr_arc = extract_data(f"{observations_path}ARC/{arc_file_name}")
    std_data, std_hdr = extract_data(f"{observations_path}SCIENCE/{std_file_name}")

    sci_data = np.flip(sci_data, axis=1)
    arc_data = np.flip(arc_data, axis=1)
    std_data = np.flip(std_data, axis=1)

    object_name = hdr_sci.get('OBJECT', 'Unknown')
    tar_RA = hdr_sci.get('TARG-RA')
    tar_DEC = hdr_sci.get('TARG-DEC')

    return sci_data, hdr_sci, arc_data, hdr_arc, std_data, std_hdr, object_name, tar_RA, tar_DEC

def visualise_raw_data(sci_data, object_name):
    plot_image(sci_data, keyword='Raw', object_name=object_name)
    spectrum_1d = extract_1D_spectrum(sci_data)
    plot_spectrum(spectrum_1d, keyword='Raw', object_name=object_name)

def prepare_output_directory(object_name):
    output_dir = Path(f"./ReducedSpectra/{object_name.replace(' ', '')}")
    output_dir.mkdir(exist_ok=True, parents=True)
    return output_dir

def save_fits_files(sci_data, hdr_sci, arc_data, hdr_arc, std_data, std_hdr, output_dir, object_name):
    science_file = output_dir / f"{object_name}_science_reduced.fits"
    arc_file = output_dir / f"{object_name}_arc_reduced.fits"
    std_file = output_dir / f"{object_name}_std_reduced.fits"

    fits.writeto(science_file, sci_data, hdr_sci, overwrite=True)
    fits.writeto(arc_file, arc_data, hdr_arc, overwrite=True)
    fits.writeto(std_file, std_data, std_hdr, overwrite=True)

    return science_file, arc_file, std_file

def create_aspired_filelists(bias_folder, flats_folder, science_file, arc_file, std_file, object_name):
    write_filelist_for_aspired(bias_dir=bias_folder, flat_dir=flats_folder,
                               science_file=science_file, arc_file=arc_file,
                               output_name=f"{object_name}_science_reduction_file.list")

    write_filelist_for_aspired(bias_dir=bias_folder, flat_dir=flats_folder,
                               science_file=std_file, arc_file=arc_file,
                               output_name=f"{object_name}_standard_reduction_file.list")

def reduce_images_with_aspired(object_name, output_directory):
    science_frame = image_reduction.ImageReduction()
    science_frame.add_filelist(f"{object_name}_science_reduction_file.list")
    science_frame.load_data()
    science_frame.reduce()
    science_frame.inspect(display=True, width=1024, height=1024)
    science_frame.save_fits(filename=str(output_directory/f"{object_name}_science_image_reduced.fits"), overwrite=True)

    standard_frame = image_reduction.ImageReduction()
    standard_frame.add_filelist(f"{object_name}_standard_reduction_file.list")
    standard_frame.load_data()
    standard_frame.reduce()
    standard_frame.inspect(display=True, width=1024, height=1024)
    standard_frame.save_fits(filename=str(output_directory/f"{object_name}_standard_image_reduced.fits"), overwrite=True)

def trim_image(data, x_min, x_max, y_min, y_max):
    return data[y_min:y_max, x_min:x_max]

def clean_cosmic_rays(data, sigclip=4.5, sigfrac=0.3, objlim=5.0):
    _, cleaned = detect_cosmics(data, sigclip=sigclip, sigfrac=sigfrac, objlim=objlim)
    return cleaned

def extract_twodspec(data, arc_data, header, n_sources=1, is_standard=False, display=True):
    twod = spectral_reduction.TwoDSpec(data, cosmicray=False, header=header)
    twod.set_properties(saxis=1, flip=False, cosmicray=False)

    trace_kwargs = dict(
        nspec=n_sources,
        nwindow=20 if is_standard else 100,
        trace_width=5,
        resample_factor=2 if is_standard else 15,
        fit_deg=2 if is_standard else 3,
        display=False if is_standard else True,
        width=1024,
        height=1024,
    )
    twod.ap_trace(**trace_kwargs)

    extract_kwargs = dict(
        apwidth=5 if is_standard else 3,
        optimal=True,
        algorithm='marsh89',
        skysep=7 if is_standard else 2,
        skywidth=10 if is_standard else 5,
        skydeg=1,
        display=True,
    )
    twod.ap_extract(**extract_kwargs)

    twod.add_arc(arc=arc_data)
    twod.extract_arc_spec()
    return twod

def wavelength_calibration(science2D, standard2D, atlas, element):
    onedspec = spectral_reduction.OneDSpec()
    onedspec.from_twodspec(science2D, stype='science')
    onedspec.from_twodspec(standard2D, stype='standard')
    onedspec.find_arc_lines(prominence=7, refine=True, display=True, stype='science+standard')

    onedspec.initialise_calibrator(stype='science+standard')
    onedspec.add_user_atlas(elements=element, wavelengths=atlas, stype='science+standard')
    onedspec.set_hough_properties(
        num_slopes=1000, xbins=100, ybins=100,
        min_wavelength=3500, max_wavelength=9000,
        range_tolerance=250, stype='science+standard'
    )
    onedspec.do_hough_transform(stype='science+standard')
    onedspec.set_ransac_properties(minimum_matches=10)
    onedspec.fit(max_tries=1000, fit_deg=3, fit_type='poly', display=True, stype='science+standard')
    onedspec.apply_wavelength_calibration(stype='science+standard')

    return onedspec

def perform_flux_calibration(onedspec, standard_name='hilt600'):
    standard_name = str.lower(standard_name.replace(' ', ''))
    onedspec.load_standard(target = standard_name)
    onedspec.inspect_standard()
    onedspec.get_sensitivity()
    onedspec.inspect_sensitivity()
    onedspec.apply_flux_calibration()

def save_final_spectra(onedspec, object_name, output_dir):
    output_dir = Path(output_dir)

    # Save full CSVs
    onedspec.save_csv(stype='science', spec_id=0, output='*',
                      filename=output_dir / object_name, overwrite=True)

    # Save just wavelength + flux
    onedspec.save_csv(stype='science', spec_id=0, output='wavelength+flux',
                      filename=output_dir / object_name, overwrite=True)

    # Merge and save
    wav = pd.read_csv(output_dir / f"{object_name}_science_wavelength.csv", skiprows=1, names=['wav'])
    flux = pd.read_csv(output_dir / f"{object_name}_science_flux.csv", skiprows=1, names=['flux', 'uflux', 'sky'])

    merged = pd.concat([wav, flux['flux']], axis=1)
    merged.to_csv(output_dir / f"{object_name}.csv", header=True, index=False)

    # SNID-style space-separated output
    snid_file = pd.concat([wav, flux['flux']], axis=1)
    snid_file.to_csv(output_dir / f"{object_name}_final.csv", sep=' ', header=True, index=False)

    return output_dir / f"{object_name}_final.csv"

def plot_final_spectrum(csv_path):
    plot_csv_final_spectrum(csv_path, scale='linear')
    plot_csv_final_spectrum(csv_path, scale='log')


pio.renderers.default = "browser"

# Set Observation Path
observations_path = "Transients Observations/0503/"
science_file_name = "Gaia3513/a6281106.fits"
arc_file_name = "a6281115.fits"
std_file_name = "HILT600/a6281069.fits"

bias_folder = f"{observations_path}BIASS"
flats_folder = f"{observations_path}FLATS"

''' 1 - Extract Data'''
# Extract Data
sci_data, hdr_sci, arc_data, hdr_arc, std_data, std_hdr, object_name, tar_RA, tar_DEC = read_in_data(observations_path, science_file_name, arc_file_name, std_file_name)

# Visualise Raw Image and Spectrum
visualise_raw_data(sci_data, object_name)

# Prepare Output Directory
output_dir = prepare_output_directory(object_name)

'''2 - Image Reduction ASPIRED'''

# Save FITS Files
science_file, arc_file, std_file = save_fits_files(
    sci_data, hdr_sci, arc_data, hdr_arc, std_data, std_hdr, output_dir, object_name
)

# Create File Lists for ASPIRED
create_aspired_filelists(bias_folder, flats_folder, science_file, arc_file, std_file, object_name.replace(" ", ""))

# Reduce with ASPIRED
reduce_images_with_aspired(object_name.replace(" ", ""), output_directory=output_dir)

'''3 - Trimming and Cosmic Ray Removal'''
# Trimming Spectra - Implement gui selection eventually
x_min = 184;
x_max = 1984;
y_min = 50; #34
y_max = 94; #94

# Read in Image Reduced Data
sci_image_reduced,_ = extract_data(output_dir/f"{object_name.replace(' ', '')}_science_image_reduced.fits")
std_image_reduced,_ = extract_data(output_dir/f"{object_name.replace(' ', '')}_standard_image_reduced.fits")

# Trim images
trimmed_sci = trim_image(sci_image_reduced, x_min, x_max, y_min, y_max)
trimmed_arc = trim_image(arc_data, x_min, x_max, y_min, y_max)
trimmed_std = trim_image(std_image_reduced, x_min, x_max, y_min, y_max)

# Cosmic Ray Removal
cleaned_sci = clean_cosmic_rays(trimmed_sci)
cleaned_arc = clean_cosmic_rays(trimmed_arc)
cleaned_std = clean_cosmic_rays(trimmed_std)

plot_image(trimmed_sci, keyword = 'Trimmed', object_name = object_name)
trimmed_1D_spectrum = extract_1D_spectrum(trimmed_sci)
plot_spectrum(trimmed_1D_spectrum, keyword = 'Trimmed', object_name = object_name)

'''4 - ASPIRED 2D Spec and Trace'''
sci2d = extract_twodspec(cleaned_sci, cleaned_arc, hdr_sci, is_standard=False)
std2d = extract_twodspec(cleaned_std, cleaned_arc, std_hdr, is_standard=True)

'''5 - Wavelength Calibration with ARC Lamp'''
atlas = [4158.59, 4200.67, 4277.53, 4609.57, 4764.86, 4879.86, 6032.127, 6965.43, 7067.21, 7272.94, 7383.9, 7503.87,
         7635.11, 7948.18, 8014.79, 8115.31, 8264.52, 8521.44, 9122.97, 9224.50] #Grating 7
element = ['CuAr']*len(atlas)
onedspec = wavelength_calibration(sci2d, std2d, atlas, element)

'''6 - Flux Calibration'''
perform_flux_calibration(onedspec, standard_name='hilt600')

'''7 - Save Final Spectrum'''
final_csv = save_final_spectra(onedspec, object_name, output_dir)

'''8 - Inspect Final Spectrum'''
plot_final_spectrum(output_dir/f"{object_name.replace(' ', '')}_final.csv")