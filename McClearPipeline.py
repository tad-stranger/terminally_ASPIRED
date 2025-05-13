from astropy.io import fits
from aspired import image_reduction, spectral_reduction
from astroscrappy import detect_cosmics
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import plotly.io as pio
import lacosmic
from pathlib import Path
import glob
import os

def extract_1D_spectrum(science_data):
    spectrum_1d = np.sum(science_data, axis = 0)
    return spectrum_1d

def extract_data(filename):
    hdul = fits.open(filename)
    header = hdul[0].header
    data = hdul[0].data
    hdul.close()
    return data, header

def plot_image(image_data, keyword):
    plt.imshow(image_data, cmap='jet', origin='lower', aspect='auto', norm=LogNorm())
    plt.colorbar()
    plt.title(f"{keyword} Spectrum of {object_name}")
    plt.xlabel("Pixel (Dispersion)")
    plt.ylabel("Spatial Axis")
    plt.show()

def plot_spectrum(spec_data, keyword):
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

def trim_and_save_folder(input_dir, output_dir, x_min, x_max, y_min, y_max):

    os.makedirs(output_dir, exist_ok=True)
    fits_files = glob.glob(f"{input_dir}/*.fits")

    for file in fits_files:
        data, header = fits.getdata(file, header=True)
        trimmed = data[y_min:y_max, x_min:x_max]
        filename = Path(file).name
        fits.writeto(os.path.join(output_dir, filename), trimmed, header, overwrite=True)

pio.renderers.default = "browser"

''' 1 - Extract Data'''

# Set Observation date path
observations_path = "Transients Observations/17042025/"

# Extract Science images
sci_data, hdr_sci = extract_data(f"{observations_path}SCIENCE/SN2024xuo/a6261089.fits")
sci_data = np.flip(sci_data, axis = 1)

# Extract Object Info
object_name = hdr_sci.get('OBJECT', 'Unkown')
tar_RA = hdr_sci.get('TARG-RA')
tar_DEC = hdr_sci.get('TARG-DEC')

# Extract ARC Images:
arc_data, hdr_arc = extract_data(f"{observations_path}ARC/a6261038.fits")
arc_data = np.flip(arc_data,axis = 1)

# Extract Standard Star (For Flux Calibration)
std_data, std_hdr = extract_data(f"{observations_path}SCIENCE/HILT600/a6261086.fits")
std_data = np.flip(std_data,axis = 1)

# Plot Raw Image for inspection
plot_image(sci_data, keyword = 'Raw')

# Extract 1D spectrum:
raw_1D_spectrum = extract_1D_spectrum(sci_data)
plot_spectrum(raw_1D_spectrum, keyword = 'Raw')



'''Image Reduction'''


'''Setting up Image Reduction'''
bias_folder = f"{observations_path}BIASS"
flats_folder = f"{observations_path}FLATS"

# Make Output Directory
output_dir = Path(f"./ReducedSpectra/{object_name.replace(' ', '')}")
output_dir.mkdir(exist_ok=True)

# File Names
science_file = output_dir / f"{object_name.replace(' ', '')}_science_reduced.fits"
arc_file = output_dir / f"{object_name.replace(' ', '')}_arc_reduced.fits"
std_file = output_dir / f"{object_name.replace(' ', '')}_std_reduced.fits"

# Save fits
fits.writeto(science_file, sci_data, hdr_sci, overwrite=True)
fits.writeto(arc_file, arc_data, hdr_arc, overwrite=True)
fits.writeto(std_file, std_data, std_hdr, overwrite=True)

# Make File list
write_filelist_for_aspired(bias_dir = bias_folder, flat_dir=flats_folder, science_file= science_file,
                           arc_file=arc_file, output_name=f"{object_name.replace(' ', '')}_science_reduction_file.list")

write_filelist_for_aspired(bias_dir = bias_folder, flat_dir=flats_folder, science_file= std_file,
                           arc_file=arc_file, output_name = f"{object_name.replace(' ', '')}_standard_reduction_file.list")

'''ASPIRED Implementation of ImageReduction'''

science_frame = image_reduction.ImageReduction()
science_frame.add_filelist(f"{object_name.replace(' ', '')}_science_reduction_file.list")
science_frame.load_data()
science_frame.reduce()
science_frame.inspect(display=True, width = 1024, height = 1024)
science_frame.save_fits(filename = f"{object_name.replace(' ', '')}_science_image_reduced.fits", overwrite=True)

standard_frame = image_reduction.ImageReduction()
standard_frame.add_filelist(f"{object_name.replace(' ', '')}_standard_reduction_file.list")
standard_frame.load_data()
standard_frame.reduce()
standard_frame.reduce()
standard_frame.inspect(display=True, width = 1024, height = 1024)
standard_frame.save_fits(filename = f"{object_name.replace(' ', '')}_standard_image_reduced.fits", overwrite=True)




'''2 - Trimming and Cosmic Ray Removal'''
# 2.1 Trimming
x_min = 184;
x_max = 1984;
y_min = 41; #34
y_max = 94; #94



science_reduced, header_sci_reduced = extract_data(f"{object_name.replace(' ', '')}_science_image_reduced.fits")
std_data, std_hdr = extract_data(f"{object_name.replace(' ', '')}_standard_image_reduced.fits")

trimmed_data = science_reduced[y_min : y_max, x_min : x_max]
trimmed_arc_data = arc_data[y_min : y_max, x_min : x_max]
trimmed_std_data = std_data[y_min : y_max, x_min : x_max]

plot_image(trimmed_data, keyword = 'Trimmed')

trimmed_1D_spectrum = extract_1D_spectrum(trimmed_data)
plot_spectrum(trimmed_1D_spectrum, keyword = 'Trimmed')


# 2.2 Remove Cosmic Rays from Data
_, science_data_cleaned = detect_cosmics(trimmed_data, sigclip=4.5, sigfrac=0.3, objlim=5.0)
plot_image(science_data_cleaned, keyword = 'Cosmic Ray Cleaned')

cosmic_cleaned_1D_spectrum = extract_1D_spectrum(science_data_cleaned)
plot_spectrum(cosmic_cleaned_1D_spectrum, keyword = 'Cosmic Ray Cleaned')

# 2.3 Remove Cosmic Rays from ARC Spectrum
_, clean_arc_data = detect_cosmics(trimmed_arc_data, sigclip=4.5, sigfrac=0.3, objlim=5.0)
plot_image(clean_arc_data, keyword = 'Cosmic Ray Cleaned')

arc_cosmic_cleaned_1D_spectrum = extract_1D_spectrum(clean_arc_data)
plot_spectrum(arc_cosmic_cleaned_1D_spectrum, keyword = 'Cosmic Ray Cleaned')

# Remove Cosmic Rays from Standard Spectrum
_, clean_std_data = detect_cosmics(trimmed_std_data, sigclip=4.5, sigfrac=0.3, objlim=5.0)


'''ASPIRED Implementation'''
reduced_sci_data = science_data_cleaned
# Science Data
n_sources = 1 # number of objects on slit, ususally 1
science2D = spectral_reduction.TwoDSpec(reduced_sci_data, cosmicray = False, header=hdr_sci) #removed header = hdr_sci
science2D.set_properties(saxis = 1, flip = False, cosmicray = False)
#science2D.apply_flip_to_arc()
science2D.ap_trace(nspec = n_sources, nwindow= 100, trace_width = 5, resample_factor= 15, fit_deg= 3, display = True, width = 1024, height = 1024)
science2D.ap_extract(apwidth = 3, optimal = True, algorithm = 'marsh89', skysep = 5, skywidth = 10, skydeg = 1, display=True)
science2D.add_arc(arc = clean_arc_data)
science2D.extract_arc_spec()

# Standard Star Data
standard2D = spectral_reduction.TwoDSpec(clean_std_data, cosmicray = False, header=std_hdr)
standard2D.set_properties(saxis = 1, flip = False, cosmicray = False)
#standard2D.apply_flip_to_arc()
standard2D.ap_trace(nspec = n_sources, nwindow= 20, trace_width = 5, resample_factor= 2, fit_deg= 3, display = False, width = 1024, height = 1024)
standard2D.ap_extract(apwidth = 3, optimal = True, algorithm = 'marsh89', skysep = 0, skywidth = 5, skydeg = 1, display=True)
standard2D.add_arc(arc = clean_arc_data)
standard2D.extract_arc_spec()


''' 3 - Wavelength Calibration with ARC Lamp'''
onedspec = spectral_reduction.OneDSpec()
onedspec.from_twodspec(science2D, stype = 'science')
onedspec.from_twodspec(standard2D, stype = 'standard')
onedspec.find_arc_lines(prominence=7, refine=True, display=True, stype='science+standard')

# Grating 6 atlas!
#atlas_full = [4259.3620, 4277.5278, 4333.5610, 4348.0640, 4400.9860, 4481.8110, 4510.7330, 4545.0520, 4579.3500, 4589.8980, 4609.5670,
         #4657.9010, 4764.8650, 4806.0210, 4847.8100, 4879.8640, 4965.0800, 5017.1630, 5062.0370, 5162.2850, 5187.7460, 5221.2729, 5495.8740, 5558.7020,
         #5606.7330, 5739.5200, 5834.2630, 5860.3100, 5888.5840, 5912.0850, 5928.8130, 6032.1270, 6043.2230, 6059.3730, 6114.9230, 6155.2388, 6172.2778, 6296.8790, 6307.6510, 6384.7170, 6416.3070, 6466.5530, 6538.1120, 6604.8530, 6677.2822,
         #6752.8330, 6766.6132, 6871.2891, 6937.6658]

#atlas = [4764.8650, 4860.0210, 4879.8640, 5912.0850, 6032.1270]

#atlas_low = [4764.8650, 4806.0210, 4847.8100, 4879.8640, 4965.0800, 5017.1630, 5062.0370, 5162.2850, 5187.7460, 5221.2729, 5495.8740, 5558.7020,
         #5606.7330, 5739.5200, 5834.2630, 5860.3100, 5888.5840, 5912.0850, 5928.8130, 6032.1270]

atlas = [4158.59, 4200.67, 4277.53, 4609.57, 4764.86, 4879.86, 6032.127, 6965.43, 7067.21, 7272.94, 7383.9, 7503.87,
         7635.11, 7948.18, 8014.79, 8115.31, 8264.52, 8521.44, 9122.97, 9224.50]


element = ['CuAr']*len(atlas)

onedspec.initialise_calibrator(stype = 'science+standard',)
onedspec.add_user_atlas(elements= element, wavelengths = atlas, stype = 'science+standard')
onedspec.set_hough_properties(num_slopes = 1000, xbins = 100, ybins = 100, min_wavelength = 3500, max_wavelength=9000, range_tolerance= 250, stype= 'science+standard')

onedspec.do_hough_transform(stype = 'science+standard')
onedspec.set_ransac_properties(minimum_matches=15)#normal = 15
onedspec.fit(max_tries = 1000, fit_deg = 3, fit_type = 'poly', display = True, stype = 'science+standard')
onedspec.apply_wavelength_calibration(stype = 'science+standard')
onedspec.inspect_reduced_spectrum(display = True, stype = 'science')
onedspec.inspect_reduced_spectrum(display = True, stype = 'standard')


'''Flux Calibration'''
onedspec.load_standard(target = 'hilt600') # Loads Literature standard star spectrum
onedspec.inspect_standard()

onedspec.get_sensitivity()
onedspec.inspect_sensitivity()

onedspec.apply_flux_calibration()

'''See the fully reduced spectrum'''
onedspec.inspect_reduced_spectrum(display = True, stype = 'science', save_fig= True, filename = object_name + '-reduced_spectrum', fig_type= 'png')
onedspec.inspect_reduced_spectrum(display = True, stype = 'standard')
