import argparse
from terminally_ASPIRED.spectral_reducer import SpectralReductionPipeline
import shutil
from importlib import resources
from pathlib import Path

def copy_defaults(target_dir: str, new_name: str = "config.json") -> Path:
    """
    Copy defaults.json from the package to target_dir under new_name.
    Returns the full path to the copied file.
    """
    target_path = Path(target_dir) / new_name
    target_path.parent.mkdir(parents=True, exist_ok=True)

    with resources.files("terminally_ASPIRED.config_files").joinpath("defaults.json").open("rb") as src:
        with open(target_path, "wb") as dst:
            shutil.copyfileobj(src, dst)

    return target_path

def make_config():
    parser = argparse.ArgumentParser(description="Create a copy of default.json config file for pipline")
    parser.add_argument("-n", "--name", type = str, default="copy.json", help = "Name of config file")
    parser.add_argument("-O", "--output-dir", type = str, default="", help = "Specify output directory name. Default is current directory" )
    args = parser.parse_args()

    config_path = copy_defaults(args.output_dir, args.name)
    print(f"Defaults copied to {config_path}")


def main():

    parser = argparse.ArgumentParser(
        description="Run ASPIRED spectral reduction pipeline for 1.9m SAAO data."
    )
    parser.add_argument("science", type = str, help = "Path to science FITS file")
    parser.add_argument("arc", type = str, help = "Path to Arc FITS file")
    parser.add_argument("standard", type = str, help = "Path to standard FITS file")
    parser.add_argument("standard_arc", type = str, help = "Path to standard arc FITS file")
    parser.add_argument("--config",
                        type = str,
                        default = "",
                        help = "Path to config JSON file (default: config_files/defaults.json)"
                        )
    parser.add_argument("-gr" ,"--grating", type = str, choices = ["7", "6", "13", "custom"],default = "7",
                        help = "Choose which grating atlas to use. Options: 7, 6, 13, custom")
    parser.add_argument("-b", "--bias", type = str, default="", help = "Path to bias directory skip on empty")
    parser.add_argument("-f", "--flat-field", type = str, default="", help = "Path to flat directory skip on empty")
    parser.add_argument("-t","--interactive-trim", action = "store_true", help = "Enables interactive trim of science images")
    parser.add_argument("--show-plots", action = "store_true", help = "Enables plotting of intermediate ASPIRED images")
    parser.add_argument("-s", "--smooth", type = int, default=1, help = "Box smoothing by n points applied to final spectrum plot")
    parser.add_argument("-v", "--verbose",action = "store_true", help = "Enables verbose mode")
    parser.add_argument("--no-warnings", action = "store_true", help = "Disables warnings")
    parser.add_argument("-O", "--output-dir",type= str, default= None, help="Specify output directory name. Default is object name from fits header")
    parser.add_argument("--show-sky", action="store_true", help = "Show sky flux on final image")
    parser.add_argument("--show-error", action = "store_true", help = "Show flux error on final image")
    args = parser.parse_args()

    # Run Pipeline
    pipeline = SpectralReductionPipeline(science_file=args.science,
                                         arc_file=args.arc,
                                         std_file=args.standard,
                                         std_arc_file=args.standard_arc,
                                         config_path=args.config,
                                         grating = args.grating,
                                         bias_path=args.bias,
                                         flat_path=args.flat_field,
                                         interactive_trim = args.interactive_trim,
                                         show_plots=args.show_plots,
                                         smooth=args.smooth,
                                         verbose=args.verbose,
                                         no_warnings=args.no_warnings,
                                         output_dir_name=args.output_dir,
                                         sky=args.show_sky,
                                         error = args.show_error)

    pipeline.run()

if __name__ == "__main__":
    main()
