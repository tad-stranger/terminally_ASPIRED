import argparse

from fontTools.svgLib.path import arc

from spectral_reducer import SpectralReductionPipeline

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
                        default = "config_files/defaults.json",
                        help = "Path to config JSON file (default: config_files/defaults.json)"
                        )
    parser.add_argument("-b", "--bias", type = str, default="", help = "Path to bias directory skip on empty")
    parser.add_argument("-f", "--flat-field", type = str, default="", help = "Path to flat directory skip on empty")
    parser.add_argument("--interactive-trim", action = "store_true", help = "Enables interactive trim of science images")
    parser.add_argument("--show-plots", action = "store_true", help = "Enables plotting of intermediate ASPIRED images")
    parser.add_argument("-s", "--smooth", type = int, default=1, help = "Box smoothing by n points applied to final spectrum plot")
    parser.add_argument("-v", "--verbose",action = "store_true", help = "Enables verbose mode")
    args = parser.parse_args()

    # Run Pipeline
    pipeline = SpectralReductionPipeline(science_file=args.science,
                                         arc_file=args.arc,
                                         std_file=args.standard,
                                         std_arc_file=args.standard_arc,
                                         config_path=args.config,
                                         bias_path=args.bias,
                                         flat_path=args.flat_field,
                                         show_plots=args.show_plots,
                                         smooth=args.smooth,
                                         verbose=args.verbose)

    if args.interactive_trim:
        pipeline.run_with_interactive_trim()
    else:
        pipeline.run()


if __name__ == "__main__":
    main()