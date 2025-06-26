import argparse
from spectral_reducer import SpectralReductionPipeline

def main():
    parser = argparse.ArgumentParser(
        description="Run ASPIRED spectral reduction pipeline for 1.9m SAAO data."
    )

    parser.add_argument("science", type = str, help = "Path to science FITS file")
    parser.add_argument("arc", type = str, help = "Path to Arc FITS file")
    parser.add_argument("standard", type = str, help = "Path to standard FITS file")
    parser.add_argument("--config",
                        type = str,
                        default = "config_files/defaults.json",
                        help = "Path to config JSON file (default: config_files/defaults.json)"
                        )
    parser.add_argument("--use_bias_flats", action = "store_true",help = "Enables bias and flat field correction")
    parser.add_argument("--bias_path", type = str, help = "Path to bias frames")
    parser.add_argument("--flat_path", type = str, help = "Path to flat frames")
    parser.add_argument("--interactive_trim", action = "store_true", help = "Enables interactive trim of science images")

    args = parser.parse_args()

    # Run Pipeline
    pipeline = SpectralReductionPipeline(science_file=args.science,
                                         arc_file=args.arc,
                                         std_file=args.standard,
                                         config_path=args.config,
                                         use_bias_flats=args.use_bias_flats,
                                         bias_path=args.bias_path,
                                         flat_path=args.flat_path)

    if args.interactive_trim:
        pipeline.run_with_interactive_trim()
    else:
        pipeline.run()


if __name__ == "__main__":
    main()
