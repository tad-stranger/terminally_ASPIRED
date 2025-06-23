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

    args = parser.parse_args()

    # Run Pipeline
    pipeline = SpectralReductionPipeline(science_file=args.science,arc_file= args.arc, std_file= args.standard, config_path = args.config)
    pipeline.run()


if __name__ == "__main__":
    main()
