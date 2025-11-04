"""
Usage: `python run_trp.py <sample_id> [--maskknowntransits]`
See docstring for trp_pipeline.run_trp for verbose explanation.

Example usage:
`python run_trp.py debug`
`python run_trp.py mystarlist.csv --writevetplot --lcpipeline unpopular`
`python run_trp.py example_hjlist.csv --maskknowntransits --writevetplot`
"""
import sys
import argparse
from trp.trp_pipeline import run_trp

def main():

    parser = argparse.ArgumentParser(description="Run TRP pipeline")

    parser.add_argument(
        "sample_id",
        help="Sample ID to process (see trp_pipeline.py)"
    )

    parser.add_argument(
        "--maskknowntransits",
        action="store_true",
        help="Search for known transits and mask them (default: False)"
    )

    parser.add_argument(
        "--writevetplot",
        action="store_true",
        help="Make a plot showing result of Prot search (default:False)"
    )

    parser.add_argument(
        "--forcerun",
        action="store_true",
        help="Overwrite pre-existing cached files (default:False)"
    )

    parser.add_argument(
        "--lcpipeline",
        default="qlp",
        help="Light curve pipeline identifier (default: qlp)"
    )

    args = parser.parse_args()

    run_trp(
        args.sample_id,
        mask_known_transits=args.maskknowntransits,
        write_vetplot=args.writevetplot,
        forcerun=args.forcerun,
        lcpipeline=args.lcpipeline
    )

if __name__ == "__main__":
    main()
