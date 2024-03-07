"""
Usage: `python run_trp.py <sample_id>`
See docstring for trp_pipeline.run_trp for verbose explanation.

(Example: `python run_trp.py debug`)
"""
import sys
from trp.trp_pipeline import run_trp

def main():
    sample_id = sys.argv[1]
    run_trp(sample_id)

if __name__ == "__main__":
    main()
