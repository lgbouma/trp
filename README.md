# tess-rotation-pipeline

### Purpose

Given a set of stars, measure their rotation periods.

### Install

Currently, the only installation method is to clone the repository and run the
`setup.py` script.
```
python setup.py develop
```

Dependencies are currently:
```
astropy==5.2.2
pandas==2.0.2
numpy==1.20.3
matplotlib==3.7.1
scipy==1.10.1
pytest==6.2.5
lightkurve==2.4.0
astrobase==0.5.3
```
which can be pip installed through something like `pip install -r
requirements.txt`.  A final dependency is the setup.py installed `cpv` package
from https://github.com/lgbouma/cpv.


### Usage

While a range of usage modes are available, a common pattern given a smaller
number of stars with known TIC8 identifiers would be to make a CSV file at
/data/targetlists, similar to `example_starlist.csv`.  Then, you can run

`python run_trp.py example_starlist.csv`

Pipeline options are currently hard-coded in, but can be changed at
`/trp/trp_pipeline.run_trp`.


### Advanced features

* [WIP] After measuring rotation periods, assess range of sinusoidal periods and
  amplitudes which could have been detected.

* Compatible with the Open Science Grid.
