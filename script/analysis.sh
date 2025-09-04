#!/bin/bash

../utility/Ensemble_avg.x
python3 ../utility/MSD_mda.py
python3 ../utility/MSD_alpha.py
python3 ../utility/MSD_alpha_window.py

python3 ../utility/RDF_mda.py
python3 ../utility/RDF_same.py