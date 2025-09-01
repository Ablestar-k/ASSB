#!/bin/bash

python3 ../utility/Ensemble_avg.x
python3 ../utility/MSD_mda.py
python3 ../utility/MSD_alpha.py
python3 ../utility/MSD_visual.py
python3 ../utility/RDF_mda.py
python3 ../utility/RDF_same.py
python3 ../utility/RDF_viwual.py



