#!/usr/bin/env bash

# full physics particle data
python download_ptcl_data.py Illutris-1 135
python download_ptcl_data.py Illutris-1 99

# dark matter only particle data
python download_ptcl_data.py Illutris-1-Dark 135
python download_ptcl_data.py Illutris-1-Dark 99
