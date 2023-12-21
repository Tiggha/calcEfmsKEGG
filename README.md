# calcEfmsKEGG
This program downloads models from KEGG database, converts KGML files to SBML and enumerates EFMs in metabolic models.

## Requirements:
    Python 3.7 or higher
    pip
    virtualenv


## mplrs: installation
    Download latest version from (http://cgm.cs.mcgill.ca/~avis/C/lrslib/archive/), use lrslib-071b
    Installation mplrs (http://cgm.cs.mcgill.ca/~avis/C/lrslib/USERGUIDE.html#mplrs), use openmpi
    Default installation: make mplrs

## Installation
    python3 -m venv env
    source env/bin/activate
    pip install -r requirements.txt

## Deactivate virtual environment
    deactivate

## Usage
The models needs to be in the 'models' directory.

    source env/bin/activate
    python3 EFM_AIO.py -f models/eco01100.sbml.xml -m efmtool
    python3 extract_reactions.py