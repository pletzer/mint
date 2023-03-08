import numpy as np
from pathlib import Path
import mint

DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')

def test_1():
    nf = mint.NcFieldRead(DATA_DIR / 'cs8_wind.nc', 'u_in_w2h')

