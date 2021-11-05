from mint import NcFieldRead
import numpy
from pathlib import Path


DATA_DIR = Path(__file__).absolute().parent.parent.parent / Path('data')

# noinspection SpellCheckingInspection
def test_read_u1():
    filename = str(DATA_DIR / Path('lfric_diag_wind.nc'))
    varname = 'u1'
    reader = NcFieldRead(filename, varname)
    

if __name__ == '__main__':
    test_read_u1()