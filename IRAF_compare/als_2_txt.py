
from os import getcwd, listdir
from os.path import join, realpath, dirname, isfile
from astropy.io import ascii
from astropy.table import hstack


"""
Convert .als files to the .txt files necessary to feed DAOMASTER
"""

mypath = realpath(join(getcwd(), dirname(__file__)))

# als_fls = []
for file in listdir(mypath):
    f = join(mypath, file)
    if isfile(f):
        if f.endswith('.als.1'):
            print(file)
            # als_fls.append(f)
            als_f = ascii.read(f, format='daophot')
            # 'CHI', 'SHARPNESS' columns are repeated in the two last positions
            # for DAOMASTER to read the values properly.
            data1 = als_f[
                'ID', 'XCENTER', 'YCENTER', 'MAG', 'MERR', 'CHI', 'SHARPNESS']
            data2 = als_f['CHI', 'SHARPNESS']
            data2.rename_column('CHI', 'CHI2')
            data2.rename_column('SHARPNESS', 'SHARPNESS2')
            data = hstack([data1, data2])
            als_n = f.split('/')[-1].replace('als.1', 'txt')
            ascii.write(
                data, als_n, format='fixed_width_no_header', delimiter=' ',
                overwrite=True)
