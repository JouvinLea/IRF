from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import Angle, SkyCoord
import numpy as np
import math
from glob import glob

# Write minimal files and observations table needed for the gammapy data store

#obs_ids = ["023523", "023526", "023559", "023592"]
#file_types = ['events', 'aeff', 'edisp', 'psf']
#listfile = glob('run*.fits')
#La commande en dessous permet d'avoir les runs number car on sait que le numero du run pour les fits de bruno se trouver entre le 5eme et 11eme caractere 
obs_ids = [int(file[6:11]) for file in glob('run*.fits')]
file_types = ['events', 'aeff', 'edisp']

# FILE TABLE

rows = []
for obs in obs_ids:
#for name in listfile:
    for filetype in file_types:
        if filetype == 'events':
#            name = "run_0"+obs+"_std_north_1b_eventlist.fits"
            name = "run_0{:06d}_std_north_1b_eventlist.fits".format(obs)
        else:
#            name = "hess_" + filetype + "_2d_" + obs + ".fits"
            name = "hess_{}_2d_{:06d}.fits".format(filetype, obs)
        #if filetype == 'events':
        #    name = "hess_events_simulated_" + obs + ".fits"
        data = dict()
        data['OBS_ID'] = obs
        data['TYPE'] = filetype
        data['NAME'] = name
        data['SIZE'] = ''
        data['MTIME'] = ''
        data['MD5'] = ''

        hdu_list = fits.open(str(name))
        hdu = hdu_list[1]
        header = hdu.header
        data['HDUNAME'] = hdu.name
        if 'HDUCLAS2' in header:
            data['HDUCLASS'] = header['HDUCLAS2']
        else:
            data['HDUCLASS'] = 'EVENTS'
        rows.append(data)

names = [
    'OBS_ID', 'TYPE',
    'NAME', 'SIZE', 'MTIME', 'MD5',
    'HDUNAME', 'HDUCLASS',
    ]
table = Table(rows=rows, names=names)
table.write("files.fits.gz", overwrite=True)



# OBS TABLE

rows = []
for obs in obs_ids:
    #filetype == 'events':
    name = "run_0{:06d}_std_north_1b_eventlist.fits".format(obs)
    hdu_list = fits.open(str(name))
    hdu = hdu_list[1]
    header = hdu.header
    data = dict()
    data['OBS_ID'] = obs
    data['ALT'] = header['ALT_PNT']
    data['AZ'] = header['AZ_PNT']
    data['N_TELS'] = header['N_TELS']
    data['MUONEFF'] = header['MUONEFF']
    data['ZEN'] = 90 - header['ALT_PNT']
    data['CosZEN'] = np.cos(data['ZEN'] * math.pi / 180.) 
    rows.append(data)

table = Table(rows=rows)
table['OBS_ID'].unit = "" 
table['ALT'].unit = "deg" 
table['AZ'] .unit = "deg" 
table['N_TELS'].unit = "" 
table['MUONEFF'].unit = "" 
table['ZEN'].unit = "deg"
table['CosZEN'].unit = "" 
table.meta['MJDREFI'] = 51544
table.meta['MJDREFF'] = 0.5
table.write("observations.fits.gz", overwrite=True)

