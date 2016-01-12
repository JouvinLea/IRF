#! /usr/bin/env python
import numpy as np
import pyfits
from scipy import interpolate
import math
from pyfits import Column
import sys
import os
import argparse
from glob import glob


def gauss(x, sigma, mean):
    f = 1 / (np.sqrt(2 * math.pi) * sigma) * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2))
    return f


"""
Extract and interpolate IRF for a specific run and config
Example of command line to run with the fits event file of the observation and the 4D numpy table matching with the config we are interested in
./Interpolation_perrun.py 'run_0054749_std_north_1b_eventlist.fits' 'IRFstd_north_1b.npz'
"""

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract and interpolate IRF for a specific run and config')
    parser.add_argument('infile', action="store", help='Input fits event file name')
    parser.add_argument('irffile', action="store", help='IRF file to be used for this run and config')
    results = parser.parse_args()
    nrun = results.infile[5:11]
    print "Building IRF files for ", results.infile
    print "Using 4D IRF file : ", results.irffile

    # Load the IRFs (Effective area and resolution (biais and sigma)) from the numpy 4D table created for each config with the script histo_array.py.
    IRF = np.load(results.irffile)
    IRFArea = IRF["TableArea"]
    IRFSigma = IRF["TableSigma"]
    IRFBiais = IRF["TableBiais"]
    enMC = IRF["enMC"]
    lnenMC = IRF["lnenMC"]
    zenMC = IRF["zenMC"]
    effMC = IRF["effMC"]
    offMC = IRF["offMC"]
    binoffMC = len(offMC)
    binEMC = len(enMC)
    binEreco = 50
    #Final dimension of the effective area table and rmf interpolated on the zenithal angle and muon efficiency of the specific observation
    bineffarea = len(offMC) * len(enMC)
    bineffresol = len(offMC) * len(enMC) * binEreco

    #Offset bin hi and low are equal
    off_low = offMC
    off_hi = offMC

    #MC energy value are in log but the width of the energy bins changes from bin to bin. Therefore the low and high value of the MC bin are determined using the middle of the MC bins.
    binlnEMC = lnenMC[1:] - lnenMC[:-1]
    # For the low edge of the MC value, for the first bin we take the width of the fisrt bin of binlnEMC
    binlnEMClow = np.insert(binlnEMC, 0, binlnEMC[0])
    # For the high edge of the MC value, for the last bin we take the width of the last bin of binlnEMC
    binlnEMCup = np.insert(binlnEMC, -1, binlnEMC[-1])
    lnEMClow = lnenMC - binlnEMClow / 2
    lnEMCup = lnenMC + binlnEMCup / 2
    E_true_low = pow(10, lnEMClow)
    E_true_up = pow(10, lnEMCup)

    # Defining Migra=Etrue/Ereco (between 0.1 and 10), logspacing
    lnEtrue_reco = np.linspace(-1, 1, binEreco)
    #the width of the migra bin is the same in log so we just take the one of the first bin
    binlnEtrue_reco = lnEtrue_reco[1] - lnEtrue_reco[0]
    lnE_true_reco_low = lnEtrue_reco - binlnEtrue_reco / 2
    lnE_true_reco_up = lnEtrue_reco + binlnEtrue_reco / 2
    Etrue_reco = pow(10, lnEtrue_reco)
    E_true_reco_low = pow(10, lnE_true_reco_low)
    E_true_reco_hi = pow(10, lnE_true_reco_up)

    AreaRun = np.zeros((binoffMC, binEMC))
    ResolRun = np.zeros((binoffMC, binEreco, binEMC))
    hdurun = pyfits.open(results.infile)
    ZenRun = 90 - hdurun[1].header["ALT_PNT"]
    EffRun = hdurun[1].header["MUONEFF"] * 100
    #For each MC true energy bin and offset, the Area and RMF are interpolated on the zenith and the muon efficiency of the run
    for (iEMC, EMC) in enumerate(enMC):
        for (ioff, off) in enumerate(offMC):
            InterArea = interpolate.interp2d(effMC, np.cos(zenMC * math.pi / 180), IRFArea[iEMC, ioff, :, :])
            InterBiais = interpolate.interp2d(effMC, np.cos(zenMC * math.pi / 180), IRFBiais[iEMC, ioff, :, :])
            InterSigma = interpolate.interp2d(effMC, np.cos(zenMC * math.pi / 180), IRFSigma[iEMC, ioff, :, :])
            AreaRun[ioff, iEMC] = InterArea(EffRun, np.cos(ZenRun * math.pi / 180))
            BiaisRun = InterBiais(EffRun, np.cos(ZenRun * math.pi / 180))
            SigmaRun = InterSigma(EffRun, np.cos(ZenRun * math.pi / 180))
            #Assume a log-normal distribution for migra
            ResolRun[ioff, :, iEMC] = gauss(lnEtrue_reco, SigmaRun, BiaisRun)
            # Normalised
            norm = np.sum(ResolRun[ioff, :, iEMC] * (E_true_reco_hi - E_true_reco_low))
            if (np.isnan(norm)):
                ResolRun[ioff, :, iEMC] = 0

            else:
                ResolRun[ioff, :, iEMC] = ResolRun[ioff, :, iEMC] / norm

    #Store the results of the AREA and RMF for the specific observation in fits files
    c1_area = Column(name='ENERG_LO', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_low))
    c2_area = Column(name='ENERG_HI', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_up))
    c3_area = Column(name='THETA_LO', format=str(binoffMC) + 'E', unit='TeV', array=np.atleast_2d(off_low))
    c4_area = Column(name='THETA_HI', format=str(binoffMC) + 'E', unit='TeV', array=np.atleast_2d(off_hi))
    c5_area = Column(name='EFFAREA', format=str(bineffarea) + 'E', unit='TeV', array=np.expand_dims(AreaRun, 0))
    c6_area = Column(name='EFFAREA_RECO', format=str(bineffarea) + 'E', unit='TeV', array=np.expand_dims(AreaRun, 0))
    tbhdu_area = pyfits.BinTableHDU.from_columns([c1_area, c2_area, c3_area, c4_area, c5_area, c6_area])
    for i in range(1, 7):
        tbhdu_area.header.comments['TTYPE' + str(i)] = 'label for field ' + str(i)
        tbhdu_area.header.comments['TFORM' + str(i)] = 'data format of field: 4-byte REAL'
        tbhdu_area.header.comments['TUNIT' + str(i)] = 'physical unit of field '

    tbhdu_area.header.set("EXTNAME", "EFFECTIVE AREA", "name of this binary table extension ")
    tbhdu_area.header.set("TDIM5", "(" + str(binEMC) + "," + str(binoffMC) + ")")
    tbhdu_area.header.set("TDIM6", "(" + str(binEMC) + "," + str(binoffMC) + ")")
    tbhdu_area.header.set("LO_THRES", "0.320357620716095", "TeV")
    tbhdu_area.header.set("HI_THRES", "31.1716079711914", "TeV")

    tbhdu_area.writeto('hess_aeff_2d_' + nrun + '.fits')

    c1_resol = Column(name='ETRUE_LO', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_low))
    c2_resol = Column(name='ETRUE_HI', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_up))
    c3_resol = Column(name='MIGRA_LO', format=str(binEreco) + 'E', unit='', array=np.atleast_2d(E_true_reco_low))
    c4_resol = Column(name='MIGRA_HI', format=str(binEreco) + 'E', unit='', array=np.atleast_2d(E_true_reco_hi))
    c5_resol = Column(name='THETA_LO', format=str(binoffMC) + 'E', unit='deg', array=np.atleast_2d(off_low))
    c6_resol = Column(name='THETA_HI', format=str(binoffMC) + 'E', unit='deg', array=np.atleast_2d(off_hi))
    c7_resol = Column(name='MATRIX ', format=str(bineffresol) + 'E', unit='TeV', array=np.expand_dims(ResolRun, 0))
    tbhdu_resol = pyfits.BinTableHDU.from_columns(
        [c1_resol, c2_resol, c3_resol, c4_resol, c5_resol, c6_resol, c7_resol])
    for i in range(1, 8):
        tbhdu_resol.header.comments['TTYPE' + str(i)] = 'label for field ' + str(i)
        tbhdu_resol.header.comments['TFORM' + str(i)] = 'data format of field: 4-byte REAL'
    tbhdu_resol.header.set("EXTNAME", "EDISP_2D", "name of this binary table extension ")
    tbhdu_resol.header.set("TDIM7", "(" + str(binEMC) + "," + str(binEreco) + "," + str(binoffMC) + ")")

    tbhdu_resol.writeto('hess_edisp_2d_' + nrun + '.fits')
