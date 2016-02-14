#! /usr/bin/env python
import numpy as np
import pyfits
from scipy import interpolate
import math
from pyfits import Column
import sys
import argparse
import os
from glob import glob


def gauss(x, sigma, mean):
    f = 1 / (np.sqrt(2 * math.pi) * sigma) * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2))
    return f


"""
Extract and interpolate IRF for a specific run and config
Example of command line to run with the fits event file of the observation and the 4D numpy table matching with the config we are interested in
./Interpolation_perrun.py 'run_0054749_std_north_1b_eventlist.fits' 'IRFstd_north_1b.npz' 'PSF_triplegauss_elm_south_stereo_save.npz' 'triplegauss'
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract and interpolate IRF for a specific run and config')
    parser.add_argument('infile', action="store", help='Input fits event file name')
    parser.add_argument('irffile', action="store", help='IRF file to be used for this run and config')
    parser.add_argument('psffile', action="store", help='PSF file to be used for this run and config')
    parser.add_argument('psf_function', action="store", help='PSF type function used to fit the PSF')
    results = parser.parse_args()
    nrun = results.infile[5:11]
    print "Building IRF files for ", results.infile
    print "Using 4D IRF file : ", results.irffile
    print "Using 4D IRF file : ", results.irffile
    print "Using 4D PSF file : ", results.psffile
    print "Using the function : ", results.psf_function, " to fit the psf"

    # Load les info sur les MCs depuis la table d'IRF ou est stocke pour toutes les nergies, zenith,offset et efficacite des MCs la valeur des la surface efficiace, du biais et sigma pour la resolution et du s1, s2, s3, A2, A3 de la tripplegauss utilisee pour fitter la psf
    PathTableIRF = "/Users/jouvin/Desktop/these/WorkGAMMAPI/IRF/CrabEventList/Crab"
    PathTablePSF = "/Users/jouvin/Desktop/these/WorkGAMMAPI/IRF/PSF"

    PSFtype = sys.argv[4]
    IRF = np.load(results.irffile)
    IRFArea = IRF["TableArea"]
    IRFSigma = IRF["TableSigma"]
    IRFBiais = IRF["TableBiais"]
    enMC = IRF["enMC"]
    lnenMC = IRF["lnenMC"]
    zenMC = IRF["zenMC"]
    effMC = IRF["effMC"]
    offMC = IRF["offMC"]
    if (PSFtype == "triplegauss"):
        PSF = np.load(results.psffile)
        PSFs1 = PSF["TableSigma1"]
        PSFs2 = PSF["TableSigma2"]
        PSFs3 = PSF["TableSigma3"]
        PSFA2 = PSF["TableA2"]
        PSFA3 = PSF["TableA3"]
        effMC_psf = effMC
        zenMC_psf = zenMC
    elif (PSFtype == "king"):
        PSF = np.load(PathTablePSF + "/PSF_king_" + coupure + ".npz")
        PSFSig = PSF["TableSig"]
        PSFGam = PSF["TableGam"]
    else:
        print "No function given for the PSF"

    binoffMC = len(offMC)
    binEMC = len(enMC)
    binEreco = 50
    bineffarea = len(offMC) * len(enMC)
    bineffresol = len(offMC) * len(enMC) * binEreco

    # reverifier qu ils ont bien ca dans leur bin PA en low edge et upper edge
    off_low = offMC
    off_hi = offMC

    # pour les extremites prendre le milieu des bin en log
    binlnEMC = lnenMC[1:] - lnenMC[:-1]
    # Pour le premier bin en energie pour defenr le edge low du bin on prend la demilargeur du premier bin
    binlnEMClow = np.insert(binlnEMC, 0, binlnEMC[0])
    # Pour le dernier bin en energie pour defenr le edge up du bin on prend la demilargeur du dernier bin
    binlnEMCup = np.insert(binlnEMC, -1, binlnEMC[-1])
    # Retrouver
    lnEMClow = lnenMC - binlnEMClow / 2
    lnEMCup = lnenMC + binlnEMCup / 2
    E_true_low = pow(10, lnEMClow)
    E_true_up = pow(10, lnEMCup)

    # Definition de Etrue/Ereco
    lnEtrue_reco = np.linspace(-1, 1, binEreco)
    # Le tableau d energie reco en log ont tous la meme largeur de bin donc on prend le premier
    binlnEtrue_reco = lnEtrue_reco[1] - lnEtrue_reco[0]
    lnE_true_reco_low = lnEtrue_reco - binlnEtrue_reco / 2
    lnE_true_reco_up = lnEtrue_reco + binlnEtrue_reco / 2
    Etrue_reco = pow(10, lnEtrue_reco)
    E_true_reco_low = pow(10, lnE_true_reco_low)
    E_true_reco_hi = pow(10, lnE_true_reco_up)

    AreaRun = np.zeros((binoffMC, binEMC))
    ResolRun = np.zeros((binoffMC, binEreco, binEMC))
    if (PSFtype == "triplegauss"):
        PSFS1Run = np.zeros((binoffMC, binEMC))
        PSFS2Run = np.zeros((binoffMC, binEMC))
        PSFS3Run = np.zeros((binoffMC, binEMC))
        PSFA2Run = np.zeros((binoffMC, binEMC))
        PSFA3Run = np.zeros((binoffMC, binEMC))
    elif (PSFtype == "king"):
        PSFSigRun = np.zeros((binoffMC, binEMC))
        PSFGamRun = np.zeros((binoffMC, binEMC))
    hdurun = pyfits.open(results.infile)
    ZenRun = 90 - hdurun[1].header["ALT_PNT"]
    EffRun = hdurun[1].header["MUONEFF"] * 100
    for (iEMC, EMC) in enumerate(enMC):
        for (ioff, off) in enumerate(offMC):
            # print ioff, " ", iEMC
            
            InterArea = interpolate.interp2d(effMC, np.cos(zenMC * math.pi / 180), IRFArea[iEMC, ioff, :, :])
            InterBiais = interpolate.interp2d(effMC, np.cos(zenMC * math.pi / 180), IRFBiais[iEMC, ioff, :, :])
            InterSigma = interpolate.interp2d(effMC, np.cos(zenMC * math.pi / 180), IRFSigma[iEMC, ioff, :, :])
            if (PSFtype == "triplegauss"):
                #interpol_param = dict(fill_value=None)
                ind_zen, ind_eff= np.where(PSFs1[iEMC, ioff, :, :] != -1)
                #Il faut qu il y ait au moins 4 eff et 4 zen pour que interp marche
                if(len(ind_zen)>4):
                    zensame=np.where(ind_zen != ind_zen[0])
                    effsame=np.where(ind_eff != ind_eff[0])
                    #Dans les minimum 4 valeurs, Il faut pas qu il y ait les memes valeurs en zenith et en eff pour toutes les simus
                    if((len(zensame[0])!=0) & (len(effsame[0])!=0)):
                        coord_eff=effMC[ind_eff]
                        coord_zen = zenMC[ind_zen]
                        InterS1 = interpolate.interp2d(coord_eff, np.cos(coord_zen * math.pi / 180), PSFs1[iEMC, ioff, ind_zen, ind_eff], fill_value=None)
                        InterS2 = interpolate.interp2d(coord_eff, np.cos(coord_zen * math.pi / 180), PSFs2[iEMC, ioff, ind_zen, ind_eff],fill_value=None)
                        InterS3 = interpolate.interp2d(coord_eff, np.cos(coord_zen * math.pi / 180), PSFs3[iEMC, ioff, ind_zen, ind_eff],fill_value=None)
                        InterA2 = interpolate.interp2d(coord_eff, np.cos(coord_zen * math.pi / 180), PSFA2[iEMC, ioff, ind_zen, ind_eff],fill_value=None)
                        InterA3 = interpolate.interp2d(coord_eff, np.cos(coord_zen * math.pi / 180), PSFA3[iEMC, ioff, ind_zen, ind_eff],fill_value=None)
            elif (PSFtype == "king"):
               if(len(ind_zen)>4):
                    if((len(zensame[0])!=0) & (len(effsame[0])!=0)): 
                        InterSig = interpolate.interp2d(coord_eff, np.cos(coord_zen * math.pi / 180), PSFSig[iEMC, ioff, ind_zen, ind_eff])
                        InterGam = interpolate.interp2d(coord_eff, np.cos(coord_zen * math.pi / 180), PSFGam[iEMC, ioff, ind_zen, ind_eff])
            AreaRun[ioff, iEMC] = InterArea(EffRun, np.cos(ZenRun * math.pi / 180))
            BiaisRun = InterBiais(EffRun, np.cos(ZenRun * math.pi / 180))
            SigmaRun = InterSigma(EffRun, np.cos(ZenRun * math.pi / 180))
            ResolRun[ioff, :, iEMC] = gauss(lnEtrue_reco, SigmaRun, BiaisRun)
            # etre sur que c est bien normalise
            norm = np.sum(ResolRun[ioff, :, iEMC] * (E_true_reco_hi - E_true_reco_low))

            if (np.isnan(norm)):
                ResolRun[ioff, :, iEMC] = 0
            else:
                ResolRun[ioff, :, iEMC] = ResolRun[ioff, :, iEMC] / norm

            if (PSFtype == "triplegauss"):
                if(len(ind_zen)>4):
                    if((len(zensame[0])!=0) & (len(effsame[0])!=0)):
                        PSFS1Run[ioff, iEMC] = InterS1(EffRun, np.cos(ZenRun * math.pi / 180))
                        PSFS2Run[ioff, iEMC] = InterS2(EffRun, np.cos(ZenRun * math.pi / 180))
                        PSFS3Run[ioff, iEMC] = InterS3(EffRun, np.cos(ZenRun * math.pi / 180))
                        PSFA2Run[ioff, iEMC] = InterA2(EffRun, np.cos(ZenRun * math.pi / 180))
                        PSFA3Run[ioff, iEMC] = InterA3(EffRun, np.cos(ZenRun * math.pi / 180))
                    else:
                        PSFS1Run[ioff, iEMC] = -1
                        PSFS2Run[ioff, iEMC] = -1
                        PSFS3Run[ioff, iEMC] = -1
                        PSFA2Run[ioff, iEMC] = -1
                        PSFA3Run[ioff, iEMC] = -1
                else:
                    PSFS1Run[ioff, iEMC] = -1
                    PSFS2Run[ioff, iEMC] = -1
                    PSFS3Run[ioff, iEMC] = -1
                    PSFA2Run[ioff, iEMC] = -1
                    PSFA3Run[ioff, iEMC] = -1
            elif (PSFtype == "king"):
                if(len(ind_zen)>4):
                    if((len(zensame[0])!=0) & (len(effsame[0])!=0)):
                        PSFSigRun[ioff, iEMC] = InterSig(EffRun, np.cos(ZenRun * math.pi / 180))
                        PSFGamRun[ioff, iEMC] = InterGam(EffRun, np.cos(ZenRun * math.pi / 180))
                    else:
                        PSFSigRun[ioff, iEMC] = -1
                        PSFGamRun[ioff, iEMC] = -1
                    

    # Ecriture des fichiers fits pour aeff, edisp et psf pour chaque observation
    # AEFF FITS FILE
    c1_area = Column(name='ENERG_LO', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_low))
    c2_area = Column(name='ENERG_HI', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_up))
    c3_area = Column(name='THETA_LO', format=str(binoffMC) + 'E', unit='deg', array=np.atleast_2d(off_low))
    c4_area = Column(name='THETA_HI', format=str(binoffMC) + 'E', unit='def', array=np.atleast_2d(off_hi))
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
    # tbhdu_area.header["EXTNAME"]='EFFECTIVE AREA'
    tbhdu_area.writeto('hess_aeff_2d_' + nrun + '.fits')

    # EDISP FITS FILE
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
    #        tbhdu_resol.header.comments['TUNIT'+str(i)]='physical unit of field '

    tbhdu_resol.header.set("EXTNAME", "EDISP_2D", "name of this binary table extension ")
    tbhdu_resol.header.set("TDIM7", "(" + str(binEMC) + "," + str(binEreco) + "," + str(binoffMC) + ")")
    # tbhdu_resol.header["EXTNAME"]='EFFECTIVE RESOL'
    tbhdu_resol.writeto('hess_edisp_2d_' + nrun + '.fits')

    # PSF FITS FILE
    c1_psf = Column(name='ENERG_LO', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_low))
    c2_psf = Column(name='ENERG_HI', format=str(binEMC) + 'E', unit='TeV', array=np.atleast_2d(E_true_up))
    c3_psf = Column(name='THETA_LO', format=str(binoffMC) + 'E', unit='deg', array=np.atleast_2d(off_low))
    c4_psf = Column(name='THETA_HI', format=str(binoffMC) + 'E', unit='deg', array=np.atleast_2d(off_hi))
    if (PSFtype == "triplegauss"):
        norm = 2 * np.pi * (PSFS1Run ** 2 + PSFA2Run * PSFS2Run ** 2 + PSFA3Run * PSFS3Run ** 2)
        norm[np.where(PSFS1Run==-1)]=-1
        c5_psf = Column(name='SIGMA_1', format=str(bineffarea) + 'E', unit='deg', array=np.expand_dims(PSFS1Run, 0))
        c6_psf = Column(name='AMPL_2', format=str(bineffarea) + 'E', unit='', array=np.expand_dims(PSFA2Run, 0))
        c7_psf = Column(name='SIGMA_2', format=str(bineffarea) + 'E', unit='deg', array=np.expand_dims(PSFS2Run, 0))
        c8_psf = Column(name='AMPL_3', format=str(bineffarea) + 'E', unit='', array=np.expand_dims(PSFA3Run, 0))
        c9_psf = Column(name='SIGMA_3', format=str(bineffarea) + 'E', unit='deg', array=np.expand_dims(PSFS3Run, 0))
        c10_psf = Column(name='SCALE', format=str(bineffarea) + 'E', unit='', array=np.expand_dims(1 / norm, 0))
        tbhdu_psf = pyfits.BinTableHDU.from_columns(
            [c1_psf, c2_psf, c3_psf, c4_psf, c5_psf, c6_psf, c7_psf, c8_psf, c9_psf, c10_psf])
    elif (PSFtype == "king"):
        c5_psf = Column(name='GAMMA', format=str(bineffarea) + 'E', unit='', array=np.expand_dims(PSFGamRun, 0))
        c5_psf = Column(name='SIGMA', format=str(bineffarea) + 'E', unit='deg', array=np.expand_dims(PSFSigRun, 0))
        tbhdu_psf = pyfits.BinTableHDU.from_columns([c1_psf, c2_psf, c3_psf, c4_psf, c5_psf])
    tbhdu_psf.header.set("EXTNAME", "PSF_2D", "name of this binary table extension ")
    tbhdu_psf.writeto('hess_psf_' + nrun + '.fits')
