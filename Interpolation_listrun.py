import numpy as np
import pyfits
from scipy import interpolate
import math
from pyfits import Column
import sys
import os
from glob import glob
def gauss(x,sigma, mean):
    f=1/(np.sqrt(2*math.pi)*sigma)*np.exp(-(x-mean)**2/(2*sigma**2))
    return f

"""
Commande a lancer pour pouvoir donner des arguments au scripts
"""
#%run Interpolation_listrun.py '/Users/jouvin/Desktop/these/WorkGAMMAPI/IRF/CrabEventList/Crab' 'std_north_1b'


PathListRun = sys.argv[1]
ListRunDirectory = glob(PathListRun+'/run*.fits')
RunNumber = [file.split('/')[-1][5:11] for file in ListRunDirectory]


#Load les info sur les MCs depuis la table d'IRF ou est stocke pour toutes les nergies, zenith,offset et efficacite des MCs la valeur des la surface efficiace, et du biais et sigma pour la resolution
PathTableIRF="/Users/jouvin/Desktop/these/WorkGAMMAPI/IRF/CrabEventList/Crab"
coupure=sys.argv[2]
IRF=np.load(PathTableIRF+"/IRF"+coupure+".npz")
IRFArea=IRF["TableArea"]
IRFSigma=IRF["TableSigma"]
IRFBiais=IRF["TableBiais"]
enMC=IRF["enMC"]
lnenMC=IRF["lnenMC"]
zenMC=IRF["zenMC"]
effMC=IRF["effMC"]
offMC=IRF["offMC"]
binoffMC=len(offMC)
binEMC=len(enMC)
binEreco=50
bineffarea=len(offMC)*len(enMC)
bineffresol=len(offMC)*len(enMC)*binEreco

#reverifier qu ils ont bien ca dans leur bin PA en low edge et upper edge
off_low=offMC
off_hi=offMC

#pour les extremites prendre le milieu des bin en log
binlnEMC=lnenMC[1:]-lnenMC[:-1]
#Pour le premier bin en energie pour defenr le edge low du bin on prend la demilargeur du premier bin
binlnEMClow=np.insert(binlnEMC,0,binlnEMC[0])
#Pour le dernier bin en energie pour defenr le edge up du bin on prend la demilargeur du dernier bin
binlnEMCup=np.insert(binlnEMC,-1,binlnEMC[-1])
#Retrouver
lnEMClow=lnenMC-binlnEMClow/2
lnEMCup=lnenMC+binlnEMCup/2
E_true_low=pow(10,lnEMClow)
E_true_up=pow(10,lnEMCup)

#Definition de Etrue/Ereco
lnEtrue_reco=np.linspace(-1,1,binEreco)
#Le tableau d energie reco en log ont tous la meme largeur de bin donc on prend le premier
binlnEtrue_reco=lnEtrue_reco[1]-lnEtrue_reco[0]
lnE_true_reco_low=lnEtrue_reco-binlnEtrue_reco/2
lnE_true_reco_up=lnEtrue_reco+binlnEtrue_reco/2
Etrue_reco=pow(10,lnEtrue_reco)
E_true_reco_low=pow(10,lnE_true_reco_low)
E_true_reco_hi=pow(10,lnE_true_reco_up)

for nrun in RunNumber:
    AreaRun=np.zeros((binoffMC,binEMC))
    ResolRun=np.zeros((binoffMC,binEreco,binEMC))
    namerun = "run_0"+nrun+"_std_north_1b_eventlist.fits"
    hdurun=pyfits.open(namerun)
    ZenRun=90-hdurun[1].header["ALT_PNT"]
    EffRun=hdurun[1].header["MUONEFF"]*100
    for (iEMC,EMC) in enumerate(enMC):
        for (ioff, off) in enumerate(offMC):
            #print ioff, " ", iEMC
            InterArea=interpolate.interp2d(effMC,np.cos(zenMC*math.pi/180),IRFArea[iEMC,ioff,:,:])
            InterBiais=interpolate.interp2d(effMC,np.cos(zenMC*math.pi/180),IRFBiais[iEMC, ioff,:,:])
            InterSigma=interpolate.interp2d(effMC,np.cos(zenMC*math.pi/180),IRFSigma[iEMC, ioff,:,:])
            AreaRun[ioff,iEMC]=InterArea(EffRun,np.cos(ZenRun*math.pi/180))
            BiaisRun=InterBiais(EffRun,np.cos(ZenRun*math.pi/180))
            SigmaRun=InterSigma(EffRun,np.cos(ZenRun*math.pi/180))
            ResolRun[ioff, : ,iEMC]=gauss(lnEtrue_reco,SigmaRun,BiaisRun)
            #etre sur que c est bien normalise
            norm=np.sum(ResolRun[ioff, : ,iEMC]*(E_true_reco_hi-E_true_reco_low))            
            print norm
            
            if(np.isnan(norm)):
                #print nrun
                #print ioff
                #print iEMC
                ResolRun[ioff, : ,iEMC]=0
                #print norm
            else:
                ResolRun[ioff, : ,iEMC]=ResolRun[ioff, : ,iEMC]/norm
                
                

    c1_area = Column(name='ENERG_LO', format=str(binEMC)+'E', unit='TeV', array=np.atleast_2d(E_true_low))
    c2_area = Column(name='ENERG_HI', format=str(binEMC)+'E', unit='TeV', array=np.atleast_2d(E_true_up))
    c3_area = Column(name='THETA_LO', format=str(binoffMC)+'E', unit='TeV', array=np.atleast_2d(off_low))
    c4_area = Column(name='THETA_HI', format=str(binoffMC)+'E', unit='TeV', array=np.atleast_2d(off_hi))
    c5_area = Column(name='EFFAREA', format=str(bineffarea)+'E', unit='TeV', array=np.expand_dims(AreaRun,0))
    c6_area = Column(name='EFFAREA_RECO', format=str(bineffarea)+'E', unit='TeV', array=np.expand_dims(AreaRun,0))
    tbhdu_area = pyfits.BinTableHDU.from_columns([c1_area,c2_area,c3_area,c4_area,c5_area,c6_area])
    for i in range(1,7):
        tbhdu_area.header.comments['TTYPE'+str(i)]='label for field '+str(i)
        tbhdu_area.header.comments['TFORM'+str(i)]='data format of field: 4-byte REAL'
        tbhdu_area.header.comments['TUNIT'+str(i)]='physical unit of field '

    tbhdu_area.header.set("EXTNAME","EFFECTIVE AREA", "name of this binary table extension ")
    tbhdu_area.header.set("TDIM5","("+str(binEMC)+","+str(binoffMC)+")")
    tbhdu_area.header.set("TDIM6","("+str(binEMC)+","+str(binoffMC)+")")
    tbhdu_area.header.set("LO_THRES","0.320357620716095","TeV")
    tbhdu_area.header.set("HI_THRES","31.1716079711914","TeV")
    #tbhdu_area.header["EXTNAME"]='EFFECTIVE AREA'
    tbhdu_area.writeto('hess_aeff_2d_'+nrun+'.fits')


    c1_resol = Column(name='ETRUE_LO', format=str(binEMC)+'E', unit='TeV', array=np.atleast_2d(E_true_low))
    c2_resol = Column(name='ETRUE_HI', format=str(binEMC)+'E', unit='TeV', array=np.atleast_2d(E_true_up))
    c3_resol = Column(name='MIGRA_LO', format=str(binEreco)+'E', unit='', array=np.atleast_2d(E_true_reco_low))
    c4_resol = Column(name='MIGRA_HI', format=str(binEreco)+'E', unit='', array=np.atleast_2d(E_true_reco_hi))
    c5_resol = Column(name='THETA_LO', format=str(binoffMC)+'E', unit='deg', array=np.atleast_2d(off_low))
    c6_resol = Column(name='THETA_HI', format=str(binoffMC)+'E', unit='deg', array=np.atleast_2d(off_hi))
    c7_resol = Column(name='MATRIX ', format=str(bineffresol)+'E', unit='TeV', array=np.expand_dims(ResolRun,0))
    tbhdu_resol = pyfits.BinTableHDU.from_columns([c1_resol,c2_resol,c3_resol,c4_resol,c5_resol,c6_resol,c7_resol])
    for i in range(1,8):
        tbhdu_resol.header.comments['TTYPE'+str(i)]='label for field '+str(i)
        tbhdu_resol.header.comments['TFORM'+str(i)]='data format of field: 4-byte REAL'
#        tbhdu_resol.header.comments['TUNIT'+str(i)]='physical unit of field '

    tbhdu_resol.header.set("EXTNAME","EDISP_2D", "name of this binary table extension ")
    tbhdu_resol.header.set("TDIM7","("+str(binEMC)+","+str(binEreco)+","+str(binoffMC)+")")
    #tbhdu_resol.header["EXTNAME"]='EFFECTIVE RESOL'
    tbhdu_resol.writeto('hess_edisp_2d_'+nrun+'.fits')

