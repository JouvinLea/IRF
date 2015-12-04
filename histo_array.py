import numpy as np
import ROOT
import root_numpy as nr
from root_numpy import hist2array
from root_numpy.testdata import get_filepath

#Values: zenithal angle, offset, efficiency and Energy used to simulate the MCs 
enMC=[0.02, 0.03,0.05,0.08, 0.125, 0.2, 0.3, 0.5, 0.8, 1.25, 2,3,5,8,12.5, 20, 30, 50, 80, 125]
lnenMC=np.log10(enMC)
zenMC=[0, 18, 26,32,37,41,46,50,53,57,60,63,67,70]
effMC=[50, 60, 70, 80, 90, 100]
offMC=[0.0, 0.5,1.0, 1.5, 2.0, 2.5]


binEMC=len(enMC)
binzen=len(zenMC)
binoff=len(offMC)
bineff=len(effMC)

#Size of the 4D table where we will stock the effective area and the resolution for each true energy, zenithal angle, offset and efficiency of the MCs
TableArea=np.zeros((binEMC,binoff,binzen,bineff))
TableBiais=np.zeros((binEMC,binoff,binzen,bineff))
TableSigma=np.zeros((binEMC,binoff,binzen,bineff))

#Il y a normalement 20 energies MCs
#Dans les .root fait avec hap pour l effective area il y a une dimension en plus dans les energies (21 bins au lieu de 20), ca descend une energie plus bas que les MCs simules. Pour la resolution, il y a 41 bin au lieu de 21 donc il faut selectionner les bins en energie du historesolution qui correspondent a l energie des MCs.
#On prend un histoarea et un histo biais donnee pour trouver les indices car ce sera les meme indices pour tous les angles zenitaux, offset et efficaicte.
HESSCONFIG="/Users/jouvin/Desktop/these/WorkGAMMAPI/IRF"
coupure="std_north_1b"
filehistoarea=HESSCONFIG+"/"+coupure+"/gFixedEnergy_paris_0-8-8-8_CamOptimal_hwtrig_eff100/CollectionArea.root"
filehistoresol=HESSCONFIG+"/"+coupure+"/gFixedEnergy_paris_0-8-8-8_CamOptimal_hwtrig_eff100/EnergyResolution.root"
name_hist_area="EffArea_67deg_2.5off_eff100_FixedE"
name_hist_biais="Resol_Biais_67deg_2.5off_eff100_FixedE"
TFileArea=ROOT.TFile(filehistoarea)
TFileResol=ROOT.TFile(filehistoresol)
histoarea=TFileArea.Get(name_hist_area)
historesol=TFileResol.Get(name_hist_biais)
#IMPORTANT :comment selctionner energie
Ehistoarea=[ histoarea.GetXaxis().GetBinLowEdge(i) for i in range(1,22) ]
Ehistoresol=[ historesol.GetXaxis().GetBinLowEdge(i) for i in range(1,42) ]
ind_area=[ np.where( Ehistoarea < i)[0][-1] for i in lnenMC ]
ind_resol=[ np.where( Ehistoresol < i)[0][-1] for i in lnenMC ]

"""
for i in range (22):
    print i, " ", histoarea.GetXaxis().GetBinLowEdge(i)," ", histoarea.GetXaxis().GetBinCenter(i), " ", histoarea.GetXaxis().GetBinUpEdge(i)," " ,histoarea.GetBinContent(i)

for i in range (42):
    print i, " ", historesol.GetXaxis().GetBinLowEdge(i)," ", historesol.GetXaxis().GetBinCenter(i), " ", historesol.GetXaxis().GetBinUpEdge(i)," " ,historesol.GetBinContent(i)
""" 
    
for (ieff,eff) in enumerate(effMC):
    print eff
    for (ioff,off) in enumerate(offMC):
        print off
        for (izen,zen) in enumerate(zenMC):
            print zen
            filehistoarea=HESSCONFIG+"/"+coupure+"/gFixedEnergy_paris_0-8-8-8_CamOptimal_hwtrig_eff"+str(eff)+"/CollectionArea.root"
            filehistoresol=HESSCONFIG+"/"+coupure+"/gFixedEnergy_paris_0-8-8-8_CamOptimal_hwtrig_eff"+str(eff)+"/EnergyResolution.root"
            TFileArea=ROOT.TFile(filehistoarea)
            TFileResol=ROOT.TFile(filehistoresol)
            if(izen==0):
                name_hist_area="EffArea_00deg_"+str(off)+"off_eff"+str(eff)+"_FixedE"
                name_hist_biais="Resol_Biais_00deg_"+str(off)+"off_eff"+str(eff)+"_FixedE"
                name_hist_sigma="Resol_Sigma_00deg_"+str(off)+"off_eff"+str(eff)+"_FixedE"
            else:    
                name_hist_area="EffArea_"+str(zen)+"deg_"+str(off)+"off_eff"+str(eff)+"_FixedE"
                name_hist_biais="Resol_Biais_"+str(zen)+"deg_"+str(off)+"off_eff"+str(eff)+"_FixedE"
                name_hist_sigma="Resol_Sigma_"+str(zen)+"deg_"+str(off)+"off_eff"+str(eff)+"_FixedE"
            histoarea=TFileArea.Get(name_hist_area)
            histobiais=TFileResol.Get(name_hist_biais)
            histosigma=TFileResol.Get(name_hist_sigma)
            AreaArray=hist2array(histoarea)[ind_area]
            AreaBiais=hist2array(histobiais)[ind_resol]
            AreaSigma=hist2array(histosigma)[ind_resol]
            TableArea[:,ioff,izen, ieff]=AreaArray
            TableBiais[:,ioff,izen, ieff]=AreaBiais
            TableSigma[:,ioff,izen, ieff]=AreaSigma

np.savez("IRF.npz",TableArea=TableArea, TableBiais=TableBiais, TableSigma=TableSigma,enMC=enMC,lnenMC=lnenMC,zenMC=zenMC, offMC=offMC, effMC=effMC)
            
