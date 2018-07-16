from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
import pdb
from astropy.table import Table
import batman
from scipy.optimize import curve_fit
import yaml
import glob

np.random.seed(3)

HSTPeriod = 96. ## minutes

class snr_sim(object):
    """ 
    Class to store the SNR simulation
    """
    def __init__(self,paramFile='params/planet_params.yaml',interactive=False):
        self.pp = yaml.load(open(paramFile))
        self.paramFile = paramFile
        self.interactive = interactive ## Interactive mode? (ie. show plots on screen?)
        self.nmForFiles = self.pp['Description']['nmForFiles']
        
        self.nOrbits = self.pp['Instrument Params']['nOrbits']
        self.gaps = self.pp['Instrument Params']['gaps']
        self.startTime = self.pp['Instrument Params']['startTime']
        
        self.transType = self.pp['Instrument Params']['transitType']
        self.make_model()
        
        self.make_snr_table()
        
        self.make_lc()
        
    def make_model(self):
        """ Makes a Transit Model """
        
        Bp = self.pp['Transit Params']
        
        BMparams = batman.TransitParams()
        BMparams.per = Bp['Period'] * 24. * 60.  #orbital period (min)
        
        if self.transType == 'secondary':
            BMparams.t_secondary = 0.
            BMparams.t0 = BMparams.per * (-0.5)
            BMparams.fp = Bp['fp']
        else:
            BMparams.t0 = 0.         #time of inferior conjunction
        
        BMparams.rp = Bp['Rp/R*']  #planet radius (in units of stellar radii)
        BMparams.a = Bp['a/R*']   #semi-major axis (in units of stellar radii)
        BMparams.inc = Bp['inc']  #orbital inclination (in degrees)
        BMparams.ecc = Bp['ecc']  #eccentricity
        BMparams.w = Bp['w']      #longitude of periastron (in degrees)
        BMparams.limb_dark = Bp['limb_dark']  #limb darkening model
        BMparams.u = Bp['u']      #limb darkening coefficients [u1, u2]
        
        cadence = self.pp['Instrument Params']['exposureCadence']
        if self.gaps == True:
            timeList = []
            for oneOrbit in np.arange(self.nOrbits):
                timeList.append(np.arange(self.startTime + oneOrbit * HSTPeriod,
                                          self.startTime + 0.5 * HSTPeriod + oneOrbit * HSTPeriod,
                                          cadence))
                self.gapText = '_gaps'
                time = np.hstack(timeList)
        else:
            self.gapText = ''
            time = np.arange(self.startTime,self.startTime + HSTPeriod * self.nOrbits,
                             cadence)
        self.BMparams = BMparams
        
        self.BMmodel = batman.TransitModel(BMparams, time,
                                           transittype=self.transType) #initializes model
        self.time = time
        
        self.xModelShow = np.linspace(np.min(time),np.max(time),1024)
        
        # def fit_fun(t,t0,rp,u1,u2):
        #     params.t0 = t0
        #     params.rp = rp
        #     #params.a = a
        #     #params.inc = inc
        #     params.u = [u1,u2]
        #     m = batman.TransitModel(params, t)
        #     return m.light_curve(params)
    
        # paramNames = ['t0 (min)','Rp/R*','u1','u2']
        # #boundsLower = [-100,0,0.1,75,-1.,-1.]
        # #boundsUpper = [ 100,1.,20,90,1.,1.]
        # boundsLower = [-100,0,-1.,-1.]
        # boundsUpper = [ 100,1.,1.,1.]
        #
        # ysims = []
        # yFits = []
        # pFits = []
        # pFitErrs = []
        #
        # t = Table()
        # t['Wave'] = np.arange(0.38,0.73,0.05)
        # t['SNR'] = 1./(3e-4 * np.ones(len(t)))
        
    def make_snr_table(self,nBinPoints=None):
        """ Makes a table of SNR values """
        instrument = self.pp['Instrument Params']['instrument']
        
        dat = ascii.read(self.pp['Instrument Params']['snrFile'])
        if 'thermal' not in dat.colnames:
            dat['thermal'] = 0 ## if thermal flux not included (in optical say)
        
        if nBinPoints is None:
            npts = self.pp['Instrument Params']['nBinPoints']
        else:
            npts = nBinPoints
        
        dat['um'] = dat['wavelength'] / 1e4
        
        nSpecPix = len(dat)
        waveArr, snrArr = [], []
        
        ## total background counts
        backTot = dat['dark_counts'] + dat['sky_counts'] + dat['thermal']
        
        for binInd in np.arange(0,len(dat),npts):
            if binInd + npts < nSpecPix - 1:
                binPts = np.arange(npts) + binInd
                waveBin = np.mean(dat['um'][binPts])
                
                fluxBin = np.sum(dat['total_counts'][binPts])
                fluxBack = np.sum(backTot[binPts])
                ## Read noise is already squared in the input file (see online table)
                readNArr = dat['readnoise'][binPts]
                
                ns = np.sqrt(fluxBin + fluxBack + np.sum(readNArr))
                snr = fluxBin / ns
                
                snrArr.append(snr)
                waveArr.append(waveBin)
                
        t = Table()
        t['Wave'] = waveArr
        t['SNR'] = snrArr
        t['frac sigma'] = 1./t['SNR']
        
        self.binTabFull = t ## full points
        ## discard points with much lower SNR
        medSNR = np.median(t['SNR'])
        goodP = t['SNR'] > 0.5 * medSNR
        
        self.binTab = t[goodP]
        self.origSNR = dat
    
    def make_lc(self):
        """ Make light curves for all wavelengths """
        ysims, yFits = [], []
        for onePt in self.binTab:
           ymodel = self.BMmodel.light_curve(self.BMparams)
           ynoise = np.random.randn(len(self.time)) * onePt['frac sigma']
           ysim = ymodel + ynoise
           ysims.append(ysim)
           
           # popt,pcov = curve_fit(fit_fun,time,ysim,
           #                       bounds=(boundsLower,boundsUpper),
           #                       p0=[0.1,0.1,0.,0.])
           #                       #p0=[0.1,0.1,6.,90.,0.,0.])
           # sigma_approx = np.sqrt(np.diag(pcov))
           # best_fitPt = fit_fun(time,*popt)
           # yModelShow = fit_fun(xModelShow,*popt)
           bm = batman.TransitModel(self.BMparams, self.xModelShow,
                                    transittype=self.transType) #initializes model
           yModelShow = bm.light_curve(self.BMparams)
           yFits.append(yModelShow)
           # pFitArr, pFitErrArr = [], []
           # for oneParam, oneErr in zip(popt,sigma_approx):
           #     pFitArr.append(oneParam)
           #     pFitErrArr.append(oneErr)
           # pFits.append(pFitArr)
           # pFitErrs.append(pFitErrArr)
        self.ysims, self.yFits = ysims, yFits
           
    def plot(self):
        """ Plots the time series """
        offset = self.pp['Configuration']['offset']
        fig, ax = plt.subplots()
        for ind,onePt in enumerate(self.binTab):
            ysim = self.ysims[ind] - offset * ind
            yfit = self.yFits[ind] - offset * ind
            ax.plot(self.time,ysim,'o',label="{:.2f} um".format(onePt['Wave']))
            
            ax.plot(self.xModelShow,yfit,label='')
        
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Normalized Flux')
        ax.legend()
        if self.interactive == True:
            fig.show()
        
        fig.savefig('plots/sim_tseries_{}.pdf'.format(self.nmForFiles))
     
# fitTab = Table()
# fitTab['Wave'] = t['Wave']
# pFits2D = np.array(pFits)
# pFitErrs2D = np.array(pFitErrs)
# for ind,oneParam in enumerate(paramNames):
#     fitTab[oneParam] = pFits2D[:,0]
#     fitTab[oneParam+' err'] = pFitErrs2D[:,0]
#
# fitTab.write('fit_results{}.csv'.format(gapText),overwrite=True)
#

#
# def time_spec():
#     fig, ax = plt.subplots()
#     ax.errorbar(t['Wave'],fitTab['t0 (min)']*60.,yerr=fitTab['t0 (min) err'] * 60.,
#                 fmt='o')
#     ax.set_xlabel('Wavelength ($\mu$m)')
#     ax.set_ylabel('Time offset (sec)')
#
#     fig.show()
#     fig.savefig('sim_t0_spectrum_wasp12{}.pdf'.format(gapText))
if __name__ == "__main__":
    paramList = glob.glob('params/*')
    for oneFile in paramList:
        s1 = snr_sim(oneFile)
        s1.plot()
    
