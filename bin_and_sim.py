from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
import pdb
from astropy.table import Table
import batman
from scipy.optimize import curve_fit
import yaml
import glob
from scipy.interpolate import interp1d

np.random.seed(3)

HSTPeriod = 96. ## minutes
tauRamp = 90. ## minutes


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
        self.include_ramp = self.pp['Systematic Params']['includeRamp']
        
        self.depth_file = self.pp['Transit Params']['depth file']
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
            expList = []
            self.exposureStarts = []
            self.exposureEnds = []
            
            for oneOrbit in np.arange(self.nOrbits):
                timesThisExposure = np.arange(self.startTime + oneOrbit * HSTPeriod,
                                             self.startTime + 0.5 * HSTPeriod + oneOrbit * HSTPeriod,
                                             cadence)
                                           
                self.exposureStarts.append(np.min(timesThisExposure))
                self.exposureEnds.append(np.max(timesThisExposure))
                
                expList.append(np.ones(len(timesThisExposure),dtype=np.int) * oneOrbit)
                
                timeList.append(timesThisExposure)
            
            self.gapText = '_gaps'
            time = np.hstack(timeList)
            self.expArray = np.hstack(expList)
            self.nExposures = self.nOrbits
            
        else:
            self.gapText = ''
            time = np.arange(self.startTime,self.startTime + HSTPeriod * self.nOrbits,
                             cadence)
            self.exposureStarts = [np.min(time)]
            self.exposureEnds = [np.max(time)]
            self.expArray = np.zeros(len(time),dtype=np.int)
            self.nExposures = 1
        
        
        self.BMparams = BMparams
        
        self.BMmodel = batman.TransitModel(BMparams, time,
                                           transittype=self.transType) #initializes model
        self.time = time
        self.t_start = np.min(self.time)
        
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
        """ Makes a table of SNR values per exposure
        """
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
    
    def make_depth_table(self):
        """ Make a table of Depth versus wavelength """
        t = Table()
        t['Wave'] = self.binTab['Wave']
        inPoints = self.BMmodel.ds < 1. ## only count w/ halfway at the limb
        outPoints = self.BMmodel.ds >= 1.
        nIn = np.sum(inPoints)
        nOut = np.sum(outPoints)
        
        tdepthErr = []
        inErr = self.binTab['frac sigma'] / np.sqrt(nIn)
        outErr = self.binTab['frac sigma'] / np.sqrt(nOut)
        totErr = np.sqrt(inErr**2 + outErr**2)
        t['Depth Err (ppm)'] = totErr * 1e6
        t.write('output/sim_depths_{}.csv'.format(self.nmForFiles),overwrite=True)
    
    def y_systematics(self,time):
        """ Systematic model to add non-IID Gaussian white noise
        This systematic model is multiplied by the theoretical transit light curve
        Parameters
        ------------
        time: numpy array
            Time in minutes
        """
        
        if self.include_ramp == True:
            ysys = np.ones_like(time,dtype=np.float) * np.nan
            
            amp = self.pp['Systematic Params']['rampAmp']
            tau = self.pp['Systematic Params']['rampTScale']
            
            for oneExp in np.arange(self.nExposures):
                ptsExposure = ((time >= self.exposureStarts[oneExp]) & 
                              (time <= self.exposureEnds[oneExp]))
                
                if oneExp >= 1:
                    useAmp = amp * 0.5
                    ## Smaller ramp effect on subsequent orbits
                else:
                    useAmp = amp
                
                tZero = self.exposureStarts[oneExp]
                ysysExposure = (1. - useAmp * np.exp(-(time - tZero)/tau))
                ysys[ptsExposure] = ysysExposure[ptsExposure]
                
        else:
            ysys = np.ones_like(time,dtype=np.float)
        
        return ysys
    
    def make_lc(self):
        """ Make light curves for all wavelengths """
        ysims, yFits = [], []
        
        if self.depth_file is not None:
            columns = ['wavelength','depth sim','depth truth','depth err']
            depthTable = ascii.read(self.depth_file,names=columns)
            f = interp1d(depthTable['wavelength'],depthTable['depth truth'])
            yInterp = f(self.binTab['Wave'])
            self.binTab['Model Depth'] = yInterp
        
        for onePt in self.binTab:
            if 'Model Depth' in self.binTab.colnames:
                if self.transType == 'secondary':
                    self.BMparams.fp = onePt['Model Depth']
                else:
                    print("not yet implemented")
                    pdb.set_trace()
            
            ymodel = self.BMmodel.light_curve(self.BMparams)
            ysys = self.y_systematics(self.time)
            ynoise = np.random.randn(len(self.time)) * onePt['frac sigma']
            ysim = ymodel * ysys + ynoise
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
            ysysModel = self.y_systematics(self.xModelShow)
            yModelShow = bm.light_curve(self.BMparams) * ysysModel
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
            waveLabel = "{:.2f} um".format(onePt['Wave'])
            ax.plot(self.time,ysim,'o',label=waveLabel)
            
            ax.plot(self.xModelShow,yfit,label='')
            
            ax.text(np.max(self.time),np.max(yfit),waveLabel,
                    verticalalignment='center',
                    horizontalalignment='left')
        
        ax.set_xlim(self.pp['Configuration']['xRange'])
        ax.set_ylim(self.pp['Configuration']['yRange'])
        ax.set_xlabel('Time (min)')
        ax.set_ylabel('Normalized Flux')
        #ax.legend()
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
        s1.make_depth_table()
    
