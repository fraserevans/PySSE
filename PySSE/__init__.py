from .pysse_module import sse_evolve
from astropy import units as u
import astropy
import numpy as np
from tqdm import tqdm

#Evolutionary phases corresponding to each type code
_ktypesLong = {
    0 :  '0: Deeply or fully convective low mass MS star',
    1 :  '1: Main Sequence star',
    2 :  '2: Hertzsprung Gap',
    3 :  '3: First Giant Branch',
    4 :  '4: Core Helium Burning',
    5 :  '5: First Asymptotic Giant Branch',
    6 :  '6: Second Asymptotic Giant Branch',
    7 :  '7: Main Sequence Naked Helium star',
    8 :  '8: Hertzsprung Gap Naked Helium star',
    9 :  '9: Giant Branch Naked Helium star',
    10 : '10: Helium White Dwarf',
    11 : '11: Carbon/Oxygen White Dwarf',
    12 : '12: Oxygen/Neon White Dwarf',
    13 : '13: Neutron Star',
    14 : '14: Black Hole',
    15 : '15: Massless Supernova'
}

_ktypesShort = {0:'lowMS', 1:'MS', 2:'HG', 3:'RGB', 4:'CHeB', 5:'AGB1', 
                6:'AGB2', 7: 'HeMS', 8:'HeHB', 9:'HeGB', 10:'HeWD',
                11:'COWD', 12:'ONeWD', 13:'NS', 14:'BH', 15:'MR'}

class SSEstar:
    '''
    Class to evolve a star using SSE (Single Star Evolution) code.
    
    Explanation of the formulae and parameters can be found in 
        https://ui.adsabs.harvard.edu/abs/2000MNRAS.315..543H/abstract

    Original SSE code found here:
        https://ui.adsabs.harvard.edu/abs/2013ascl.soft03015H/abstract
    
    Parameters:
    ----------

    M: float or astropy mass quantity
        Initial mass of the star in solar masses.
        Assumed to be in units of Msun if supplied as float.
    tevo: float or astropy time quantity
        Time for which to evolve the star.
        Assumed to be in units of Myr if supplied as float.

    Optional Parameters:
    -------------------

    Z: float (default = 0.02)
        Initial TOTAL metallicity of the star.
    neta: float (default = 0.5)
        Reimers mass-loss coefficient.
    bwind: float (default = 0.0)
        Binary enhanced mass loss parameter.
    hewind: float (default = 0.5)
        Helium star mass loss parameter.
    sigma: float of astropy speed quantity (default = 190.0 km/s)
        Dispersion of the Maxwellion supernova kick velocity distribution.
    ifflag: int (default = 0)
        If 1, uses WD IFMR of HPE, 1995, MNRAS, 272, 800
    wdflag: int (default = 1)
        If 1, uses modified-Mestel cooling for WDs
    bhflag: int (default = 0)
        If 1, allows velocity kick at BH formation
    nsflag: int (default = 1)
        If 1, takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407
    mxns: float or astropy mass quantity (default = 3.0 Msun)
        Maximum neutron star mass.
    pts1 : float (default = 0.05)
        Reciprocal of number of timesteps used for MS phase
    pts2 : float (default = 0.01)
        Reciprocal of number of timesteps used for GB/CHeB/AGB/HeGB/phases
    pts3 : float (default = 0.02)
        Reciprocal of number of timesteps used for HG/HeMS phases
    returnall: bool (default = False)
        If True, returns all the timesteps of the star evolution.
        If False, returns only the final state of the star.

    Attributes:
    ----------
    time: astropy time quantity
        Age of the star since zero age main sequence (ZAMS).
    type: astropy dimensionless quantity
        Type of the star at the current time.
    Mo: astropy mass quantity
        Initial mass of the star
    mass: astropy mass quantity
        Mass of the star at the current time
    Lum: astropy luminosity quantity
        Luminosity of the star at the current time in solar luminosities
    radius: astropy length quantity
        Radius of the star at the current time
    Teff: astropy temperature quantity
        Effective temperature of the star at the current time
    mcore: astropy mass quantity
        Mass of the core of the star at the current time
    Rcore: astropy length quantity
        Radius of the core of the star at the current time
    menv: astropy mass quantity
        Mass of the envelope of the star at the current time
    Renv: astropy length quantity
        Radius of the envelope of the star at the current time
    epoch: astropy time quantity
        Epoch of the star at the current time
    spin: astropy speed quantity
        Spin of the star at the current time
    strtype: str
        String representation of the type of the star at the current time
    '''

    def __init__(self, M, tevo, Z=0.02, \
                 neta=0.5, bwind=0., hewind=0.5, sigma=190., 
                 ifflag = 0, wdflag = 1, bhflag = 0, nsflag = 1, 
                 mxns = 3.0, pts1 = 0.05, pts2 = 0.01, pts3 = 0.02,
                 returnall=False, **kwargs):

        attrs = {'type': u.dimensionless_unscaled, 'Mo': u.Msun, 
                'mass': u.Msun, 'Lum': u.Lsun, 'radius': u.Rsun, 
                'Teff': u.K, 'mcore': u.Msun, 'Rcore': u.Rsun, 'menv': u.Msun, 'Renv': u.Rsun, 'epoch': u.Myr, 'spin': u.km/u.s}

        def broadcast_args(**kwargs):
            #If any provided argument is a numpy array,
            #broadcast all other arguments to the same length

            # Determine the target length based on the FIRST array argument
            target_length = None
            for value in kwargs.values():
                if isinstance(value, list):
                    # If the value is a list, convert it to a numpy array
                    value = np.array(value)
                if isinstance(value, np.ndarray):
                    target_length_tmp = len(value)
                    
                    if target_length is None:
                        target_length = target_length_tmp
                    elif target_length != target_length_tmp:
                        raise ValueError("Array arguments must be same length")
            
            # If no array arguments are found, return the original kwargs
            if target_length is None:
                return kwargs

            # Broadcast scalars to arrays of the same length
            result = {}
            for key, value in kwargs.items():
                if isinstance(value, np.ndarray):
                    result[key] = value
                else:
                    result[key] = np.full(target_length, value)
            
            return result

        #Check if any of the arguments are astropy quantities
        #If so, convert them to the appropriate units for the SSE code
        if isinstance(M, astropy.units.Quantity):
            M = M.to(u.Msun).value
        if isinstance(tevo, astropy.units.Quantity):
            tevo = tevo.to(u.Myr).value
        if isinstance(mxns, astropy.units.Quantity):
            mxns = mxns.to(u.Msun).value
        if isinstance(sigma, astropy.units.Quantity):
            sigma = sigma.to(u.km/u.s).value

        #Check if any of the arguments are numpy arrays
        #If so, convert all other arguments to arrays of the same length
        args = broadcast_args(M=M, tevo=tevo, Z=Z, neta=neta, bwind=bwind,
                              hewind=hewind, sigma=sigma, ifflag=ifflag,
                              wdflag=wdflag, bhflag=bhflag, nsflag=nsflag,
                              mxns=mxns, pts1=pts1, pts2=pts2, pts3=pts3)

        if isinstance(args['M'], np.ndarray):

            #Initialize phases as an empty list of length len(args['M'])

            self.time, self.type, self.Mo, self.mass, self.Lum, \
                self.radius, self.Teff, self.mcore, self.Rcore, self.menv, \
                self.Renv, self.epoch, self.spin, self.strtype, self.phases \
                = [ [None] * len(args['M'])  for _ in range(15)]

            #Loop over each star and evolve it
            for i in tqdm(range(len(args['M'])), 
                        desc='Evolving stars', unit='star'):

                #Call the SSE code for each star.
                #The outputs are stored in a 2D array outs.
                #The timestamps where the evolutionary phases change
                #are stored in the arrays ts and phases.
                outs, phases, ts = sse_evolve(
                                        args['M'][i], args['Z'][i], 
                                        args['tevo'][i], args['neta'][i], 
                                        args['bwind'][i], args['hewind'][i], 
                                        args['sigma'][i], args['ifflag'][i], 
                                        args['wdflag'][i], args['bhflag'][i], 
                                        args['nsflag'][i], args['mxns'][i], 
                                        args['pts1'][i], args['pts2'][i], 
                                        args['pts3'][i])

                #Outputs will have many trailing zeros to accommodate
                #for the maximum number of timesteps since declaring arrays of 
                #ambiguous size is not possible in Fortran77.
                #Need to trim off some trailing zeros, therefore.
                self.time[i] = np.trim_zeros(outs[:, 0], 'b')*u.Myr

                #Assign the outputs to the attributes with appropriate units
                for j, name in enumerate(attrs.keys()):
                    if returnall:
                        #Outputs at all timesteps are recorded
                        getattr(self, name)[i] = outs[:len(self.time[i]),j+1] \
                                                        * attrs[name]
                    else:
                        #Outputs at the final timestep are recorded
                        getattr(self, name)[i] = outs[len(self.time[i])-1,j+1]\
                                                        * attrs[name]
                
                #Assign the string type of the star
                if returnall:
                    
                    self.strtype[i] = [_ktypesShort[k] 
                                    for k in self.type[i].value]
                else:                       
                    self.time[i] = self.time[i][-1]
                    self.strtype[i] = _ktypesShort[self.type[i].value]
                
                #Assign the phases and timestamps to the star as a length-N
                #list of nested length-2 lists of [phase, timestamp].
                #Again need to trim off some trailing zeros.
                self.phases[i] = np.array(list(zip(phases, ts)))
                self.phases[i] = self.phases[i][
                                    np.any(self.phases[i] != 0, axis=1)]

        else:
            #If the arguments are not numpy arrays, evolve a single star
            
            #Call the SSE code for the star.
            outs, phases, ts = sse_evolve(M, Z, tevo, neta, bwind, hewind,
                                          sigma, ifflag, wdflag, bhflag,
                                          nsflag, mxns, pts1, pts2, pts3)

            #Same as above, but can now assign the entire attributes at once

            self.time = np.trim_zeros(outs[:, 0], 'b')*u.Myr
            
            for i, name in enumerate(attrs.keys()):
                if returnall:
                    setattr(self, name, outs[:len(self.time), i+1] \
                            * attrs[name])
                else:
                    setattr(self, name, outs[len(self.time)-1, i+1] \
                            * attrs[name])

            if not returnall:
                self.time = self.time[-1]
                self.strtype = _ktypesShort[int(self.type)]
            else:
                self.strtype = [_ktypesShort[k] 
                                    for k in self.type.value]

            self.phases = np.array(list(zip(phases, ts)))
            self.phases = self.phases[ np.any(self.phases != 0, axis=1)]

        #If a single star is being evolved or only the final state
        #of an array of stars is being returned, convert the attributes
        #to simple astropy quantity arrays
        if (not isinstance(args['M'], np.ndarray) ) or (not returnall):

            for name in attrs.keys():
                setattr(self,name, u.Quantity(getattr(self, name)))

            self.time = u.Quantity(self.time)
            self.strtype = np.array(self.strtype)
        else:
            #Otherwise, attributes remain as lists of astropy quantity arrays
            pass
    def printtypes(self,k=np.arange(0, 16)):
        '''
        Helper function to descibe the types of stars in the SSE code.
        
        Parameters:
        ----------
        k: int or array-like (default = np.arange(0, 10))
            Index or indices of the types to show.
        '''
        if isinstance(k, int):
            print(_ktypesLong[k])
        else:
            for i in k:
                print(_ktypesLong[i])
