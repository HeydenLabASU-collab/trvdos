# %%
import ctypes as ct
import MDAnalysis as mda
import MDA_unwrap_PBC as unwrap
from scipy.fft import fft, ifft, dct, idct
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# %%
# ctype data structures for rotatable bonds
class t_rotBond(ct.Structure):
    '''datatype describing a rotatable bonds and its satellites'''
    _fields_ = (("bondAtomIndices", ct.c_int32*2),
                ("sat0AtomIndices", ct.c_int32*4),
                ("nSat0", ct.c_int32),
                ("sat3AtomIndices", ct.c_int32*4),
                ("nSat3", ct.c_int32),
                ("inertia", ct.c_double),
                ("logInertia", ct.c_double),
                ("wOmegaBuffer", ct.POINTER(ct.c_double)))

class t_rotBondSet(ct.Structure):
    '''datatype describing a list of rotatable bonds (e.g., for a single residue)'''
    _fields_ = (("rotBonds", ct.POINTER(t_rotBond)),
                ("nRotBonds", ct.c_int32))

class t_rotBondSets(ct.Structure):
    '''datatype describing a list of rotatable bond sets (e.g., for each residue in a selection)'''
    _fields_ = (("rotBondSets", ct.POINTER(t_rotBondSet)),
                ("nRotBondSets", ct.c_int32))

# %%
#load the shared library with C routines
clib = ct.cdll.LoadLibrary("vdos.so")

# %%
#define argument types of function 'inertiaTensorAngularMomentumLabFrame' in imported library 'clib'
clib.inertiaTensorAngularMomentumLabFrame.argtypes = [
    #number of atoms of one residue
    ct.c_int32,
    #list of masses in one residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #atomic coordinates of one residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #atomic velocities of one residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #inertia Tensor of one residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #angular momentum of one residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
]
#define return type of function 'inertiaTensorAngularMomentumLabFrame' in imported library 'clib'
clib.inertiaTensorAngularMomentumLabFrame.restype = ct.c_int32

#define argument types of function 'applyAxes' in imported library 'clib'
clib.applyAxes.argtypes = [
    #list of moments of inertia, e.g., eigenvalues of inertia tensor
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    #list of rotational axes (transposed rmatrix), e.g., eigenvectors of inertia tensor
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #angular momentum in lab frame
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    #angular momentum in molecular frame (to be computed)
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    #weighted rotational velocity (sqrt(I) x omega) in molecular frame (to be computed)
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    #number of atoms of one residue
    ct.c_int32,
    #atomic coordinates of one residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #atomic velocities of one residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')
]
#define return type of function 'applyAxes' in imported library 'clib'
clib.applyAxes.restype = ct.c_int32

#define argument types of function 'getRotBonds' in imported library 'clib'
clib.getRotBonds.argtypes = [
    #sets of rotatable bonds for each residue
    ct.POINTER(t_rotBondSets),
    #number of residues
    ct.c_int32,
    #list of atom indices for dihedrals
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='C_CONTIGUOUS'),
    #number of dihedrals
    ct.c_int32,
    #pair of atom indices indicating the first and last atom index of each residue
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='C_CONTIGUOUS'),
    #number of correlation times to allocate arrays
    ct.c_int32
]
#define return type of function 'getRotBonds' in imported library 'clib'
clib.getRotBonds.restype = ct.c_int32

#define argument types of function 'analyzeRotBondsInResidue' in imported library 'clib'
clib.analyzeRotBondsInResidue.argtypes = [
    #set of rotatable bonds for one residue
    ct.POINTER(t_rotBondSet),
    #list of atom masses in residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #list of atom positions in residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #list of atom velocities in residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #time step
    ct.c_int32,
    #number of correlation times
    ct.c_int32
]
#define return type of function 'analyzeRotBondsInResidue' in imported library 'clib'
clib.analyzeRotBondsInResidue.restype = ct.c_int32

#define argument types of function 'corrRotBonds' in imported library 'clib'
clib.corrRotBonds.argtypes = [
    #sets of rotatable bonds for each residue
    ct.POINTER(t_rotBondSets),
    #array of correlation functions for each residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #start, i.e., index of reference time in circular array
    ct.c_int32,
    #number of correlation times
    ct.c_int32
]
#define return type of function 'corrRotBonds' in imported library 'clib'
clib.corrRotBonds.restype = ct.c_int32

#define argument types of function 'corrAtomsRes' in imported library 'clib'
clib.corrAtomsRes.argtypes = [
    #number of atoms in residue
    ct.c_int32,
    #residue index
    ct.c_int32,
    #number of residues
    ct.c_int32,
    #array of atomic massed for one residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #circular buffer of atomic velocites for one resiude
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='C_CONTIGUOUS'),
    #array for correlation function for one residue
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #start, i.e., index of reference time in circular array
    ct.c_int32,
    #number of correlation times
    ct.c_int32
]
#define return type of function 'corrAtomsRes' in imported library 'clib'
clib.corrAtomsRes.restype = ct.c_int32

# %%
class vdos:
    def __init__(self,sel,nCorr):
        self.sel = sel
        self.nCorr = nCorr
        self.nRes = sel.residues.n_residues
        # initialize lists and arrays
        #
        #Note: atMassLists & atVelBuffers are lists of numpy arrays
        #      => they are not a numpy arrays themselves
        self.atMassLists = [] # list of numpy arrays of atomic masses for each residue
        self.atVelBuffers = [] # list of numpy array buffers of atomic velocities for each residue
        for res in self.sel.residues:
            #Note: np.newaxis is used to make the array 2D
            self.atMassLists.append(res.atoms.masses[:,np.newaxis].astype('float64'))
            self.atVelBuffers.append(np.zeros((nCorr, len(res.atoms), 3), dtype = np.float64))
        # numpy array of residue masses from MDAnalysis
        self.resMassList = sel.residues.masses
        # numpy array buffers for COM position, COM velocity, angular momentum, and (sorted) moments of inertia
        self.COMposBuffer = np.zeros((nCorr, self.nRes, 3), dtype = np.float64)
        self.COMvelBuffer = np.zeros((nCorr, self.nRes, 3), dtype = np.float64)
        self.lastEvecs = np.zeros((self.nRes, 3, 3), dtype = np.float64)
        self.wOmegaBuffer = np.zeros((nCorr, self.nRes, 3), dtype = np.float64)
        self.inertia = np.zeros((self.nRes, 3), dtype = np.float64)
        self.logInertia = np.zeros((self.nRes, 3), dtype = np.float64)
        # time and frequency axes for correlation functions and VDoS
        self.tau = np.zeros(nCorr, dtype = np.float64)
        self.wavenumber = np.zeros(nCorr, dtype = np.float64)
        # numpy arrays for correlation functions and VDoS
        # -> rigid body translation
        self.trCorr = np.zeros((nCorr, self.nRes, 3), dtype = np.float64)
        self.trVDoS = np.zeros((nCorr, self.nRes, 3), dtype = np.float64)
        # -> rigid body rotation
        self.rotCorr = np.zeros((nCorr, self.nRes, 3), dtype = np.float64)
        self.rotVDoS = np.zeros((nCorr, self.nRes, 3), dtype = np.float64)
        # -> rotatable bonds
        self.rotBondCorr = np.zeros((nCorr, self.nRes), dtype = np.float64)
        self.rotBondVDoS = np.zeros((nCorr, self.nRes), dtype = np.float64)
        # -> total VDoS
        self.totalCorr = np.zeros((nCorr, self.nRes), dtype = np.float64)
        self.totalVDoS = np.zeros((nCorr, self.nRes), dtype = np.float64)
        # initialize counter for normalization
        self.corrCnt = np.zeros(self.nRes, dtype = int)
        # flag for normalization
        self.normalized = 0
        # initialize rotatable bonds for selection (see definition of getRotBonds)
        self.rotBondSets = self.getRotBonds()
        
    def processStep(self,tStep,time):
        '''process a single time step of the simulation'''
        idx = tStep % self.nCorr
        if tStep < self.nCorr:
            self.tau[tStep] = time
        pos = []
        vel = []
        r = 0
        #compute center of mass position and velocity
        for res in self.sel.residues:
            pos.append(res.atoms.positions.astype('float64'))
            vel.append(res.atoms.velocities.astype('float64'))
            #store all atomic velocities for current residue in buffer for total corr and VDoS
            self.atVelBuffers[r][idx] = vel[-1]
            #self.COMposBuffer[idx,r] = res.atoms.center_of_mass() # => too slow & no equivalent for velocity
            self.COMposBuffer[idx,r] = np.sum(self.atMassLists[r] * pos[-1] , axis = 0) / self.resMassList[r]
            self.COMvelBuffer[idx,r] = np.sum(self.atMassLists[r] * vel[-1], axis = 0) / self.resMassList[r]
            r += 1
        r = 0
        #remove center of mass position and velocity before computing angular momentum
        for res in self.sel.residues:
            pos[r] = pos[r] - self.COMposBuffer[idx,r]
            vel[r] = vel[r] - self.COMvelBuffer[idx,r]
            r += 1
        #compute weighted rotational velocity (sqrt(I) x omega) and accumulate inertia I for each residue
        #automatically removes rigid body rotation from velocities
        self.wOmegaBuffer[idx] = self.wOmega(tStep,pos,vel)
        #compute angular momenta for rotatable bonds
        if self.rotBondSets.nRotBondSets > 0:
            for r in range(self.nRes):
                error = clib.analyzeRotBondsInResidue(
                    ct.pointer(self.rotBondSets.rotBondSets[r]),
                    self.atMassLists[r],
                    pos[r],
                    vel[r],
                    tStep,
                    self.nCorr
                )
        #if sufficient data is available in buffers, compute correlation functions
        if tStep >= self.nCorr - 1:
            self.calcCorr(tStep+1)
            
    def calcCorr(self,start):
        '''compute correlation functions for all data in buffers'''
        #compute time correlation function for COM translation and rotation (for each residue)
        for i in range(self.nCorr):
            j = start % self.nCorr
            k = (j + i) % self.nCorr
            self.trCorr[i] += self.COMvelBuffer[j] * self.COMvelBuffer[k]
            self.rotCorr[i] += self.wOmegaBuffer[j] * self.wOmegaBuffer[k]
        #compute time correlation function for rotatable bonds (for each residue, averaged over rotatable bonds)
        error = clib.corrRotBonds(
            ct.pointer(self.rotBondSets),
            self.rotBondCorr,
            start,
            self.nCorr)
        #compute time correlation function for all atoms (for each residue, averaged over atoms)
        #for i in range(self.nCorr):
            # j = start % self.nCorr
            # k = (j + i) % self.nCorr
            # for r in range(self.nRes):
            #     self.totalCorr[i][r] += np.sum(self.atMassLists[r] * self.atVelBuffers[r][j] * self.atVelBuffers[r][k])
        for r in range(self.nRes):
            error = clib.corrAtomsRes(
                len(self.atMassLists[r]),
                r,
                self.nRes,
                self.atMassLists[r],
                self.atVelBuffers[r],
                self.totalCorr,
                start,
                self.nCorr
            )
        self.corrCnt += 1

    def normalize(self):
        '''normalize correlation functions by number of data points'''
        if self.normalized == 0:
            for i in range(self.nRes):
                self.trCorr[:,i] *= self.resMassList[i] / self.corrCnt[i]
                self.rotCorr[:,i] /= self.corrCnt[i]
                self.inertia[i] /= self.corrCnt[i]
                self.logInertia[i] /= self.corrCnt[i]
                self.rotBondCorr[:,i] /= self.corrCnt[i]
                self.totalCorr[:,i] /= self.corrCnt[i]
            for i in range(self.rotBondSets.nRotBondSets):
                for j in range(self.rotBondSets.rotBondSets[i].nRotBonds):
                    self.rotBondSets.rotBondSets[i].rotBonds[j].inertia /= self.corrCnt[i]
                    self.rotBondSets.rotBondSets[i].rotBonds[j].logInertia /= self.corrCnt[i]
            self.normalized = 1
    
    def calcVDOS(self):
        '''compute vibrational density of states from time correlation functions'''
        period = (self.tau[1] - self.tau[0]) * (2 * self.nCorr - 1)
        wn0 = (1.0 / period) * 33.35641
        self.wavenumber = np.arange(0,self.nCorr) * wn0
        tmp1 = np.zeros(2 * self.nCorr - 1, dtype = np.float64)
        tmp2 = np.zeros(2 * self.nCorr - 1, dtype = np.float64)
        tmp3 = np.zeros(2 * self.nCorr - 1, dtype = np.float64)
        tmp4 = np.zeros(2 * self.nCorr - 1, dtype = np.float64)
        tmp5 = np.zeros(2 * self.nCorr - 1, dtype = np.float64)
        tmp6 = np.zeros(2 * self.nCorr - 1, dtype = np.float64)
        tmp7 = np.zeros(2 * self.nCorr - 1, dtype = np.float64)
        tmp8 = np.zeros(2 * self.nCorr - 1, dtype = np.float64)
        for i in range(self.nRes):
            for j in range(self.nCorr):
                tmp1[j] = self.trCorr[j][i][0]
                tmp2[j] = self.trCorr[j][i][1]
                tmp3[j] = self.trCorr[j][i][2]
                tmp4[j] = self.rotCorr[j][i][0]
                tmp5[j] = self.rotCorr[j][i][1]
                tmp6[j] = self.rotCorr[j][i][2]
                tmp7[j] = self.rotBondCorr[j][i]
                tmp8[j] = self.totalCorr[j][i]
            for j in range(1,self.nCorr):
                k = 2 * self.nCorr - j - 1
                tmp1[k] = tmp1[j]
                tmp2[k] = tmp2[j]
                tmp3[k] = tmp3[j]
                tmp4[k] = tmp4[j]
                tmp5[k] = tmp5[j]
                tmp6[k] = tmp6[j]
                tmp7[k] = tmp7[j]
                tmp8[k] = tmp8[j]
            tmp1 = fft(tmp1)
            tmp2 = fft(tmp2)
            tmp3 = fft(tmp3)
            tmp4 = fft(tmp4)
            tmp5 = fft(tmp5)
            tmp6 = fft(tmp6)
            tmp7 = fft(tmp7)
            tmp8 = fft(tmp8)
            for j in range(self.nCorr):
                self.trVDoS[j][i][0]   = tmp1[j].real
                self.trVDoS[j][i][1]   = tmp2[j].real
                self.trVDoS[j][i][2]   = tmp3[j].real
                self.rotVDoS[j][i][0]  = tmp4[j].real
                self.rotVDoS[j][i][1]  = tmp5[j].real
                self.rotVDoS[j][i][2]  = tmp6[j].real
                self.rotBondVDoS[j][i] = tmp7[j].real
                self.totalVDoS[j][i]   = tmp8[j].real

    # def inertiaTensor(self,pos):
    #     '''compute inertia tensor for each residue in selection'''
    #     I = np.zeros((self.nRes,3,3), dtype = np.float64)
    #     for r in range(self.nRes):
    #         masses = (self.atMassLists[r])[:,0]
    #         crd = (pos[r])
    #         I[r,0,0] = np.sum(masses * (crd[:,1]**2 + crd[:,2]**2))
    #         I[r,1,1] = np.sum(masses * (crd[:,0]**2 + crd[:,2]**2))
    #         I[r,2,2] = np.sum(masses * (crd[:,0]**2 + crd[:,1]**2))
    #         I[r,0,1] = -np.sum(masses * crd[:,0] * crd[:,1])
    #         I[r,0,2] = -np.sum(masses * crd[:,0] * crd[:,2])
    #         I[r,1,2] = -np.sum(masses * crd[:,1] * crd[:,2])
    #         I[r,1,0] = I[r,0,1]
    #         I[r,2,0] = I[r,0,2]
    #         I[r,2,1] = I[r,1,2]
    #     return I

    # def angMomLabFrame(self,pos,vel):
    #     '''compute angular momentum for each residue in selection'''
    #     #Note: atMassLists is a list of numpy arrays
    #     #      => atMassLists itself is not a numpy array
    #     #Note: pos and vel are lists of numpy arrays
    #     #      => pos and vel themselves are not numpy arrays
    #     #Note: center of mass coordinates are expected to have been substracted from pos
    #     #Note: center of mass velocities are expected to have been substracted from vel
    #     L = np.zeros((self.nRes,3), dtype = np.float64)
    #     for r in range(self.nRes):
    #         L[r] += np.sum(self.atMassLists[r] * np.cross(pos[r],vel[r]), axis = 0)
    #     return L

    # def rotAxes(self,I):
    #     '''compute principal axes of rotation for a given inertia tensor'''
    #     evals,evecsTmp = np.linalg.eig(I)
    #     #sort eigenvalues and eigenvectors
    #     #Purpose: ordering of the rotational axis needs to be consistent between time steps for correlation functions
    #     #Question: --> is this necessary or is the ordering of generated eigenvalues consistent???
    #     #              e.g., small to large eigenvalues: need to check docs for np.linalg.eig
    #     idx = evals.argsort()[::-1]   
    #     evals = evals[idx]
    #     evecsTmp = evecsTmp[:,idx]
    #     evecs = np.zeros((3,3), dtype = np.float64)
    #     evecs[0] = evecsTmp[0]
    #     evecs[1] = evecsTmp[1]
    #     evecs[2] = evecsTmp[2]
    #     return evals, evecs

    def rotAxes(self,I):
        '''compute principal axes of rotation for a given inertia tensor'''
        #use this version for python implementation of clib.applyAxes
        evals,evecs = np.linalg.eig(I)
        #transpose eigenvectors to make sure that evecs[0] is the first principal axis, etc.
        evecs=evecs.T
        #sort eigenvalues and eigenvectors
        #Purpose: ordering of the rotational axis needs to be consistent between time steps for correlation functions
        #Question: --> is this necessary or is the ordering of generated eigenvalues consistent???
        #              e.g., small to large eigenvalues: need to check docs for np.linalg.eig
        idx = evals.argsort()[::-1]   
        evals = evals[idx]
        evecs = evecs[idx]
        # #use this version for C implementation of clib.applyAxes
        # evals,evecsTmp = np.linalg.eig(I)
        # #sort eigenvalues and eigenvectors
        # #Purpose: ordering of the rotational axis needs to be consistent between time steps for correlation functions
        # #Question: --> is this necessary or is the ordering of generated eigenvalues consistent???
        # #              e.g., small to large eigenvalues: need to check docs for np.linalg.eig
        # idx = evals.argsort()[::-1]   
        # evecsTmp = evecsTmp[:,idx]
        # evecs = np.zeros((3,3), dtype = np.float64)
        # evecs[0] = evecsTmp[0]
        # evecs[1] = evecsTmp[1]
        # evecs[2] = evecsTmp[2]
        return evals, evecs

    def angMomMolFrame(self,L,evecs):
        '''rotate moment of inertia into molecular frame (defined by principal axes of rotation, i.e., evecs of inertia tensor)'''
        return np.dot(evecs,L)

    def omegaToLabFrame(self,omegaMol,evecs):
        '''rotate rotational velocity back into lab frame'''
        return np.dot(evecs.T,omegaMol)

    def wOmega(self,tStep,pos,vel):
        '''compute weighted rotational vemocities (sqrt(I) x omega) and moments of inertia in molecular frame for each residue in selection'''
        # I = self.inertiaTensor(pos)
        # L = self.angMomLabFrame(pos,vel)
        I = np.zeros((self.nRes,3,3), dtype = np.float64)
        L = np.zeros((self.nRes,3), dtype = np.float64)
        r = 0
        for res in self.sel.residues:
            error = clib.inertiaTensorAngularMomentumLabFrame(
                len(res.atoms),
                self.atMassLists[r],
                pos[r],
                vel[r],
                I[r],
                L[r]
            )
            r += 1
        Lmol = np.zeros((self.nRes,3), dtype = np.float64)
        wOmegaMol= np.zeros((self.nRes,3), dtype = np.float64)
        inertia = np.zeros((self.nRes,3), dtype = np.float64)
        for r in range(self.nRes):
            #eigenvalus are moments of inertia in molecular frame
            inertia[r],evecs = self.rotAxes(I[r])
            if tStep>0:
                if np.dot(evecs[0],self.lastEvecs[r][0]) < 0:
                    evecs[0] = -evecs[0]
                if np.dot(evecs[1],self.lastEvecs[r][1]) < 0:
                    evecs[1] = -evecs[1]
                if np.dot(evecs[2],self.lastEvecs[r][2]) < 0:
                    evecs[2] = -evecs[2]
            self.lastEvecs[r] = evecs
            #C implementation requires changes in self.rotAxes and is slower
            # error = clib.applyAxes(inertia[r],evecs,L[r],Lmol[r],wOmegaMol[r],len(res.atoms),pos[r],vel[r])
            #BEGIN section that can be replaced by clib.applyAxes
            #rotate angular momentum to molecular frame
            Lmol[r] = self.angMomMolFrame(L[r],evecs)
            #weighted rotational velocity in molecular frame: wOmega = sqrt(I) x omega
            wOmegaMol[r] = Lmol[r] / np.sqrt(inertia[r]).real
            #rotational velocity in radians: omega = L/I
            omegaMol = Lmol[r] / inertia[r]
            #rotate
            omegaLab = self.omegaToLabFrame(omegaMol,evecs)
            #subtract rigid body rotation from velocities
            vel[r] -= np.cross(omegaLab,pos[r])
            #END section that can be replaced by clib.applyAxes
            # # if np.any(evals <= 0.0):
            # #     print("Error: non-positive eigenvalue in inertia tensor")
            # #     print("Eigenvalues: ",evals)
            # # if np.any(isinstance(evals,complex)):
            # #     print("Error: complex eigenvalue in inertia tensor")
            # #     print("Eigenvalues: ",evals)
        self.inertia += inertia
        self.logInertia += np.log(inertia)
        return wOmegaMol

    def getRotBonds(self):
        '''construct list of rotatable bonds for selection'''
        dihed = self.sel.intra_dihedrals.indices
        resAtomRangeList = np.zeros((len(self.sel.residues),2), dtype = np.int32)
        i = 0
        for res in self.sel.residues:
            resAtomRangeList[i][0] = np.amin(res.atoms.indices)
            resAtomRangeList[i][1] = np.amax(res.atoms.indices)
            i += 1
        rotBondSets = t_rotBondSets()
        error = clib.getRotBonds(
            ct.pointer(rotBondSets),
            ct.c_int(len(self.sel.residues)),
            dihed,
            ct.c_int(len(dihed)),
            resAtomRangeList,
            self.nCorr
        )
        if error != 0:
            print(f'ERROR reported by \'getRotBonds\' function\n')
        return rotBondSets

# %%
TOPOL = "/Users/mheyden/Dropbox (ASU)/ASU-Research/POPC/run-NVE.tpr"
TRAJ = "/Users/mheyden/Dropbox (ASU)/ASU-Research/POPC/run-NVE_test.trr"
u = mda.Universe(TOPOL,TRAJ)
trees = unwrap.unwrap.buildTrees(u)
sel = u.select_atoms("resname POPC")
vdos = vdos(sel,200)

# %%
tStep = 0
for ts in tqdm(u.trajectory):
    u.trajectory.ts._pos = unwrap.unwrap.unwrap(u,trees)
    vdos.processStep(tStep,ts.time)
    tStep += 1

# %%
vdos.normalize()
vdos.calcVDOS()

# %%
averMass = np.sum(vdos.resMassList)/vdos.nRes
averInertia = np.sum(vdos.inertia, axis=0)/vdos.nRes
averLogInertia = np.sum(vdos.logInertia, axis=0)/vdos.nRes

averRotBondInertia = 0.0
averRotBondLogInertia = 0.0
averRotBondCount = 0
count = 0
for r in range(vdos.rotBondSets.nRotBondSets):
    averRotBondCount += vdos.rotBondSets.rotBondSets[r].nRotBonds
    for b in range(vdos.rotBondSets.rotBondSets[r].nRotBonds):
        averRotBondInertia += vdos.rotBondSets.rotBondSets[r].rotBonds[b].inertia
        averRotBondLogInertia += vdos.rotBondSets.rotBondSets[r].rotBonds[b].logInertia
        count += 1
averRotBondCount /= vdos.nRes
averRotBondInertia /= count
averRotBondLogInertia /= count

outFile = open("residueProperties.dat","w")
outFile.write("#Average Mass (g/mol):\n%20.6e\n" % averMass)
outFile.write("#Average Inertia (g/mol*A^2):\n%20.6e %20.6e %20.6e\n" % (averInertia[0],averInertia[1],averInertia[2]))
outFile.write("#Average Log(Inertia):\n%20.6e %20.6e %20.6e\n" % (averLogInertia[0],averLogInertia[1],averLogInertia[2]))
outFile.write("#Average Rotatable Bond Inertia (g/mol*A^2): (%d)\n%20.6e\n" % (averRotBondCount,averRotBondInertia))
outFile.write("#Average Rotatable Bond Log(Inertia): (%d)\n%20.6e\n" % (averRotBondCount,averRotBondLogInertia))
outFile.close()

# %%
averTransCorr = np.sum(vdos.trCorr, axis = 1) / vdos.nRes
averRotCorr = np.sum(vdos.rotCorr, axis = 1) / vdos.nRes
averRotBondCorr = np.sum(vdos.rotBondCorr, axis = 1) / vdos.nRes
averTotalCorr = np.sum(vdos.totalCorr, axis = 1) / vdos.nRes

outFile = open("corr.dat","w")
outFile.write("%-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s\n" % ("#Time (ps)","Translation x","Translation y","Translation z","Rotation x","Rotation y","Rotation z","Rotatable Bonds","Total"))
for i in range(vdos.nCorr):
    outFile.write("%-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e\n" % (vdos.tau[i],averTransCorr[i][0],averTransCorr[i][1],averTransCorr[i][2],averRotCorr[i][0],averRotCorr[i][1],averRotCorr[i][2],averRotBondCorr[i],averTotalCorr[i]))
outFile.close()

# %%
averTransVDoS = np.sum(vdos.trVDoS, axis = 1) / vdos.nRes
averRotVDoS = np.sum(vdos.rotVDoS, axis = 1) / vdos.nRes
averRotBondVDoS = np.sum(vdos.rotBondVDoS, axis = 1) / vdos.nRes
averTotalVDoS = np.sum(vdos.totalVDoS, axis = 1) / vdos.nRes

outFile = open("vdos.dat","w")
outFile.write("%-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s\n" % ("#Wavenumber (cm^-1)","Translation x","Translation y","Translation z","Rotation x","Rotation y","Rotation z","Rotatable Bonds","Total"))
for i in range(vdos.nCorr):
    outFile.write("%-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e\n" % (vdos.wavenumber[i],averTransVDoS[i][0],averTransVDoS[i][1],averTransVDoS[i][2],averRotVDoS[i][0],averRotVDoS[i][1],averRotVDoS[i][2],averRotBondVDoS[i],averTotalVDoS[i]))
outFile.close()
