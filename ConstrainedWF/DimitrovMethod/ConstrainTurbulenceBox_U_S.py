from IPython import display
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import os 
import scipy.io
from scipy import stats
import scipy.io as sio
#from spectrum import nextpow2
#from MannTensor import MannTensor
from MannTensor_lib import MannTensor,TrapezoidalSum2D
import cmath
from scipy.interpolate import interp1d

def nextpow2(x):
    """returns the smallest power of two that is greater than or equal to the
    absolute value of x.

    This function is useful for optimizing FFT operations, which are
    most efficient when sequence length is an exact power of two.

    :Example:

    .. doctest::

        >>> from spectrum import nextpow2
        >>> x = [255, 256, 257]
        >>> nextpow2(x)
        array([8, 8, 9])

    """
    res = np.ceil(np.log2(x))
    return res.astype('int')  #we want integer values only but ceil gives float



def ConstrainTurbulenceBox_U_S(Input):

    
    ## ASSIGN INPUT VARIABLES
    Udatafile = Input['Udatafile']
    Constraints = Input['Constraints']
    
    ## ASSIGN INPUT VARIABLES

    if not('CalculateSpectrum' in Input.keys()):    
        CalculateSpectrum = 1
    else:
        CalculateSpectrum = Input['CalculateSpectrum']
    
    if CalculateSpectrum == 0:
        CrossSpectrumFilenames = Input['CrossSpectrumFilenames']
    
    if not('CrossSpectrumComponent1' in Input.keys()):
        CrossSpectrumComponent1 = 1
        CrossSpectrumComponent2 = 1
    else:
        CrossSpectrumComponent1 = Input['CrossSpectrumComponent1']
        CrossSpectrumComponent2 = Input['CrossSpectrumComponent2']
    
    Gamma = Input['Gamma']
    L = Input['L']
    if not('alphaepsilon' in Input.keys()):
        alphaepsilon = 1
    else:
        alphaepsilon = Input['alphaepsilon']
    
    if not('SpectrumSteps' in Input.keys()):
        kdepth = 120
    else:
        kdepth = Input['SpectrumSteps']
    
    if not('PrintDetailedOutput' in Input.keys()):
        PrintOutput = False # 999 use this to inspect the code
    else:
        PrintOutput = Input['PrintDetailedOutput']
    
    if not('MaxIter' in Input.keys()):
        MaxIter = 0
    
    Nx = Input['Nx']
    Ny = Input['Ny']
    Nz = Input['Nz']
    T = Input['T']
    Umean = Input['Umean']
    BoxWidth = Input['BoxWidth']
    BoxHeight = Input['BoxHeight']
    
    if not('SaveToFile' in Input.keys()):
        SaveToFile = 0
    else:
        SaveToFile = Input['SaveToFile']
        if ~('OutputFileName' in Input.keys()):
            OutputFileName = ('%s_c.bin' %Udatafile[0:-4])
        else:
            OutputFileName = Input['OutputFileName']
            
    
    ##BASIC GEOMETRY
    
    # dx = (Umean*T)/Nx
    dy = BoxWidth/Ny
    dz = BoxHeight/Nz
    
    ## SPECTRAL INPUT, CROSS-CORRELATION STRUCTURE
    # One main thing to check is [1] and [0] indexing.. differences between 
    # matlab and python 
    
    fs = Nx/T
    dt = 1/fs
    
    k1range = [np.max([np.log10(1e-3/L)+np.spacing(2), np.log10(np.pi/(Nx*Umean*dt))]), np.min([np.log10(1e3/L)-np.spacing(2), np.log10(np.pi/(Umean*dt))])]
    # k1range = [log10(pi/(Nx*Umean*dt)), log10(pi/(Umean*dt))]
    # krange = [log10(1e-4/L) log10(pi/(Umean*dt))]
    krange = [np.min([np.log10(1e-3/L), np.log10(np.pi/(Umean*T))]), np.max([np.log10(1e3/L), np.log10(np.pi/(Umean*dt))])]
    
    Nmin = 2**nextpow2(2*np.pi/(10**k1range[0]*Umean*T/Nx))
    Nfactor = Nmin/Nx
    Nx_s = Nfactor*Nx
    # T_s = Nfactor*T
    
    L0 = np.arange(1,Nx_s/2,1)
    f0 = L0/(dt*Nx_s)
    tx = dt*np.arange(0,(Nx),1)
    # tx_s = dt*(0:(Nx_s-1))
    
    ksim = 2*np.pi*f0/Umean
    
    k1points = 64
    # kdepth = 240
    kmin = ksim[0]
    kmax = ksim[-1:]
    kvarrange = [np.log10(kmin)-np.spacing(2), np.log10(kmax)+np.spacing(2)]    
    
    k1 = 10**np.linspace(kvarrange[0],kvarrange[1],k1points)
    
    if CalculateSpectrum > 0:
        print('Calculating Mann spectrum.. ')
    
        if MaxIter == 0:
            Sarray = np.zeros(shape=(len(k1),Ny,Nz))
            Scomponent = 10*CrossSpectrumComponent1 + CrossSpectrumComponent2
            for ik in range(0,len(k1)):   # range(0,len(k1)):range(0,len(k1)):        
                if PrintOutput:
                    print('ik=%s ' %ik)
                
                klog0 = np.linspace(krange[0],krange[-1:],kdepth)
                k0p = 10**klog0        
                k0n = -np.flip(10**klog0,axis=0)
                k2_p = np.concatenate((k0n,k0p), axis=0)
                k3_p = np.concatenate((k0n,k0p), axis=0)
                k2grid,k3grid=np.meshgrid(k2_p,k3_p)
                k1_in=k1[ik]*np.ones(k2grid.shape)
                k2_in=np.tile(k2_p,(len(k3_p),1))
                k3_in=np.tile(k3_p.T,(len(k2_p),1)).T
                PsiK = MannTensor(k1_in,k2_in,k3_in,Gamma,L,alphaepsilon,Scomponent)
                # ok till here ! PsiK is verified 
                for iY in range(0,Ny):
                    for iZ in range(0,Nz):
                        dyi = (iY+1-1)*dy
                        dzi = (iZ+1-1)*dz
                        # check PsiIJ is complex number 
                        PsiIJ = PsiK*np.exp(  cmath.sqrt(-1)*(k2grid*dyi + k3grid*dzi))
                        ksum = TrapezoidalSum2D(PsiIJ,k2_p,k3_p)
                        # correct until here" 
                        Sarray[ik,iY,iZ] = ksum
    
            Psi0 = Sarray[:,0,0]
            k1data = k1
            
        else:
            tol = 1e-4
            reltol = 5e-3
            Sarray = np.zeros(shape=(len(k1),Ny,Nz))
            Scomponent = 10*CrossSpectrumComponent1 + CrossSpectrumComponent2
    
            for ik in range(0,2):   # range(0,len(k1)):
                if PrintOutput:
                    print('ik = %s ' %(ik))
    
                klog0 = np.linspace(krange[0],krange[-1:],kdepth)
                k0p = 10**klog0        
                k0n = -np.flip(10**klog0)
                k2_p = np.concatenate((k0n,k0p), axis=0)
                k3_p = np.concatenate((k0n,k0p), axis=0)
                k2grid,k3grid=np.meshgrid(k2_p,k3_p)      
                PsiK = MannTensor(k1_in,k2_in,k3_in,Gamma,L,alphaepsilon,Scomponent)
                for iY in range(0,Ny):
                    for iZ in range(0,Nz):
                        dyi = (iY+1-1)*dy
                        dzi = (iZ+1-1)*dz
                        PsiIJ = PsiK*np.exp(  cmath.sqrt(-1)*(k2grid*dyi + k3grid*dzi))
                        ksum = TrapezoidalSum2D(PsiIJ,k2_p,k3_p)
                        Sarray[ik,iY,iZ] = ksum
    
    
                F1 = Sarray[ik,0,0]
    
                go = 1
                passcount = 0
        #         MaxIter = 1
    # ------------------------- 
    # code the remaining code ...  99999
    # -----------------------  
    
    #            while go == 1:
    #                passcount = passcount + 1
    #                F0 = F1     
    #                klog1 = np.sort([klog0, klog0(1:end-1) + diff(klog0)./2])
    #                k1p = 10.^klog1
    #                k1n = -fliplr(10.^klog1)
    #                k2 = [k1n k1p]
    #                k3 = [k1n k1p]
    #                [k2grid,k3grid] = meshgrid(k2,k3)        
    #                PsiK = MannTensor(k1(ik)*ones(size(k2grid)),ones(length(k3),1)*k2,k3'*ones(1,length(k2)),Gamma,L,alphaepsilon,Scomponent)
    #                for iY = 1:Ny
    #                    for iZ = 1:Nz
    #                        dyi = (iY-1)*dy
    #                        dzi = (iZ-1)*dz
    #                        PsiIJ = PsiK.*exp( sqrt(-1)*(k2grid.*dyi + k3grid.*dzi))   
    #                        ksum = TrapezoidalSum2D(PsiIJ,k2,k3)
    #                        Sarray(ik,iY,iZ) = ksum
    #                    end
    #                end
    #
    #                F1 = Sarray(ik,1,1)
    #                res = F1 - F0
    #                resrel = res/F1
    #                if PrintOutput
    #                    disp(['Absolute residual = ' num2str(res)])
    #                    disp(['Relative residual = ' num2str(resrel)])
    #                end
    #
    #                if abs(res) <= tol
    #                    if PrintOutput
    #                        disp(['Converged (abstol) after ' num2str(passcount) ' iterations'])
    #                    end
    #                    go = 0
    #                elseif abs(resrel) <= reltol
    #                    if PrintOutput
    #                        disp(['Converged (reltol) after ' num2str(passcount) ' iterations'])
    #                    end
    #                    go = 0
    #                else
    #                    if passcount >= MaxIter
    #                        if PrintOutput
    #                            disp(['Max number of iterations exceeded: ' num2str(passcount)])
    #                        end
    #                        go = 0
    #                    end        
    #                end
    #
    #                klog0 = klog1
    #            end
    #        end
    #
    #        Psi0 = Sarray(:,1,1)
    #        k1data = k1
    #    end
    #else
    #    Psi0 = MannCrossSpectrum_T(Gamma,L,alphaepsilon,krange,k1,CrossSpectrumComponent1,CrossSpectrumComponent2,0,0)
    #    k1data = 10.^linspace(log10(0.5e-4),log10(2e2),120)
    #end    
    
    PsiR = 10**(interp1d(np.log10(k1),np.log10(Psi0))(np.log10(ksim)))
    PsiR[np.isnan(PsiR)]=0
    
    # varU = 2*trapz(ksim,PsiR)
    MuS = np.mean(PsiR)
    
    
    # xrangeTwoside = (-Nx_s/2)*dx:dx:((Nx_s-1)/2)*dx
    # xrangeOneside = 0:dx:(Nx-1)*dx
    
    print('Calculating correlations.. ')
    
    ROneside = np.zeros(shape=(Nx,Ny,Nz))
    # RTwoside = zeros(length(xrangeTwoside),Ny,Nz)
    NormalizeRatio = 1
    
    for iY in range(0,Ny):
    #     iY
        for iZ in range(0,Nz):
            dyi = (iY+1-1)*dy
            dzi = (iZ+1-1)*dz
            if CalculateSpectrum == 1:
                SiRe = np.real(Sarray[:,iY,iZ])
            #else: # davide 9999 to code 
                #if (CrossSpectrumComponent1 == 1) and (CrossSpectrumComponent2 == 1):
                #    SiRe = load([CrossSpectrumFilenames num2str(Gamma) '_L' num2str(L) '_ae1_dy' num2str(dyi) '_dz' num2str(dzi) '_Re.txt'])
                #else
                #    SiRe = load([CrossSpectrumFilenames num2str(Gamma) '_L' num2str(L) '_ae1_dy'...
                #        num2str(dyi) '_dz' num2str(dzi) '_' num2str(CrossSpectrumComponent1)...
                #        num2str(CrossSpectrumComponent2) '_Re.txt'])
                #end
            #end
    #         SiIm = load([CrossSpectrumFilenames num2str(Gamma) '_L' num2str(L) '_ae1_dy' num2str(dyi) '_dz' num2str(dzi) '_Im.txt'])
    #         Sii = SiRe + sqrt(-1)*SiIm
            Sii = SiRe
            SiiR = alphaepsilon*(10**(interp1d(np.log10(k1data),np.log10(Sii))(np.log10(ksim))))  
            SiiR[np.isnan(SiiR)]=0
          
    #         Rii = (1/mean(SiiR))*abs(ifft([0, SiiR(2:end), 0, fliplr(conj(SiiR(2:end)))]))
    #         Rii = (1/mean(SiiR))*real(ifft([SiiR, fliplr(conj(SiiR))]))
            Rii = (1/MuS)*np.real(np.fft.ifft(np.concatenate((SiiR,np.flip(np.conj(SiiR),axis=0))) ))
            if (dyi == 0) and (dzi == 0):
                NormalizeRatio = 1/Rii[0]
            
            Rii = Rii*NormalizeRatio
    #         RiiTwoside = [Rii(Nx_s/2+1:end) Rii(1:Nx_s/2)]
            RiiOneside = Rii[0:Nx]
    #         RTwoside(:,iY,iZ) = RiiTwoside
            ROneside[:,iY,iZ] = RiiOneside
    
    
    ## GRID DEFINITION
    
    Yvalues = dy/2 + np.arange(0,(Ny),1)*dy
    Zvalues = dz/2 + np.arange(0,(Nz),1)*dz  
    
    ## LOAD SOURCE TURBULENCE VECTOR
    print('Loading turbulence box.. ')
    #s_nr=40
    fn = r'./%s' %Input['Udatafile']
    #u, v, w = [np.fromfile(fn % uvw, np.dtype('<f'), -1).reshape(Nx , Ny*Nz) for uvw in ['u', 'v', 'w']  ]
    #u, v, w = [np.fromfile(fn % uvw, np.dtype('<f'), -1) for uvw in ['u', 'v', 'w']  ]
    u = np.fromfile(fn, np.dtype('<f'), -1)
    Udata=u
    
    
    #Ufile = fopen(Udatafile)
    #Udata = fread(Ufile,'single')
    #fclose(Ufile)
    
    ## RESHAPE FROM TURBULENCE VECTOR TO TURBULENCE BOX
    
    Unorm = np.zeros(shape=(Nx,Ny,Nz))
    varUdata = np.var(Udata)
    muUdata = np.mean(Udata)
    stdUdata = np.sqrt(varUdata)
    
    for xloc in range(0,Nx):
    #     xloc
        for yloc in range(0,Ny):
            dataindex = Ny*Nz*(xloc) + Ny*(yloc)
            UperZ = Udata[dataindex:dataindex+Nz]
            Unorm[xloc,yloc,:] = (UperZ - muUdata)/stdUdata
    
    Constraints=np.asarray(Constraints)
    ConstraintValuesNorm = (Constraints[:,3] - muUdata)/stdUdata
    
    ## CONSTRAINT LOCATION 
    
    Clocx = interp1d(tx,np.arange(0,Nx,1),'nearest')(Constraints[:,0])
    Clocy = interp1d(Yvalues,np.arange(0,Ny,1),'nearest',fill_value='extrap')(Constraints[:,1])
    Clocz = interp1d(Zvalues,np.arange(0,Nz,1),'nearest',fill_value='extrap')(Constraints[:,2])
    
    Clocx=Clocx.astype(int)
    Clocy=Clocy.astype(int)
    Clocz=Clocz.astype(int)
    
    # Eliminate overlapping constraints
    ClocA = [Clocx, Clocy, Clocz]
    ClocA=np.vstack((Clocx,Clocy,Clocz)).T
    ClocA=ClocA.astype(int)
    # here is the mistake: 
    df1 = pd.DataFrame(ClocA, columns=['a','b','c'])
    ClocIndex=df1.sort_values(['a', 'b','c'], ascending=[True, True,True]).index
    Constraints = Constraints[ClocIndex]
    Clocx = Clocx[ClocIndex]
    Clocy = Clocy[ClocIndex]
    Clocz = Clocz[ClocIndex]
    ConstraintValuesNorm = ConstraintValuesNorm[ClocIndex]
    Nconstraints = len(Constraints[:,0])
    
    # Eliminate constraints too close to each other
    Xconstr = Constraints[:,0]*Umean
    Xdist = np.tile(Xconstr,(Nconstraints,1)).T - np.tile(Xconstr,(Nconstraints,1))
    Ydist = np.tile(Constraints[:,1],(Nconstraints,1)).T - np.tile(Constraints[:,1],(Nconstraints,1))
    Zdist = np.tile(Constraints[:,2],(Nconstraints,1)).T - np.tile(Constraints[:,2],(Nconstraints,1))
    
    Rdist = np.sqrt(Xdist**2 + Ydist**2 + Zdist**2)
    ValidDistIndex = np.ones((Rdist[:,0]).shape, dtype=bool)
    for i in range(0,Nconstraints):
        Rdisti = Rdist[0:i,i]
        Rdistindex = [(Rdisti[k] > 0 and Rdisti[k] < L/10) for k in range(0,len(Rdisti))] #np.find(Rdisti > 0 and Rdisti < L/10) davide 9999
        indices=[indices for indices, x in enumerate(Rdistindex) if x]
        ValidDistIndex[indices]= False
    
    ValidDistIndex=ValidDistIndex[:]
    Constraints = Constraints[ValidDistIndex]
    Clocx = Clocx[ValidDistIndex]
    Clocy = Clocy[ValidDistIndex]
    Clocz = Clocz[ValidDistIndex]
    ConstraintValuesNorm = ConstraintValuesNorm[ValidDistIndex]
    Nconstraints = len(Constraints[:,0])
    
    # here we are: 
    
    CpointsY = np.unique(Clocy)
    CpointsZ = np.unique(Clocz)
    
    CvectY = np.reshape(np.tile(CpointsY,(len(CpointsZ),1)).T,len(CpointsY)*len(CpointsZ),1)
    CvectZ = np.reshape(np.tile(CpointsZ,(len(CpointsY),1)),len(CpointsY)*len(CpointsZ),1)
    
    NCp = np.zeros((CvectY).shape)
    for iV in range(0,len(CvectY)):
        NCp[iV] = np.sum( np.logical_and(Clocy == CvectY[iV], Clocz == CvectZ[iV]))
    
    CPlocy = CvectY[NCp>0]
    CPlocz = CvectZ[NCp>0]
    # davide 99999 we need to scale it at some point 
    CPlocy_index=CPlocy-1
    CPlocz_index=CPlocz-1
    Cplane = np.vstack((np.asarray([Yvalues[k] for k in CPlocy_index.astype(int)]), np.asarray([Yvalues[k] for k in CPlocz_index.astype(int)]))).T
    NCpoints = len(Cplane[:,0])
    
    
    ## COVARIANCE STRUCTURE
    
    print('Assembling the covariance matrix..')
    
    # FIND COVARIANCES WITHIN THE SOURCE TURBULENCE BOX
    
    CorrCMann = np.zeros(shape=(Nconstraints,Nconstraints))
    CpointIndex = np.zeros(shape=(Nconstraints,1))
    
    for jCp in range(0,NCpoints):
        CpointIndex[ np.logical_and(Clocy == CPlocy[jCp], Clocz == CPlocz[jCp] )  ]= jCp
    
    CpointIndex=CpointIndex[:,0]
    # tic
    for iC in range(0,Nconstraints):
        if PrintOutput:
            if np.remainder(iC,100) == 0:
                print('iC = %s' %(iC))
    
    #     ti = Constraints(iC,1)
        xloci = Clocx[iC]
        yloci = Clocy[iC]
        zloci = Clocz[iC]
        xlocCij = np.abs(xloci-Clocx) #+ 1
        ylocCij = np.abs(yloci-Clocy) #+ 1
        zlocCij = np.abs(zloci-Clocz) #+ 1
        # as integer 
        xlocCij=xlocCij.astype(int)
        ylocCij=ylocCij.astype(int)
        zlocCij=zlocCij.astype(int)
    
        Corrij = np.zeros(shape=(Nconstraints,1))
        for jC in range(0,Nconstraints):
    #        Corrij = ROneside[xlocCij,ylocCij,zlocCij]
            Corrij[jC] = ROneside[xlocCij[jC],ylocCij[jC],zlocCij[jC]]
        
        CorrCMann[:,iC] = Corrij[:,0]
        CorrCMann[iC,:] = CorrCMann[:,iC].T
    
    ## APPLY CONSTRAINTS
    print('Calculating constrained field..')
    UconstrainedNorm = np.zeros(shape=(Nx,Ny,Nz))
    Ucontemporaneous = np.zeros(shape=(Nconstraints,1))
    # as integer 
    Clocx=Clocx.astype(int)
    Clocy=Clocy.astype(int)
    Clocz=Clocz.astype(int)
    
    for iC in range(0,Nconstraints):
        Ucontemporaneous[iC] = Unorm[Clocx[iC],Clocy[iC],Clocz[iC]]
    
    
    # 999 davide issue in A\B solution 
    #x1 = np.linalg.lstsq(CorrCMann, (ConstraintValuesNorm - Ucontemporaneous))[0]
    #x2 = np.linalg.solve(CorrCMann, (ConstraintValuesNorm - Ucontemporaneous))
    #CorrCMann.dot(x1)
    CConst=np.linalg.lstsq(CorrCMann,(ConstraintValuesNorm - Ucontemporaneous[:,0]))[0]
    #CConst_as_per_matlab=np.linalg.lstsq(CorrCMann_matlab,(ConstraintValuesNorm_matlab - Ucontemporaneous_matlab[:,0]))[0]
    CConst1=np.linalg.solve(CorrCMann,(ConstraintValuesNorm - Ucontemporaneous[:,0]))

    #aa=CorrCMann.dot(CConst)
    #"CConst = CorrCMann\(ConstraintValuesNorm - Ucontemporaneous)
    # cond(CorrCMann)
    #del CorrCMann
    
    
    if Nconstraints <= 18000:
        YZcount = 0
    
        for yloc in range(0,Ny):
            for zloc in range(0,Nz):
                xloc = np.arange(0,Nx,1) # 999 
                dyCov = np.abs(yloc - Clocy) +1
                dzCov = np.abs(zloc - Clocz) +1
                YZcount = YZcount + 1
                #if PrintOutput:
                   # print('YZcount = %s' %(YZcount))    
    
                Corrxyz = np.zeros(shape=(Nx,Nconstraints))
                for iC in range(0,Nconstraints):
                    dxCov = np.abs(xloc - Clocx[iC])+1
        #         Ri = Rvect(DXvect == dxCov(iC) & DYvect == dyCov(iC) & DZvect == dzCov(iC))
        # davide 9999, here is the error, in the selection of the elements for some array 
                    Corrxyz[:,iC] = ROneside[dxCov,dyCov[iC],dzCov[iC]]
            
            # davide 9999 to check also thisn one
            # CConst is slighty different 
            # Corrxyz is slightly different 
                CtermYZ = Corrxyz.dot(CConst)
                UconstrainedNorm[:,yloc,zloc] = Unorm[:,yloc,zloc] + CtermYZ
    
    #     YZcount = 0
    # 
    #     for yloc = 1:Ny
    #         for zloc = 1:Nz
    # #             yloc = 5
    # #             zloc = 17
    # #     for yloc = 16:17
    # #         for zloc = 16:17
    #             YZcount = YZcount + 1
    #             disp(['YZcount = ' num2str(YZcount)])
    #             Corrxyz = zeros(Nx,Nconstraints)
    #             for Cloc = 1:NCpoints
    #                 dyloc = abs(CPlocy(Cloc) - yloc) + 1
    #                 dzloc = abs(CPlocz(Cloc) - zloc) + 1               
    #                 CorrMannYZ = ROneside(:,dyloc,dzloc)
    #                 for xloc = 1:Nx
    #                     dxCov = abs(xloc-Clocx) + 1
    #                     Corrxyz(xloc,CpointIndex==Cloc) = CorrMannYZ(dxCov(CpointIndex==Cloc))
    #                 end
    #             end
    # 
    #             CtermYZ = Corrxyz*CConst
    #             UconstrainedNorm(:,yloc,zloc) = Unorm(:,yloc,zloc) + CtermYZ
    #         end
    #     end
        
   # else:
#        # davide 99999 to check it out too 
#        # try a simulation where this loops is econturerd. 
#        for xloc in range(0,Nx):
#        #     if mod(xloc,100) == 0
#            print('xloc = %s' %(xloc))
#        #     end
#            dxCov = np.abs(xloc-Clocx) + 1
#    
#            for yloc in range(0,Ny):
#                for zloc in range(0,Nz):
#                    Corrxyz = np.zeros(shape=(1,Nconstraints))
#                    for Cloc in range(0,NCpoints):
#                        dyloc = np.abs(CPlocy[Cloc] - yloc)
#                        dzloc = np.abs(CPlocz[Cloc] - zloc)
#                        # as integer 
#                        dyloc=dyloc.astype(int)
#                        dzloc=dzloc.astype(int)
#                        dxCov.dxCov.astype(int)
#        #                     CorrMannYZ = ROneside(:,dyloc,dzloc)             
#        #                     Corrxyz(CpointIndex==Cloc) = CorrMannYZ(dxCov(CpointIndex==Cloc))
#                        Corrxyz[0,CpointIndex==Cloc] = ROneside[dxCov[CpointIndex==Cloc],dyloc,dzloc]
#                    UconstrainedNorm[xloc,yloc,zloc] = Unorm[xloc,yloc,zloc] + Corrxyz[0,:]*CConst
#    
    
    #del Unorm
    #del Corrxyz
    # error in here
    Uconstrained = UconstrainedNorm*stdUdata + muUdata
    
    ## SAVE OUTPUT FILE
    if SaveToFile == 1:
        Uvect = np.zeros(shape=(Nx*Ny*Nz,1))
        for xloc in range(0,Nx):
        #     xloc
            for yloc in range(0,Ny):
                dataindex = Ny*Nz*(xloc) + Ny*(yloc)
                UperZ = Uconstrained[xloc,yloc,:]
                Uvect[dataindex:dataindex+Nz,0] = UperZ                        
        
        Uvect=Uvect.astype('single')
        F = open('%s' %Input['OutputFileName'],'wb') 
        Uvect[:,0].tofile(F)
        F.close()
    
    print('DONE!!')

