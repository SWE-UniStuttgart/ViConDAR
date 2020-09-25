from scipy.special import factorial
import math 
import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad
import numpy as np
import scipy

# HypergeometricFunction(1/3,17/6,4/3,-(kL.^(-2)),1,50) ) )
#a=1/3 b=17/6 c=4/3 z=-(kL**(-2)) method=1 accuracy=50 

def HypergeometricFunction(a,b,c,z,method,accuracy):
    
    if not(method):
        method = 1
    
    if method == 1:
        if not(accuracy):
            npoints = 100
        else:
            npoints = accuracy
        
        
        n = np.arange(0,npoints+1,1)
        nfact = scipy.misc.factorial(n)
        Ga= math.gamma(a)
        Gb = math.gamma(b)
        Gc = math.gamma(c)
        Gcab = math.gamma(c - a - b)
        Gca = math.gamma(c - a)
        Gcb = math.gamma(c - b)    
        Gan = np.array([math.gamma(a + k) for k in n])
        Gbn = np.array([math.gamma(b + k) for k in n])
        Gcn = np.array([math.gamma(c + k) for k in n])
        ai = a
        bi = c - b
        ci = c
        Gai = math.gamma(ai)
        Gbi = math.gamma(bi)
        Gci = math.gamma(ci)
        Gain = np.array([math.gamma(ai + k) for k in n]) 
        Gbin = np.array([math.gamma(bi + k) for k in n])
        Gcin = np.array([math.gamma(ci + k) for k in n])
             
        zsize = z.shape    
        F = np.zeros(shape=(zsize[0],zsize[1]))
        
        for i in range(0,zsize[0]):
            for j in range(0,zsize[1]):
                           
                if (z[i,j]) == 1:
                    F[i,j] = Gc*Gcab/(Gca*Gcb)
                elif (z[i,j] >= -0.5 and z[i,j] < 1) or (z[i,j] > 1 and z[i,j] < 2):
                    zi = z[i,j]
                    zn = zi**n
    
                    Gprod = (Gan*Gbn*zn)/(Gcn*nfact)
                    # here we are !!! 
                    Gprod = Gprod[np.logical_and((~np.isnan(Gprod)), (~np.isinf(Gprod)))]
                    F[i,j] = (Gc/(Ga*Gb))*np.sum(Gprod)
                else:
                    zi = z[i,j]/(z[i,j] - 1)
                    q0 = (1 - z[i,j])**(-a)
                    zn = zi**n
    
                    Gprod = (Gain*Gbin*zn)/(Gcin*nfact)
                    Gprod = Gprod[np.logical_and((~np.isnan(Gprod)), (~np.isinf(Gprod)))]
                    F[i,j] = q0*(Gci/(Gai*Gbi))*np.sum(Gprod)            
    
                
    #             if abs(z[i,j]) == 1
    #                 Gc = gamma(c)
    #                 Gcab = gamma(c - a - b)
    #                 Gca = gamma(c - a)
    #                 Gcb = gamma(c - b)
    #                 F[i,j] = Gc.*Gcab./(Gca.*Gcb)
    #             elseif (z[i,j] >= -0.5 && z[i,j] < 1) || (z[i,j] > 1 && z[i,j] < 2)
    #                 zi = z[i,j]
    #                 n = 0:1:npoints
    #                 Gan = gamma(a + n)
    #                 Gbn = gamma(b + n)
    #                 Gcn = gamma(c + n)
    #                 zn = zi.^n
    #                 nfact = factorial(n)
    # 
    #                 Ga = gamma(a)
    #                 Gb = gamma(b)
    #                 Gc = gamma(c)
    #                 Gprod = (Gan.*Gbn.*zn)./(Gcn.*nfact)
    #                 Gprod = Gprod(~isnan(Gprod))
    #                 F[i,j] = (Gc/(Ga*Gb))*sum(Gprod)
    #             else
    #                 ai = a
    #                 bi = c - b
    #                 ci = c
    #                 zi = z[i,j]/(z[i,j] - 1)
    #                 q0 = (1 - z[i,j])^(-a)
    # 
    #                 n = 0:1:npoints
    #                 Gan = gamma(ai + n)
    #                 Gbn = gamma(bi + n)
    #                 Gcn = gamma(ci + n)
    #                 zn = zi.^n
    #                 nfact = factorial(n)
    # 
    #                 Ga = gamma(ai)
    #                 Gb = gamma(bi)
    #                 Gc = gamma(ci)
    #                 Gprod = (Gan.*Gbn.*zn)./(Gcn.*nfact)
    #                 Gprod = Gprod(~isnan(Gprod))
    #                 F[i,j] = q0*(Gc/(Ga*Gb))*sum(Gprod)            
    # 
    
        
    elif method == 2:
        
        if not(accuracy):
            npoints = 100
        else:
            npoints = accuracy
        
        
        F = np.zeros(shape=(1,len(z)))
        for i in range(1,len(z)):
            zi = z[i]
            n = np.arange(0,npoints+1,1)
            Gan = [math.gamma(a + k) for k in n]
            Gbn = [math.gamma(b + k) for k in n]
            Gcn = [math.gamma(c + k) for k in n]
            zn = zi**n
            nfact = scipy.misc.factorial(n)
    
            Ga = math.gamma(a)
            Gb = math.gamma(b)
            Gc = math.gamma(c)
            Gprod = (Gan*Gbn*zn)/(Gcn*nfact)
            Gprod = Gprod[~np.isnan(Gprod)]
            F[i] = (Gc/(Ga*Gb))*np.sum(Gprod)
        
    elif method == 3:
        integralfunc = lambda t: ( t**(b-1) ) * ( (1 - t)**(c - b - 1) )* ( (1 - t*z)**(-a) )
        #fff = lambda x: x**2+c
        integralvalue = integrate.quad(integralfunc, 0, 1)
    #    integral(integralfunc,0,1,'ArrayValued',true,'AbsTol',1e-12,'RelTol',1e-12)
    #     integralvalue = zeros(1,length(z))
    #     for i = 1:length(z)
    #         integralfunci = @(t) ( t.^(b-1) ) .* ( (1 - t).^(c - b - 1) ) .* ( (1 - t.*z(i)).^(-a) )
    #         integralvalue(i) = integral(integralfunci,0,1,'AbsTol',1e-2,'RelTol',1e-10)
    #     end
    
        F = (math.gamma(c)/math.gamma(b))/math.gamma(c-b)*integralvalue
        
    elif method == 4:
        
        if not(accuracy):
            npoints = 100
        else:
            npoints = accuracy
        
        F = np.zeros(shape=(1,len(z)))
        
        for i in range(1,len(z)):
            if np.abs(z[i]) == 1:
                Gc = math.gamma(c)
                Gcab = math.gamma(c - a - b)
                Gca = math.gamma(c - a)
                Gcb = math.gamma(c - b)
                F[i] = Gc*Gcab/(Gca*Gcb)
            else:
                zi = z[0,i]
                n = np.arange(0,npoints+1,1)
                qa = np.concatenate(([1], a + n[1:] - 1))
                qb = np.concatenate(([1], b + n[1:] - 1))
                qc = np.concatenate(([1], c + n[1:] - 1))
                
                qsuma = np.cumsum(np.log(qa))
                qsumb = np.cumsum(np.log(qb))
                qsumc = np.cumsum(np.log(qc))           
                neven = 2*np.floor(n/2)
                oddindex = np.remainder(n,2)
                znvect = (neven/2)*np.log(zi**2)
                nvect = [0, np.cumsum(np.log(n[1:]))]        
                allsum = qsuma + qsumb - qsumc + znvect - nvect
                zminusvect = zi*oddindex
                expsum = np.exp(allsum)*zminusvect
                F[i] = np.sum(expsum)
    
    return F

# ------------
#gamma_par=3 
#L=23 
#alphaepsilon=0.1 
#ElementChoice=1
# temp 
#k2_p = np.concatenate((k0n,k0p), axis=0)
#k3_p = np.concatenate((k0n,k0p), axis=0)
#k1=k1[ik]*np.ones(k2grid.shape),
#k2=np.tile(k2_p,(len(k3_p),1))
#k3=np.tile(k3_p.T,(len(k2_p),1)).T
#ElementChoice=11
#MannTensor(k1,k2,k3,gamma_par,L,alphaepsilon,ElementChoice)

def MannTensor(k1,k2,k3,gamma_par,L,alphaepsilon,ElementChoice):
    
    # once this work, use pandas dataframe to speed-up to process?
    
    #function [Phi,Phi11,Phi22,Phi33] = MannTensor(k1,k2,k3,gamma_par,L,alphaepsilon,ElementChoice)
# 99999    
#    if nargin == 6:
#        ElementChoice = 'All'
    
    k = np.sqrt( k1**2 + k2**2 + k3**2)
    kL = k*L
    # correct till here 9999 
    Beta = gamma_par/( (kL**(2/3)) * np.sqrt( HypergeometricFunction(1/3,17/6,4/3,-(kL**(-2)),1,50) ) )
    # Beta1 = gamma_par./exp( (2/3)*log(kL) + 0.5*log(HypergeometricFunction(1/3,17/6,4/3,-(kL.^(-2)),1,50)))
    
    k30 = k3 + Beta*k1
    k0 = np.sqrt(k1**2 + k2**2 + k30**2)
    kL0 = k0*L
    # correct till hre 
    C1 = Beta * (k1**2) * (k0**2 - 2*(k30**2) + Beta * k1 * k30) / ( (k**2)*(k1**2 + k2**2) )
    # slightly different in here k1 to 2. check 99999 
    C2 = (k2*(k0**2)/( (k1**2 + k2**2)**(3/2) ))*np.arctan2( Beta*k1*np.sqrt( k1**2 + k2**2),( (k0**2) - k30*k1*Beta ) )
    
    Ek0 = alphaepsilon*(L**(5/3))*(kL0**4)/((1 + kL0**2)**(17/6))
    
    zeta1 = C1 - (k2/k1)*C2
    zeta2 = (k2/k1)*C1 + C2
    
    
    # to do: implemente a switcher instead of the if,else conditions.
    if ElementChoice == 11:
        Phi = (Ek0/(4*np.pi*(k0**4)))*( k0**2 - k1**2 - 2*k1*k30*zeta1 + (k1**2 + k2**2)*(zeta1**2))
    elif ElementChoice == 22:
        Phi = (Ek0/(4*np.pi*(k0**4)))*( k0**2 - k2**2 - 2*k2*k30*zeta2 + (k1**2 + k2**2)*(zeta2**2))
    elif ElementChoice == 33:
        Phi = (Ek0/(4*np.pi*(k**4)))*(k1**2 + k2**2)
    elif ElementChoice ==12:
        Phi = (Ek0/(4*np.pi*(k0**4)))*( -k1*k2 - k1*k30*zeta2 - k2*k30*zeta1 + (k1**2 + k2**2)*zeta1*zeta2)
    elif ElementChoice == 13:
        Phi = (Ek0/(4*np.pi*(k0**2)*(k**2)))*( -k1*k30 + (k1**2 + k2**2)*zeta1)
    elif ElementChoice == 23:
        Phi = (Ek0/(4*np.pi*(k0**2)*(k**2)))*( -k2*k30 + (k1**2 + k**2)*zeta2)
    else: 
        Phi11 = (Ek0/(4*np.pi*(k0**4)))*( k0**2 - k1**2 - 2*k1*k30*zeta1 + (k1**2 + k2**2)*(zeta1**2))
        Phi22 = (Ek0/(4*np.pi*(k0**4)))*( k0**2 - k2**2 - 2*k2*k30*zeta2 + (k1**2 + k2**2)*(zeta2**2))
        Phi33 = (Ek0/(4*np.pi*(k**4)))*(k1**2 + k2**2)
    
        Phi12 = (Ek0/(4*np.pi*(k0**4)))*( -k1*k2 - k1*k30*zeta2 - k2*k30*zeta1 + (k1**2 + k2**2)*zeta1*zeta2)
        Phi13 = (Ek0/(4*np.pi*(k0**2)*(k**2)))*( -k1*k30 + (k1**2 + k2**2)*zeta1)
        Phi23 = (Ek0/(4*np.pi*(k0**2)*(k**2)))*( -k2*k30 + (k1**2 + k2**2)*zeta2)
    # 99999 to do 
   #     if nargout == 1:
   #         Phi = [Phi11, Phi22, Phi33, Phi12, Phi13, Phi23]                
   #     else:
   #         Phi = []
    return Phi #,Phi11,Phi22,Phi33


def TrapezoidalSum2D(f,x,y):
    xa = x[0:-1]
    xb = x[1:]
    ya = y[0:-1] 
    yb = y[1:] 
    dx = np.ones(shape=(len(xa),1))*(xb - xa)
    dx = np.tile((xb - xa),(len(xa),1))
    dy = np.tile((yb - ya),(len(ya),1)).T
    #dy = (yb - ya).T*np.ones(shape(1,len(ya)))
    darea = dx*dy
    fa = f[0:-1,0:-1]
    fb = f[0:-1,1:]
    fc = f[1:,0:-1]
    fd = f[1:,1:] 
    Int = np.sum( np.sum( darea*(fa + fb + fc + fd)/4 ) );
    return Int