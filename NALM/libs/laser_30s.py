import numpy as np
from numpy import exp,abs
from numpy.fft import fft,ifft
from numpy import sum as summation

def coupler(uu,coupling):
    u1 = coupling**0.5*uu
    u2 = 1j*(1-coupling)**0.5*uu
    return u1,u2

def amplifier(uu, omega, halfstep1, dt, dz, nz, gamma, deltawg,Pp0,Psatp,Psats,alpha1p,alpha1s,alphap,alphas,mu,twindow):
	nt = len(uu)
	Ttwo=2/deltawg
	u1 = uu
	halfe = exp(halfstep1*dz/2)
	omegasq =1- (omega*omega)*(Ttwo*Ttwo)

	Pp=Pp0
	for iz in range(nz):
		u_abs = abs(u1)
		u_abs *= u_abs        
		Pen = summation(u_abs)*dt
		Ppz = ((alpha1p*Pp/Psatp+alpha1p*(Pen/twindow)/(mu*Psats))/(1+(Pen/twindow)/Psats+Pp/Psatp))*Pp-alphap*Pp;
    
    #%-----------------------------
    
    #%g1=g0/(1+Pen/Penergy);
		g1=(alpha1s*mu*Pp/Psatp + alpha1s*(Pen/twindow)/Psats)/(1+(Pen/twindow)/Psats +Pp/Psatp)-alphas
		g=g1*omegasq
    #%g=g1./(1+(w.^2)*(Ttwo^2));
    
    #%------------------------------

		#halfstep=halfstep1+g/2    
		halfstep = halfe*exp(g*dz/4)
    #%-------------------------------        
		uhalf = ifft(halfstep*(fft(u1)))
		uv = uhalf*exp(1j*gamma*(u_abs )*dz)
		ufft = halfstep*(fft(uv)) 
		uv = ifft(ufft)
		u1 = uv
		Pp=Pp +Ppz*dz
	#print ('Pen is', Pen)
	return u1

def fiber(u0,halfstep,dz,nz,gamma):
	halfstep = exp(halfstep*dz/2)

	u1 = u0
	Gamma = 1j*gamma*dz
	for iz in range(nz):
		u_abs = abs(u1)
		#u_abs *= u_abs
		ufft = fft(u1)
		uhalf = ifft(halfstep*ufft)
		uv = uhalf*exp(Gamma*(u_abs*u_abs))
		uv = fft(uv)
		ufft = halfstep*uv
		uv = ifft(ufft)
		u1 = uv


	u1=uv
	return u1
	
def combiner(U1,U2,coupling):

	
	
	return coupling**0.5*U1 + 1j*(1-coupling)**0.5*U2
	
def outcoupler(Uin,loss):
	

	Uref=((1-loss)**0.5)*Uin
	Uout=(loss**0.5)*Uin
	return Uref, Uout

def amplifier2(c0,omega,halfstep1,dt,dz,nz,gamma,deltawg,Pp0,Psatp,Psats,alpha1p,alpha1s,alphap,alphas,mu,twindow):
	nt = len(c0)
	Ttwo=2/deltawg
	u1 = c0
	omegasq =1- (omega*omega)*(Ttwo*Ttwo)
	Pp=Pp0
	halfe = exp(halfstep1*dz/2)
	for _ in range(nz):
		u_abs = abs(u1)
		u_abs *= u_abs
		Pen=summation(u_abs)*dt
		#Pen=sum(u_abs*u_abs)*dt
		Ppz=-((alpha1p*Pp/Psatp+alpha1p*(Pen/twindow)/(mu*Psats))/(1+(Pen/twindow)/Psats+Pp/Psatp))*Pp+alphap*Pp
    
    #%-----------------------------
    
    #%g1=g0/(1+Pen/Penergy);
		g1 =(alpha1s*mu*Pp/Psatp + alpha1s*(Pen/twindow)/Psats)/(1+(Pen/twindow)/Psats +Pp/Psatp)-alphas
		g=g1*omegasq
    #%g=g1./(1+(w.^2)*(Ttwo^2));
    
    #%------------------------------

		#halfstep=halfstep1+g/2    
		halfstep = halfe*exp(g*dz/4)#exp(halfstep*dz/2)
		temp = halfstep*(fft(u1))
    #%-------------------------------        
		uhalf = ifft(temp)
		#uv = uhalf*exp(1j*gamma*(u_abs*u_abs )*dz)
		uv = uhalf*exp(1j*gamma*(u_abs)*dz)        
		ufft = halfstep*(fft(uv))
		uv = ifft(ufft)
		u1 = uv
		Pp=Pp +Ppz*dz
	return u1
