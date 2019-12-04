"""
===========================
Hammersley vs Random 
===========================

In this example we compare the tiling efficiency of just choosing random 
numbers for :math:`\\theta` and :math:`\\phi` tiling angles compared to using Hammersley Points
"""

import numpy as np
from scipy import fftpack
from matplotlib import pyplot as plt
import sys


def return_Hammersley_points_array( l, n, p):
    """
    l is the power of x to go out to p^m
    n is the maximun number of points
    p is the order of the Hammersley point, 1,2,3,4,... etc
    
    returns
    --------
    np.array of double
    
    """
    
    vvv = np.zeros(n)
    
    for m in range(n):
        m1=1*m
        if p == 1:
            vvv[m] =  m1/n
        else:        
            v = 0.0
           
            for j in range(l,-1,-1):
                num = m1//p**j
                
                if num > 0:
                    m1 -= num*p**j
                    v  += num / ( p ** (j+1) )
                    
            vvv[m]=v
                
    return(vvv)


def omega_cs( theta, phi, iso_cs=0.0, asymm_cs=100, eta_cs=1.0):
        
    return (iso_cs +0.5* asymm_cs*(3.0 * (np.cos(theta)**2) -1.0 - eta_cs*(np.sin(theta)**2)*np.cos( 2.0 * phi ))), np.sin(theta)

if __name__ == "__main__":
    
    # Define CSA powder pattern 
    
    # Principal components of the chemical shift shielding tensor
    
    s_zz = -120.0
    s_yy = -50.0
    s_xx =  100.0
    
    ncols = 256
    nrows = 256
    
    sw = 500.0
    dt = 1/sw
    time_axis = np.arange(ncols)*dt
    lb = 2.0
    axis_Hz = np.linspace( -sw/2., sw/2., ncols )
    
    # Check for Haeberlens convention
    
    iso_cs =(s_zz+s_yy+s_xx)/3.
    
    if abs(s_zz-iso_cs) >= abs(s_xx-iso_cs) and abs(s_xx-iso_cs) >= abs(s_yy-iso_cs):
        h_zz = s_zz
        h_yy = s_yy
        h_xx = s_xx
        
    elif abs(s_zz-iso_cs) < abs(s_xx-iso_cs) and abs(s_xx-iso_cs) >= abs(s_yy-iso_cs):
        h_zz = s_xx
        h_yy = s_yy
        h_xx = s_zz

    else:
        print("problem with assignment of cs tensors")
        sys.exit()
    
    asymm_cs = h_zz-iso_cs
    eta_cs   = (h_xx-h_yy)/(h_zz-iso_cs)
    
    # Calculate Hammersley Points and Powder pattern
    N_particles = 2**17
     
    theta = return_Hammersley_points_array(22, N_particles, 2) 
    phi   = return_Hammersley_points_array(22, N_particles, 3) 
    
    omega_ham, solid_angle_ham = omega_cs(theta*np.pi,2*np.pi*phi, eta_cs=eta_cs, iso_cs=iso_cs, asymm_cs=asymm_cs)
    
    # Calculate random Points 
    
    pts = np.random.rand(2,N_particles)
    
    theta = pts[0]*np.pi
    phi   = pts[1]*2.0*np.pi
    
    omega_rand, solid_angle_rand = omega_cs( pts[0]*np.pi,pts[1]*2.0*np.pi, eta_cs=eta_cs, iso_cs=iso_cs, asymm_cs=asymm_cs)
    
    # Propogate FID
    
    ccc_rand = np.cos( np.outer( 2*np.pi*omega_rand, time_axis ))
    sss_rand = np.sin( np.outer( 2*np.pi*omega_rand, time_axis ))
    
    ccc_ham = np.cos( np.outer( 2*np.pi*omega_ham, time_axis ))
    sss_ham = np.sin( np.outer( 2*np.pi*omega_ham, time_axis ))
    
    ## Apply solid angle weighting
    
    ccc_rand =  solid_angle_rand.reshape((N_particles,1)) *ccc_rand
    sss_rand =  solid_angle_rand.reshape((N_particles,1)) *sss_rand

    ccc_rand = ccc_rand.sum( axis=0)
    sss_rand = sss_rand.sum( axis=0)

    ccc_ham =  solid_angle_ham.reshape((N_particles,1)) *ccc_ham
    sss_ham =  solid_angle_ham.reshape((N_particles,1)) *sss_ham

    ccc_ham = ccc_ham.sum( axis=0)
    sss_ham = sss_ham.sum( axis=0)
    
    # create complex array, apply exponential linebroadeining, FFT
    
    fid_rand = ccc_rand*np.exp(-lb*time_axis) + 1J * sss_rand*np.exp(-lb*time_axis)
    spec_rand = fftpack.fft(fid_rand)
    spec_rand = fftpack.fftshift(spec_rand)
    
    fid_ham = ccc_ham*np.exp(-lb*time_axis) + 1J * sss_ham*np.exp(-lb*time_axis)
    spec_ham = fftpack.fft(fid_ham)
    spec_ham = fftpack.fftshift(spec_ham)
    
    
    # display Spectra side by side
    
    fig, axs = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    
    ax = axs.flatten()
    
    ax[0].plot(axis_Hz,(spec_rand.real/(spec_rand.real).max())  );
    ax[1].plot(axis_Hz,(spec_ham.real/(spec_ham.real).max())  );
    l,r = ax[0].get_xlim()
    ax[0].set_xlim(r,l)
    ax[1].set_xlim(r,l)

    ax[0].spines['top'].set_visible(False)
    ax[0].spines['left'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    
    ax[0].set_yticks([])
    ax[1].set_yticks([])
    
    ax[0].set_xlabel( 'Hz', fontsize=14 )
    ax[1].set_xlabel( 'Hz', fontsize=14 )
    
    ax[0].set_title(str(N_particles) +" Random Pts", fontsize=14)
    ax[1].set_title(str(N_particles) +" Hammersley Pts", fontsize=14)
    
    plt.show()