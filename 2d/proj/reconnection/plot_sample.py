# For batch execution on a headless machine, run as
# $ MPLBACKEND=agg python plot_sample.py

#load modules
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def test_param(fn):
    #parameters
    print_param = {
        'nx' : 'number of grid in x',
        'ny' : 'number of grid in y',
        'np' : 'maximum number of particles in y',
        'c' : 'speed of light',
        'r' : 'mass',
        'q' : 'charge',
        'wpi' : 'ion plasma frequency',
        'wpe' : 'proper electron plasma frequency',
        'wgi' : 'ion gyro frequency',
        'wge' : 'electron gyro frequency',
        'vti' : 'ion thermal velocity',
        'vte' : 'electron thermal velocity',
        'vai' : 'ion Alfven velocity',
        'vae' : 'electron Alfven velocity',
        'delx' : 'grid size',
        'delt' : 'time step',
        'n0' : 'number of particle / cell',
        'ls' : 'electron skin depth',
    }

    # read all parameters
    with h5py.File(fn, 'r') as f:
        param = dict(f.attrs)

    # some additional parameters
    param['ls']   = param['c']/param['wpe']

    print('*** print parameters ***')
    for key, desc in print_param.items():
        print('- {:40s} : {:}'.format(desc, param[key]))

    return param


def test_moment(fn, it, param, batch=True):
    # read moment data
    with h5py.File(fn, 'r') as dat:
        uf   = dat['uf'][()]
        den  = dat['den'][()]
        vel  = dat['vel'][()]
        temp = dat['temp'][()]
        bx   = uf[...,0]
        by   = uf[...,1]
        bz   = uf[...,2]
        ex   = uf[...,3]
        ey   = uf[...,4]
        ez   = uf[...,5]
        vx   = vel[...,0]
        vy   = vel[...,1]
        vz   = vel[...,2]
        txx  = temp[...,0]
        tyy  = temp[...,1]
        tzz  = temp[...,2]

    # plot moment
    nx    = param['nx']
    ny    = param['ny']
    n0    = param['n0']
    dt    = param['delt']
    dx    = param['delx']
    ls    = param['ls']
    c     = param['c']
    mi    = param['r'][0]
    me    = param['r'][1]
    wpe   = param['wpe']
    wpi   = param['wpi']
    wge   = param['wge']
    wgi   = param['wgi']
    b0    = np.sqrt(4.0*np.pi*n0*mi*param['vti']**2)
    v0    = wge/wpe*c*np.sqrt(me/mi)
    lsize = 14
    tsize = 14
    pad   = 0.1
    ls    = ls*np.sqrt(mi/me)
    xmin  = -0.5*nx*dx/ls
    xmax  = +0.5*nx*dx/ls
    ymin  = 0
    ymax  = ny*dx/ls
    bmax  = np.max(np.abs(bz)/b0)
    bmin  = -bmax
    vmax  = 3.0
    vmin  = -3.0
    nmax  = np.floor(np.max(den/n0))
    nmin  = np.floor(np.min(den/n0))

    plt.figure(figsize=(15,6))
    plt.subplots_adjust(wspace=0.4,hspace=0.3)
    cmap1 = plt.cm.seismic
    cmap2 = plt.cm.hot
    x = np.linspace(xmin,xmax,nx)
    y = np.linspace(ymin,ymax,ny)

    plt.subplot(141)
    plt.title(r'$N_e/N_0$'.format(wgi*it*dt) ,fontsize = lsize,y=1.01)
    plt.streamplot(x,y,bx,by)
    plt.imshow(den[1,:,:]/n0,extent=[xmin,xmax,ymax,ymin],vmin=nmin,vmax=nmax,cmap=cmap2)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    plt.xticks(fontsize=tsize)
    plt.yticks(fontsize=tsize)
    plt.xlabel(r'$x /(c/\omega_{pi})$',fontsize=lsize)
    plt.ylabel(r'$y /(c/\omega_{pi})$',fontsize=lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='2%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(nmin,nmax,6))
    cbar.ax.tick_params(labelsize=tsize)

    plt.subplot(142)
    plt.title(r'$V_{yi}/V_A$',fontsize = lsize,y=1.01)
    plt.streamplot(x,y,bx,by)
    plt.imshow(vy[0,:,:]/v0,extent=[xmin,xmax,ymax,ymin],vmin=vmin,vmax=vmax,cmap=cmap1)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    plt.xticks(fontsize=tsize)
    plt.yticks(fontsize=tsize)
    plt.xlabel(r'$x /(c/\omega_{pi})$',fontsize=lsize)
#    plt.ylabel(r'$y /(c/\omega_{pi})$',fontsize=lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='2%',pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(vmin,vmax,6))
    cbar.ax.tick_params(labelsize=tsize)

    plt.subplot(143)
    plt.title(r'$V_{ye}/V_A$' ,fontsize = lsize,y=1.01)
    plt.streamplot(x,y,bx,by)
    plt.imshow(vy[1,:,:]/v0,extent=[xmin,xmax,ymax,ymin],vmin=vmin,vmax=vmax,cmap=cmap1)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    plt.xticks(fontsize=tsize)
    plt.yticks(fontsize=tsize)
    plt.xlabel(r'$x /(c/\omega_{pi})$',fontsize=lsize)
#    plt.ylabel(r'$y /(c/\omega_{pi})$',fontsize=lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='2%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(vmin,vmax,6))
    cbar.ax.tick_params(labelsize=tsize)

    plt.subplot(144)
    plt.title(r'$B_z/B_0$',fontsize = lsize,y=1.01)
    plt.streamplot(x,y,bx,by)
    plt.imshow(bz/b0,extent=[xmin,xmax,ymax,ymin],vmin=bmin,vmax=bmax,cmap=cmap1)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    plt.xticks(fontsize=tsize)
    plt.yticks(fontsize=tsize)
    plt.xlabel(r'$x /(c/\omega_{pi})$',fontsize=lsize)
#    plt.ylabel(r'$y /(c/\omega_{pi})$',fontsize=lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='2%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(bmin,bmax,6))
    cbar.ax.tick_params(labelsize=tsize)

    if batch:
        plt.savefig('moment.png')


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        datadir  = sys.argv[1] + '/'
    else:
        datadir = './'
    it   = 5000
    fn0  = datadir+'init_param.h5'
    fn1  = datadir+'{:07d}_mom.h5'.format(it)

    param = test_param(fn0)
    test_moment(fn1, it, param)

    if mpl.get_backend() != 'agg':
        plt.show()
