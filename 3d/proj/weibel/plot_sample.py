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
        den0 = dat['den'][()]
        vel  = dat['vel'][()]
        temp = dat['temp'][()]
        bx   = uf[0,:,:,0]
        by   = uf[0,:,:,1]
        bz   = uf[0,:,:,2]
        ex   = uf[0,:,:,3]
        ey   = uf[0,:,:,4]
        ez   = uf[0,:,:,5]
        den  = den0[:,0,:,:]
        vx   = vel[:,0,:,:,0]
        vy   = vel[:,0,:,:,1]
        vz   = vel[:,0,:,:,2]
        txx  = temp[:,0,:,:,0]
        tyy  = temp[:,0,:,:,1]
        tzz  = temp[:,0,:,:,2]

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
    lsize = 16
    tsize = 16
    pad   = 0.1
    xmin  = 0
    xmax  = nx*dx/ls
    ymin  = 0
    ymax  = ny*dx/ls
    bmax  = np.floor(np.max([np.max(np.abs(bx)/b0),np.max(np.abs(by)/b0)]))
    bmin  = -bmax
    nmax  = np.floor(np.max(den/n0))
    nmin  = np.floor(np.min(den/n0))

    plt.figure(figsize=(12,12))
    plt.subplots_adjust(wspace=0.4,hspace=0.3)
    cmap1 = plt.cm.seismic
    cmap2 = plt.cm.hot

    plt.subplot(221)
    plt.title(r'$\omega _{pe}t = %5d$' % (wpe*it*dt) ,fontsize = lsize,y=1.05)
    plt.imshow(den[1,:,:]/n0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=nmin,vmax=nmax,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    plt.xticks(fontsize=tsize)
    plt.yticks(fontsize=tsize)
    plt.xlabel(r'$x /(c/\omega_{pe})$',fontsize=lsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='2%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(nmin,nmax,6))
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$N_e/N_0$',fontsize=lsize)

    plt.subplot(222)
    plt.title(r'$\omega _{pe}t = %5d$' % (wpe*it*dt) ,fontsize = lsize,y=1.05)
    plt.imshow(den[0,:,:]/n0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=nmin,vmax=nmax,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    plt.xticks(fontsize=tsize)
    plt.yticks(fontsize=tsize)
    plt.xlabel(r'$x /(c/\omega_{pe})$',fontsize=lsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize=lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='2%',pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(nmin,nmax,6))
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$N_i/N_0$',fontsize=lsize)

    plt.subplot(223)
    plt.title(r'$\omega _{pe}t = %5d$' % (wpe*it*dt) ,fontsize = lsize,y=1.05)
    plt.imshow(bx/b0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=bmin,vmax=bmax,cmap=cmap1,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    plt.xticks(fontsize=tsize)
    plt.yticks(fontsize=tsize)
    plt.xlabel(r'$x /(c/\omega_{pe})$',fontsize=lsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='2%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(bmin,bmax,6))
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$B_x/B_0$',fontsize=lsize)

    plt.subplot(224)
    plt.title(r'$\omega _{pe}t = %5d$' % (wpe*it*dt) ,fontsize = lsize,y=1.05)
    plt.imshow(by/b0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=bmin,vmax=bmax,cmap=cmap1,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    plt.xticks(fontsize=tsize)
    plt.yticks(fontsize=tsize)
    plt.xlabel(r'$x /(c/\omega_{pe})$',fontsize=lsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='2%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(bmin,bmax,6))
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$B_y/B_0$',fontsize=lsize)

    if batch:
        plt.savefig('moment.png')


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        datadir  = sys.argv[1] + '/'
    else:
        datadir = './'
    it   = 300
    fn0  = datadir+'init_param.h5'
    fn1  = datadir+'{:07d}_mom.h5'.format(it)

    param = test_param(fn0)
    test_moment(fn1, it, param)

    if mpl.get_backend() != 'agg':
        plt.show()
