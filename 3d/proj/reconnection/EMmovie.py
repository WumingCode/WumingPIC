import numpy as np
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import ArtistAnimation
plt.rcParams["font.size"]=14
plt.rcParams["lines.linewidth"]=2.5

datadir = './'
iz = 10

#----- Parameter reading -----
def test_param(fn):
    #parameters
    print_param = {
        'nx' : 'number of grid in x',
        'ny' : 'number of grid in y',
        'nz' : 'number of grid in z',
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

fn = datadir+'init_param.h5'
param = test_param(fn)

nx    = param['nx']
ny    = param['ny']
nz    = param['nz']
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
#lsize = 14
#tsize = 14
#pad   = 0.1
ls    = ls*np.sqrt(mi/me)
xmin  = -0.5*nx*dx/ls
xmax  = +0.5*nx*dx/ls
ymin  = 0
ymax  = ny*dx/ls
zmin  = 0
zmax  = nz*dx/ls
#bmax  = np.max(np.abs(bz)/b0)
#bmin  = -bmax
#vmax  = 3.0
#vmin  = -3.0
#nmax  = np.floor(np.max(den/n0))
#nmin  = np.floor(np.min(den/n0))

cm1 = plt.cm.seismic
cm2 = plt.cm.hot
#plt.subplots_adjust(wspace=0.4,hspace=0.3)
x = np.linspace(xmin,xmax,nx)
y = np.linspace(ymin,ymax,ny)

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(8, 8), sharex=True, sharey=True)
axs[0, 0].set_ylabel(r'$y /(c/\omega_{pi})$')
axs[0, 0].set_title('$B_x$')
axs[0, 1].set_title('$B_y$')
axs[0, 2].set_title('$B_z$')
axs[1, 0].set_xlabel(r'$x /(c/\omega_{pi})$')
axs[1, 0].set_ylabel(r'$y /(c/\omega_{pi})$')
axs[1, 0].set_title('$E_x$')
axs[1, 1].set_xlabel(r'$x /(c/\omega_{pi})$')
axs[1, 1].set_title('$E_y$')
axs[1, 2].set_xlabel(r'$x /(c/\omega_{pi})$')
axs[1, 2].set_title('$E_z$')

#----- Making movie -----
frames = []
nt = 20
dt_mom = 500
for it in range(nt):
    t = (it+1)*dt_mom
    file = datadir+'{:07d}_mom.h5'.format(t)
    with h5py.File(file, 'r') as dat:
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

    frame = []
    im0 = axs[0, 0].imshow(bx[iz, :, :]/b0, extent=[xmin,xmax,ymax,ymin], cmap=cm1, animated=True)
    im1 = axs[0, 1].imshow(by[iz, :, :]/b0, extent=[xmin,xmax,ymax,ymin], cmap=cm1, animated=True)
    im2 = axs[0, 2].imshow(bz[iz, :, :]/b0, extent=[xmin,xmax,ymax,ymin], cmap=cm1, animated=True)
    im3 = axs[1, 0].imshow(ex[iz, :, :]/b0, extent=[xmin,xmax,ymax,ymin], cmap=cm1, animated=True)
    im4 = axs[1, 1].imshow(ey[iz, :, :]/b0, extent=[xmin,xmax,ymax,ymin], cmap=cm1, animated=True)
    im5 = axs[1, 2].imshow(ez[iz, :, :]/b0, extent=[xmin,xmax,ymax,ymin], cmap=cm1, animated=True)
    frame.append(im0)
    frame.append(im1)
    frame.append(im2)
    frame.append(im3)
    frame.append(im4)
    frame.append(im5)
    title = axs[0, 1].text(0., -2., f'it={t}', ha='center')
    frame.append(title)
    frames.append(frame)
    
anim  = ArtistAnimation(fig, frames, interval=500)
anim.save("./EMmovie.mp4")