# Import librairies
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox

# Plot RGB images
def plot_rgb(gal, bands=[5,6,7], ax=None, band_first=True, zoom=1.5, shifts=None, clip = False):
    '''
    Return plot in rgb of image
    '''
    if ax is None:
        ax = plt.subplot()
    if band_first:
        tr = [1,2,0]
    else:
        tr = [0,1,2]
    
    imsize = float(gal.shape[1])/2
    if clip==True:
        ax.imshow(np.clip(gal[:,:,:].transpose(tr)[:,:,bands], a_min=0.0, a_max=1.), extent=(-imsize,imsize,-imsize,imsize), origin='lower left')
    else:
        # Here, in the extent kwargs, the bottom is put at the top of image and vice versa, to mimic the basic behavior of matplotlib if imshow in only one band
        ax.imshow(gal[:,:,:].transpose(tr)[:,:,bands], extent=(-imsize,imsize,imsize,-imsize), origin='lower left')
    if shifts is not None:
        for (x,y) in shifts:
            if x !=0 and y!=0: # No scatter if there is no shift (i.e. no galaxy)
                # As a consequence from previous comment, need to put -y here.
                ax.scatter(x, -y,  marker='+', c='r')
    ax.set_xlim(-imsize/zoom,imsize/zoom)
    ax.set_ylim(-imsize/zoom,imsize/zoom)
    ax.axis('off')



def scatter_galaxies(image, shift, pixel_scale, stamp_size, scatter = 'numbers', blendedness = None, ax=None):
    """
    Return plot on single band and scatter number on each galaxies

    Parameters:
    ----------
    image: single band image
    shift: list of shift (output from blended images generation function)
    pixel_scale: pixel scale of the bandpass filter used to plot the image
    stamp_size: size of the stamp
    """
    ax.imshow(image)

    if scatter == 'numbers':
        for k in range (len(shift)):
            ax.scatter((stamp_size) + shift[k][0]/pixel_scale, (stamp_size) + shift[k][1]/pixel_scale, s = 50 ,c='red', marker="${}$".format(k))
    elif scatter == 'blendedness':
        for k in range (len(blendedness)):
            ax.scatter((stamp_size) + shift[k][0]/pixel_scale, (stamp_size) + shift[k][1]/pixel_scale, s = 500 ,c='red', marker="${0:.2f}$".format(blendedness[k]))



def mean_var(x,y,bins):
    """
    Return mean and variance in each bins of the histogram
    """
    n,_ = np.histogram(x,bins=bins, weights=None)
    ny,_ = np.histogram(x,bins=bins, weights=y)
    mean_y = ny/n
    ny2,_ = np.histogram(x,bins=bins, weights=y**2)
    var_y = (ny2/n - mean_y**2)/n
    
    return (mean_y, var_y)


def createCircularMask(h, w, center=None, radius=None):
    '''
    Return circular mask around the center of the galaxy
    '''
    if center is None: # use the middle of the image
        center = [int(w/2), int(h/2)]
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask




# To plot corner plot of latente space
def plot_corner_latent(z, lim=3, nbins=25, show_title=True):
    """
    Make a corner plot of standard gaussian distributed latent variables.
    Parameters
    ----------
    z : latent variables, array of size (n_samples, latent_dim)
    lim : int, optional
        [description], by default 3
    nbins : int, optional
        [description], by default 25
    show_title : bool, optional
        [description], by default True
    Example
    -------
    z = np.random.normal(size=(1000,8))
    plot_corner_latent(z)
    """

    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import Grid
    
    latent_dim = z.shape[1]
    
    fig = plt.figure(figsize=(latent_dim*2,latent_dim*2))
    grid = Grid(fig, rect=111, nrows_ncols=(latent_dim,latent_dim), axes_pad=0.25, label_mode='L', share_y=False)
    
    colors = mpl.cm.jet(np.linspace(0,1,latent_dim))
    bins = nbins
    
    for i in range(latent_dim):
        for j in range(latent_dim):
            ax = grid[i*latent_dim+j]
            if i == j :
                n,_,_ = ax.hist(z[:,i], bins=bins, normed=True, color=colors[i])
                ax.set_yticks([])
                if show_title:
                    ax.set_title('$z_{}$'.format(i))
            if i > j :
                ax.hist2d(z[:,j], z[:,i], bins=bins, cmap=mpl.cm.gray)
            if i < j :
                ax.axis('off')
    
    plt.tight_layout()
    plt.xlim(-lim,lim)