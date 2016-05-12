import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from matplotlib import rc
rc('text', usetex=True)
from scipy.stats import nanmean


def main():

    path = "./"
    totISRFfile    = "total212Labs.dat"
    Lsun = 3.846e26 # Watts
    # Total energy absorbed in the new component
    # Derived from ../../SKIRTrun/models/testHeating/MappingsHeating/plotSEDs.py
    Lnew = 176495776.676 # in Lsun

    print 'reading data... '
    input = np.loadtxt(path+totISRFfile,skiprows=1)
    ID       = input[:,0]
    volume   = input[:,1]
    density  = input[:,2]
    massFrac = input[:,3]
    density_new = input[:,5]
    x   = input[:,6]
    y   = input[:,7]
    tot = input[:,9] / Lsun
    old = input[:,10] / Lsun
    yng = input[:,11] / Lsun
    new = input[:,12] / Lsun
    
    energy_new = volume*density_new*Lnew
    r = np.array([np.sqrt(x[i]**2 + y[i]**2) for i in range(len(ID))])
    F_abs_yng = (yng+new+energy_new)/(old+yng+new+energy_new)

    print 'Diffusely Absorbed energy from all stars:                 ', np.sum(tot)
    print 'Diffusely Absorbed energy from old stellar populations:   ', np.sum(old)
    print 'Diffusely Absorbed energy from young stellar populations: ', np.sum(yng)
    print 'Diffusely Absorbed energy from new stellar populations:   ', np.sum(new)
    print 'Internal  Absorbed energy from new stellar populations:   ', np.sum(energy_new)

    print ID[10000], tot[10000], old[10000], yng[10000], new[10000], energy_new[10000]
    
    print 'plotting... '
    fig  = plt.figure(figsize=(15,5))
    locplot=[[0.05,0.12,0.28,0.85],[0.38,0.12,0.28,0.85],[0.71,0.12,0.28,0.85]]

    fig_a = plt.axes(locplot[0])
    fig_a.set_ylabel(r'$L^\mathrm{abs}_\mathrm{tot} \, [L_\odot]$',fontsize=18)
    fig_a.set_xlabel(r'$L^\mathrm{abs}_\mathrm{evolved}+L^\mathrm{abs}_\mathrm{young}+L^\mathrm{abs}_\mathrm{new} \, [L_\odot]$',fontsize=18)
    fig_a.scatter(old+yng+new,tot,color='green',marker='.')
    fig_a.plot([1e0,1e8],[1e0,1e8],'r-')
    #fig_a.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
    #fig_a.set_xlim(4.9,8.9)
    #fig_a.set_ylim(-10,90)
    fig_a.set_xscale('log')
    fig_a.set_yscale('log')
    fig_a.set_xlim(1e0,1e6)
    fig_a.set_ylim(1e0,1e6)

    fig_b = plt.axes(locplot[1])
    fig_b.set_ylabel('Mass weighted frequency',fontsize=18)
    fig_b.set_xlabel(r'$\mathcal{F}_\mathrm{unev.}^\mathrm{abs}$',fontsize=18)
    fig_b.hist(F_abs_yng,bins=50, normed=True, weights=massFrac, label = 'young', color='green' )
    med = np.median(F_abs_yng)
    print med
    fig_b.plot([med,med],[0,8],'r-')
    #fig_b.hist(old/(old+yng),bins=50, normed=True, weights=massFrac, label = 'old' )
    fig_b.set_xlim(0.0,0.65)
    #fig_b.legend()

    rBins_F, FBins_r = getRadBins(r/1000.,F_abs_yng,1 , massFrac)
    rBins_F[rBins_F > 25] = np.nan
    

    # Estimate the 2D histogram
    nbins = 200
    H, xedges, yedges = np.histogram2d(r/1000.,F_abs_yng,bins=nbins, normed=True, weights=massFrac)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
    
    # Plot 2D histogram using pcolor
    
    fig_c = plt.axes(locplot[2])
    fig_c.set_ylabel('$\mathcal{F}_\mathrm{unev.}^\mathrm{abs}$',fontsize=18)
    fig_c.set_xlabel('R (kpc)',fontsize=18)
    #fig_c.hexbin(r/1000.,F_abs_yng,gridsize=150,bins='log',cmap=plt.cm.autumn, mincnt=1,linewidths=0)
    fig_c.pcolormesh(xedges, yedges, Hmasked)
    fig_c.plot(rBins_F,FBins_r, 'k-',linewidth=2)
    fig_c.plot(rBins_F,FBins_r, 'w-',linewidth=1)

    fig_c.errorbar(1.7,0.88,xerr=1.4, color='k')
    fig_c.text(1.8,0.90,'Bulge', ha='center')
    fig_c.errorbar(11.,0.88,xerr=2.75, color='k')
    fig_c.text(11.,0.90,'main SF ring', ha='center')
    fig_c.errorbar(16.,0.88,xerr=1, color='k')
    fig_c.text(15.,0.90,r'$2^\mathrm{nd}$ SF ring', ha='left')
    
    fig_c.set_ylim(0.0,1.0)
    fig.savefig(path+"plot212Labs.png",format='png', dpi=600)

# Plot figures separetly
#
#    fig  = plt.figure(figsize=(5,5))
#    fig_a = plt.axes([0.12,0.11,0.87,0.87])
#    fig_a.set_ylabel('LabsTot',fontsize=18)
#    fig_a.set_xlabel(r'LabsOld$+$LabsYoung',fontsize=18)
#    fig_a.plot(old+yng,tot,'b.')
#    fig_a.plot([10,1e6],[10,1e6],'r-')
#    #fig_a.scatter(x,y, c=z, s=10,cmap=plt.get_cmap('autumn') ,edgecolor='')
#    #fig_a.set_xlim(4.9,8.9)
#    #fig_a.set_ylim(-10,90)
#    fig_a.set_xscale('log')
#    fig_a.set_yscale('log')
#    
#    fig.savefig(path+"plot_ISRF.png",format='png')
#    
#    fig  = plt.figure(figsize=(5,5))
#    fig_a = plt.axes([0.11,0.12,0.86,0.86])
#    fig_a.set_ylabel('Mass weighted frequency',fontsize=18)
#    fig_a.set_xlabel(r'$F_\mathrm{young}^\mathrm{abs}$',fontsize=18)
#    fig_a.hist(F_abs_yng,bins=50, normed=True, weights=massFrac, label = 'young', color='green' )
#    med = np.median(F_abs_yng)
#    fig_a.plot([med,med],[0,7],'r-')
#    #fig_a.hist(old/(old+yng),bins=50, normed=True, weights=massFrac, label = 'old' )
#    fig_a.set_xlim(0.0,0.65)
#    #fig_a.legend()
#    fig.savefig(path+"plothist_Fyoung.png",format='png')
#    
#    
#    fig  = plt.figure(figsize=(5,5))
#    fig_a = plt.axes([0.14,0.12,0.85,0.85])
#    fig_a.set_ylabel('$F_\mathrm{young}^\mathrm{abs}$',fontsize=18)
#    fig_a.set_xlabel('R (kpc)',fontsize=18)
#    #fig_a.plot(r, yng/(old+yng) , c='b',markersize=1,alpha=0.1)
#    fig_a.hexbin(r/1000.,F_abs_yng,gridsize=150,bins='log',cmap=plt.cm.autumn, mincnt=1,linewidths=0)
#    fig_a.set_ylim(0.0,0.65)
#    fig.savefig(path+"plotrad_Fyoung.pdf",format='pdf')

def getRadBins(xarr, yarr, binstep, weights):
    
    idx = np.argsort(xarr)
    sortx = xarr[idx]
    sorty = yarr[idx]
    sortweights = weights[idx]
    
    avx = np.array([])
    avy = np.array([])
    
    i = 0
    counter = 0
    while i*binstep < sortx[-1]:
        n=0
        sumy = 0
        sumweight = 0
        while sortx[counter] < (i+1)*binstep and counter < len(sortx)-1:
            sumy += sorty[counter]*sortweights[counter]
            sumweight += sortweights[counter]
            n += 1
            counter += 1
        avx = np.append(avx,i*binstep+0.5)
        avy = np.append(avy,sumy/sumweight)
        i += 1
    return avx, avy




if __name__ == '__main__':
    main()