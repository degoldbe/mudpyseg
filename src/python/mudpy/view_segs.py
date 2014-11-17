
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
cdict = {'red': ((0., 1, 1),
                 (0.05, 1, 1),
                 (0.11, 0, 0),
                 (0.66, 1, 1),
                 (0.89, 1, 1),
                 (1, 0.5, 0.5)),
         'green': ((0., 1, 1),
                   (0.05, 1, 1),
                   (0.11, 0, 0),
                   (0.375, 1, 1),
                   (0.64, 1, 1),
                   (0.91, 0, 0),
                   (1, 0, 0)),
         'blue': ((0., 1, 1),
                  (0.05, 1, 1),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65, 0, 0),
                  (1, 0, 0))}

whitejet = matplotlib.colors.LinearSegmentedColormap('whitejet',cdict,256)

def tile_slip(rupt,nstrike,ndip):
    '''
    Quick and dirty plot of a .rupt file
    '''
    
    from numpy import genfromtxt,unique,where,zeros
    import matplotlib.pyplot as plt
    
    f=genfromtxt(rupt)
    num=f[:,0]
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    ss=zeros(len(unum))
    ds=zeros(len(unum))
    for k in range(len(unum)):
        i=where(unum[k]==num)
        ss[k]=all_ss[i].sum()
        ds[k]=all_ds[i].sum()
    #Sum them
    slip=(ss**2+ds**2)**0.5
    #Get unit rake vector
    rakess=ss/slip
    rakeds=ds/slip
    #Get indices for plot
    istrike=zeros(nstrike*ndip)
    idip=zeros(nstrike*ndip)
    k=0
    for i in range(ndip):
         for j in range(nstrike):
             istrike[k]=nstrike-j
             idip[k]=ndip-i
             k+=1           
    #Plot
    plt.figure()
    plt.scatter(istrike,idip,marker='o',c=slip,s=250,cmap=whitejet)
    plt.ylabel('Along-dip index')
    plt.xlabel('Along-strike index')
    cb=plt.colorbar()
    cb.set_label('Slip (m)')
    plt.axis('equal')
    plt.xlim(istrike.min()-1,istrike.max()+1)
    plt.ylim(idip.min()-1,idip.max()+1)
    plt.quiver(istrike,idip,rakess,rakeds,color='green',width=0.002)
    plt.grid()
    plt.title(rupt)
    plt.show()


def tile_slip_segs(rupt,nstrike,ndip):
    '''
    Quick and dirty plot of a .rupt file with multiple fault segments
    '''
    
    from numpy import genfromtxt,unique,where,zeros
    import matplotlib.pyplot as plt
    
    f=genfromtxt(rupt)
    num=f[:,0]
    all_ss=f[:,8]
    all_ds=f[:,9]
    #Now parse for multiple rupture speeds
    unum=unique(num)
    ss=zeros(len(unum))
    ds=zeros(len(unum))
    for k in range(len(unum)):
        i=where(unum[k]==num)
        ss[k]=all_ss[i].sum()
        ds[k]=all_ds[i].sum()
    #Sum them
    slip=(ss**2+ds**2)**0.5
    #Get unit rake vector
    rakess=ss/slip
    rakeds=ds/slip
    #Get indices for plot
    for r in range(len(nstrike)):
        istrike=zeros(nstrike[r]*ndip)
        idip=zeros(nstrike[r]*ndip)
        k=0
        for i in range(ndip):
            for j in range(nstrike[r]):
                istrike[k]=nstrike[r]-j
                idip[k]=ndip-i
                k+=1           
        #Plot
        plt.figure()
        plt.scatter(istrike,idip,marker='o',c=slip,s=250,cmap=whitejet)
        plt.ylabel('Along-dip index')
        plt.xlabel('Along-strike index')
    	cb=plt.colorbar()
        cb.set_label('Slip (m)')
        plt.axis('equal')
        plt.xlim(istrike.min()-1,istrike.max()+1)
        plt.ylim(idip.min()-1,idip.max()+1)
        plt.quiver(istrike,idip,rakess,rakeds,color='green',width=0.002)
        plt.grid()
        plt.title(rupt)
    plt.show()
