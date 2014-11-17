def laplace_stencil(ifault,nstrike,ndip,bounds):
    '''
    Find the index of the subfaults that make the laplacian stencil of fault number ifault
    It assumes all boundaries are initally locked. After assigning stencil values it parses
    the variable 'bounds' and makes corrections if any boundaries were requested to be
    'free'.
    
    Usage:
        stencil=laplace_stencil(ifault,nstrike,ndip)
        
    IN:
        ifault: subfault index number
        nstrike: number of along-strike fault segments
        ndip: number of along dip subfaults
        bounds: A tuple with 4 strings corresponding to the boundary conditions requested by
            the user on the edges fo the fault model. The ordering is top,bototm,left and right edges.
            Possible values for each element of the tuple are 'free' for a free boundary condition
            and 'locked' for a locked one. For example a bounds tuple with 3 locked edges and the top
            edge free would be bounds=('free', 'locked', 'locked', 'locked')
    OUT:
        stencil: indices of the subfaults that contribute to the laplacian
        values: nuemrical values of the stencil
        
    DEG CHANGES:
        line 38: changed to nstrike.any<3 to account for all segments
    '''
    
    from numpy import array
    
    #Get boundary conditions
    top=bounds[0]
    bottom=bounds[1]
    left=bounds[2]
    right=bounds[3]
    #Create stencil
    row=ifault/nstrike #Row number corresponding to this subfault
    column=ifault-(nstrike*row)
    if nstrike.any<3 or ndip<3:
        print "ERROR: The fault model is too small for Laplacian regualrization. You need a minimum of 3 rows and 3 columns in the model."
        return False,False
    if row==0 and column==0: #Top right corner
        stencil=array([ifault,ifault+1,ifault+nstrike])
        values=array([-4,1,1])
        if top.lower()=='free':
            values[2]+=1
        if right.lower=='free':
            values[1]+=1
        return stencil,values
    if row==0 and column==(nstrike-1): #Top left corner
        stencil=array([ifault,ifault-1,ifault+nstrike])
        values=array([-4,1,1])
        if top.lower()=='free':
            values[2]+=1
        if left.lower=='free':
            values[1]+=1
        return stencil,values
    if row==(ndip-1) and column==0: #Bottom right corner
        stencil=array([ifault,ifault+1,ifault-nstrike])
        values=([-4,1,1])
        if bottom.lower()=='free':
            values[2]+=1
        if right.lower=='free':
            values[1]+=1
        return stencil,values
    if row==(ndip-1) and column==(nstrike-1): #Bottom left corner
        stencil=array([ifault,ifault-1,ifault-nstrike])
        values=array([-4,1,1,])
        if bottom.lower()=='free':
            values[2]+=1
        if left.lower=='free':
            values[1]+=1
        return stencil,values
    if row==0: #Top edge, NOT the corner
        stencil=array([ifault,ifault+1,ifault-1,ifault+nstrike])
        values=array([-4,1,1,1])
        if top.lower()=='free':
            values[3]+=1
        return stencil,values
    if row==(ndip-1): #Bottom edge, NOT the corner
        stencil=array([ifault,ifault-1,ifault+1,ifault-nstrike])
        values=array([-4,1,1,1])
        if bottom.lower()=='free':
            values[3]+=1
        return stencil,values
    if column==0: #Right edge, NOT the corner
        stencil=array([ifault,ifault-nstrike,ifault+nstrike,ifault+1])
        values=array([-4,1,1,1])
        if right.lower()=='free':
            values[3]+=1
        return stencil,values
    if column==(nstrike-1): #left edge, NOT the corner
        stencil=array([ifault,ifault-nstrike,ifault+nstrike,ifault-1])
        values=array([-4,1,1,1])
        if left.lower()=='free':
            values[3]+=1
        return stencil,values
    else: #Somewhere in the middle
        stencil=array([ifault,ifault-1,ifault+1,ifault-nstrike,ifault+nstrike])
        values=array([-4,1,1,1,1])
        return stencil,values


def getLs(home,project_name,fault_name,nfaults,num_windows,bounds):
    '''
    Make spatial regularization matrix based on finite difference Lapalce operator.
    This routine will request adjustments depending on the boundary conditions requested
    on the edges of the fault model.
    
    IN:
        home: Home directory
        project_name: Name of the problem
        fault_name: Name of fault description file
        nfaults: Total number of faults in the model
        num_windows: Number of rupture windows
        bounds: A tuple with 4 strings corresponding to the boundary conditions requested by
            the user on the edges fo the fault model. The ordering is top,bototm,left and right edges.
            Possible values for each element of the tuple are 'free' for a free boundary condition
            and 'locked' for a locked one. For example a bounds tuple with 3 locked edges and the top
            edge free would be bounds=('free', 'locked', 'locked', 'locked')

    OUT:
        Lout: The regularization matrix
    '''
    
    from numpy import loadtxt,zeros,tile
    
    #Load source
    source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
    N=len(source) #No. of subfaults
    nstrike=nfaults[0]
    ndip=nfaults[1]
    #Initalize
    L=zeros((2*N,2*N))
    #Which L am I building?
    print 'Making discrete Laplace operator regularization matrix...'
    for kfault in range(N):#Loop over faults and fill regularization matrix
        stencil,values=laplace_stencil(kfault,nstrike,ndip,bounds)
        #Add strike slip branches of stencil
        L[2*kfault,2*stencil]=values
        #Add dip slip branches of stencil
        L[2*kfault+1,2*stencil+1]=values
    if num_windows==1: #Only one rupture speed
        Lout=L 
    else: #Multiple rupture speeds, smooth total moment laplacian
        Lout=L
        Lout=tile(Lout,(1,num_windows))/num_windows
        #for k in range(num_windows-1):
        #    Lout=block_diag(Lout,L)
    return Lout

#def getLs(home,project_name,fault_name,nfaults,num_windows,bounds):
#    '''
#    Make spatial regularization matrix based on finite difference Lapalce operator.
#    This routine will request adjustments depending on the boundary conditions requested
#    on the edges of the fault model.
#    
#    IN:
#        home: Home directory
#        project_name: Name of the problem
#        fault_name: Name of fault description file
#        nfaults: Total number of faults in the model
#        num_windows: Number of rupture windows
#        bounds: A tuple with 4 strings corresponding to the boundary conditions requested by
#            the user on the edges fo the fault model. The ordering is top,bototm,left and right edges.
#            Possible values for each element of the tuple are 'free' for a free boundary condition
#            and 'locked' for a locked one. For example a bounds tuple with 3 locked edges and the top
#            edge free would be bounds=('free', 'locked', 'locked', 'locked')
#
#    OUT:
#        Lout: The regularization matrix
#    
#    DEG CHANGES:
#        line 145: loop through laplace_stencil for each fault segment, output,
#        then rearrarange matrix to correct dimensions using block_diag
#    '''
#    
#    from numpy import loadtxt,zeros,tile,delete
#    from scipy.linalg import block_diag
#    
#    #Load source
#    source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
#    #N=len(source) #No. of subfaults
#    nstrike=nfaults[0]
#    ndip=nfaults[1]
#    #Initalize
#    #L=zeros((2*N,2*N))
#    Pout=[]
#    #Which L am I building?
#    print 'Making discrete Laplace operator regularization matrix...'
#    for i in range(len(nstrike)):
#        N=nstrike[i]*ndip
#        L=zeros((2*N,2*N))
#        for kfault in range(N):#Loop over faults and fill regularization matrix
#            stencil,values=laplace_stencil(kfault,nstrike[i],ndip,bounds)
#            #Add strike slip branches of stencil
#            L[2*kfault,2*stencil]=values
#            #Add dip slip branches of stencil
#            L[2*kfault+1,2*stencil+1]=values
#        if num_windows==1: #Only one rupture speed
#            Lout=L 
#        else: #Multiple rupture speeds, smooth total moment laplacian
#            Lout=L
#            Lout=tile(Lout,(1,num_windows))/num_windows
#            Pout=block_diag(Pout,Lout)
#            #for k in range(num_windows-1):
#            #    Lout=block_diag(Lout,L)
#    Pout=delete(Pout,(0),axis=0)
#    return Pout
    
    
def tile_slip(rupt,nstrike,ndip):
    '''
    Quick and dirty plot of a .rupt file
    DEG Changes:
        line 192: change definition to assign plotting indices to work for fault
        segments
    '''
    
    from numpy import genfromtxt,unique,where,zeros
    import numpy as np
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
    istrike=np.zeros(sum(nstrike)*ndip)
    idip=np.ones(sum(nstrike)*ndip)
    pdip=[]
    pstr=[]
    for r in range(len(nstrike)):
        for i in range(ndip):
            idip=np.ones(nstrike[r])*ndip-i
            pdip=np.append(pdip,idip)
    for r in range(len(nstrike)):      
        istr=np.array(range(sum(nstrike[r:(len(nstrike))])-nstrike[r]+1,sum(nstrike[r:(len(nstrike))])+1))
        for l in range(ndip):
            pstr=np.append(pstr,istr[::-1])                    
    #Plot
    plt.figure()
    plt.scatter(pstr,pdip,marker='o',c=slip,s=250,cmap=whitejet)
    plt.ylabel('Along-dip index')
    plt.xlabel('Along-strike index')
    cb=plt.colorbar()
    cb.set_label('Slip (m)')
    plt.axis('equal')
    plt.xlim(pstr.min()-1,pstr.max()+1)
    plt.ylim(pdip.min()-1,pdip.max()+1)
    plt.quiver(pstr,pdip,rakess,rakeds,color='green',width=0.002)
    plt.grid()
    plt.title(rupt)
    plt.show()
    
    
    #Now make synthetics for source/station pairs
def make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,static,tsunami,beta,
                    hot_start,coord_type,time_epi,nfaults):
    '''
    This routine will take the impulse response (GFs) and pass it into the routine that will
    convovle them with the source time function according to each subfaults strike and dip.
    The result fo this computation is a time series dubbed a "synthetic"
    
    IN:
        home: Home directory
        project_name: Name fo the problem
        station_file: File with coordinates of stations
        fault_name: Name of fault file
        model_Name: Name of Earth structure model file
        integrate: =0 if you want output to be velocity, =1 if you want output to be displacement
        static: =0 if computing full waveforms, =1 if computing only the static field
        hot_start: =k if you want to start computations at k-th subfault, =0 to compute all
        coord_type: =0 if problem is in cartesian coordinates, =1 if problem is in lat/lon
        
    OUT:
        Nothing
        
    DEG Changes:
        line 271: begin loop through subfaults to determine which segment each
        subfault is in, to force use of appropriate beta value
    '''
    from mudpy import green
    import datetime
    from numpy import loadtxt,append
    from string import rjust
    import gc
    
    green_path=home+project_name+'/GFs/'
    station_file=home+project_name+'/data/station_info/'+station_file
    fault_file=home+project_name+'/data/model_info/'+fault_name
    logpath=home+project_name+'/logs/'
    #Time for log file
    now=datetime.datetime.now()
    now=now.strftime('%b-%d-%H%M')
    #First read fault model file
    source=loadtxt(fault_file,ndmin=2)
    #Now compute synthetics please, one sub fault at a time
    #for k in range(hot_start,source.shape[0]):
    nstrike=nfaults[0]
    ndip=nfaults[1]
    strike=[]
    for i in range(ndip):
        strike=append(strike,nstrike)
    for k in range(hot_start,source.shape[0]):
        subfault=rjust(str(k+1),4,'0')
        j=1
        while int(subfault) > sum(strike[0:j]):
            j+=1
        rem=(j-1)%(len(nstrike))
        log=green.run_syn(home,project_name,source[k,:],station_file,green_path,model_name,integrate,static,tsunami,
                subfault,coord_type,time_epi,beta[rem])
        f=open(logpath+'make_synth.'+now+'.log','a')
        f.write(log)
        f.close()
        gc.collect()
        
def rot2ds(sol,beta,nfaults,num_windows):
    '''
    Reverses the operation described in function ds2rot()
    DEG Changes:
        line 306: begin loop to determine which segment a subfault is in to force
        use of appropriate beta
        line 316: in order to use correct beta for individual subfaults, calculate
        ssds one at a time and append, rather than calculating entire matrix at once.
    '''
    from numpy import array,deg2rad,cos,sin,arange,vstack,zeros,append
    
    nstrike=nfaults[0]
    ndip=nfaults[1]
    #Split into strike-slip and dip-slip
    iss=arange(0,len(sol),2)
    ids=arange(1,len(sol),2)
    if len(iss)==1:
        ssrot=sol[0]
        dsrot=sol[1]
    else:
        ssrot=sol[iss]
        dsrot=sol[ids]
    #Rotate
    beta=deg2rad(beta)
    strike=[]
    for i in range(ndip*num_windows):
        strike=append(strike,nstrike)
    ss=[]
    ds=[]
    for k in range(len(ssrot)):
        j=1
        while k>sum(strike[0:j]):
            j+=1
        rem=(j-1)%(len(nstrike))
        ssds_a=array([[cos(beta[rem]),-sin(beta[rem])],[sin(beta[rem]),cos(beta[rem])]]).dot(vstack((ssrot[k].transpose(),dsrot[k].transpose())))
        ss=append(ss,ssds_a[0])
        ds=append(ds,ssds_a[1])
    #Re-insert in output vector
    out=zeros(sol.shape)
    out[iss,0]=ss[:]
    out[ids,0]=ds[:]
    return out
    
    