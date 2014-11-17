def getLs_segs(home,project_name,fault_name,nfaults,num_windows,bounds,segs):
    '''
    Make spatial regularization matrix based on finite difference Lapalce operator.
    This routine will request adjustments depending on the boundary conditions requested
    on the edges of the fault model.
    
    IN:
        home: Home directory
        project_name: Name of the problem
        fault_name: Name of fault description file
        ndip: Number of subfaults down dip
        nstrike_segs: numpy.array of nstrikes for each fault segment
        num_windows: Number of rupture windows
        bounds: A tuple with 4 strings corresponding to the boundary conditions requested by
            the user on the edges fo the fault model. The ordering is top,bototm,left and right edges.
            Possible values for each element of the tuple are 'free' for a free boundary condition
            and 'locked' for a locked one. For example a bounds tuple with 3 locked edges and the top
            edge free would be bounds=('free', 'locked', 'locked', 'locked')

    OUT:
        Pout: The regularization matrix
    '''
    
    from numpy import loadtxt,zeros,tile,delete
    from mudpy.inverse import laplace_stencil
    from scipy.linalg import block_diag
    
    #Load source
    source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
    M=len(source) #No. of subfaults
    nstrike_segs=nfaults[0]
    ndip=nfaults[1]
    #Initalize
    P=zeros((2*M,2*M))
    Pout=[]
    #Which L am I building?
    print 'Making discrete Laplace operator regularization matrix...'
    for i in range(len(nstrike_segs)):
        print 'Working on segment of length ' + nstrike_segs[i].astype(str)
        N = nstrike_segs[i]*ndip
        L=zeros((2*N,2*N))
        for kfault in range(N):#Loop over faults and fill regularization matrix
            stencil,values=laplace_stencil(kfault,nstrike_segs[i],ndip,bounds,segs)
            #Add strike slip branches of stencil
            L[2*kfault,2*stencil]=values
            #Add dip slip branches of stencil
            L[2*kfault+1,2*stencil+1]=values
        print 'Creating L'
        if num_windows==1: #Only one rupture speed
            Lout=L 
        else: #Multiple rupture speeds, smooth total moment laplacian
            Lout=L
            Lout=tile(Lout,(1,num_windows))/num_windows
            #for k in range(num_windows-1):
            #    Lout=block_diag(Lout,L)
            Pout=block_diag(Pout,Lout)
    Pout=delete(Pout,(0),axis=0)
    return Pout
    
  
# I'm pretty sure this following function is not necessary, let's just save it here for a moment in case      
def get_ABIC(G,GTG,sol,d,lambda_s,lambda_t,Ls,LsLs,Lt,LtLt,nfaults,segs):
    '''
    Compute Akaike's Bayesian information criterion, for details see Ide et al. (1996)
    in BSSA, specifically equation 33.
    
    IN:
        G: GFs matrix
        sol Solution vector from inversion
        d: Data vector
        lambda_s: Spatial regularization parameter
        lambda_t: Temporal regularization parameter
        Ls: Spatial regularization matrix
        Lt:Temporal reularization matrix
        Ls_rank: Rank of Ls (#eigenvalues>0)
        Lt_rank: Rank of Lt
    OUT:
        ABIC: Akaike's Bayesian information criterion
    
    '''
    
    from numpy import log
    from numpy.linalg  import norm,slogdet
    
    #Data points
    N=d.size
    #Model parameters
    M=sol.size
    #nfaults, segments
    nstrike_segs = nfaults[0]
    #Off you go, compute it
    if lambda_t==0: #There is only one contraint (no temporal regularization)
        if segs==0:
            s=norm(d-G.dot(sol))**2+(lambda_s**2)*norm(Ls.dot(sol))**2
            a1=N*log(s)
            a2=M*log(lambda_s**2)
            sq,a3=slogdet(GTG+(lambda_s**2)*LsLs)
            #Add 'em up
            ABIC=a1-a2+a3
        elif segs==1:
            for i in range(len(nstrike_segs)):
                s0=norm(d-G.dot(sol))**2+(lambda_s**2)*norm(Ls[i].dot(sol))**2
                a10=N*log(s0)
                a20=M*log(lambda_s**2)
                sq0,a30=slogdet(GTG+(lambda_s**2)*LsLs[i])
                ABIC0=a10-a20+a30
                s.append(s0)
                a1.append(a10)
                a2.append(a20)
                sq.append(sq0)
                a3.append(a30)
            #Add 'em up
            ABIC.append(ABIC0)
        return ABIC
    else: #There is a double regularization, use Fukahata et al. definition
        print '... computing 2d-ABIC'
        if segs==0:
            s=(norm(d-G.dot(sol))**2)+((lambda_s**2)*(norm(Ls.dot(sol))**2))+((lambda_t**2)*(norm(Lt.dot(sol))**2))
            a1=N*log(s)
            sq,a2=slogdet((lambda_s**2)*LsLs+(lambda_t**2)*LtLt)
            sq,a3=slogdet(GTG+(lambda_s**2)*LsLs+(lambda_t**2)*LtLt)
            #Add 'em up
            ABIC=a1-a2+a3
        elif segs==1:
            for i in range(len(nstrike_segs)):
                s0=(norm(d-G.dot(sol))**2)+((lambda_s**2)*(norm(Ls[i].dot(sol))**2))+((lambda_t**2)*(norm(Lt.dot(sol))**2))
                a10=N*log(s0)
                sq0,a20=slogdet((lambda_s**2)*LsLs[i]+(lambda_t**2)*LtLt)
                sq0,a30=slogdet(GTG+(lambda_s**2)*LsLs[i]+(lambda_t**2)*LtLt)
                ABIC0=a10-a20+a30
                s.append(s0)
                a1.append(a10)
                a2.append(a20)
                sq.append(sq0)
                a3.append(a30)
            #Add 'em up
            ABIC.append(ABIC0)
        return ABIC

        