from numpy import savetxt,genfromtxt,zeros,array
from mudpy import inverse as inv

home='/Users/degoldbe/Documents/SOPAC/slip_inversion/'
project_name='El_Mayor'
fault_name='elmayor_3124_merge.fault'
nstrike=array([11,6,9,3]) ; ndip=4 ; nfaults=(nstrike,ndip)
#nstrike=array([29]) ; ndip=4 ; nfaults=(nstrike,ndip)
num_windows=1
top='free' ; bottom='locked' ; left='locked' ; right='locked' #'locked' or 'free'
bounds=(top,bottom,left,right)

ruptfile=genfromtxt(u'/Users/degoldbe/Documents/SOPAC/slip_inversion/El_Mayor_Files/segsmoothing_test.rupt',skip_header=1,dtype='f8')
#lon=genfromtxt(u'/Users/degoldbe/Documents/SOPAC/slip_inversion/El_Mayor/data/model_info/elmayor_3124_segtest.fault',usecols=[2],dtype='f8')
#lat=genfromtxt(u'/Users/degoldbe/Documents/SOPAC/slip_inversion/El_Mayor/data/model_info/elmayor_3124_segtest.fault',usecols=[3],dtype='f8')
ss=genfromtxt(u'/Users/degoldbe/Documents/SOPAC/slip_inversion/El_Mayor_Files/segsmoothing_test.rupt',usecols=[8],skip_header=1,dtype='f8')
ds=genfromtxt(u'/Users/degoldbe/Documents/SOPAC/slip_inversion/El_Mayor_Files/segsmoothing_test.rupt',usecols=[9],skip_header=1,dtype='f8')
m=zeros(2*len(ss))
i1=range(0,2*len(ss),2)
i2=range(1,2*len(ss),2)
m[i1]=ss
m[i2]=ds
Ls=inv.getLs(home,project_name,fault_name,nfaults,num_windows,bounds)
mL=Ls.dot(m)
j1=range(0,len(m),2)
j2=range(1,len(m),2)
ssout=mL[j1]
dsout=mL[j2]
ruptout=ruptfile[0:116,:]
ruptout[:,8]=ssout
ruptout[:,9]=dsout
#ruptout[:,2]=lon
#ruptout[:,3]=lat
savetxt('/Users/degoldbe/Documents/SOPAC/slip_inversion/El_Mayor_Files/mL.rupt',ruptout)