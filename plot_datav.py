from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

#add location of data vector file for plotting
plotfile = "plots/hires_omnu0_v16smooth_A1_0.70_A2_-1.36_BTA_1.00_eta1_-1.70_eta2_-2.50.png"

# New files:
#cosmosisfile = "cosmosis/Feb2020/Lucas_hires_omnu0_v16smooth_sim_TATT_A1_0.70_A2_-1.36_BTA_1.00_eta1_-1.70_eta2_-2.50.fits"
datavfile1 = "datav_elaheh/xi_nz_v0.16_smoothed_A1_0.70_A2_-1.36_BTA_1.00_eta1_-1.70_eta2_-2.50_Om_nu_0.0"
datavfile2 = "datav_elaheh/xi_nz_v0.16_smoothed_A1_0.70_A2_-1.36_BTA_1.00_eta1_-1.70_eta2_-2.50_Om_nu_0.0"
#datavfile2 = "cosmolike_baseline/xi_final_test_nz_2_26+lens_nz_A1_0.70_A2_-1.36_BTA_1.00_eta1_-1.70_eta2_-2.50"
# v1 N(z)'s files here:
#cosmosisfile = "cosmosis/Feb2020/final_test_nz_2_26_sim_TATT_A1_0.70_A2_-1.36_BTA_1.00_eta1_-1.70_eta2_-2.50.fits"
#datavfile2 = "datav_v16/xi_nz_v0.16_smoothed"

# d1 = np.zeros(900)
# f = fits.open(cosmosisfile)
# #xip
# d1[0:200] = f[2].data['value']
# #xim
# d1[200:400] = f[3].data['value']
# #gammat
# d1[400:800] = f[4].data['value']
# #w(theta)m autobins
# d1[800:] = f[5].data['value']

d1 = np.genfromtxt(datavfile1)[:,1]
d2 = np.genfromtxt(datavfile2)[:,1]

#use this covariance for redmagic lens sample
covfile = "./cov_y3/cov_cosmolike_y3_rm.txt"

ndata = d1.shape[0]
#file with Y1 scale cuts
m = np.genfromtxt("mask_Y1_mcal")[:,1]

ind1 = np.where(m)
ind0 = np.where(m-1.0)
data = np.genfromtxt(covfile)
cov =np.zeros((ndata,ndata))
for i in range(0,data.shape[0]):
	cov[int(data[i,0]),int(data[i,1])] = data[i,2]
	cov[int(data[i,1]),int(data[i,0])] = data[i,2]
	if (int(data[i,0])-int(data[i,1])):
		cov[int(data[i,0]),int(data[i,1])]*= m[int(data[i,0])]*m[int(data[i,1])]  	
		cov[int(data[i,1]),int(data[i,0])]*= m[int(data[i,0])]*m[int(data[i,1])]  	
#s now contains sqrt(cov[i,i])
s = np.sqrt(np.diag(cov))
ind = np.arange(0,ndata)
print(d2/d1*m)
print("chi2 calculated using Y1 scale cuts")
chi =0.0
inv = LA.inv(cov)
for i in range(0,ndata):
	for j in range(0,ndata):
		chi +=(d1[i]-d2[i])*inv[i,j]*(d1[j]-d2[j])*m[i]*m[j]
print("3x2pt: Delta chi2 = %f" %(chi))

chi =0.0
inv = LA.inv(cov[400:,400:])
for i in range(0,500):
	for j in range(0,500):
		chi +=(d1[400+i]-d2[400+i])*inv[i,j]*(d1[400+j]-d2[400+j])*m[400+i]*m[400+j]
print("2x2pt: Delta chi2 = %f" %(chi))

chi =0.0
inv = LA.inv(cov[800:,800:])
for i in range(0,100):
	for j in range(0,100):
		chi +=(d1[800+i]-d2[800+i])*inv[i,j]*(d1[800+j]-d2[800+j])*m[800+i]*m[800+j]
print("wtheta: Delta chi2 = %f" %(chi))

nxip = 200
nxim = 200
nw = 100
nggl = ndata - nw - nxip -nxim
s = np.sqrt(np.diag(cov))

plt.figure(figsize=(8,8), dpi=400)
fs = 18
plt.subplot(4,2,1)
plt.yscale('log')
plt.ylim(2.e-7,1.2e-4)
plt.xlim(0,nxip-1)
#plt.title(r'$\xi_+$')
plt.ylabel(r'$\xi_+$', fontsize = fs)
plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.25)
plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)


plt.subplot(4,2,2)
plt.ylim(-0.1,0.1)
plt.plot([0,1000],[0,0],linestyle ='--',color='k')
plt.xlim(0,nxip-1)
plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
plt.errorbar(ind,d1*0,s/d1,marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)

plt.subplot(4,2,3)
plt.yscale('log')
plt.ylim(2.e-7,6.e-5)
plt.xlim(nxip,nxip+nxim-1)
#plt.title(r'$\xi_-$')
plt.ylabel(r'$\xi_-$', fontsize = fs)
plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.25)
plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)

plt.subplot(4,2,4)
plt.ylim(-0.1,0.1)
plt.plot([0,1000],[0,0],linestyle ='--',color='k')
plt.xlim(nxip,nxip+nxim-1)
plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
plt.errorbar(ind,d1*0,s/d1,marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)

plt.subplot(4,2,5)

plt.yscale('log')
plt.ylim(2.e-6,2.5e-3)
plt.xlim(nxip+nxim,nxip+nxim+nggl-1)
#plt.title(r'$\gamma_t$')
plt.ylabel(r'$\gamma_t$', fontsize = fs)
plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.2)
plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)
#plt.plot(ind,d3,linestyle = '-')

plt.subplot(4,2,6)
plt.ylim(-0.1,0.1)
plt.plot([0,1000],[0,0],linestyle ='--',color='k')
plt.xlim(nxip+nxim,nxip+nxim+nggl-1)
#plt.title(r'$\gamma_t$')
plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
plt.errorbar(ind,d1*0,s/d1,marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)


plt.subplot(4,2,7)
plt.yscale('log')
plt.ylim(1.e-4,0.6)
plt.xlim(nxip+nxim+nggl,ndata)
#plt.title(r'$w$')
plt.ylabel(r'$w$', fontsize = fs)
plt.xlabel(r'bin number', fontsize = fs)
plt.errorbar(ind,d1,s,marker='o', color='k',linestyle = '',markersize = 0.5,alpha = 0.4)
plt.plot(ind,d1,marker='o', color='r',linestyle = '',markersize = 1.5)
#plt.plot(ind,d3,linestyle = '-')


plt.subplot(4,2,8)
plt.ylim(-0.1,0.1)
plt.plot([0,1000],[0,0],linestyle ='--',color='k')
plt.xlim(nxip+nxim+nggl,ndata)
#plt.title(r'$w$')
plt.xlabel(r'bin number', fontsize = 18)
plt.ylabel(r'(d2-d1)/d2', fontsize = fs)
plt.errorbar(ind,d1*0,s/d1,marker='o', color='k',linestyle = '',markersize = 0.0,alpha = 0.1)
plt.plot(ind[ind0],(d2[ind0]-d1[ind0])/d2[ind0],marker='x', color='k',linestyle = '',markersize = 1.0)
plt.plot(ind[ind1],(d2[ind1]-d1[ind1])/d2[ind1],marker='o', color='r',linestyle = '',markersize = 1.0)

plt.tight_layout()
plt.savefig(plotfile,dpi=400)

