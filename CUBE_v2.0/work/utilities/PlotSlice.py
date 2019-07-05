import sys
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

nc=512

zs='1.000'
dim='_x'

fig=plt.figure(1,figsize=(8,8))
axes=plt.subplot2grid((100,100),(0,0),colspan=100,rowspan=100)

fn='./'+zs+'nuslice'+dim+'.bin'
data=np.fromfile(fn,dtype='float32')
data=np.reshape(data,(nc,nc))#,order='F')
print('max,min,mean',np.max(data),np.min(data),np.mean(data))
data=np.log10(2+data)

im=axes.imshow(data,cmap='magma',origin='lower left',interpolation='bilinear')

fn='./'+zs+'nuslice'+dim+'.pdf'
print("Saving to file: "+fn)
plt.savefig(fn,bbox_inches='tight')#,rasterized=True)
plt.close(1)
