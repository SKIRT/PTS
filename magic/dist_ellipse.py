
def dist_ellipse(n,xc,yc,ratio,pa=0):
  
    ang = np.radians(pa)
    cosang = np.cos(ang)
    sinang = np.sin(ang)
    nx=n[0]
    ny=n[1]
x=np.arange(-xc,nx-xc)
y=np.arange(-yc,ny-yc)
im=np.empty(n)
xcosang = x*cosang
xsinang = x*sinang
for i in range(0,ny):
xtemp = xcosang + y[i]*sinang
ytemp = -xsinang + y[i]*cosang
im[i,:] = np.sqrt( (xtemp*ratio)**2 + ytemp**2 )
return im