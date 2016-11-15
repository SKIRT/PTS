#import aplpy
import pyparsing
#import pyregion
import pyfits
import matplotlib.pyplot as pyplot
import matplotlib as mpl

path = 'FinalRun/'

fig = pyplot.figure(figsize=(9,10))
fig.text(0.385,0.97,"Offset from centre (degrees)",color='black',size='16',weight='bold')
fig.text(0.02,0.615,"Offset from centre (degrees)",color='black',size='16',weight='bold',rotation='vertical')

def standard_setup(sp):
  sp.set_frame_color('black')
  sp.set_tick_labels_font(size='10')
  sp.set_axis_labels_font(size='12')
    #sp.set_tick_labels_format(xformat='hh:mm',yformat='dd:mm')
  sp.set_xaxis_coord_type('scalar')
  sp.set_yaxis_coord_type('scalar')
  sp.set_tick_color('black')
  sp.recenter(x=0.0, y=0.0,width=3.,height=0.6)
  sp.set_tick_xspacing(0.4)
  sp.set_tick_yspacing(0.25)
  sp.set_system_latex(True)
  sp.tick_labels.hide()
  sp.axis_labels.hide()


plotloc = [[0.09,0.825,0.3,0.115],[0.39,0.825,0.3,0.115],[0.69,0.825,0.3,0.115],
           [0.09,0.710,0.3,0.115],[0.39,0.710,0.3,0.115],[0.69,0.710,0.3,0.115],
           [0.09,0.595,0.3,0.115],[0.39,0.595,0.3,0.115],[0.69,0.595,0.3,0.115],
           [0.09,0.480,0.3,0.115],[0.39,0.480,0.3,0.115],[0.69,0.480,0.3,0.115],
           [0.09,0.365,0.3,0.115],[0.39,0.365,0.3,0.115],[0.69,0.365,0.3,0.115],
           [0.09,0.250,0.3,0.115],[0.39,0.250,0.3,0.115],[0.69,0.250,0.3,0.115],
           [0.09,0.135,0.3,0.115],[0.39,0.135,0.3,0.115],[0.69,0.135,0.3,0.115],
           [0.09,0.020,0.3,0.115]]
# First row

f1 = aplpy.FITSFigure(path+'maps/plotimFUVJy.fits', figure=fig, subplot=plotloc[0])
standard_setup(f1)
#f1.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f1.show_colorscale(vmin=0, vmax=0.00025, cmap='hot')
f1.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
#f1.axis_labels.show_y()
f1.tick_labels.set_xposition('top')
f1.tick_labels.show()


f2 = aplpy.FITSFigure(path+'maps/plotimNUVJy.fits', figure=fig, subplot=plotloc[1])
standard_setup(f2)
#f2.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f2.show_colorscale(vmin=0, vmax=0.0004, cmap='hot')
f2.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
f2.tick_labels.set_xposition('top')
f2.tick_labels.show_x()

f3 = aplpy.FITSFigure(path+'maps/plotimuJy.fits', figure=fig, subplot=plotloc[2])
standard_setup(f3)
#f3.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f3.show_colorscale(vmin=0, vmax=0.004, cmap='hot')
f3.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
f3.tick_labels.set_xposition('top')
f3.tick_labels.show_x()

# Next rows

f4 = aplpy.FITSFigure(path+'maps/plotimgJy.fits', figure=fig, subplot=plotloc[3])
standard_setup(f4)
#f4.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f4.show_colorscale(vmin=0, vmax=0.015, cmap='hot')
f4.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
#f4.axis_labels.show_y()
f4.tick_labels.show_y()


f5 = aplpy.FITSFigure(path+'maps/plotimrJy.fits', figure=fig, subplot=plotloc[4])
standard_setup(f5)
#f5.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f5.show_colorscale(vmin=0, vmax=0.045, cmap='hot')
f5.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')

f6 = aplpy.FITSFigure(path+'maps/plotimiJy.fits', figure=fig, subplot=plotloc[5])
standard_setup(f6)
#f6.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f6.show_colorscale(vmin=0, vmax=0.05, cmap='hot')
f6.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')

f7 = aplpy.FITSFigure(path+'maps/plotimzJy.fits', figure=fig, subplot=plotloc[6])
standard_setup(f7)
#f7.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f7.show_colorscale(vmin=0, vmax=0.07, cmap='hot')
f7.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
#f7.axis_labels.show_y()
f7.tick_labels.show_y()


f8 = aplpy.FITSFigure(path+'maps/plotimW1Jy.fits', figure=fig, subplot=plotloc[7])
standard_setup(f8)
#f8.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f8.show_colorscale(vmin=0, vmax=0.075, cmap='hot')
f8.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')

f9 = aplpy.FITSFigure(path+'maps/plotim3.6Jy.fits', figure=fig, subplot=plotloc[8])
standard_setup(f9)
#f9.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f9.show_colorscale(vmin=0, vmax=0.075, cmap='hot')
f9.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')

f10 = aplpy.FITSFigure(path+'maps/plotim4.5Jy.fits', figure=fig, subplot=plotloc[9])
standard_setup(f10)
#f10.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f10.show_colorscale(vmin=0, vmax=0.055, cmap='hot')
f10.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
#f10.axis_labels.show_y()
f10.tick_labels.show_y()


f11 = aplpy.FITSFigure(path+'maps/plotimW2Jy.fits', figure=fig, subplot=plotloc[10])
standard_setup(f11)
#f11.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f11.show_colorscale(vmin=0, vmax=0.055, cmap='hot')
f11.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')

f12 = aplpy.FITSFigure(path+'maps/plotim5.8Jy.fits', figure=fig, subplot=plotloc[11])
standard_setup(f12)
#f12.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f12.show_colorscale(vmin=0, vmax=0.075, cmap='hot')
f12.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')

f13 = aplpy.FITSFigure(path+'maps/plotim8Jy.fits', figure=fig, subplot=plotloc[12])
standard_setup(f13)
#f13.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f13.show_colorscale(vmin=0, vmax=0.08, cmap='hot')
f13.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
#f13.axis_labels.show_y()
f13.tick_labels.show_y()

f14 = aplpy.FITSFigure(path+'maps/plotimW3Jy.fits', figure=fig, subplot=plotloc[13])
standard_setup(f14)
#f14.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f14.show_colorscale(vmin=0, vmax=0.06, cmap='hot')
f14.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')

f15 = aplpy.FITSFigure(path+'maps/plotimW4Jy.fits', figure=fig, subplot=plotloc[14])
standard_setup(f15)
#f15.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f15.show_colorscale(vmin=0, vmax=0.05, cmap='hot')
f15.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')

f16 = aplpy.FITSFigure(path+'maps/plotim24Jy.fits', figure=fig, subplot=plotloc[15])
standard_setup(f16)
#f16.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f16.show_colorscale(vmin=0, vmax=0.035, cmap='hot')
f16.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
#f16.axis_labels.show_y()
f16.tick_labels.show_y()


f17 = aplpy.FITSFigure(path+'maps/plotim70Jy.fits', figure=fig, subplot=plotloc[16])
standard_setup(f17)
#f17.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f17.show_colorscale(vmin=0, vmax=0.35, cmap='hot')
f17.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')

f18 = aplpy.FITSFigure(path+'maps/plotim100Jy.fits', figure=fig, subplot=plotloc[17])
standard_setup(f18)
#f18.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f18.show_colorscale(vmin=0, vmax=0.7, cmap='hot')
f18.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')

f19 = aplpy.FITSFigure(path+'maps/plotim160Jy.fits', figure=fig, subplot=plotloc[18])
standard_setup(f19)
#f19.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f19.show_colorscale(vmin=0, vmax=1.3, cmap='hot')
f19.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
#f19.axis_labels.show_y()
f19.tick_labels.show_y()


f20 = aplpy.FITSFigure(path+'maps/plotim250Jy.fits', figure=fig, subplot=plotloc[19])
standard_setup(f20)
#f20.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f20.show_colorscale(vmin=0, vmax=1.3, cmap='hot')
f20.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
#f20.axis_labels.show_x()
f20.tick_labels.show_x()

f21 = aplpy.FITSFigure(path+'maps/plotim350Jy.fits', figure=fig, subplot=plotloc[20])
standard_setup(f21)
#f21.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f21.show_colorscale(vmin=0, vmax=0.6, cmap='hot')
f21.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
#f21.axis_labels.show_x()
f21.tick_labels.show_x()

f22 = aplpy.FITSFigure(path+'maps/plotim500Jy.fits', figure=fig, subplot=plotloc[21])
standard_setup(f22)
#f22.show_colorscale(pmax=99.25, pmin=0.50,  cmap='hot')
f22.show_colorscale(vmin=0, vmax=0.3, cmap='hot')
f22.show_beam(major=0.01, minor=0.01, angle=0,fill=True,color='white')
#f22.axis_labels.show()
f22.tick_labels.show()


# Add a colourbar

axisf3 = fig.add_axes([0.45,0.07,0.5,0.02])
cmapf3 = mpl.cm.hot
normf3 = mpl.colors.Normalize(vmin=0, vmax=1)
cbf3 = mpl.colorbar.ColorbarBase(axisf3, cmap=cmapf3, norm=normf3, orientation='horizontal')
cbf3.set_label('Flux (arbitrary units)')

# Add labels

fig.text(0.38,0.915,"FUV",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.68,0.915,"NUV",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.98,0.915,"u",color='white',size='14',weight='bold', horizontalalignment='right')

fig.text(0.38,0.8,"g",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.68,0.8,"r",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.98,0.8,"i",color='white',size='14',weight='bold', horizontalalignment='right')

fig.text(0.38,0.685,"z",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.68,0.685,"W1",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.98,0.685,"3.6",color='white',size='14',weight='bold', horizontalalignment='right')

fig.text(0.38,0.57,"4.5",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.68,0.57,"W2",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.98,0.57,"5.8",color='white',size='14',weight='bold', horizontalalignment='right')

fig.text(0.38,0.455,"8",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.68,0.455,"W3",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.98,0.455,"W4",color='white',size='14',weight='bold', horizontalalignment='right')

fig.text(0.38,0.34,"24",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.68,0.34,"70",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.98,0.34,"100",color='white',size='14',weight='bold', horizontalalignment='right')

fig.text(0.38,0.225,"160",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.68,0.225,"250",color='white',size='14',weight='bold', horizontalalignment='right')
fig.text(0.98,0.225,"350",color='white',size='14',weight='bold', horizontalalignment='right')

fig.text(0.38,0.11,"500",color='white',size='14',weight='bold', horizontalalignment='right')

fig.patch.set_facecolor('#3f3f3f')
fig.canvas.draw()
#fig.savefig("ngc891_figure.eps", dpi=300)
fig.savefig(path+"dataset.eps", dpi=600)
#fig.savefig("ngc891_figure.png", dpi=300)
pyplot.show()