import yt

# Loading data
ds = yt.load('test/2D_amr_stable0014.dat')

#############################################################
# 1D slice plot
# create a slice from (0,0,0) to (0,0,1) using 1000 points
lpl = yt.LinePlot(ds, 'velocity_1',(1,0), (5,0), 1000)

#Save 1D slice plot
lpl.save('plots/')


#############################################################
# 2D slice plot

#make a 2D plot
p = yt.plot_2d(ds, "m1")
p.set_log("m1", False)

# overplot the grid
p.annotate_grids()

# save the 2D plot
p.save('plots/')
