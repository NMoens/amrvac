import yt

# Loading data
ds = yt.load('output_br/bw_3d0003.dat')
#############################################################
# 1D slice plot
# create a slice from (0,0,0) to (0,0,1) using 1000 points
lpl = yt.LinePlot(ds, 'density',(0,0,0), (0,0,0.5), 1000)

#Save 1D slice plot
lpl.save('slice1D')

#############################################################
# 1D slice plot in multiple anglesi
# Loading data
ds = yt.load('output_br/bw_3d0003.dat')

# create an empty array of lines
lines = []

# Add slices to your array
lines.append(yt.LineBuffer(ds, (0,0,0), (0,0.0,0.5),100, label='$n_z$'))
lines.append(yt.LineBuffer(ds, (0,0,0), (0,0.5,0.5),100, label='$n_y+n_z$'))
lines.append(yt.LineBuffer(ds, (0,0,0), (0.5,0.5,0.5),100, label='$n_x+n_y+n_z$'))

# plot lines
plot = yt.LinePlot.from_lines(ds, 'density', lines)

# add legend()
plot.annotate_legend('density')

#Save 1D slice plot
plot.save('multi_slice1D')

#############################################################
# # 1D slice plot on multiple times
# # Loading data at different timesteps
# es = [yt.load('output/bw_3d0000.dat'), yt.load('output/bw_3d0010.dat'), yt.load('output/bw_3d0020.dat'), yt.load('output/bw_3d0030.dat')]
#
# # create an empty array of profiles
# profiles = []
# labels = []
#
# # Add slices to your array
# for ds in es:
#     #profiles.append(yt.LineBuffer(ds, (0,0,0), (0,0.0,0.5),100, label='ds.get_time()'))
#     line = yt.LineBuffer(ds, (0,0,0), (0,0.0,0.5),100, label='t = ' + str(ds.current_time))
#     plot = yt.LinePlot.from_lines(ds,'density',[line])
#
#
# # plot lines
# #plot = yt.LinePlot.from_lines(es,'density',profiles)
#
# # add legend()
# plot.annotate_legend('density')
#
# #Save 1D slice plot
# plot.save('multi_t_slice1D')



#############################################################
# 2D slice plot
# Loading data
ds = yt.load('output_br/bw_3d0003.dat')
# make a 2D plot of the density
slc = yt.SlicePlot(ds, 'z', 'density')

# overplot the grid
slc.annotate_grids()

# save the 2D plot
slc.save('slice2D')
