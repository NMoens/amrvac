from amrvac_tools.datfiles.reading import amrvac_reader
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors

import os, glob
import tkinter as tk
from tkinter import filedialog
from mpl_toolkits.axes_grid1 import make_axes_locatable
import asyncio as syn

# Current possibilities:
# press right
# press left
# press up
# press down
# press shift right
# press shift left
# press shift up
# press shift down
# press t
# press l
# press d
# press u
# press m
# press i
# press shift i






file_counter = 0
var_counter = 0
loglin = False
lineout = False
rel_diff = False
cgs_units = False
my_var = False



# class variables:
#
#     def g_rad(self,filepath):
#
#         return g_rad
#
#     def g_grav(self,filepath):
#
#         return g_grav
#
#     def Gamma(self,filepath):
#
#         return Gamma
#
#     def Trad(self,filepath):
#
#         return Trad
#
#     def Tgas(self,filepath):
#
#         return Tgas


class Plotter:
    def __init__(self, filepath,variable):
        self.fig = plt.figure()
        self.cid = self.fig.canvas.mpl_connect('key_press_event', self.press)
        self.plotter(filepath,variable)
        plt.show()

    def press(self, event):
        print('press', event.key)
        global file_counter
        global var_counter
        global loglin
        global lineout
        global rel_diff
        global cgs_units
        global my_var

        if event.key == 'right':
            if (file_counter + 1) < len(files):
                file_counter = file_counter + 1
            else:
                file_counter = 0
            print('file counter:  ',file_counter)
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'left':
            if (file_counter - 1) >= 0:
                file_counter = file_counter - 1
            else:
                file_counter = len(files) -1
            print('file counter:  ',file_counter)
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'down':
            if (var_counter + 1) < len(variables):
                var_counter = var_counter + 1
            else:
                var_counter = 0
            print('variable counter:  ',var_counter)
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'up':
            if (var_counter - 1) >= 0:
                var_counter = var_counter - 1
            else:
                var_counter = len(variables) -1
            print('variable counter:  ',var_counter)
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'shift+left':
            file_counter = 0
            print('going to the first variable')
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'shift+right':
            file_counter = len(files) - 1
            print('going to the last variable')
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'shift+up':
            var_counter = 0
            print('going to the first variable')
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'shift+down':
            var_counter = len(variables) - 1
            print('going to the last variable')
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 't':
            print('Toggle between linear and logarithmic colorbar')
            if loglin:
                loglin = False
            else:
                 loglin= True
            print('logarithmic scale:  ',loglin)
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'l':
            print('Toggle between lineout and pseudocolor')
            if lineout:
                lineout = False
            else:
                lineout = True
            print('Lineout plot:   ',lineout)
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'd':
            print('Toggle between relative difference with respect to first snapshot')
            if rel_diff:
                rel_diff = False
            else:
                rel_diff = True
            print('relative difference plot:   ',rel_diff)
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'u':
            print('Toggle between dimensionless and cgs units')
            if cgs_units:
                cgs_units = False
            else:
                cgs_units = True
            print('cgs units:   ',cgs_units)
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'm':
            print('Toggle standard variables and my defined variable')
            if my_var:
                my_var = False
            else:
                my_var = True
            print('plotting my variable:   ',my_var)
            self.plotter(files[file_counter],'my_variable1')

        if event.key == 'i':
            self.fig.canvas.mpl_disconnect(self.cid)
            print('Go to variable:')
            # gotovar = str(input("Enter variable name"))
            self.cid = self.fig.canvas.mpl_connect('key_press_event', self.press)
            try:
                var_counter = variables.index(gotovar)
                self.plotter(files[file_counter],variables[var_counter])
            except:
                print("Variable not recognised in list")
                self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'I':
            print('Go to snap:')

            gotosnap = int(input("Enter index of snap"))
            try:
                file_counter = gotosnap
                self.plotter(files[file_counter],variables[var_counter])
            except:
                print("snap number not recognised in list")
                self.plotter(files[file_counter],variables[var_counter])


    def units(self, filepath, varname):
        ds = amrvac_reader.load_file(filepath)

        ds.units.set_units(unit_length=69599000000.0, unit_numberdensity=729723637293892.50, unit_temperature=189666.48683552662)

        if varname == 'rho':
            return ds.units.unit_density
        elif varname == 'v1' or varname == 'v2':
            return ds.units.unit_velocity
        elif varname == 'm1'  or varname == 'm2':
            return ds.units.unit_density*ds.units.unit_velocity
        elif varname == 'e' or varname == 'r_e' or varname == 'p':
            return ds.units.unit_pressure
        elif varname == 'Edd11' or varname == 'Edd12' or varname == 'Edd21' or varname == 'Edd22':
            return ds.units.unit_pressure
        elif varname == 'Kappa':
            return 1.0/(ds.units.unit_density*ds.units.unit_length)
        elif varname == 'F1' or varname == 'F2':
            return ds.units.unit_velocity*ds.units.unit_pressure
        elif varname == 'Lambda' or varname == 'gamma' or varname == 'fld_R' or varname == 'D':
            return 1.e0
        else:
            print('Varname not known/defined')
            return 1.e0

    def calc_myvariable(self,filepath):
        ds = amrvac_reader.load_file(filepath)
        ad = ds.load_all_data()
        x,y = ds.get_coordinate_arrays()
        data = (ad['e']-(ad['m1']**2 + ad['m2']**2)/(0.5*ad['rho']))*(1.666667 -1.0)
        return data

    def plotter(self, filepath, varname):

        ds = amrvac_reader.load_file(filepath)
        bounds_x, bounds_y = ds.get_bounds()
        time = ds.get_time()
        ad = ds.load_all_data()
        x,y = ds.get_coordinate_arrays()
        nx = len(x)
        ny = len(y)

        if varname == 'my_variable1':
            data = self.calc_myvariable(filepath)[0]
        else:
            data = ad[varname]

        if rel_diff:
            ds0 = amrvac_reader.load_file(files[0])
            ad0 = ds0.load_all_data()
            if varname == 'my_variable1':
                data0 = self.calc_myvariable(ad0)
            else:
                data0 = ad0[varname]

            data = abs((data-data0)/data0)



        if cgs_units:
            time = time*ds.units.unit_time
            if not rel_diff:
                data = data*self.units(filepath,varname)



        plt.clf()

        ax = self.fig.add_subplot(111)

        if lineout:
            if loglin and data.any() > 0:
                for i in range(nx):
                    ax.loglog(y,data[i],'b-')
                ax.loglog(y,np.mean(data,axis=0),'ro')

            else:
                for i in range(nx):
                    ax.plot(y,data[i],'b-')
                ax.plot(y,np.mean(data,axis=0),'ro')

        else:
            norm = None
            if loglin:
                norm = matplotlib.colors.LogNorm()

            im = ax.imshow(np.rot90(data), extent=[*bounds_x, *bounds_y], norm = norm)
            ax.set_aspect('equal')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            self.fig.colorbar(im, cax=cax)

        fname = os.path.split(filepath)[1]
        plt.title(varname + '   at T=' + str(round(time,4)) )
        plt.draw()




if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()

    dir_path = filedialog.askdirectory()
    files = list(glob.glob(os.path.join(dir_path, 'const*.dat')))
    files.sort()

    ds = amrvac_reader.load_file(files[0])
    variables = ds.get_varnames()
    #> I HAD TO DEFINE THIS FUNCTION IN NIELS' TOOLS


    Plot = Plotter(files[file_counter],variables[var_counter])
