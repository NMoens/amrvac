from amrvac_tools.datfiles.reading import amrvac_reader
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors

import os, glob
import tkinter as tk
from tkinter import filedialog
from mpl_toolkits.axes_grid1 import make_axes_locatable
import asyncio as syn

from ipywidgets import widgets, interactive

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
# press a
# press o
# press pgd
# press pgu
################
# press i
# press shift i

unit_length=69599000000.0
unit_numberdensity=729723637293892.50
unit_temperature=189666.48683552662

c_light = 2.99792458e10
a_rad = 7.5657e-15
G_grav = 6.67191e-8
k_b = 1.380648520e-16
m_p = 1.6737236e-24
mu = 0.6
M_sun = 1.9891000e33
M_star = 1e1*M_sun

hd_gamma = 1.66667


file_counter = 0
var_counter = 0
loglin = False
lineout = False
rel_diff = False
cgs_units = False
my_var = False
fix_axes = False
prim_cons = False
plot_orig = False



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
        global fix_axes
        global prim_cons
        global plot_orig

        if event.key == 'right':
            if (file_counter + 1) < len(files):
                file_counter = file_counter + 1
            else:
                file_counter = 0
            print('file counter:  ',file_counter)
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'pagedown':
            if (file_counter + 20) < len(files):
                file_counter = file_counter + 20
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

        if event.key == 'pageup':
            if (file_counter - 20) >= 0:
                file_counter = file_counter - 20
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

        if event.key == 'a':
            print('Toggle fixed axes vs automatic axes')
            if fix_axes:
                fix_axes = False
            else:
                fix_axes = True
            print('Fixing axes:   ',fix_axes)
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'c':
            print('Toggle primitive vs conservative')
            if prim_cons:
                prim_cons = False
            else:
                prim_cons = True
            print('Primitive variables:   ',prim_cons)
            self.plotter(files[file_counter],variables[var_counter])

        if event.key == 'o':
            print('Overplot t=0 profile')
            if plot_orig:
                plot_orig = False
            else:
                plot_orig = True
            print('Plotting t=0 profile:   ',plot_orig)
            self.plotter(files[file_counter],variables[var_counter])

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

        ds.units.set_units(unit_length=unit_length, unit_numberdensity=unit_numberdensity, unit_temperature=unit_temperature)

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
        elif varname =='M_dot':
            return ds.units.unit_density*ds.units.unit_length**3/ds.units.unit_time*(365.25*24*60*60)/M_sun
        elif varname =='T_gas' or varname == 'T_rad':
                return ds.units.unit_temperature
        elif varname == 'Lambda' or varname == 'gamma' or varname == 'fld_R' or varname == 'D':
            return 1.e0
        else:
            print('Varname not known/defined')
            return 1.e0

    def calc_myvariable(self,filepath,varname):
        ds = amrvac_reader.load_file(filepath)
        ds.units.set_units(unit_length=unit_length, unit_numberdensity=unit_numberdensity, unit_temperature=unit_temperature)
        ad = ds.load_all_data()
        x,y = ds.get_coordinate_arrays()
        nx = len(x)
        ny = len(y)

        r = []
        for i in range(nx):
            r.append(y)
        r = np.array(r)

        #> Define your variables here:
        if varname == 'a2':
            p = (hd_gamma - 1)*(ad['e'] - (0.5*(ad['m1']**2+ad['m2']**2)/ad['rho']))
            a2 = p/ad['rho']
            return a2

        if varname == 'g_rad':
            grad = ad['Kappa']*ad['F2']/c_light*ds.units.unit_velocity
            return grad

        if varname == 'g_grav':
            ggrav = G_grav*M_star/(r*ds.units.unit_length)**2\
            *(ds.units.unit_time**2/ds.units.unit_length)
            return ggrav

        if varname == 'Gamma':
            grad = ad['Kappa']*ad['F2']/(c_light/ds.units.unit_velocity)
            ggrav = G_grav*M_star/(r*ds.units.unit_length)**2\
            *(ds.units.unit_time**2/ds.units.unit_length)
            Gamma = grad/ggrav
            return Gamma

        if varname == 'M_dot':
            M_dot = 4.*np.pi*ad['m2']*r**2.
            return M_dot

        if varname == 'T_rad':
            T_rad = (ad['r_e']*ds.units.unit_pressure/a_rad)**0.25/ds.units.unit_temperature
            return T_rad

        if varname == 'T_gas':
            p = (hd_gamma - 1)*(ad['e'] - (0.5*(ad['m1']**2+ad['m2']**2)/ad['rho']))
            T_gas = m_p*mu/k_b*(p*ds.units.unit_pressure/ad['rho']/ds.units.unit_density)/ds.units.unit_temperature
            return T_gas

        if varname == 'Av_rho':
            Av_rho = ad['int_r']/ad['int_dt']
            return Av_rho

        if varname == 'Av_v':
            Av_v = ad['int_v']/ad['int_dt']
            return Av_v

        if varname == 'Av_re':
            Av_re = ad['int_re']/ad['int_dt']
            return Av_re

        if varname == 'Gamma_weighted':
            grad = ad['Kappa']*ad['F2']/(c_light/ds.units.unit_velocity)
            ggrav = G_grav*M_star/(r*ds.units.unit_length)**2\
            *(ds.units.unit_time**2/ds.units.unit_length)
            Gamma = grad/ggrav
            Gamma_w = Gamma*1./3./ad['Lambda']
            return Gamma_w

        if varname == 'tau':
            dy = (max(y) - min(y))/ny
            dtau = ad['Kappa']*ad['rho']*dy
            tau = np.zeros(np.shape(dtau))
            nx = np.shape(dtau)[0]
            ny = np.shape(dtau)[1]
            for j in range(nx):
                tau[j][ny-1] = dtau[j][ny-1]
                for i in range(ny-2,-1,-1):
                    tau[j][i] = tau[j][i+1] + dtau[j][i]
            return tau

        else:
            print('variable not defined')


    def plotter(self, filepath, varname):

        #Read data from current step
        ds = amrvac_reader.load_file(filepath)
        bounds_x, bounds_y = ds.get_bounds()
        time = ds.get_time()
        ad = ds.load_all_data()
        x,y = ds.get_coordinate_arrays()
        nx = len(x)
        ny = len(y)

        #Read data from iniial step
        ds0 = amrvac_reader.load_file(files[0])
        ad0 = ds0.load_all_data()


        #For time-integrated data
        if varname == 'Av_rho':
            data = self.calc_myvariable(filepath, varname)
            data0 = ad0['rho']
        elif varname == 'Av_v':
            data = self.calc_myvariable(filepath, varname)
            data0 = ad0['m2']/ad0['rho']
        elif varname == 'Av_re':
            data = self.calc_myvariable(filepath, varname)
            data0 = ad0['r_e']
        #For selfdefined variables, read/calculate data for current and initial snap
        elif varname in extra_vars:
            data = self.calc_myvariable(filepath, varname)
            data0 = self.calc_myvariable(files[0], varname)
        #For original variables: rho, m, e, Er,....
        else:
            data = ad[varname]
            data0 = ad0[varname]



        if rel_diff:
            data = abs((data-data0)/data0)


        #Convert to pirmitive variables if necessary
        if prim_cons:
            if varname == 'm1':
                varname = 'v1'
                data = data/ad['rho']
                data0 = data0/ad0['rho']
            if varname == 'm2':
                varname = 'v2'
                data = data/ad['rho']
                data0 = data0/ad0['rho']
            if varname == 'e':
                varname = 'p'
                data = (hd_gamma - 1)*(data - (0.5*(ad['m1']**2+ad['m2']**2)/ad['rho']))
                data0 = (hd_gamma - 1)*(data0 - (0.5*(ad0['m1']**2+ad0['m2']**2)/ad0['rho']))

        #Multiply data with correct units
        if cgs_units:
            if not rel_diff:
                data = data*self.units(filepath,varname)
                data0 = data0*self.units(filepath,varname)


        #Clear figure
        plt.clf()

        ax = self.fig.add_subplot(111)

        #Plot the 1D profiles
        if lineout:
            ax.grid(which='both')
            #Loglog plot
            if loglin and data.any() > 0:
                for i in range(nx):
                    ax.loglog(y,data[i],'b-')
                if plot_orig:
                    ax.loglog(y,np.mean(data0,axis=0),'k-')
                ax.loglog(y,np.mean(data,axis=0),'ro')
            #Linear plot
            else:
                for i in range(nx):
                    ax.plot(y,data[i],'b-')
                if plot_orig:
                    ax.plot(y,np.mean(data0,axis=0),'k-')
                ax.plot(y,np.mean(data,axis=0),'ro')


            #Fix the ylim axes
            if fix_axes:
                d_lim = 0.05*(data0.max() - data0.min())
                ax.set_ylim(data0.min() - d_lim , data0.max() + d_lim)
            else:
                d_lim = 0.05*(np.mean(data,axis=0).max() - np.mean(data,axis=0).min())
                ax.set_ylim(np.mean(data,axis=0).min() - d_lim ,np.mean(data,axis=0).max() + d_lim)



        #Plot the pseudocolor using imshow
        else:
            norm = None
            #logarithmic colormap
            if loglin:
                #Fix colorbar to value in initial snap
                if fix_axes:
                    norm = matplotlib.colors.LogNorm(vmin=data0.min(), vmax=data0.max())
                else:
                    norm = matplotlib.colors.LogNorm()

            #Linear colormap
            else:
                #Fix colorbar to value in initial snap
                if fix_axes:
                    norm = matplotlib.colors.Normalize(vmin=data0.min(), vmax=data0.max())


            im = ax.imshow(np.rot90(data), extent=[*bounds_x, *bounds_y], norm = norm)
            ax.set_aspect('equal')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            self.fig.colorbar(im, cax=cax)

        #Set title including variable name and timestamp
        fname = os.path.split(filepath)[1]


        if cgs_units:
            plt.title(varname + '   at T=' + str(round(time/ds.units.unit_time,4)) )
        else:
            plt.title(varname + '   at T=' + str(round(time,4)) )
        plt.draw()




if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()

    dir_path = filedialog.askdirectory()
    files = list(glob.glob(os.path.join(dir_path, 'const*.dat')))
    files.sort()

    ds = amrvac_reader.load_file(files[0])
    orig_variables = ds.get_varnames()
        #> I HAD TO DEFINE THIS FUNCTION IN NIELS' TOOLS
    extra_vars = ['M_dot', 'g_rad', 'Gamma', 'T_gas', 'Av_rho',  'Av_v',  'Av_re','tau']
    delete_vars = []
    variables = ds.get_varnames() #.extend(extra_vars)
    variables.extend(extra_vars)

    #Get rid of the following variables
    # variables.remove('Edd11')
    variables.remove('Edd12')
    variables.remove('Edd21')
    variables.remove('Edd22')

    variables.remove('int_r')
    variables.remove('int_v')
    variables.remove('int_re')
    variables.remove('int_dt')

    Plot = Plotter(files[file_counter],variables[var_counter])
