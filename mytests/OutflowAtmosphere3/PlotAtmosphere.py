from amrvac_tools.datfiles.reading import amrvac_reader
import numpy as np







def PlotSnapshot(file):
    ds = amrvac_reader.load_file(file)
    ds.get_info()

    ad = ds.load_all_data()
    rho = ad['lambda']

    bounds_x, bounds_y = ds.get_bounds()

    x,y = ds.get_coordinate_arrays()

    print(rho)

    print(x)

    return

PlotSnapshot('BigSimComplete_pert/const0000.dat')
