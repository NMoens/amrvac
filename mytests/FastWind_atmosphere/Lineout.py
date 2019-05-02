import sys
import os
sys.path.append("/lhome/nicolasm/Visit/2.13.2/linux-x86_64/lib/site-packages")
from visit import *
Launch()


local = os.getcwd()
OpenDatabase(local + "/output/const0000.vtu")

DeleteAllPlots()
AddPlot("Pseudocolor","rho")
DrawPlots()
Lineout((0.0,0.0), (0.0,1.0))





# AddOperator("Lineout")
#
# C = CurveAttributes()
# C.CurveColor('r')
# SetCurveOptions(C)
#
# L = LineoutAttributes()
# L.point1 = (0.0,0.0)
# L.point2 = (0.0,1.0)
# SetOperatorOptions(L)
#
#
# DrawPlots()
