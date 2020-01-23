# AGS_OPT_2D

The integration of Algebraic Graphic Statics and Layout Optimization. 
The algorithm produces form diagram (i.e. a valid strut-and-tie model), force diagram and stress fields, which accounts for the boundary conditions and the design boundary of the reiforced concrete block.

The algorithm is written in Python 3.73 (https://www.python.org), where the bsaic implementations use NumPy (http://www.numpy.org) and SciPy (http://www.scipy.org) packages, the linear optimization is performed with the CVXPY package (https://www.cvxpy.org),
which  solves  convex  optimization  problems.   Also,  the  data structures such as meshes and graphs are created or modified
using COMPAS (http://compas-dev.github.io).

The Sour
