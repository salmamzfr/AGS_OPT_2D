# AGS_OPT_2D

The integration of Algebraic Graphic Statics and Layout Optimization. 
The algorithm produces form diagram (i.e. a valid strut-and-tie model), force diagram and stress fields, which accounts for the boundary conditions and the design boundary of the reiforced concrete block.

The algorithm is developed by Salma Moaffari at ETH Zurich. It is written in Python 3.73 (https://www.python.org), where the bsaic implementations use NumPy (http://www.numpy.org) and SciPy (http://www.scipy.org) packages, the linear optimization is performed with the CVXPY package (https://www.cvxpy.org), which  solves  convex  optimization  problems.   Also,  the  data structures such as meshes and graphs are created or modified using COMPAS (http://compas-dev.github.io).

The "Source" folder includes the main scripts. "AGS_OPT_2D" runs the algorithm after loading a specified mesh from the "MeshObjects" folder. The pre-defined meshes in "MeshObjects" are generated in Rhino (https://www.rhino3d.com) and exported as .obj file. The insctruction on how to run the algorithms is included in "AGS_OPT_2D" file. 

The scripts in "Rhino" folder help in re-prodcution of the results (i.e. form, force and stress fields) in Rhino, after running the main script "AGS_OPT_2D". 

To cite the work, please use below BibTeX:

@misc{ags_opt,<br/>
author = {Mozaffari, Salma},<br/>
title = {AGS_OPT: The integration of Algebraic Graphic Statics and Layout Optimization.},<br/>
year = {2020},<br/>
url = {https://github.com/salmamzfr/AGS_OPT_2D},<br/>
}

In case of questions, please contact mozaffari@ibk.baug.ethz.ch
