# AGS-OPT
<a href="https://zenodo.org/badge/latestdoi/235868917"><img src="https://zenodo.org/badge/235868917.svg" alt="DOI"></a>

## The integration of Algebraic Graphic Statics and Layout Optimization

The algorithm produces a form diagram(i.e., a valid strut-and-tie model), force diagram, and stress fields, which account for the boundary conditions and the design boundary of the reinforced concrete block.

The algorithm is developed by Salma Moaffari at ETH Zurich. It is written in Python 3.7(https://www.python.org), where the basic implementations use NumPy(http://www.numpy.org) and SciPy(http://www.scipy.org) packages, the linear optimization is performed with the CVXPY package(https://www.cvxpy.org), which solves convex optimization problems. Also, the data structures such as meshes and graphs are created or modified using COMPAS 15.2(http://compas-dev.github.io).

The "Source" folder includes the main scripts. "AGS_OPT_2D" runs the algorithm after loading a specified mesh from the "MeshObjects" folder. The pre-defined meshes in "MeshObjects" are generated in Rhino (https://www.rhino3d.com) and exported as a *.obj* file. The user can also define any mesh object, load it to the "AGS_OPT_2D", and modify the boundary conditions to run an arbitrary problem. The instruction on how to run the algorithm is included in the "AGS_OPT_2D" file. 

The scripts in the "Rhino" folder help in the re-production of results (i.e., form, force, and stress fields) in Rhino, after running the main script "AGS_OPT_2D". 

To cite the work, please use below BibTeX:

@misc{ags-opt,<br/>
  author = {Mozaffari, Salma},<br/>
  title = {AGS-OPT: The integration of Algebraic Graphic Statics and Layout Optimization},<br/>
  year = {2020},<br/>
  version = {1.0.0},<br/>
  doi = {10.5281/zenodo.3628672},<br/>
  url = { https://github.com/salmamzfr/AGS_OPT_2D },<br/>
}

In case of questions, please contact mozaffari@ibk.baug.ethz.ch
