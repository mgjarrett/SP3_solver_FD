This code solves a 3D finite difference diffusion or SP3 problem with CMFD 
acceleration. The solution is obtained using superLU decomposition, or GMRES 
if the system is too large for superLU.

The number of meshes, solution method, and numerical method are input options. 
Other chagnes require simple modifications to the source code, such as the 
dimension of the problem or the cross sections.
