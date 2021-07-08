1. New class derived from pointPatchFields to calculate patch displacement in dynamicMeshes. It is conceived for airfoil or lifting surfaces deflection.


2. Custom pimpleFoam-based solver (pimpleRemeshFoam) that checks mesh quality and calls remeshing routines in Pointwise if needed


(These features are not fully functional)
The point displacement class uses a PID controller to reach a specified setPoint among the following variables:
- deflection angle
- target lift force in any given direction
- target lift force coefficient in any given direction

Contributors:
- Sócrates Fernández
