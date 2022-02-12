## Work in progress
This repository contains the code in a raw state that allowed me to complete my MSc thesis, but there are plenty of details to correct before I can call it a "release". This README will be updated when it is ready.

## Summary of contents

1. New class PIDangularDisplacement derived from fixeValuePointPatchField<vector> to calculate patch displacement in dynamicMeshes. It is conceived for airfoil or lifting surfaces automatic deflection through rotation around a point using a PID controller


2. A solver called pimpleRemeshFoam that checks mesh quality, writes and stops the simulation if bad quality is encountered. To be used along with *pimpleRemesh.py* as process handler for remeshing. Pointwise is the default remesher but the call can be customised in the file *parameters.py*


The point displacement class uses a PID controller to reach a specified setPoint among the following variables:
- deflection angle
- target lift force in any given direction
- target lift force coefficient in any given direction

Contributors:
- Sócrates Fernández
- Tom-Robin Teschner [Thesis supervisor]: (automated pointwise airfoil mesher: meshUnstructured.glf)
