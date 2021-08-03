caseName="testPimpleRemesh"
solver="pimpleRemeshFoam"
parallel = True
crescent = False
rotateCenter = [0.5828,0,0] # only needed for flap? Yes for the moment
meshType = "airfoil"  # do I need meshtype + AIRFOILREMESH or is it enough with REMESHCOMMAND? Yes bc of rotationCenter
AIRFOILREMESH="pointwise -b meshUnstructuredRotate.glf parameters.dat"
FLAPREMESH="pointwise -b rotateFlap.glf parameters.dat"
testing=True