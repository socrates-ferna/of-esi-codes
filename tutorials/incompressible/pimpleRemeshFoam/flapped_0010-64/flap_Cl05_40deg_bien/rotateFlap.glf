# Load Pointwise Glyph package and Tk
package require PWI_Glyph 2.4
pw::Application reset

set scriptDir [file dirname [info script]]
puts $scriptDir
# Source external procedure script
source [file join $scriptDir "myprocs.glf"]


set parameterfile "parameters.dat"

# read the parameters from the parameter file
puts "Reading parameter file from: $parameterfile ..."
set file [open $parameterfile r]
while {![eof $file]} {
    set part [split [gets $file] "="]
    set parameters([string trimright [lindex $part 0]]) [string trimleft [lindex $part 1]]
}
close $file

# Open the templated mesh
Openpw $scriptDir

# Rotate the flap connectors. Domain should reinitialize automatically
rotateFlap

# outputting mesh to disk
puts "Writing mesh to disk ..."
writeMeshToDisk
puts "Finished with mesh generation!"