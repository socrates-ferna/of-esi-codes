# AIRFOIL GENERATION PROCEDURE
# -----------------------------------------------
proc airfoilGen {} {

# AIRFOIL INPUTS
# -----------------------------------------------
# m = maximum camber 
# p = maximum camber location 
# t = maximum thickness
set m [expr {[string index $::parameters(naca) 0]/100.0}]  
set p [expr {[string index $::parameters(naca) 1]/10.0}] 
set a [string index $::parameters(naca) 2]
set b [string index $::parameters(naca) 3]
set c "$a$b"
scan $c %d c
set t [expr {$c/100.0}]

# GENERATE AIRFOIL COORDINATES
# -----------------------------------------------
# Initialize Arrays
set x {}
set xu {}
set xl {}
set yu {}
set yl {}
set yc {0}
set yt {}

# Airfoil step size
set ds 0.001

# Check if airfoil is symmetric or cambered
if {$m == 0 && $p == 0 || $m == 0 || $p == 0} {set symm 1} else {set symm 0}

# Get x coordinates
for {set i 0} {$i < [expr {1+$ds}]} {set i [expr {$i+$ds}]} {lappend x $i}

# Calculate mean camber line and thickness distribution
foreach xx $x {

	# Mean camber line definition for symmetric geometry
	if {$symm == 1} {lappend yc 0}

	# Mean camber line definition for cambered geometry
	if {$symm == 0 && $xx <= $p} {
		lappend yc [expr {($m/($p**2))*(2*$p*$xx-$xx**2)}]
	} elseif {$symm == 0 && $xx > $p} {
		lappend yc [expr {($m/((1-$p)**2)*(1-2*$p+2*$p*$xx-$xx**2))}]
	}

	# Thickness distribution
	lappend yt [expr {($t/0.20)*(0.29690*sqrt($xx)-0.12600*$xx- \
	                  0.35160*$xx**2+0.28430*$xx**3-0.10150*$xx**4)}]

	# Theta
	set dy [expr {[lindex $yc end] - [lindex $yc end-1]}]
	set th [expr {atan($dy/$ds)}]

	# Upper x and y coordinates
	lappend xu [expr {$xx-[lindex $yt end]*sin($th)}]
	lappend yu [expr {[lindex $yc end]+[lindex $yt end]*cos($th)}]

	# Lower x and y coordinates
	lappend xl [expr {$xx+[lindex $yt end]*sin($th)}]
	lappend yl [expr {[lindex $yc end]-[lindex $yt end]*cos($th)}]

}

# GENERATE AIRFOIL GEOMETRY
# -----------------------------------------------
# Create upper airfoil surface
set airUpper [pw::Application begin Create]
set airUpperPts [pw::SegmentSpline create]

for {set i 0} {$i < [llength $x]} {incr i} {
	$airUpperPts addPoint [list [lindex $xu $i] [lindex $yu $i] 0]
}

set airUpperCurve [pw::Curve create]
$airUpperCurve addSegment $airUpperPts
$airUpper end

# Create lower airfoil surface
set airLower [pw::Application begin Create]
set airLowerPts [pw::SegmentSpline create]

for {set i 0} {$i < [llength $x]} {incr i} {
	$airLowerPts addPoint [list [lindex $xl $i] [lindex $yl $i] 0]
}

set airLowerCurve [pw::Curve create]
$airLowerCurve addSegment $airLowerPts
$airLower end

# Create flat trailing edge
set airTrail [pw::Application begin Create]
set airTrailPts [pw::SegmentSpline create]
$airTrailPts addPoint [list [lindex $xu end] [lindex $yu end] 0]
$airTrailPts addPoint [list [lindex $xl end] [lindex $yl end] 0]
set airTrailCurve [pw::Curve create]
$airTrailCurve addSegment $airTrailPts
$airTrail end

# Zoom to airfoil
pw::Display resetView

}

# BOUNDARY LAYER MESH GENERATION PROCEDURE
# -----------------------------------------------
proc airfoilMesh {} {

# CONNECTOR CREATION, DIMENSIONING, AND SPACING
# -----------------------------------------------
# Get all database entities
set dbEnts [pw::Database getAll]

# Get the curve length of all db curves
foreach db $dbEnts {
    lappend crvLength [$db getLength 1.0]
}

# Find trailing edge from minimum curve length
if {[lindex $crvLength 0] < [lindex $crvLength 1]} {
    set min 0
} else {
    set min 1
}

if {[lindex $crvLength $min] < [lindex $crvLength 2]} {
    set min $min
} else {
    set min 2
}

set dbTe [lindex $dbEnts $min]

# Get upper and lower surfaces
foreach db $dbEnts {
    if {$db != $dbTe} {
        lappend upperLower $db
    }
}

# Find y values at 50 percent length of upper and lower surfaces
set y1 [lindex [[lindex $upperLower 0] getXYZ -arc 0.5] 1]
set y2 [lindex [[lindex $upperLower 1] getXYZ -arc 0.5] 1]

# Determine upper and lower surface db entities
if {$y1 < $y2} {
    set dbLower [lindex $upperLower 0]
    set dbUpper [lindex $upperLower 1]
} else {
    set dbLower [lindex $upperLower 1]
    set dbUpper [lindex $upperLower 0]
}

# Create connectors on database entities
set upperSurfCon [pw::Connector createOnDatabase -type Unstructured $dbUpper]
set lowerSurfCon [pw::Connector createOnDatabase -type Unstructured $dbLower]
set trailSurfCon [pw::Connector createOnDatabase -type Unstructured $dbTe]
$upperSurfCon setName airfoilUpper
$lowerSurfCon setName airfoilLower
$trailSurfCon setName airfoilTrailing
set cons "$upperSurfCon $lowerSurfCon $trailSurfCon"

# Calculate main airfoil connector dimensions
foreach con $cons {lappend conLen [$con getLength -arc 1]}
set upperSurfConLen [lindex $conLen 0]
set lowerSurfConLen [lindex $conLen 1]
set trailSurfConLen [lindex $conLen 2]
set conDim [expr int($::parameters(numpts))]

# Dimension upper and lower airfoil surface connectors
$upperSurfCon setDimension $conDim
$lowerSurfCon setDimension $conDim

# Dimension trailing edge airfoil connector
set teDim [expr int($trailSurfConLen/(100*$::parameters(initds)))+1]
if {$teDim < 11} {
  set teDim 11
} elseif {$teDim > 25} {
  set teDim 25
}
$trailSurfCon setDimension $teDim

# Set leading and trailing edge connector spacings
# set ltDs [expr 10*$initds]
set ltDs $::parameters(telespacing)

set upperSurfConDis [$upperSurfCon getDistribution 1]
set lowerSurfConDis [$lowerSurfCon getDistribution 1]
set trailSurfConDis [$trailSurfCon getDistribution 1]

$upperSurfConDis setBeginSpacing $ltDs
$upperSurfConDis setEndSpacing $ltDs
$lowerSurfConDis setBeginSpacing $ltDs
$lowerSurfConDis setEndSpacing $ltDs
$trailSurfConDis setBeginSpacing $ltDs
$trailSurfConDis setEndSpacing $ltDs

# rotate connectors from user specified
set _CN(1) [pw::GridEntity getByName airfoilUpper]
set _CN(2) [pw::GridEntity getByName airfoilLower]
set _CN(3) [pw::GridEntity getByName airfoilTrailing]
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(1) $_CN(2) $_CN(3)]]
  pw::Entity transform [pwu::Transform rotation -anchor {0.25 0 0} {0 0 1} $::parameters(rotateAngle)] [$_TMP(mode_1) getEntities]
$_TMP(mode_1) end
unset _TMP(mode_1)


# create north domain boundary
set northBC [pw::Application begin Create]
set northSegment [pw::SegmentSpline create]
$northSegment addPoint [list 0 $::parameters(domainHeight) 0]
$northSegment addPoint [list $::parameters(domainLength) $::parameters(domainHeight) 0]
set northConnector [pw::Connector create]
$northConnector addSegment $northSegment
unset northSegment
$northConnector setDimensionFromSpacing -resetDistribution $::parameters(farfieldspacing)
$northConnector setName north
$northBC end
unset northBC

# create east domain boundary
set eastBC [pw::Application begin Create]
set eastSegment [pw::SegmentSpline create]
$eastSegment addPoint [list $::parameters(domainLength) $::parameters(domainHeight) 0]
$eastSegment addPoint [list $::parameters(domainLength) -$::parameters(domainHeight) 0]
set eastConnector [pw::Connector create]
$eastConnector addSegment $eastSegment
unset eastSegment
$eastConnector setDimensionFromSpacing -resetDistribution $::parameters(farfieldspacing)
$eastConnector setName east
$eastBC end
unset eastBC

# create south domain boundary
set southBC [pw::Application begin Create]
set southSegment [pw::SegmentSpline create]
$southSegment addPoint [list $::parameters(domainLength) -$::parameters(domainHeight) 0]
$southSegment addPoint [list 0 -$::parameters(domainHeight) 0]
set southConnector [pw::Connector create]
$southConnector addSegment $southSegment
unset southSegment
$southConnector setDimensionFromSpacing -resetDistribution $::parameters(farfieldspacing)
$southConnector setName south
$southBC end
unset southBC

# create west domain boundary
set westBC [pw::Application begin Create]
set westSegment [pw::SegmentCircle create]
$westSegment addPoint [list 0 -$::parameters(domainHeight) 0]
$westSegment addPoint [list 0 0 0]
$westSegment setEndAngle 180 {0 0 1}
$westSegment flipArc
set westConnector [pw::Connector create]
$westConnector addSegment $westSegment
$westConnector setDimensionFromSpacing -resetDistribution $::parameters(farfieldspacing)
$westConnector setName west
unset westSegment
$westBC end
unset westBC

# generate unstructured domain around airfoil
pw::Application setGridPreference Unstructured

# generate unstructured domain
set temp(mode1) [pw::Application begin Create]
set temp(edge1) [pw::Edge create]
$temp(edge1) addConnector $northConnector
$temp(edge1) addConnector $eastConnector
$temp(edge1) addConnector $southConnector
$temp(edge1) addConnector $westConnector
set temp(edge2) [pw::Edge create]
$temp(edge2) addConnector $lowerSurfCon
$temp(edge2) addConnector $trailSurfCon
$temp(edge2) addConnector $upperSurfCon
set domain [pw::DomainUnstructured create]
$domain addEdge $temp(edge1)
$domain addEdge $temp(edge2)
$domain setName airfoil-mesh
unset temp(edge2)
unset temp(edge1)
$temp(mode1) end
unset temp(mode1)

# create cylinder source for mesh refinement
set temp(mode1) [pw::Application begin Create]
set source [pw::SourceShape create]
$source cylinder -radius 0.5 -topRadius 3 -length [expr {1.25 * $::parameters(domainLength)}]
$source setTransform [list 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1]
$source setPivot Base
$source setSectionMinimum 0
$source setSectionMaximum 360
$source setSidesType Plane
$source setBaseType Plane
$source setTopType Plane
$source setEnclosingEntities {}
$temp(mode1) end
unset temp(mode1)

# rotate cylinder source into correct position
set temp(mode1) [pw::Application begin Modify [list $source]]
pw::Entity transform [pwu::Transform rotation -anchor {0 0 0} {0 1 0} 90] [$temp(mode1) getEntities]
pw::Entity transform [pwu::Transform rotation -anchor {0 0 0} {0 0 1} $::parameters(aoa)] [$temp(mode1) getEntities]
$temp(mode1) end
unset temp(mode1)

# set source distribution around cylinder
set temp(mode1) [pw::Application begin Modify [list $source]]
$source setBeginSpacing $::parameters(startspacing)
$source setEndSpacing $::parameters(endspacing)
$temp(mode1) end
unset temp(mode1)

# refine the mesh based on sources
set temp(mode1) [pw::Application begin UnstructuredSolver [list $domain]]
set temp(solve1) [pw::TRexCondition create]
$temp(solve1) setName farfield
unset temp(solve1)
set temp(solve1) [pw::TRexCondition getByName farfield]
$temp(solve1) apply [list [list $domain $northConnector Same] [list $domain $eastConnector Same] [list $domain $southConnector Same]]
$temp(solve1) setConditionType Match
$temp(solve1) setAdaptation On

# add T-REX for boundary layer inflation cells
set temp(solve2) [pw::TRexCondition create]
$temp(solve2) setName airfoil
unset temp(solve2)
set temp(solve2) [pw::TRexCondition getByName airfoil]
$temp(solve2) apply [list [list $domain $upperSurfCon Opposite] [list $domain $lowerSurfCon Same] [list $domain $trailSurfCon Opposite]]
$temp(solve2) setConditionType Wall
$temp(solve2) setValue $::parameters(initds)
$domain setUnstructuredSolverAttribute TRexMaximumLayers $::parameters(maxlayers)
$domain setUnstructuredSolverAttribute TRexFullLayers $::parameters(fulllayers)
$domain setUnstructuredSolverAttribute TRexGrowthRate $::parameters(cellgr)
$domain setUnstructuredSolverAttribute TRexPushAttributes True
$domain setUnstructuredSolverAttribute TRexCellType TriangleQuad
$domain setUnstructuredSolverAttribute TRexIsotropicHeight 0.5

# re-initialise the 2D domain based on mesh refinement and T-REX
$temp(mode1) run Initialize

$temp(mode1) end
unset temp(solve1)
unset temp(solve2)
unset temp(mode1)

# print mesh statistics
set domain [pw::GridEntity getByName "airfoil-mesh"]
set cells [$domain getCellCount]
puts "Generated domain with $cells cells"

# set boundary conditions
createBoundaryConditions

# Reset view
pw::Display resetView

}

# PROCEDURE TO GENERATE BOUNDARY CONDITIONS
# -----------------------------------------------
proc createBoundaryConditions {} {
    # solver attribute
    pw::Application setCAESolver OpenFOAM 2
    pw::Application setCAESolverAttribute SideBCExport Single
    pw::Application setCAESolverAttribute Thickness 1

    # get connectors
    set domain [pw::GridEntity getByName "airfoil-mesh"]
    set au [pw::GridEntity getByName "airfoilUpper"]
    set al [pw::GridEntity getByName "airfoilLower"]
    set te [pw::GridEntity getByName "airfoilTrailing"]
    set north [pw::GridEntity getByName "north"]
    set south [pw::GridEntity getByName "south"]
    set east [pw::GridEntity getByName "east"]
    set west [pw::GridEntity getByName "west"]

    set temp(BC1) [pw::BoundaryCondition create]
    $temp(BC1) setName airfoil
    $temp(BC1) setPhysicalType -usage CAE wall
    $temp(BC1) apply [list [list $domain $au] [list $domain $al] [list $domain $te]]

    set temp(BC2) [pw::BoundaryCondition create]
    $temp(BC2) setName freestream
    $temp(BC2) setPhysicalType -usage CAE patch
    $temp(BC2) apply [list [list $domain $north] [list $domain $south] [list $domain $east] [list $domain $west]]
}

# PROCEDURE OUTPUTING CASE AS OPENFOAM FILE
# -----------------------------------------------
proc writeMeshToDisk {} {
    set pathToWrite $::parameters(polyMeshDir)
    file mkdir $pathToWrite

    set domain [pw::GridEntity getByName "airfoil-mesh"]
    set temp [pw::Application begin CaeExport [pw::Entity sort [list $domain]]]
    $temp initialize -strict -type CAE $pathToWrite
    $temp verify
    $temp write
    $temp end
    unset temp
}

# -----------------------------------------------
# MAIN SCRIP
# -----------------------------------------------

# Load Pointwise Glyph package and Tk
package require PWI_Glyph 2.4
pw::Application reset

# ensure we know where the parameter file is located
if {$argc != 1} {
    puts "Parameter file needs to be specified via the command line, exiting ..."
    return
}
set parameterfile [lindex $argv 0]

# read the parameters from the parameter file
puts "Reading paremter file from: $parameterfile ..."
set file [open $parameterfile r]
while {![eof $file]} {
    set part [split [gets $file] "="]
    set parameters([string trimright [lindex $part 0]]) [string trimleft [lindex $part 1]]
}
close $file

# generating the airfoil geometry
puts "Generating NACA $::parameters(naca) airfoil ..."
airfoilGen

# generating the mesh around the airfoil
puts "Generating mesh for NACA $::parameters(naca) airfoil ..."
airfoilMesh

# outputting mesh to disk
puts "Writing mesh to disk ..."
writeMeshToDisk
puts "Finished with mesh generation!"