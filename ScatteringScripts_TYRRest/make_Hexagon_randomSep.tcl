
proc randomHex { minr maxr } {
     set sep [expr {int($maxr - $minr)}] 
     set comdist [expr { $minr + rand()*$sep}]
     set hexheight [expr { $comdist*sqrt(3)/2}]
     set halfhex [expr {$comdist/2}] 
     return [list $comdist $hexheight $halfhex] 
}

proc getMoveVector { hexgeom position } {
    if { $position == 1 } {
        puts "Moving MT to position 1"
        set hexvec [list [expr {1.0*[lindex $hexgeom 2]}] [expr {1.0*[lindex $hexgeom 1]}] 0.0];
    } elseif { $position == 2 } {
        puts "Moving MT to position 2"
        set hexvec [list [lindex $hexgeom 0] 0.0 0.0];
    } elseif { $position == 3 } {
        puts "Moving MT to position 3"
        set hexvec [list [expr {1.0*[lindex $hexgeom 2]}] [expr {-1.0*[lindex $hexgeom 1]}] 0.0];
    } elseif { $position == 4 } {
        puts "Moving MT to position 4"
        set hexvec [list [expr {-1.0*[lindex $hexgeom 2]}] [expr {-1.0*[lindex $hexgeom 1]}] 0.0];
    } elseif { $position == 5 } {
        puts "Moving MT to position 5"
        set hexvec [list [expr {-1.0*[lindex $hexgeom 0]}] 0.0 0.0];
    } elseif { $position == 6 } {
        puts "Moving MT to position 1"
        set hexvec [list [expr {-1.0*[lindex $hexgeom 2]}] [expr {1.0*[lindex $hexgeom 1]}] 0.0];  
    } else  {
        puts "position should be including and between 1 to 6"
    }
    return $hexvec

}

set psf [lindex $argv 0]
set pdb [lindex $argv 1]

mol new $psf 
mol addfile $pdb
mol addfile $pdb
mol addfile $pdb
mol addfile $pdb
mol addfile $pdb
mol addfile $pdb
mol addfile $pdb

set m1 [atomselect top "all" frame 1]
set m2 [atomselect top "all" frame 2]
set m3 [atomselect top "all" frame 3]
set m4 [atomselect top "all" frame 4]
set m5 [atomselect top "all" frame 5]
set m6 [atomselect top "all" frame 6]

##move to position 1
set hex1 [randomHex 312 332.0]
set hexvector1 [getMoveVector $hex1 1]   
$m1 moveby $hexvector1

set hex1 [randomHex 312 332.0]
set hexvector2 [getMoveVector $hex1 2]
$m2 moveby $hexvector2

set hex1 [randomHex 312 332.0]
set hexvector3 [getMoveVector $hex1 3]
$m3 moveby $hexvector3

set hex1 [randomHex 312.0 332.0]
set hexvector4 [getMoveVector $hex1 4]
$m4 moveby $hexvector4

set hex1 [randomHex 312 332.0]
set hexvector5 [getMoveVector $hex1 5]
$m5 moveby $hexvector5

set hex1 [randomHex 312 332.0]
set hexvector6 [getMoveVector $hex1 6]
$m6 moveby $hexvector6

animate write pdb 3jalfull_Hexagon.pdb beg 0 end 6 waitfor all sel [atomselect top "all"] top
animate write dcd 3jalfull_Hexagon.dcd beg 0 end 6 waitfor all sel [atomselect top "all"] top
quit
