
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

proc randShiftCC1 { fnum min_dist max_dist } {
    set cc1segs [ atomselect top "chain C and resid 1 and name N" frame $fnum ]
    set mtb_atoms [atomselect top "chain A B" frame $fnum]
    set commtb [measure center $mtb_atoms]
    foreach cc1id [$cc1segs get segid] {
        
	set shiftR [expr { $min_dist + rand()*$max_dist }]
        puts "Shifting $cc1id by $shiftR"
        set ccmove [atomselect top "segname $cc1id" frame $fnum]
        set cogmove [measure center $ccmove]
        set vecmove [vecnorm [vecsub $cogmove $commtb]]
        $ccmove moveby [vecscale $shiftR $vecmove] 
    }
} 


set psf [lindex $argv 0]
set pdb [lindex $argv 1]
set shiftR [lindex $argv 2]

mol new $psf 
mol addfile $pdb
mol addfile $pdb
mol addfile $pdb
mol addfile $pdb
mol addfile $pdb
mol addfile $pdb
mol addfile $pdb

randShiftCC1 0 $shiftR 10.0
set m0 [atomselect top "all" frame 0]
$m0 writepdb ../3jalfull_gmoltype_RandConf_ShiftCC1.pdb
 
set m1 [atomselect top "all" frame 1]
randShiftCC1 1 $shiftR 10.0
set m2 [atomselect top "all" frame 2]
randShiftCC1 2 $shiftR 10.0
set m3 [atomselect top "all" frame 3]
randShiftCC1 3 $shiftR 10.0
set m4 [atomselect top "all" frame 4]
randShiftCC1 4 $shiftR 10.0
set m5 [atomselect top "all" frame 5]
randShiftCC1 5 $shiftR 10.0
set m6 [atomselect top "all" frame 6]
randShiftCC1 6 $shiftR 10.0

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
