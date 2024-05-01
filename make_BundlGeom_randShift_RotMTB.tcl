
proc randomHex { minr maxr } {
     set sep [expr {int($maxr - $minr)}] 
     set comdist [expr { $minr + rand()*$sep}]
     set hexheight [expr { $comdist*sqrt(3)/2}]
     set halfhex [expr {$comdist/2}] 
     return [list $comdist $hexheight $halfhex] 
}

proc getMoveVector { hexgeom position } {
    
    if { $position == 0 } {
        puts "Not moving MT"
        set hexvec [list 0.0 0.0 0.0];
    } elseif { $position == 1 } {
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

proc randShiftCC1 { fnum max_dist } {
    set cc1segs [ atomselect top "chain C and resid 1 and name N" frame $fnum ]
    set mtb_atoms [atomselect top "chain A B" frame $fnum]
    set commtb [measure center $mtb_atoms]
    foreach cc1id [$cc1segs get segid] {
        
	set shiftR [expr { 60.0 + rand()*$max_dist }]
        puts "Shifting $cc1id by $shiftR"
        set ccmove [atomselect top "segname $cc1id" frame $fnum]
        set cogmove [measure center $ccmove]
        set vecmove [vecnorm [vecsub $cogmove $commtb]]
        $ccmove moveby [vecscale $shiftR $vecmove] 
    }
} 

proc randRotMTB { mtbsel } {

    set commtb [measure center $mtbsel weight mass]
    set rotangle [expr rand()*2*3.14159]
    set rotZmat [transaxis z $rotangle rad] 
    $mtbsel moveby [vecscale -1.0 $commtb]
    $mtbsel move $rotZmat
    $mtbsel moveby $commtb
}

set psf [lindex $argv 0]
set pdb [lindex $argv 1]
set comdsep [lindex $argv 2]
set geombndl [lindex $argv 3]

#set dRcom 10.0
set minRcom [expr $comdsep - 12.5]
set maxRcom [expr $comdsep + 7.0]

if { $geombndl == "Linear" } {
    set positionlist { 0 2 5 }      
} elseif { $geombndl == "Triangle" } {
    set positionlist { 0 1 2 }
} elseif { $geombndl == "Rectangle" } {
    set positionlist { 1 3 4 6 }
} elseif { $geombndl == "Diamond" } {
    set positionlist { 0 1 2 3 }
} elseif { $geombndl == "Trapezoid" } {
    set positionlist { 0 1 2 5 6 }
} elseif { $geombndl == "Hexagon" } {
    set positionlist { 0 1 2 3 4 5 6 }
} else {
    echo "invalid bundle geometry entered; defaulting to Hexagon"
    set geombndl "Hexagon" 
    set positionlist { 1 2 3 4 5 6 }
}

puts $geombndl
puts $positionlist
mol new $psf 
for { set ploc 0 } {$ploc < [llength $positionlist] } {incr ploc} {
    mol addfile $pdb
    puts $ploc
    set mtbatsel($ploc) [atomselect top "all" frame $ploc] 
    randRotMTB $mtbatsel($ploc)
    set hexd [randomHex $minRcom $maxRcom] 
    puts [lindex $positionlist $ploc]
    set hexvector [getMoveVector $hexd [lindex $positionlist $ploc]]
    puts $hexvector
    $mtbatsel($ploc) moveby $hexvector       
}

puts [llength $positionlist] 
molinfo top get numframes
animate write pdb 3jalfull_${geombndl}.pdb beg 0 waitfor all sel [atomselect top "all"] top
animate write dcd 3jalfull_${geombndl}.dcd beg 0 waitfor all sel [atomselect top "all"] top
quit
