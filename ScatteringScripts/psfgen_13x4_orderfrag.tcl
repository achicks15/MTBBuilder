
package require psfgen 

proc loop_PF_prot { pfmax pref chid attype rn } {
set pfnum 1 

while {$pfnum <= $pfmax} {
    set sname $attype$chid$pfnum
    puts $sname
    set dirloc "/lustre/or-scratch/cades-bsd/9cq/CC1/3jal_Microtubule/PF$pfnum"
    segment $sname { pdb "$dirloc/3jal_$pref$chid.pdb"
    first None
    last None
    }
    patch NONE $sname:1
    patch NONE $sname:$rn
    coordpdb "$dirloc/3jal_$pref$chid.pdb" $sname
    
    set pfnum [expr {$pfnum +1}]     
    }
}

proc loop_PF_lig { pfmax pref chid attype } {
set pfnum 1

while {$pfnum <= $pfmax} {
    set sname $attype$chid$pfnum
    puts $sname
    set dirloc "/lustre/or-scratch/cades-bsd/9cq/CC1/3jal_Microtubule/PF$pfnum"
    segment $sname { pdb "$dirloc/3jal_$pref$chid.pdb"
    first None
    last None
    }
    coordpdb "$dirloc/3jal_$pref$chid.pdb" $sname

    set pfnum [expr {$pfnum +1}]
    }
} 

proc align2PF { pfn chidtarg cc1loc chalign } { 

     ## 
     mol new /lustre/or-scratch/cades-bsd/9cq/CC1/3jal_Microtubule/PF${pfn}/3jal_PRO${chidtarg}.pdb 
     set targetid [molinfo top]
     mol new $cc1loc
     set moveid [molinfo top]
     set cc1_move [atomselect $moveid "chain Z"]
     set align_target [atomselect $targetid "name CA"]
     ## doesn't matter where the alignment goes as long as its consistent with the mapping I think  
     set align_move [atomselect $moveid "chain ${chalign} and name CA"]
     set M [measure fit $align_move $align_target]
     
     $cc1_move move $M
     $cc1_move writepdb ./cc1nt_PROZ.pdb 
     mol delete $targetid
     mol delete $moveid
}

proc loop_PF_CC1NT { pfmax pref chid attype cc1fname } {
set pfnum 1

while {$pfnum <= $pfmax} {
    set sname $attype$chid$pfnum
    puts $sname
    set dirloc "."
    ## alignment stuff 
    align2PF $pfnum $chid $cc1fname "Q"
    ## psfgen stuff
    pdbalias residue HIS HSD
    pdbalias atom ILE CD1 CD
    segment $sname { pdb "$dirloc/cc1nt_PROZ.pdb"
    first None
    last None
    }
    patch NONE $sname:1
    patch NONE $sname:135
    coordpdb "$dirloc/cc1nt_PROZ.pdb" $sname
      
    set pfnum [expr {$pfnum +1}]
    }

}

#psfgen_logfile psfprep/topology.log

topology /lustre/or-scratch/cades-bsd/9cq/CC1/3jal_Microtubule/charmm-gui-3jal/toppar/top_all36_prot.rtf
topology /lustre/or-scratch/cades-bsd/9cq/CC1/3jal_Microtubule/charmm-gui-3jal/toppar/top_all36_na.rtf
topology /lustre/or-scratch/cades-bsd/9cq/CC1/3jal_Microtubule/charmm-gui-3jal/toppar/top_all36_cgenff.rtf
topology /lustre/or-scratch/cades-bsd/9cq/CC1/3jal_Microtubule/charmm-gui-3jal/toppar/top_all36_carb.rtf
topology /lustre/or-scratch/cades-bsd/9cq/CC1/3jal_Microtubule/charmm-gui-3jal/toppar/toppar_water_ions.str
topology /lustre/or-scratch/cades-bsd/9cq/CC1/3jal_Microtubule/charmm-gui-3jal/toppar/toppar_all36_na_nad_ppi.str

#psfgen_logfile close

set cc1ntfile [lindex $argv 0]
#set restrtype [lindex $argv 1]
#set outdir [lindex $argv 1]

#set protchlist [list  "A" "B" "C" "D" "E" "F" "G" "H" ]
set gtpchlist [list "I" "J" ]
# "K" "L" "M" "N" "O" "P" ]
set mgchlist [list "Q" "R" ]
#"S" "T" "U" "V" "W" "X" ]
#set resnum [list "451" "445" "451" "445" "451" "445" "451" "445"] 

set atublist [list "A" "C" "E" "G" ]
set atubresn 451
set btublist [list "B" "D" "F" "H" ]
set btubresn 445


#set pfnum 1
#set pfmax 13
#set dirloc "../PF$pfnum"

#psfgen_logfile psfprep/structure_prep.log 

foreach pch $atublist {
    puts $pch
    loop_PF_prot 13 "PRO" $pch "AT" $atubresn   
}

foreach pch $btublist {
    puts $pch
    loop_PF_prot 13 "PRO" $pch "BT" $btubresn
}

foreach pch [list "B" "D" "F" "H"] {
     loop_PF_CC1NT 13 "PRO" $pch "C1" "${cc1ntfile}" 
}
foreach gtpch $gtpchlist {
    loop_PF_lig 13 "GTP" $gtpch "GTP" 
}

foreach mgch $mgchlist {
    loop_PF_lig 13 "MG" $mgch "MG"
}


guesscoord 
writepsf 3jalfull_13x1_gmoltype_1Xcc1nt.psf
writepdb 3jalfull_13x1_gmoltype_1Xcc1nt.pdb 

#psfgen_logfile close 

#package require topotools

#set topoloc "/lustre/or-scratch/cades-bsd/9cq/CC1/3jal_Microtubule/charmm-gui-3jal/toppar"
#set topolist [list ${topoloc}/par_all36m_prot.prm ${topoloc}/par_all36_na.prm ${topoloc}/par_all36_cgenff.prm ${topoloc}/toppar_water_ions.str ${topoloc}/toppar_all36_na_nad_ppi.str]

#mol new 3jalfull_13x1_gmoltype_1Xcc1nt.psf
#mol addfile 3jalfull_13x1_gmoltype_1Xcc1nt.pdb

#topo writegmxtop ./topol_3jal_13x1Full_molgroup_1Xcc1nt.top $topolist

quit 
