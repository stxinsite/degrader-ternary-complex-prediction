# --------------------------------------------------------
# User INPUTs (THESE ARE HREMD SIMULATIONS)
mol new 0ps.pdb 
mol addfile cluster.mol-center-lig.fit.xtc waitfor all 
# --------------------------------------------------------

# --------------------------------------------------------
# Topology and trajectory of CRL (these are CRL-VHL data)
mol new /bgfs01/insite02/asghar.razavi/ub_zone/vhl/relax_run_run2_align.pdb
mol addfile /bgfs01/insite02/asghar.razavi/ub_zone/vhlfinal_align_on_vhl_all_proper.trr waitfor all
# --------------------------------------------------------

# --------------------------------------------------------
# Aligning Inputs on CRL based on VHL

set v1 [atomselect 1 "chain J and name CA and resid 1474 to 1560 1598 to 1611"]

# The below line needs adjustment for VHL residues. Should match these residues: SER68 to PRO154 and PRO192 to ARG205
set v0 [atomselect 0 "chain B and name CA and resid 9 to 95 133 to 146"]


set v0num [${v0} num]
set v1num [${v1} num]
puts "${v0num}"
puts "${v1num}"
set all0 [atomselect 0 all]
mol top 0
set numframes [molinfo 0 get numframes]
puts "number of frames: ${numframes}"
for {set frame 0} {${frame} < ${numframes}} {incr frame} {
    puts "$frame"
    ${all0} frame ${frame}
    ${v0} frame ${frame}
    set tmat [measure fit ${v0} ${v1}]
    ${all0} move ${tmat}
}
# --------------------------------------------------------

# --------------------------------------------------------
# Writing coordinates of Lys residues
animate write trr acbi1.trr sel [atomselect 0 "noh and chain A and resname LYS"] beg 0 end -1 skip 1 0
animate write pdb acbi1.pdb sel [atomselect 0 "noh and chain A and resname LYS"] beg 0 end 0 skip 1 0

# Writing coordinates of C atom of ubiquitin C-terminus
animate write pdb ub.pdb sel [atomselect 1 "chain G and resid 1347 and name C and resname GLY"] beg 0 end -1 skip 1 1
# --------------------------------------------------------
