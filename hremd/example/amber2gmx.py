#!/usr/bin/python

import parmed as pmd

# convert AMBER topology to GROMACS format
amber = pmd.load_file('complex.parm7', 'complex.rst7')
# Save as GROMACS topology and GRO files
amber.save('complex.top')
amber.save('complex.gro')
