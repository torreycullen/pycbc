#!/usr/bin/env python

# Copyright (C) 2015 Ian W. Harry, Patricia Schmidt
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""
This program adds an additional column "f_ref" to an existing sim_inspiral xml table. This additional parameter is relevant for NR injections.
"""


import sys
#import numpy as np
from optparse import OptionParser
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import ilwd
from glue.ligolw import utils
from glue.ligolw.utils import process as ligolw_process
import lal
import numpy as np

class ContentHandler(ligolw.LIGOLWContentHandler):
    pass
lsctables.use_in(ContentHandler)

cols = lsctables.SimInspiralTable.validcolumns

def fill_missing_columns(sngl):
    for entry in cols.keys():
        if not(hasattr(sngl,entry)):
            if cols[entry] in ['real_4','real_8']:
                setattr(sngl,entry,0.)
            elif cols[entry] == 'int_4s':
                setattr(sngl,entry,0)
            elif cols[entry] == 'lstring':
                setattr(sngl,entry,'')
            elif entry == 'process_id':
                sngl.process_id = ilwd.ilwdchar("sim_inspiral:process_id:0")
            elif entry == 'simulation_id':
                sngl.event_id = ilwd.ilwdchar("sim_inspiral:simulation_id:0")
            else:
                print >> sys.stderr, "Column %s not recognized" %(entry)
                raise ValueError

parser = OptionParser(
    usage   = "%prog [OPTIONS]",
    description = "Add fref column to inspinj output" )

parser.add_option("-V", "--verbose", action="store_true", help="print extra debugging information", default=False )
parser.add_option("-i", "--input-file", action="store", type="string",  help="Input file with sngl_inspiral xml table")
parser.add_option("-o", "--output-file", action="store", type="string",  help="Output file name")
(opts,args) = parser.parse_args()

if not opts.input_file:
    print >> sys.stderr, "--input-file must be supplied"
    sys.exit(1)
if not opts.output_file:
    print >> sys.stderr, "--output-file must be supplied"
    sys.exit(1)

oldxml = utils.load_filename(opts.input_file, gz=opts.input_file.endswith('gz'),
                             contenthandler=ContentHandler)
oldSimTable = table.get_table(oldxml,"sim_inspiral")

#hdfxml = utils.load_filename(opts.catalog_file, gz=opts.catalog_file.endswith('gz'),
#contenthandler=ContentHandler)

#hdfSimTable = table.get_table(hdfxml, "sim_inspiral")
data_O1extreme_final = np.loadtxt('/home/torrey.cullen/lambdas/lambdas.txt')
fit = np.poly1d(np.polyfit(np.arange(1,2,.01),data_O1extreme_final, 10))

for key in cols.keys():
    if key not in oldSimTable.columnnames:
        oldSimTable.appendColumn(key)

for sim in oldSimTable:
    sim.lambda1 = fit(sim.mass1)
    sim.lambda2 = fit(sim.mass2)
    #print sim.mass1,sim.lambda1
    #print sim.mass2,sim.lambda2
utils.write_filename(oldxml, opts.output_file,
                     gz=opts.output_file.endswith('gz'))
