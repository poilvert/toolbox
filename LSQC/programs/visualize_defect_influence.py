#!/usr/bin/env python

"""
   ===============================================
   POILVERT Nicolas
   Department of Materials Science and Engineering
   Massachusetts Institute of Technology
   Cambridge, USA
   ===============================================

 :Program:    visualize_defect_influence.py
 :Purpose:    plot the 'on-site' matrix element to
              see the influence of the defect onto
              the Wannier Functions
 :Input(s):   name of the *_tran_info.dat file
 :Output(s):  3D plot with Mayavi and optionally
              a file with all the informations
"""

# import needed standard modules

import os
import sys
import time
import optparse

# import needed modules from the "utilities" directory

sys.path.append(os.path.join(os.path.dirname(__file__),'../utilities'))
from constraints import constraints
from modules import header, read_tran_info

# main() function

def main():

  # command-line options
  parser = optparse.OptionParser()
  parser.add_option('-f', '--file',     
                    dest="tran_info_file",
                    help="name of the *_tran_info.dat file")
  parser.add_option('-p', '--print',
                    dest="output_deviations",
                    action="store_true",
                    help="prints the deviations in a file")
  options, remainder = parser.parse_args()
  # if the program is called without options it prints the help menu
  if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

  # Printing some information with date and time
  print header
  print " program visualize_defect_influence.py"
  print " program started on %s at %s" %(time.strftime("%A, %d %B %Y",time.localtime()),
                                         time.strftime("%H:%M:%S",time.localtime()))
  print ""

  # making sure the user properly tuned the constraints function
  print " Did you modify the 'constraints.py' module"
  print " to fit your system ?"
  constraints_modified = 'dummy'
  while constraints_modified not in ['yes','no']:
    constraints_modified = str(raw_input("\n -> type 'yes' or 'no' : "))
  if constraints_modified=='no':
    print ""
    sys.exit(1)
  print ""

  # loading the data in the *_tran_info.dat file
  (wf_total,
   wf_per_pl,
   cells_per_pl,
   cond_dir,
   second_dir,
   wannier_functions,
   atom_number,
   atomic_coordinates) = read_tran_info(options.tran_info_file)

  # calling the program's main routine
  main_routine(wf_total,
               wf_per_pl,
               cells_per_pl,
               cond_dir,
               wannier_functions,
               atom_number,
               atomic_coordinates,
               options.output_deviations)
  print ""

  return

def main_routine(wf_total,
                 wf_per_pl,
                 cells_per_pl,
                 cond_dir,
                 wannier_functions,
                 atom_number,
                 atomic_coordinates,
                 print_output):

  # making sure that a working version of mayavi is installed
  try:
    from enthought.mayavi import mlab
  except:
    print ""
    print " Error when trying to import the mlab module from mayavi"
    print " Please install mayavi2 to use this program"
    print " Exiting the program..."
    print ""
    sys.exit(1)
  # making sure that numpy is installed
  try:
    import numpy
  except:
    print ""
    print " Error when trying to import numpy"
    print " Please install numpy to use this program"
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # number of Wannier Functions per unit cell
  num_of_wf_in_cell = int(wf_per_pl/cells_per_pl)

  # removing the Wannier Functions belonging to the defect
  wf_index_final = []
  x_final = []
  y_final = []
  z_final = []
  mat_ele_final = []
  for i in xrange(len(wannier_functions)):
    if constraints(wannier_functions[i].x,wannier_functions[i].y,wannier_functions[i].z):
      wf_index_final.append(wannier_functions[i].index)
      x_final.append(wannier_functions[i].x)
      y_final.append(wannier_functions[i].y)
      z_final.append(wannier_functions[i].z)
      mat_ele_final.append(wannier_functions[i].onsite)

  # taking into account the direction of conduction
  if cond_dir=='x':
    first_coord  = y_final
    second_coord = z_final
  elif cond_dir=='y':
    first_coord  = z_final
    second_coord = x_final
  else:
    first_coord  = x_final
    second_coord = y_final

  # trying to find the matrix element deviations taking as a reference
  # the left-most unit cell
  deviation = []
  for i in xrange(num_of_wf_in_cell):              # left-most unit cell is the reference
    deviation.append(0.0)
  for i in xrange(num_of_wf_in_cell,len(x_final)): # rest of the Wannier Functions
    closest_wf = []
    closest_me = []
    for j in xrange(num_of_wf_in_cell):
      if distance(first_coord[j],second_coord[j],
                  first_coord[i],second_coord[i]) <= 0.5 and abs(mat_ele_final[i]-mat_ele_final[j]) <= 1.5:
        closest_wf.append(distance(first_coord[j],second_coord[j],first_coord[i],second_coord[i]))
        closest_me.append(mat_ele_final[i]-mat_ele_final[j])
    if len(closest_wf)==0:                         # if no closest WF, then deviation defaults to 0.0
      deviation.append(0.0)
    else:                                          # looking for the closest WF both in geometric distance
                                                   # and matrix element "distance"
      if numpy.argmin(closest_wf)==numpy.argmin(numpy.abs(closest_me)):
        deviation.append(closest_me[numpy.argmin(numpy.abs(closest_me))])
      else:                                        # if we don't find the  above WF, then we take the
                                                   # deviation corresponding to the closest WF in terms
                                                   # of geometric distance
        deviation.append(closest_me[numpy.argmin(closest_wf)])

  # if asked, prints the deviations in a file
  if print_output:
    out_file = open('deviations.dat','w')
    out_file.write(' # Index         x           y           z     deviation \n')
    for i in xrange(len(wf_index_final)):
      out_file.write('%6i  %10.6f  %10.6f  %10.6f    %10.3e \n' %(wf_index_final[i],
                                                                  x_final[i],
                                                                  y_final[i],
                                                                  z_final[i],
                                                                  deviation[i]))
    out_file.close()

  # now we will plot the Wannier Function centres with a ball the size of which
  # is depending on the magnitude of the deviation
  figure_hot  = mlab.points3d(numpy.array(x_final),
                              numpy.array(y_final),
                              numpy.array(z_final),
                              numpy.array(numpy.abs(deviation)),
                              colormap='Reds',
                              scale_factor=1.0,
                              opacity=0.70)
  x_positions = []
  y_positions = []
  z_positions = []
  for i in xrange(len(atomic_coordinates)):
    x_positions.append(atomic_coordinates[i].x)
    y_positions.append(atomic_coordinates[i].y)
    z_positions.append(atomic_coordinates[i].z)
  # atomic positions in black-white-grey
  figure_atoms = mlab.points3d(numpy.array(x_positions),
                               numpy.array(y_positions),
                               numpy.array(z_positions),
                               colormap='black-white',
                               scale_factor=0.5,
                               opacity=0.30)
  
  mlab.show()

  return

def distance(x1,y1,x2,y2):

  from math import sqrt
  return sqrt((x1-x2)**2+(y1-y2)**2)

if __name__=="__main__":
  main()

