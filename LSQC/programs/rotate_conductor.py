#!/usr/bin/env python

"""
   ===============================================
   POILVERT Nicolas
   Department of Materials Science and Engineering
   Massachusetts Institute of Technology
   Cambridge, USA
   ===============================================

 :Program:    rotate_conductor.py
 :Purpose:    the purpose is to take an Hamiltonian matrix
              expressed in a certain ordered Wannier Function
              basis, and then rotate those Wannier Functions
              about the direction of conduction by a certain angle.
              Then the program recomputes the Hamiltonian matrix
              in this "rotated" basis.
 :Input(s):   a string corresponding to the name of the \*_tran_info.dat
              file
              a string corresponding to the name of the Hamiltonian
              matrix file
              a float corresponding to the angle of rotation
 :Output(s):  a file containing the "rotated" matrix in Wannier90
              format
              a file corresponding to the new \*_tran_info.dat
"""

# import needed standard modules

import os
import sys
import time
try:
  import numpy
except:
  print ""
  print " Error in rotate_conductor :"
  print " Needed package 'Numpy' could not be found."
  print " Exiting the program..."
  print ""
  sys.exit(1)
import optparse

# import needed modules from the "utilities" directory

sys.path.append(os.path.join(os.path.dirname(__file__),'../utilities'))
from modules import header, read_tran_info, read_htC,    \
barycenter, wannier_function, print_wf, sort, write_htC, \
write_tran_info, atom

# main() function

def main():

  # command-line options
  parser = optparse.OptionParser()
  parser.add_option('-f', '--file',     
                    dest="tran_info_file",
                    help="name of the *_tran_info.dat file")
  parser.add_option('-m', '--matrix',
                    dest="hamiltonian_matrix",
                    help="name of the Hamiltonian matrix (*_htC.dat)")
  parser.add_option('-a', '--angle',
                    dest="angle",
                    type="float",
                    default=0.0,
                    help="angle of rotation in radian")
  parser.add_option('-d', '--delta',
                    dest="delta",
                    type="float",
                    default=0.15,
                    help="distance threshold  for  WF grouping.\n"+
                         "use the same value as the one in the\n"+
                         "Wannier90 master input file")
  options, remainder = parser.parse_args()
  # if the program is called without options it prints the help menu
  if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

  # Printing some information with date and time
  print header
  print " program rotate_conductor.py"
  print " program started on %s at %s" %(time.strftime("%A, %d %B %Y",time.localtime()),
                                         time.strftime("%H:%M:%S",time.localtime()))
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

  # loading the Hamiltonian matrix
  hamiltonian_matrix = read_htC(options.hamiltonian_matrix)

  # calling the program's main routine
  t_start = time.time()
  main_routine(wf_total,
               wf_per_pl,
               cells_per_pl,
               cond_dir,
               second_dir,
               wannier_functions,
               atom_number,
               atomic_coordinates,
               hamiltonian_matrix,
               options.angle,
               options.delta)
  t_stop = time.time()

  # printing some timing informations
  print ""
  elapsed_time = t_stop - t_start
  if elapsed_time < 1.0:
    print " Execution time : %6.2f ms" %(elapsed_time*1000)
  else:
    print " Execution time : %6.2f s " %(elapsed_time)
  print ""

  return

def main_routine(wf_total,
                 wf_per_pl,
                 cells_per_pl,
                 cond_dir,
                 second_dir,
                 wannier_functions,
                 atom_number,
                 atomic_coordinates,
                 matrix,
                 angle,
                 delta):

  # first we separate the left-most principal layer, the conductor region
  # and the right-most principal layer
  PL1       = wannier_functions[:wf_per_pl]
  conductor = wannier_functions[wf_per_pl:-wf_per_pl]
  PL4       = wannier_functions[-wf_per_pl:]

  # now we compute the unit vector for the axis of rotation
  point_A = barycenter(PL1)
  point_B = barycenter(PL4)
  unit_vector = (1.0/numpy.linalg.norm(point_B-point_A))*(point_B-point_A)
  if cond_dir=='x':
    reference = numpy.array([1.0,0.0,0.0])
  elif cond_dir=='y':
    reference = numpy.array([0.0,1.0,0.0])
  else:
    reference = numpy.array([0.0,0.0,1.0])
  if numpy.dot(unit_vector,reference)<=0.95:
    print ""
    print " Warning in rotate_conductor :"
    print " The  computed  unit  vector for the axis of"
    print " rotation seems to be off from the reference"
    print " unit vector : %s" %(cond_dir)
    print " The dot product of unit vector with %s is : %3.2f" %(cond_dir,numpy.dot(unit_vector,reference))
    print ""
    user_answer = ""
    while user_answer not in ["yes", "no"]:
      user_answer = str(raw_input(" is this OK ? ('yes' or 'no') "))
    if user_answer=='no':
      print " Exiting the program..."
      print ""
      sys.exit(1)

  # we search the groups of WF with similar centres. To do this
  # a group is identified by the index of the first WF in the group
  # and the group size
  gp_dict = {}
  gp_bool = numpy.zeros(len(conductor),dtype='bool')
  for i in xrange(len(conductor)):
    if gp_bool[i]==False:
      gp_bool[i] = True
      gp_size    = 1
      # for each wf we search for similar centres just around that wf
      for j in xrange(max(i-10,0),min(i+11,len(conductor))):
        if gp_bool[j]==False:
          if same_center(conductor[i],conductor[j]):
            gp_bool[j] = True
            gp_size   += 1
      # if gp_size is at least 2 we add this group to the dictionnary
      if gp_size >= 2:
        gp_dict[conductor[i].index] = gp_size

  # now we need to remove the redondant WF from the conductor
  not_to_add = []
  for wf in conductor:
    if wf.index in gp_dict.keys():
      for j in xrange(1,gp_dict[wf.index]):
        not_to_add.append(wf.index+j)
  reduced_conductor = []
  for wf in conductor:
    if wf.index not in not_to_add:
      reduced_conductor.append(wf)

  # with the reduced_conductor list, we can rotate the centers
  reduced_rotated_conductor = rotate_wf(reduced_conductor,angle,unit_vector,point_A)

  # we then sort the rotated Wannier Functions
  sorted_reduced_conductor = sort(reduced_rotated_conductor,
                                  cond_dir,
                                  second_dir,
                                  delta)

  # we re-insert the WF with similar centers
  sorted_conductor = []
  for wf in sorted_reduced_conductor:
    if wf.index not in gp_dict.keys():
      sorted_conductor.append(wf)
    else:
      sorted_conductor.append(wf)
      for j in xrange(1,gp_dict[wf.index]):
        new_wf = wannier_function()
        pos_in_list = wf.index+j-conductor[0].index
        new_wf.setindex(conductor[pos_in_list].index)
        new_wf.setx(wf.x)
        new_wf.sety(wf.y)
        new_wf.setz(wf.z)
        new_wf.setonsite(conductor[pos_in_list].onsite)
        sorted_conductor.append(new_wf)

  # we now find the mapping between the original order of WF
  # and the final order after sorting
  final_order = []
  reference   = conductor[0].index
  for wf in sorted_conductor:
    final_order.append(wf.index-reference)

  # now we can reconstruct the matrix
  final_matrix = numpy.zeros((matrix.shape[0],matrix.shape[1]),dtype='float')
  for i in xrange(matrix.shape[0]):
    for j in xrange(matrix.shape[1]):
      final_matrix[i,j] = matrix[final_order[i],final_order[j]]
  write_htC('rotated_matrix_htC.dat',final_matrix)

  # at last we re-generate the *_tran_info.dat file
  final_wf_list = PL1 + sorted_conductor + PL4
  final_atomic_coordinates = rotate_atoms(atomic_coordinates,angle,unit_vector,point_A)
  write_tran_info(wf_total,
                  wf_per_pl,
                  cells_per_pl,
                  cond_dir,
                  second_dir,
                  final_wf_list,
                  atom_number,
                  final_atomic_coordinates)

  return

def same_center(wf1,wf2):
  """
   Takes 2 Wannier Function objects and gives back a boolean
   depending on whether the Wannier Function centres are almost
   identical (up to an uncertainty delta)
  """

  from math import sqrt

  # parameter that determines whether 2 Wannier Functions
  # have identical centres
  delta = 0.01

  # computing the Euclidian distance between the two Wannier
  # Function centers
  x1 = wf1.x
  x2 = wf2.x
  y1 = wf1.y
  y2 = wf2.y
  z1 = wf1.z
  z2 = wf2.z
  distance = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
  answer = False
  if distance <= delta:
    answer = True

  return answer

def rotate_wf(wf_list,theta,axis,ref_pt):
  """
   This function takes a Wannier Function list, and
   rotates the coordinates of the WF centres by an angle
   theta about the axis 'axis'.
   It returns a new list of Wannier Functions with the
   rotated coordinates.
  """

  import numpy
  from numpy import dot as dot
  from numpy import cross as cross
  from math import cos as cos
  from math import sin as sin

  wf_list_final = []
  for wf in wf_list:
    new_wf = wannier_function()

    # first we compute the relative vector
    v_rel = numpy.array([wf.x,wf.y,wf.z]) - ref_pt

    # then we rotate
    v_rot = dot(v_rel,axis)*axis+cos(theta)*(v_rel-dot(v_rel,axis)*axis)+ \
            sin(theta)*cross(axis,v_rel)

    # then we reintroduce the absolute coordinates
    v_final = v_rot + ref_pt

    # now we update the wf coordinates
    new_wf.setindex(wf.index)
    new_wf.setx(v_final[0])
    new_wf.sety(v_final[1])
    new_wf.setz(v_final[2])
    new_wf.setonsite(wf.onsite)

    # adding the new wf to the final list
    wf_list_final.append(new_wf)

  return wf_list_final

def rotate_atoms(atom_list,theta,axis,ref_pt):
  """
   This function takes a list of atom objects, and
   rotates the coordinates of those atoms by an angle
   theta about the axis 'axis'.
   It returns a new list of atom objects with the
   rotated coordinates.
  """

  import numpy
  from numpy import dot as dot
  from numpy import cross as cross
  from math import cos as cos
  from math import sin as sin

  new_atom_list = []
  for atome in atom_list:
    new_atom = atom()

    # first we compute the relative vector
    v_rel = numpy.array([atome.x,atome.y,atome.z]) - ref_pt

    # then we rotate
    v_rot = dot(v_rel,axis)*axis+cos(theta)*(v_rel-dot(v_rel,axis)*axis)+ \
            sin(theta)*cross(axis,v_rel)

    # then we reintroduce the absolute coordinates
    v_final = v_rot + ref_pt

    # now we update the atom coordinates
    new_atom.setsymbol(atome.symbol)
    new_atom.setx(v_final[0])
    new_atom.sety(v_final[1])
    new_atom.setz(v_final[2])

    # adding the new wf to the final list
    new_atom_list.append(new_atom)

  return new_atom_list

if __name__=="__main__":
  main()

