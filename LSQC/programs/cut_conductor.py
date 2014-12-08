#!/usr/bin/env python

"""
   ===============================================
   POILVERT Nicolas
   Department of Materials Science and Engineering
   Massachusetts Institute of Technology
   Cambridge, USA
   ===============================================

 :Program:    cut_conductor.py
 :Purpose:    this program takes an Hamiltonian matrix
              and "cuts" some lead pristine unit cells
              on both side of the conductor region
 :Input(s):   a string corresponding to the name of the
              *_tran_info.dat file
              a string corresponding to the name of the
              matrix file (usually ending with *_htC.dat)
              optionally a string corresponding to the
              name of the atomic position file
 :Output(s):  one or more directories with a hamiltonian
              matrix and possibly an atomic position file
"""

# import needed standard modules

import os
import sys
import time
try:
  import numpy
except:
  print ""
  print " Error in cut_conductor :"
  print " Needed package 'Numpy' could not be found."
  print " Exiting the program..."
  print ""
  sys.exit(1)
import optparse

# import needed modules from the "utilities" directory

sys.path.append(os.path.join(os.path.dirname(__file__),'../utilities'))
from modules import header, read_tran_info,   \
sort_positions, write_htC, translate_to_zero, \
write_xyz, read_htC

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
  parser.add_option('-p', '--positions',
                    dest="atomic_positions",
                    action="store_true",
                    help="whether of not the atomic positions are\n"+
                         "also printed")
  parser.add_option('-a', '--all',
                    dest="all_permutations",
                    action="store_true",
                    help="if activated, the program prints all the\n"+
                         "hamiltonian matrices resulting from  all\n"+
                         "possible cutting schemes")
  options, remainder = parser.parse_args()
  # if the program is called without options it prints the help menu
  if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

  # Printing some information with date and time
  print header
  print " program cut_conductor.py"
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
   atomic_coordinates_unsorted) = read_tran_info(options.tran_info_file)

  # loading the Hamiltonian matrix
  hamiltonian_matrix = read_htC(options.hamiltonian_matrix)

  # sort the atomic coordinates in the direction of conduction
  atomic_coordinates = sort_positions(atomic_coordinates_unsorted,cond_dir)

  # calling the program's main routine
  t_start = time.time()
  main_routine(wf_total,
               wf_per_pl,
               cells_per_pl,
               cond_dir,
               wannier_functions,
               atom_number,
               atomic_coordinates,
               hamiltonian_matrix,
               options.atomic_positions,
               options.all_permutations)
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
                 wannier_functions,
                 atom_number,
                 atomic_coordinates,
                 hamiltonian_matrix,
                 print_atomic_coordinates,
                 print_all_permutations):

  from math import floor
  import os

  # asking the user how many unit cells on the left are to be removed
  try:
    left_cells  = int(raw_input("\n Number of left  cells to remove in the conductor : "))
    right_cells = int(raw_input("\n Number of right cells to remove in the conductor : "))
    print ""
  except:
    print ""
    print " Error in cut_conductor :"
    print " Could not transform the number of cells into"
    print " an integer"
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # checking that some of the information in the *_tran_info.dat file
  # is actually compatible with some of the information about the matrix
  if wf_total-2*wf_per_pl!=hamiltonian_matrix.shape[0]:
    print ""
    print " Error in cut_conductor :"
    print " the number of Wannier Functions in the conductor"
    print " region is incompatible."
    print " In file *_tran_info.dat, this number is : %6i" %(wf_total-2*wf_per_pl)
    print " In the matrix file, this number is      : %6i" %(hamiltonian_matrix.shape[0])
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # checking whether the number of cells are not too big
  num_wann_per_cell = int(wf_per_pl/cells_per_pl)
  max_cell_number = int(floor(hamiltonian_matrix.shape[0]/(2.0*num_wann_per_cell)))
  user_continue = "dummy"
  if left_cells > max_cell_number or right_cells > max_cell_number:
    print ""
    print " Warning in cut_conductor :"
    print " The number of unit cells to remove seems too big."
    print " The program expects at most %i to be removed on " %(max_cell_number)
    print " one side (for systems with the defect 'in the middle')"
    print " Do you still want to continue ?"
    while user_continue not in ['yes','no']:
      user_continue = str(raw_input("\n -> type 'yes' or 'no' : "))
  if user_continue=='no':
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # asking the user the number of atoms per unit cell
  if print_atomic_coordinates:
    try:
      atoms_per_cell = int(raw_input(" How many atoms in the lead unit cell : "))
      print ""
    except:
      print ""
      print " Error in cut_conductor :"
      print " Could not transform the number of atoms per cell"
      print " into an integer"
      print " Exiting the program..."
      print ""
      sys.exit(1)

  # constructing the matrix and the atomic position files
  if not print_all_permutations:
    # creating a directory for storing the reduced matrix and
    # possibly the atomic coodinates
    os.mkdir("./%i_left_%i_right" %(left_cells,right_cells))
    # first we construct the reduced matrix
    final_matrix = hamiltonian_matrix[left_cells*num_wann_per_cell:hamiltonian_matrix.shape[0]-right_cells*num_wann_per_cell,
                                      left_cells*num_wann_per_cell:hamiltonian_matrix.shape[0]-right_cells*num_wann_per_cell]
    write_htC("./%i_left_%i_right/%i_left_%i_right_htC.dat" %(left_cells,
                                                              right_cells,
                                                              left_cells,
                                                              right_cells),final_matrix)
    # then optionnaly the atomic positions
    if print_atomic_coordinates:
      # computing the size of the reduced structure
      if cond_dir=='x':
        conductor_size = atomic_coordinates[-(cells_per_pl+right_cells)*atoms_per_cell].x - \
                         atomic_coordinates[(cells_per_pl+left_cells)*atoms_per_cell].x
      elif cond_dir=='y':
        conductor_size = atomic_coordinates[-(cells_per_pl+right_cells)*atoms_per_cell].y - \
                         atomic_coordinates[(cells_per_pl+left_cells)*atoms_per_cell].y
      else:
        conductor_size = atomic_coordinates[-(cells_per_pl+right_cells)*atoms_per_cell].z - \
                         atomic_coordinates[(cells_per_pl+left_cells)*atoms_per_cell].z
      # selecting the atoms in the reduced structure
      final_positions = atomic_coordinates[(cells_per_pl+left_cells)*atoms_per_cell:len(atomic_coordinates)-(cells_per_pl+right_cells)*atoms_per_cell]
      # translating the coordinates in the direction of conduction such
      # that the left-most atom is at 0
      final_positions = translate_to_zero(final_positions,cond_dir)
      # writing the positions to an xyz file
      write_xyz("./%i_left_%i_right/%i_left_%i_right.xyz" %(left_cells,
                                                            right_cells,
                                                            left_cells,
                                                            right_cells),final_positions,conductor_size)
  else:
    # here we do the same thing as above but for all possible cutting schemes
    for i in xrange(left_cells+1):
      for j in xrange(right_cells+1):
        # creating a directory for storing the reduced matrix and
        # possibly the atomic coodinates
        os.mkdir("./%i_left_%i_right" %(i,j))
        # first we construct the reduced matrix
        final_matrix = hamiltonian_matrix[i*num_wann_per_cell:hamiltonian_matrix.shape[0]-j*num_wann_per_cell,
                                          i*num_wann_per_cell:hamiltonian_matrix.shape[0]-j*num_wann_per_cell]
        write_htC("./%i_left_%i_right/%i_left_%i_right_htC.dat" %(i,j,i,j),final_matrix)
        # then optionnaly the atomic positions
        if print_atomic_coordinates:
          # computing the size of the reduced structure
          if cond_dir=='x':
            conductor_size = atomic_coordinates[-(cells_per_pl+j)*atoms_per_cell].x - \
                             atomic_coordinates[(cells_per_pl+i)*atoms_per_cell].x
          elif cond_dir=='y':
            conductor_size = atomic_coordinates[-(cells_per_pl+j)*atoms_per_cell].y - \
                             atomic_coordinates[(cells_per_pl+i)*atoms_per_cell].y
          else:
            conductor_size = atomic_coordinates[-(cells_per_pl+j)*atoms_per_cell].z - \
                             atomic_coordinates[(cells_per_pl+i)*atoms_per_cell].z
          # selecting the atoms in the reduced structure
          final_positions = atomic_coordinates[(cells_per_pl+i)*atoms_per_cell:len(atomic_coordinates)-(cells_per_pl+j)*atoms_per_cell]
          # translating the coordinates in the direction of conduction such
          # that the left-most atom is at 0
          final_positions = translate_to_zero(final_positions,cond_dir)
          # writing the positions to an xyz file
          write_xyz("./%i_left_%i_right/%i_left_%i_right.xyz" %(i,j,i,j),final_positions,conductor_size)

  # At last a simple Warning to the user for checking the atomic position files
  if print_atomic_coordinates: 
    print ""
    print " WARNING :"
    print " In constructing the atomic position files, this program uses a"
    print " very simple sorting algorithm. Please visualize the structures"
    print " yourself in order to double check them."

  return

if __name__=="__main__":
  main()

