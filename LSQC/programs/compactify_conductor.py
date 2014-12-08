#!/usr/bin/env python

"""
   ===============================================
   POILVERT Nicolas
   Department of Materials Science and Engineering
   Massachusetts Institute of Technology
   Cambridge, USA
   ===============================================

 :Program:    compactify_conductor.py
 :Purpose:    the purpose of the program is to take
              a conductor hamiltonian matrix, expressed
              in a Wannier Function basis, and reconstruct
              it in such a way that the Wannier Functions
              of the defect are gathered in the 'center'
              of the (re-ordered) Wannier Function basis.
 :Input(s):   this program needs the *_tran_info.dat file
              and a Hamiltonian matrix *_htC.dat both given
              by the Wannier90 code
 :Output(s):  the program produces another Hamiltonian matrix
              *_htC_compactified.dat
              a new \*_tran_info.dat file
"""

# import needed standard modules

import sys
import time
try:
  import numpy
except:
  print ""
  print " Error in compactify_conductor :"
  print " Needed package 'Numpy' could not be found."
  print " Exiting the program..."
  print ""
  sys.exit(1)
import optparse

# import needed modules from the "utilities" directory

sys.path.append(os.path.join(os.path.dirname(__file__),'../utilities'))
from constraints import constraints
from modules import header, read_tran_info, read_htC, write_htC, \
write_tran_info

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
  options, remainder = parser.parse_args()
  # if the program is called without options it prints the help menu
  if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

  # Printing some information with date and time
  print header
  print " program compactify_conductor.py"
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
               hamiltonian_matrix)
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
                 hamiltonian_matrix):

  # making sure the user properly tuned the constraints function
  print ""
  print " Did you modify the 'constraints.py' module"
  print " to fit your system ?"
  constraints_modified = 'dummy'
  while constraints_modified not in ['yes','no']:
    constraints_modified = str(raw_input("\n -> type 'yes' or 'no' : "))
  if constraints_modified=='no':
    print ""
    sys.exit(1)
  print ""

  # checking that some of the information in the *_tran_info.dat file
  # is actually compatible with some of the information about the matrix
  if wf_total-2*wf_per_pl!=hamiltonian_matrix.shape[0]:
    print ""
    print " Error in compactify_conductor :"
    print " the number of Wannier Functions in the conductor"
    print " region is incompatible."
    print " In file *_tran_info.dat, this number is : %6i" %(wf_total-2*wf_per_pl)
    print " In the matrix file, this number is      : %6i" %(hamiltonian_matrix.shape[0])
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # creating two lists of indices. One for the backbone WF and one for
  # the defect WF. Those lists contain the indices of the WF
  backbone_list = []
  defect_list   = []
  for i in xrange(wf_per_pl,wf_total-wf_per_pl):  # this loop restricts to the WF in the
                                                  # conductor region
    if constraints(wannier_functions[i].x,wannier_functions[i].y,wannier_functions[i].z):
      backbone_list.append(i-wf_per_pl)
    else:
      defect_list.append(i-wf_per_pl)

  # now we will re-order the Wannier Function basis in such a way that
  # the defect Wannier Functions will be "in the middle"
  middle_point = len(backbone_list)/2
  final_list = backbone_list[:middle_point] + defect_list + backbone_list[middle_point:]

  # all we have to do now is re-construct the hamiltonian matrix with
  # this newly ordered basis
  reordered_hamiltonian_matrix = numpy.zeros((hamiltonian_matrix.shape[0],
                                              hamiltonian_matrix.shape[1]),dtype='float')
  for row in xrange(len(final_list)):
    for col in xrange(len(final_list)):
      reordered_hamiltonian_matrix[row,col] = hamiltonian_matrix[final_list[row],final_list[col]]

  # writing the new hamiltonian to a Wannier90 compatible file
  write_htC("compactified_htC.dat",reordered_hamiltonian_matrix)

  # generating the new *_tran_info.dat file
  PL1       = wannier_functions[:wf_per_pl]
  conductor = wannier_functions[wf_per_pl:-wf_per_pl]
  PL4       = wannier_functions[-wf_per_pl:]
  final_conductor = []
  for index in final_list:
    final_conductor.append(conductor[index])
  final_wf_list = PL1 + final_conductor + PL4
  
  write_tran_info(wf_total,
                  wf_per_pl,
                  cells_per_pl,
                  cond_dir,
                  second_dir,
                  final_wf_list,
                  atom_number,
                  atomic_coordinates)

  return

if __name__=="__main__":
  main()

