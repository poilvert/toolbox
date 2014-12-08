#!/usr/bin/env python

"""
   ===============================================
   POILVERT Nicolas
   Department of Materials Science and Engineering
   Massachusetts Institute of Technology
   Cambridge, USA
   ===============================================

 :Program:    system_builder.py
 :Purpose:    the purpose of this program is to assemble
              many different defect Hamiltonian matrices
              into one large system. The detailed structure
              of that system is entirely controlled by the
              user.
 :Input(s):   a command line option '-c' or '-r' to know
              whether a custom made or a random system is to
              be build.
              a string to indicate the full path to the directory
              that contains the defects.
              the rest is given by the user interactively
 :Output(s):  a directory containing the Hamiltonian matrix
              of the user-defined system and optionnaly the
              atomic coordinates of that system as well.

.. Note::     
    the directories containing the defects should
    contain at most 2 files : one for the Hamiltonian
    matrix and one for the atomic positions
"""

# import needed standard modules

import os
import sys
import time
try:
  import numpy
except:
  print ""
  print " Error in system_builder :"
  print " Needed package 'Numpy' could not be found."
  print " Exiting the program..."
  print ""
  sys.exit(1)
import random
import optparse

# import needed modules from the "utilities" directory

sys.path.append(os.path.join(os.path.dirname(__file__),'../utilities'))
from modules import header, read_htC, read_xyz, \
write_htC, read_lead

# main() function

def main():

  # command-line options
  parser = optparse.OptionParser()
  parser.add_option('-c', '--custom',
                    dest="custom",
                    action="store_true",
                    help="if selected, the program builds a custom\n"+
                         "structure from the user input")
  parser.add_option('-r', '--random',
                    dest="random",
                    action="store_true",
                    help="if selected, the program builds a random\n"+
                         "structure from the user input")
  parser.add_option('-p', '--path',
                    dest="path",
                    default="./",
                    help="full path to the directory containing the\n"+
                         "defect(s). Default is : './'")
  options, remainder = parser.parse_args()
  # if the program is called without options it prints the help menu
  if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

  # some checks
  if options.custom and options.random:
    print ""
    print " Error in system_builder :"
    print " You selected both a custom and random structure!"
    print " You can only use one of these options."
    print " Exiting the program..."
    print ""
    sys.exit(1)
  if not options.custom and not options.random:
    print ""
    print " Error in system_builder :"
    print " You did not select any type of system to build!"
    print " Options are : '-c' for 'custom' system"
    print "               '-r' for 'random' system"
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # Printing some information with date and time
  print header
  print " program system_builder.py"
  print " program started on %s at %s" %(time.strftime("%A, %d %B %Y",time.localtime()),
                                         time.strftime("%H:%M:%S",time.localtime()))
  print ""

  # now extracting the defect directory names
  try:
    defect_list = os.listdir(options.path)
  except:
    print ""
    print " Error in system_builder :"
    print " Could not read the content of directory : %s" %(options.path)
    print " Is this a directory ?"
    print " Exiting the program..."
    print ""
    sys.exit(1)
  items_to_remove = []
  for item in defect_list:
    # impossible item ends or beginnings
    if item.startswith(".") or item.startswith("README") or item.endswith("~") or item.endswith(".txt"):
      items_to_remove.append(item)
    # removing the content of the large scale transport directory
    # from the list of potential defects
    elif item in ['Tran_info_file_format.txt',
                  'compactify_conductor',
                  'compactify_conductor.py',
                  'compactify_conductor.pyc',
                  'constraints',
                  'constraints.py',
                  'constraints.pyc',
                  'cut_conductor',
                  'cut_conductor.py',
                  'cut_conductor.pyc',
                  'modules',
                  'modules.py',
                  'modules.pyc',
                  'progressbar',
                  'progressbar.py',
                  'progressbar.pyc',
                  'system_builder',
                  'system_builder.py',
                  'system_builder.pyc',
                  'template',
                  'template.py',
                  'template.pyc',
                  'visualize_defect_influence',
                  'visualize_defect_influence.py',
                  'visualize_defect_influence.pyc',
                  'rotate_conductor',
                  'rotate_conductor.py',
                  'rotate_conductor.pyc']:
      items_to_remove.append(item)
    # if the item is not a directory, then one removes it from the list
    # of potential defects
    elif not os.path.isdir(options.path+"/"+item):
      items_to_remove.append(item)
    # also if the item is a directory but contains more than 2 files
    # (so 3 files or more), then one removes it from the list of potential defects
    elif os.path.isdir(options.path+"/"+item) and len(os.listdir(options.path+"/"+item))>2:
      items_to_remove.append(item)
  for not_defect in items_to_remove:
    defect_list.remove(not_defect)

  # checking that each file contains either one file ending in '_htC.dat'
  # or two files ending in '_htC.dat' and '.xyz'. Also determining whether
  # all the position files are present.
  all_positions_here = True
  defect_to_remove = []
  for defect in defect_list:
    # if there is no file, we drop the defect from the list
    if len(os.listdir(options.path+"/"+defect))==0:
      print ""
      print " Warning in system_builder :"
      print " no file detected in %s" %(defect)
      print " removing from the list of possible defects."
      print ""
      defect_to_remove.append(defect)
    # we're expecting only a hamiltonian matrix ending in '_htC.dat'
    elif len(os.listdir(options.path+"/"+defect))==1:
      filename = os.listdir(options.path+"/"+defect)[0]
      if not filename.endswith('_htC.dat'):
        print ""
        print " Warning in system_builder :"
        print " the file in %s does not end with '_htC.dat'" %(defect)
        print " removing from the list of possible defects."
        print ""
        defect_to_remove.append(defect)
      else:
        all_positions_here = False
    # we're expecting a hamiltonian matrix ending in '_htC.dat' and
    # a position file ending in '.xyz'
    elif len(os.listdir(options.path+"/"+defect))==2:
      files = os.listdir(options.path+"/"+defect)
      filename0 = files[0]
      filename1 = files[1]
      if not filename0.endswith('_htC.dat') and not filename1.endswith('_htC.dat'):
        print ""
        print " Warning in system_builder :"
        print " there seems to be a missing hamiltonian matrix"
        print " in %s" %(defect)
        print " removing from the list of possible defects."
        print ""
        defect_to_remove.append(defect)
      elif not filename0.endswith('.xyz') and not filename1.endswith('.xyz'):
        all_positions_here = False
  if len(defect_to_remove)>0:
    for defect in defect_to_remove:
      defect_list.remove(defect)

  # print the list of final defects to the screen
  num_defect_per_line = 2 # number of defect names to print on one line
  counter = 1
  defect_list_string = ""
  for defect in defect_list:
    if counter <= num_defect_per_line:
      defect_list_string += "%50s" %(defect)
      counter += 1
    else:
      counter = 2
      defect_list_string += "\n"
      defect_list_string += "%50s" %(defect)
  print " The program found the following defects :"
  print defect_list_string

  # asking the user which directory corresponds to the lead
  principal_layer = ""
  while principal_layer not in defect_list:
    principal_layer = str(raw_input("\n which defect represents the lead : "))
  defect_list.remove(principal_layer)

  # if all atomic position files are here, we ask the user whether
  # the atomic positions of the whole structure should be built
  print_xyz = False
  if all_positions_here:
    user_answer = ""
    while user_answer not in ['yes','no']:
      user_answer = str(raw_input("\n do you wish to print the atomic positions\n for the whole structure ('yes' or 'no') : "))
    if user_answer=='yes':
      print_xyz = True
    else:
      print_xyz = False

  # now it is time to build the list of defects making up the system
  system_defect_list = []
  if options.random:
    try:
      defect_number = int(raw_input("\n how many defects in the system : "))
    except:
      print ""
      print " Error in system_builder :"
      print " the program could not convert the number of defects"
      print " into an integer"
      print " Exiting the program..."
      print ""
      sys.exit(1)
    for i in xrange(defect_number):
      system_defect_list.append(random.choice(defect_list))
  elif options.custom:
    user_answer = ""
    print ""
    print " Enter the list of defects one by one (type 'stop' to finish the list)"
    while user_answer!="stop":
      user_answer = str(raw_input(" - defect name : "))
      if user_answer!="stop":
        if user_answer not in defect_list:
          print "   '%s' is not in the list of defects" %(user_answer,)
        else:
          system_defect_list.append(user_answer)

  # calling the routine to build the Hamiltonian matrix
  if options.random:
    system_type = 'random'
  else:
    system_type = 'custom'
  t_matrix_start = time.time()
  (num_cell_in_pl,buffer_num) = build_hamiltonian(system_defect_list,
                                                  principal_layer,
                                                  options.path,
                                                  system_type)
  t_matrix_stop  = time.time()

  # calling the routine to build the atomic positions
  if print_xyz:
    t_positions_start = time.time()
    build_structure(system_defect_list,
                    principal_layer,
                    num_cell_in_pl,
                    buffer_num,
                    options.path,
                    system_type)
    t_positions_stop  = time.time()

  # printing some timing informations
  print ""
  matrix_elapsed_time = t_matrix_stop - t_matrix_start
  if matrix_elapsed_time < 1.0:
    print " Time to build matrix    : %10.2f ms" %(matrix_elapsed_time*1000)
  else:
    print " Time to build matrix    : %10.2f s " %(matrix_elapsed_time)
  print ""
  if print_xyz:
    positions_elapsed_time = t_positions_stop - t_positions_start
    if positions_elapsed_time < 1.0:
      print " Time to build positions : %10.2f ms" %(positions_elapsed_time*1000)
    else:
      print " Time to build positions : %10.2f s " %(positions_elapsed_time)
    print ""

  return

def build_hamiltonian(defect_list,principal_layer,path,system_type):

  # finding out the list of distinct defects
  distinct_defects = []
  for defect in defect_list:
    if defect not in distinct_defects:
      distinct_defects.append(defect)

  # now we will load the matrices of the needed distinct defects
  # and gather them in a dictionnary. Also we will extract the matrix
  # sizes
  hamiltonian_matrices     = {}
  hamiltonian_matrix_sizes = {}
  for defect in distinct_defects:
    filename = ""
    for string in os.listdir(path+"/"+defect):
      if '_htC.dat' in string:
        filename = string
    conductor_matrix = read_htC(path+"/"+defect+"/"+filename)
    hamiltonian_matrices[defect]     = conductor_matrix
    hamiltonian_matrix_sizes[defect] = conductor_matrix.shape[0]

  # loading the H00 and H01 matrices of the lead
  filename = ""
  for string in os.listdir(path+"/"+principal_layer):
    if '_htC.dat' in string:
      filename = string
  (h00_matrix,h01_matrix) = read_lead(path+"/"+principal_layer+"/"+filename)
  h00_size = h00_matrix.shape[0]
  h01_size = h01_matrix.shape[0]
  hamiltonian_matrices[principal_layer]     = (h00_matrix,h01_matrix)
  hamiltonian_matrix_sizes[principal_layer] = (h00_size,h01_size)

  # determining the nbandc parameter
  H_sizes = []
  for defect in hamiltonian_matrix_sizes.keys():
    if defect!=principal_layer:
      H_sizes.append(hamiltonian_matrix_sizes[defect]+1)
  H_sizes.append(h00_size+h01_size)
  nbandc = max(H_sizes)

  # asking the user for the number of unit cells per principal layer
  try:
    num_cell_in_pl = int( raw_input('\n how many unit cells in one principal layer : ') )
  except:
    print ""
    print " Error in system_builder :"
    print " the program could not convert the number of unit"
    print " cells into an integer"
    print " Exiting the program..."
    print ""
    sys.exit(1)
  num_wf_in_cell = h00_size / num_cell_in_pl # a small check
  if num_cell_in_pl*num_wf_in_cell!=h00_size:
    print ""
    print " Error in system_builder :"
    print " The number of unit cells in a principal layer "
    print " multiplied by the number of Wannier Functions "
    print " per cell does not add up to the number of WFs "
    print " in  a principal layer. You probably entered a "
    print " wrong number of unit cells.                   "
    print " Exiting the program...                        "
    print ""
    sys.exit(1)

  # asking the user about the number of buffer unit cells to
  # "squeeze in" between the defects in the system
  try:
    buffer_num = int( raw_input('\n how many buffer unit cells between defects : ') )
  except:
    print ""
    print " Error in system_builder :"
    print " the program could not convert the number of buffer"
    print " unit cells into an integer"
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # size of the final Hamiltonian matrix
  total_size = 2*h00_size
  for defect in defect_list:
    if defect==principal_layer:
      total_size += h00_size
    else:
      total_size += hamiltonian_matrix_sizes[defect]
  total_size += len(defect_list)*num_wf_in_cell*buffer_num

  ### it is now time to build the final matrix ###
  from math import ceil

  # starting with a zero matrix with the right size
  dis_matrix = numpy.zeros((total_size,total_size),dtype='float')

  # updating the defect list by adding one principal layer at
  # the beginning of the structure and one at the end. This ensures
  # a clean connection with the leads.
  defect_list = [principal_layer]+defect_list+[principal_layer]

  # Main loop over the defect index
  offset = 0 # this parameter is used in the loop as a reference point
  for i in xrange(len(defect_list)):
    # Starting by inserting the defect matrices with their "left"
    # padding. Special attention is to be taken for the very first
    # and the very last defects that correspond to the beginning and
    # terminating lead principal layers
    if i==0 or i==len(defect_list)-1:
      dis_matrix[offset:offset+h00_size,offset:offset+h00_size] = hamiltonian_matrices[principal_layer][0]
      offset += h00_size
    else:
      # inserting the padding matrix first
      if defect_list[i]==principal_layer:
        def_size   = h00_size
        def_matrix = hamiltonian_matrices[principal_layer][0]
        num_cell   = int( buffer_num + num_cell_in_pl )
      else:
        def_size   = hamiltonian_matrix_sizes[defect_list[i]]
        def_matrix = hamiltonian_matrices[defect_list[i]]
        num_cell   = int( buffer_num + float(def_size)/(2.0*num_wf_in_cell) )
      if buffer_num >= 1:
        num_pl       = int( ceil(float(num_cell) / float(num_cell_in_pl)) )
        dummy_matrix = build_pristine( num_pl,
                                       hamiltonian_matrices[principal_layer][0],
                                       hamiltonian_matrices[principal_layer][1] )
        dis_matrix[offset:offset+num_cell*num_wf_in_cell,offset:offset+num_cell*num_wf_in_cell] = \
        dummy_matrix[:num_cell*num_wf_in_cell,:num_cell*num_wf_in_cell]
      # then inserting the defect matrix itself
      dis_matrix[offset+buffer_num*num_wf_in_cell:offset+buffer_num*num_wf_in_cell+def_size,
                 offset+buffer_num*num_wf_in_cell:offset+buffer_num*num_wf_in_cell+def_size] = def_matrix
      offset += buffer_num*num_wf_in_cell+def_size
    # Now we have to connect the defects together (off-diagonal elements)
    # We exclude the last defect since this one is not connected with any
    # other defect "to the right"
    if i <= len(defect_list)-2:
      if defect_list[i]==principal_layer:
        def_size_a   = h00_size
        def_matrix_a = hamiltonian_matrices[principal_layer][0]
      else:
        def_size_a   = hamiltonian_matrix_sizes[defect_list[i]]
        def_matrix_a = hamiltonian_matrices[defect_list[i]]
      if defect_list[i+1]==principal_layer:
        def_size_b   = h00_size
        def_matrix_b = hamiltonian_matrices[principal_layer][0]
      else:
        def_size_b   = hamiltonian_matrix_sizes[defect_list[i+1]]
        def_matrix_b = hamiltonian_matrices[defect_list[i+1]]
      # if both the current defect and the next one are principal layers, we 
      # insert simply H01 as a connection matrix
      if defect_list[i]==principal_layer and defect_list[i+1]==principal_layer:  
        dis_matrix[offset-h00_size:offset,offset:offset+h00_size] = hamiltonian_matrices[principal_layer][1]
        dis_matrix[offset:offset+h00_size,offset-h00_size:offset] = numpy.transpose(hamiltonian_matrices[principal_layer][1])
      # case where the current defect is a principal layer
      elif defect_list[i]==principal_layer:
        num_cell_a = num_cell_in_pl
        num_cell_b = int( buffer_num + float(def_size_b)/(2.0*num_wf_in_cell) )
        if num_cell_b==0:
          num_cell_b = int( buffer_num + ceil(float(def_size_b)/(2.0*num_wf_in_cell)) )
        num_cell_c = min(num_cell_a,num_cell_b)
        if num_cell_c >= num_cell_in_pl:
          dis_matrix[offset-h00_size:offset,offset:offset+h00_size] = hamiltonian_matrices[principal_layer][1]
          dis_matrix[offset:offset+h00_size,offset-h00_size:offset] = numpy.transpose(hamiltonian_matrices[principal_layer][1])
        else:
          dis_matrix[offset-num_cell_c*num_wf_in_cell:offset,offset:offset+num_cell_c*num_wf_in_cell] = \
                          hamiltonian_matrices[principal_layer][1][h00_size-num_cell_c*num_wf_in_cell:,:num_cell_c*num_wf_in_cell]
          dis_matrix[offset:offset+num_cell_c*num_wf_in_cell,offset-num_cell_c*num_wf_in_cell:offset] = \
                          numpy.transpose(hamiltonian_matrices[principal_layer][1][h00_size-num_cell_c*num_wf_in_cell:,:num_cell_c*num_wf_in_cell])
      # case where the next defect is a principal layer
      elif defect_list[i+1]==principal_layer:
        num_cell_b = num_cell_in_pl + buffer_num
        num_cell_a = int( float(def_size_a)/(2.0*num_wf_in_cell) )
        if num_cell_a==0:
          num_cell_a = int( ceil(float(def_size_a)/(2.0*num_wf_in_cell)) )
        num_cell_c = min(num_cell_a,num_cell_b)
        if num_cell_c >= num_cell_in_pl:
          dis_matrix[offset-h00_size:offset,offset:offset+h00_size] = hamiltonian_matrices[principal_layer][1]
          dis_matrix[offset:offset+h00_size,offset-h00_size:offset] = numpy.transpose(hamiltonian_matrices[principal_layer][1])
        else:
          dis_matrix[offset-num_cell_c*num_wf_in_cell:offset,offset:offset+num_cell_c*num_wf_in_cell] = \
                          hamiltonian_matrices[principal_layer][1][h00_size-num_cell_c*num_wf_in_cell:,:num_cell_c*num_wf_in_cell]
          dis_matrix[offset:offset+num_cell_c*num_wf_in_cell,offset-num_cell_c*num_wf_in_cell:offset] = \
                          numpy.transpose(hamiltonian_matrices[principal_layer][1][h00_size-num_cell_c*num_wf_in_cell:,:num_cell_c*num_wf_in_cell])
      # case where both the current and the next defect are not principal layers
      else:
        num_cell_a = int( float(def_size_a)/(2.0*num_wf_in_cell) )
        if num_cell_a==0:
          num_cell_a = int( ceil(float(def_size_a)/(2.0*num_wf_in_cell)) )
        num_cell_b = int( buffer_num + float(def_size_b)/(2.0*num_wf_in_cell) )
        if num_cell_b==0:
          num_cell_b = int( buffer_num + ceil(float(def_size_b)/(2.0*num_wf_in_cell)) )
        num_cell_c = min(num_cell_a,num_cell_b)
        if num_cell_c >= num_cell_in_pl:
          dis_matrix[offset-h00_size:offset,offset:offset+h00_size] = hamiltonian_matrices[principal_layer][1]
          dis_matrix[offset:offset+h00_size,offset-h00_size:offset] = numpy.transpose(hamiltonian_matrices[principal_layer][1])
        else:
          dis_matrix[offset-num_cell_c*num_wf_in_cell:offset,offset:offset+num_cell_c*num_wf_in_cell] = \
                          hamiltonian_matrices[principal_layer][1][h00_size-num_cell_c*num_wf_in_cell:,:num_cell_c*num_wf_in_cell]
          dis_matrix[offset:offset+num_cell_c*num_wf_in_cell,offset-num_cell_c*num_wf_in_cell:offset] = \
                          numpy.transpose(hamiltonian_matrices[principal_layer][1][h00_size-num_cell_c*num_wf_in_cell:,:num_cell_c*num_wf_in_cell])

  # writing the system's Hamiltonian matrix to file
  write_htC("%s_system_htC.dat" %("./"+system_type),dis_matrix)

  # creating the final directory
  os.system("rm -rf %s_system/" %(system_type))
  os.mkdir("%s_system" %(system_type))
  os.system("mv %s_system_htC.dat %s_system/." %(system_type,system_type))

  # adding the H_L, H_LC and H_CR matrices
  hl_file = open("./%s_system/%s_system_htL.dat" %(system_type,system_type),'w')
  hl_file.write(' # Left lead matrix with H00 and H01 \n')
  hl_file.write(' %6i \n' %(h00_size,))
  for j in xrange(h00_size):
    for i in xrange(h00_size):
      hl_file.write(' %10.6f' %(hamiltonian_matrices[principal_layer][0][i,j],))
    hl_file.write(' \n')
  hl_file.write(' %6i \n' %(h01_size,))
  for j in xrange(h01_size):
    for i in xrange(h01_size):
      hl_file.write(' %10.6f' %(hamiltonian_matrices[principal_layer][1][i,j],))
    hl_file.write(' \n')
  hl_file.close()
  hlc_file = open("./%s_system/%s_system_htLC.dat" %(system_type,system_type),'w')
  hlc_file.write(' # Left lead - Conductor matrix H_LC \n')
  hlc_file.write(' %6i  %6i \n' %(h01_size,h01_size))
  for j in xrange(h01_size):
    for i in xrange(h01_size):
      hlc_file.write(' %10.6f' %(hamiltonian_matrices[principal_layer][1][i,j],))
    hlc_file.write(' \n')
  hlc_file.close()
  hcr_file = open("./%s_system/%s_system_htCR.dat" %(system_type,system_type),'w')
  hcr_file.write(' # Conductor - Right lead matrix H_CR \n')
  hcr_file.write(' %6i  %6i \n' %(h01_size,h01_size))
  for j in xrange(h01_size):
    for i in xrange(h01_size):
      hcr_file.write(' %10.6f' %(hamiltonian_matrices[principal_layer][1][i,j],))
    hcr_file.write(' \n')
  hcr_file.close()

  # At last the Wannier90 master input file
  w90_input = open("./%s_system/%s_system.win" %(system_type,system_type),'w')
  w90_input.write("# %s system with %6i defect(s) \n" %(system_type,len(defect_list)-2,))
  w90_input.write("# written by 'system_builder.py' on %s at %s \n" %(time.strftime("%A, %d %B %Y",time.localtime()),
                                                                      time.strftime("%H:%M:%S",time.localtime())))
  w90_input.write("\n")
  w90_input.write("transport          = .true. \n")
  w90_input.write("transport_mode     =  lcr   \n")
  w90_input.write("tran_read_ht       = .true. \n")
  w90_input.write("tran_use_same_lead = .true. \n")
  w90_input.write("tran_win_min       = -3.0   \n")
  w90_input.write("tran_win_max       =  2.0   \n")
  w90_input.write("tran_energy_step   =  0.01  \n")
  w90_input.write("tran_num_ll        = %6i \n" %(h00_size,))
  w90_input.write("tran_num_rr        = %6i \n" %(h00_size,))
  w90_input.write("tran_num_cc        = %6i \n" %(total_size,))
  w90_input.write("tran_num_lc        = %6i \n" %(h00_size,))
  w90_input.write("tran_num_cr        = %6i \n" %(h00_size,))
  w90_input.write("tran_num_bandc     = %6i \n" %(nbandc,))
  w90_input.close()

  return (num_cell_in_pl,buffer_num)

def build_structure(defect_list,pl,num_cell_in_pl,buffer_num,path,system_type):

  # Asking the user the direction of conduction
  cond_dir = ''
  while cond_dir not in ['x','y','z']:
    cond_dir = str( raw_input("\n what is the direction of conduction ('x', 'y' or 'z') : ") )
    print ""

  # List of distinct defects
  distinct_defects = []
  for defect in defect_list:
    if defect not in distinct_defects:
      distinct_defects.append(defect)
  if pl not in distinct_defects:
    distinct_defects.append(pl)

  # Loading the atomic positions and fragment sizes
  positions = {}
  atom_number = {}
  length = {}
  for defect in distinct_defects:
    filename = ''
    for string in os.listdir(path+"/"+defect):
      if '.xyz' in string:
        filename = string
    temp_file = open(path+"/"+defect+"/"+filename,'r')
    lines = temp_file.readlines()
    atom_number[defect] = int(   lines[0].split()[0] )
    length[defect]      = float( lines[1].split()[0] )
    position_list = []
    for i in xrange(2,len(lines)):
      symbol  = lines[i].split()[0]
      x_coord = float( lines[i].split()[1] )
      y_coord = float( lines[i].split()[2] )
      z_coord = float( lines[i].split()[3] )
      position_list.append([symbol,x_coord,y_coord,z_coord])
    positions[defect] = position_list
    temp_file.close()

  # Finding out the number of atoms per unit cell
  num_atom_in_cell = atom_number[pl] / num_cell_in_pl
  if num_atom_in_cell*num_cell_in_pl!=atom_number[pl]:
    print ""
    print " Error in system_builder :"
    print " The total number of atoms in a principal layer "
    print " is not divisible  by the number  of unit cells "
    print " in a principal layer. Something is wrong !     "
    print " Dropping the construction of the position file."
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # Building the positions for the buffer only if
  # the buffer is at least one unit cell long
  if buffer_num >= 1:
    buffer           = []
    buffer_length    = 0.0
    buffer_num_atoms = 0
    unit_cell_size   = length[pl]/float(num_cell_in_pl)
    for i in xrange(buffer_num):
      for j in xrange(num_atom_in_cell):
        symbol  = positions[pl][j][0]
        x_coord = positions[pl][j][1]
        y_coord = positions[pl][j][2]
        z_coord = positions[pl][j][3]
        if cond_dir=='x':
          x_coord += buffer_length
        elif cond_dir=='y':
          y_coord += buffer_length
        else:
          z_coord += buffer_length
        buffer.append([symbol,x_coord,y_coord,z_coord])
      buffer_length    += length[pl]/float(num_cell_in_pl)
      buffer_num_atoms += num_atom_in_cell
    positions['buffer']   = buffer
    length['buffer']      = buffer_length
    atom_number['buffer'] = buffer_num_atoms

  # Building the structure
  final_structure = open('./%s_system/%s_system.xyz' %(system_type,system_type),'w')
  # computing the structure length
  structure_length  = 0.0
  structure_length += 2*length[pl]
  if buffer_num >= 1:
    structure_length += len(defect_list)*length['buffer']
  for defect in defect_list:
    structure_length += length[defect]
  # computing the total number of atoms in the structure
  structure_num_atoms  = 0
  structure_num_atoms += 2*atom_number[pl]
  if buffer_num >= 1:
    structure_num_atoms += len(defect_list)*atom_number['buffer']
  for defect in defect_list:
    structure_num_atoms += atom_number[defect]
  # inserting the above informations into the xyz file
  final_structure.write('%6i \n' %(structure_num_atoms,))
  final_structure.write('%10.6f # total length of the system in Angstroms \n' %(structure_length,))
  # First the left-most principal layer
  length_counter = 0.0
  for i in xrange(len(positions[pl])):
    symbol  = positions[pl][i][0]
    x_coord = positions[pl][i][1]
    y_coord = positions[pl][i][2]
    z_coord = positions[pl][i][3]
    if cond_dir=='x':
      x_coord += length_counter
    elif cond_dir=='y':
      y_coord += length_counter
    else:
      z_coord += length_counter
    final_structure.write('%s  %10.6f  %10.6f  %10.6f \n' %(symbol,x_coord,y_coord,z_coord))
  length_counter += length[pl]
  # Now the defects with their buffer
  for defect in defect_list:
    if buffer_num >= 1:
      for i in xrange(len(positions['buffer'])):
        symbol  = positions['buffer'][i][0]
        x_coord = positions['buffer'][i][1]
        y_coord = positions['buffer'][i][2]
        z_coord = positions['buffer'][i][3]
        if cond_dir=='x':
          x_coord += length_counter
        elif cond_dir=='y':
          y_coord += length_counter
        else:
          z_coord += length_counter
        final_structure.write('%s  %10.6f  %10.6f  %10.6f \n' %(symbol,x_coord,y_coord,z_coord))
        length_counter += length['buffer']
    for i in xrange(len(positions[defect])):
      symbol  = positions[defect][i][0]
      x_coord = positions[defect][i][1]
      y_coord = positions[defect][i][2]
      z_coord = positions[defect][i][3]
      if cond_dir=='x':
        x_coord += length_counter
      elif cond_dir=='y':
        y_coord += length_counter
      else:
        z_coord += length_counter
      final_structure.write('%s  %10.6f  %10.6f  %10.6f \n' %(symbol,x_coord,y_coord,z_coord))
    length_counter += length[defect]
  # finally the right-most principal layer
  for i in xrange(len(positions[pl])):
    symbol  = positions[pl][i][0]
    x_coord = positions[pl][i][1]
    y_coord = positions[pl][i][2]
    z_coord = positions[pl][i][3]
    if cond_dir=='x':
      x_coord += length_counter
    elif cond_dir=='y':
      y_coord += length_counter
    else:
      z_coord += length_counter
    final_structure.write('%s  %10.6f  %10.6f  %10.6f \n' %(symbol,x_coord,y_coord,z_coord))
  length_counter += length[pl]
  final_structure.close()

  return

def build_pristine(num_pl,H00,H01):

  import numpy

  dummy_array = numpy.zeros((num_pl*H00.shape[0],num_pl*H00.shape[0]),dtype='float')
  offset = 0
  for i in xrange(num_pl):
    dummy_array[offset:offset+H00.shape[0],offset:offset+H00.shape[0]] = H00
    if i <= num_pl-2:
      dummy_array[offset:offset+H00.shape[0],offset+H00.shape[0]:offset+2*H00.shape[0]] = H01
      dummy_array[offset+H00.shape[0]:offset+2*H00.shape[0],offset:offset+H00.shape[0]] = numpy.transpose(H01)
    offset += H00.shape[0]
  return dummy_array

if __name__=="__main__":
  main()

