#!/usr/bin/env python
# modules.py

"""
   **Module information** :

   :Program:   Module containing  all the classes, functions
               and parameters needed by the main programs
   :Purpose:   Prevent cluttering of the other programs
   :Input(s):  none
   :Output(s): none

"""

### Header bar to display in all program's output ###

header = """
          LARGE SCALE QUANTUM CONDUCTANCE PACKAGE
          Nicolas Poilvert
"""

### Default widgets for the progress bar ###

from progressbar import Percentage, Bar, ETA

widgets = [' Progress : ', Percentage(), ' ',
             Bar(marker="#",left="[",right="]"),' ',
             ETA()]

### Function designed to read a *_tran_info.dat file ###
### and load its content                             ###

import os
import re
import sys
from os import path
from progressbar import ProgressBar

def read_tran_info(file):
  """
   This function  is designed  to read a \*_tran_info.dat file
   and extract its content

   The format of the file is as follows:
     - first line contains the total number of Wannier Funtions
       in the system
     - second line contains the total number of WF in one principal layer
     - third line contains the number of lead unit cells in one
       principal layer
     - forth line contains the direction of conduction "x", "y"
       or "z" followed by the second direction of sorting
     - the next lines contain:
         WF index i, x\ :sub:`i`\ , y\ :sub:`i`\ , z\ :sub:`i`\ , <w\ :sub:`i`\ \| H \|w\ :sub:`i`\ >
     - then comes the total number of atoms in the system
     - at last the remaining lines contain:
         Atomic symbol j, x\ :sub:`j`\ , y\ :sub:`j`\ , z\ :sub:`j`\ 

   x\ :sub:`i`\  is the x coordinate of the i-th WF

   y\ :sub:`j`\  is the y coordinate of the j-th atom

   w\ :sub:`i`\  is the i-th WF

   <w\ :sub:`i`\ \| H \|w\ :sub:`i`\ >  is  the "on-site" matrix element associated
   with w\ :sub:`i`\ 

  :Input(s):   a string corresponding to the name of the \*_tran_info.dat file

  :Output(s):  an integer representing the total number of Wannier Functions
               in the system

               an integer representing the number of Wannier Functions in one
               principal layer

               an integer representing the number of unit cells in one principal
               layer

               a string giving the direction of conduction ("x", "y" or "z")

               a string giving the second direction of sorting ("x", "y" or "z")

               a list of Wannier Function objects with the proper attributes as
               defined in the *wannier_function* class

               an integer representing the total number of atoms in the system

               a list of atom objects with the proper attributes as defined in
               the *atom* class

  """

  # checking that the file actually exists
  if not path.isfile(file):
    print ""
    print " Error in read_tran_info :"
    print " %s is either not in the current working directory" %(file)
    print " or is not a file"
    print " Exiting the program..."
    print ""
    sys.exit(1)
 
  file_unit = open(file,'r')
  lines     = file_unit.readlines()
  file_unit.close()

  # creating a progress bar
  loop_max = len(lines)
  load_widgets = [' Reading %s : ' %(file), Percentage(), ' ',
                  Bar(marker="#",left="[",right="]"),' ',
                  ETA()]
  pgb = ProgressBar(widgets=load_widgets,maxval=loop_max,term_width=100).start()

  # regular expressions for integers and floats
  integer_re = r'\d+'
  float_re   = r'[+-]?[0-9]*\.[0-9]*[eE]?[+-]?[0-9]*'
  #float_re = r'[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
  atom_re    = r'(\w+)\s+('+float_re+r')\s+('+float_re+r')\s+('+float_re+r')'
  wf_re      = r'('+integer_re+r')\s+('+float_re+r')\s+('+float_re+r')\s+('+float_re+r')\s+('+float_re+r')'

  # extracting the total number of WF on the first line
  try:
    match = re.search(r'('+integer_re+r')',lines[0])
    if match:
      wf_total = int(match.group(1))
  except:
    print ""
    print " Error in read_tran_info :"
    print " Could not extract the total number of Wannier"
    print " Functions from file : %s" %(file)
    print " Exiting the program..."
    print ""
    sys.exit(1)
  pgb.update(1)

  # extracting the total number of WF in a principal layer
  try:
    match = re.search(r'('+integer_re+r')',lines[1])
    if match:
      wf_per_pl = int(match.group(1))
  except:
    print ""
    print " Error in read_tran_info :"
    print " Could not extract the total number of Wannier"
    print " Functions per principal layer from file : %s" %(file)
    print " Exiting the program..."
    print ""
    sys.exit(1)
  pgb.update(2)

  # extracting the number of lead unit cells in a principal layer
  try:
    match = re.search(r'('+integer_re+r')',lines[2])
    if match:
      cells_per_pl = int(match.group(1))
  except:
    print ""
    print " Error in read_tran_info :"
    print " Could not extract the total number of lead"
    print " unit cells per principal layer from file : %s" %(file)
    print " Exiting the program..."
    print ""
    sys.exit(1)
  pgb.update(3)

  # sanity check (wf_per_pl must be divisible by cells_per_pl)
  if wf_per_pl%cells_per_pl!=0:
    print ""
    print " Error in read_tran_info :"
    print " the number of Wannier Functions per principal layer"
    print " must be divisible by the number of cells per principal"
    print " layer."
    print " But here, the number of WF per PL is : %6i" %(wf_per_pl)
    print " and the number of cells per PL is    : %6i" %(cells_per_pl)
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # extracting the direction of conduction and second direction
  # of sorting
  try:
    match = re.search(r'(\w+)\s+(\w+)',lines[3])
    if match:
      cond_dir   = str(match.group(1))
      second_dir = str(match.group(2))
  except:
    print ""
    print " Error in read_tran_info :"
    print " Could not extract the direction of conduction and/or"
    print " the second direction of sorting"
    print " from file : %s" %(file)
    print " Exiting the program..."
    print ""
    sys.exit(1)
  if cond_dir not in ['x','y','z'] or second_dir not in ['x','y','z']:
    print ""
    print " Error in read_tran_info :"
    print " The direction of conduction is %s" %(cond_dir)
    print " The second direction of sorting is %s" %(second_dir)
    print " The expected directions are : 'x', 'y' or 'z'"
    print ""
    sys.exit(1)
  pgb.update(4)

  # extracting the Wannier Function information
  wannier_functions = []
  wf_iterator = 4
  try:
    for wf_counter in xrange(4,4+wf_total):
      match = re.search(wf_re,lines[wf_counter])
      if match:
        wf_to_add = wannier_function()
        wf_to_add.setindex(int(match.group(1)))
        wf_to_add.setx(float(match.group(2)))
        wf_to_add.sety(float(match.group(3)))
        wf_to_add.setz(float(match.group(4)))
        wf_to_add.setonsite(float(match.group(5)))
        wannier_functions.append(wf_to_add)
      pgb.update(wf_counter+1)
      wf_iterator += 1
  except:
    print ""
    print " Error in read_tran_info :"
    print " Could not extract Wannier Function information"
    print " from file   : %s" %(file)
    print " line number : %i" %(wf_iterator+1)
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # extracting the total number of atoms
  try:
    match = re.search(r'('+integer_re+r')',lines[4+wf_total])
    if match:
      atom_number = int(match.group(1))
  except:
    print ""
    print " Error in read_tran_info :"
    print " Could not extract the total number of atoms"
    print " from file : %s" %(file)
    print " Exiting the program..."
    print ""
    sys.exit(1)
  pgb.update(5+wf_total)

  # extracting the atomic coordinates
  atomic_coordinates = []
  at_iterator = 5+wf_total
  try:
    for at_counter in xrange(5+wf_total,len(lines)):
      match = re.search(atom_re,lines[at_counter])
      if match:
        atom_to_add = atom()
        atom_to_add.setsymbol(str(match.group(1)))
        atom_to_add.setx(float(match.group(2)))
        atom_to_add.sety(float(match.group(3)))
        atom_to_add.setz(float(match.group(4)))
        atomic_coordinates.append(atom_to_add)
      pgb.update(at_counter+1)
      at_iterator += 1
  except:
    print ""
    print " Error in read_tran_info :"
    print " Could not extract the atomic coordinates"
    print " from file   : %s" %(file)
    print " line number : %i" %(at_iterator+1)
    print " Exiting the program..."
    print ""
    sys.exit(1)

  pgb.finish()

  return (wf_total,
          wf_per_pl,
          cells_per_pl,
          cond_dir,
          second_dir,
          wannier_functions,
          atom_number,
          atomic_coordinates)

### function writing a *_tran_info.dat file in proper format ###

def write_tran_info(wf_num,
                    wf_num_per_pl,
                    cell_per_pl,
                    cond_dir,
                    second_dir,
                    wf_list,
                    atom_num,
                    atom_list):
  """
   This function is designed to write a \*_tran_info.dat file from a list of data

   The newly created file is called `new_tran_info.dat`

   :Input(s):  an integer for the total number of Wannier Functions in the system

               an integer for the number of Wannier Functions in one principal layer

               an integer for the number of unit cells in one principal layer

               a string giving the direction of conduction (which is also the first
               direction of sorting)

               a string giving the second direction of sorting

               a list of Wannier Function objects

               an integer giving the number of atoms in the system

               a list of atom objects

   :Output(s): a file in a format similar to the \*_tran_info.dat files given by Wannier90

  """

  # creating the ouput file
  outfile = open('./new_tran_info.dat','w')
  outfile.write('%6i\n' %(wf_num,))
  outfile.write('%6i\n' %(wf_num_per_pl,))
  outfile.write('%6i\n' %(cell_per_pl,))
  outfile.write('     %s %s\n' %(cond_dir,second_dir))
  iterator = 1
  # here we forget about the original Wannier Function indices because the order
  # given by wf_list is going to be the new order of the Wannier Function basis
  # so we simply set the Wannier Function indices to be an increasing sequence of
  # integers starting at 1
  for wf in wf_list:
    outfile.write('%6i   %14.8f   %14.8f   %14.8f   %14.8f   \n' %(iterator,wf.x,wf.y,wf.z,wf.onsite))
    iterator += 1
  outfile.write('%6i\n' %(atom_num,))
  for atom in atom_list:
    outfile.write('%5s   %14.8f   %14.8f   %14.8f   \n' %(atom.symbol,atom.x,atom.y,atom.z))
  outfile.close()

  return

### function designed to read a Hamiltonian matrix ###
### as given by the Wannier90 code                 ###

def read_htC(matrix_file):
  """
   This function reads a (conductor) Hamiltonian matrix given by the
   Wannier90 code.

   The name of the file generally ends with `*_htC.dat`

   The format of this file is as follows :

    - the first line contains comments
    - the second line contains the size N of the (square) N*N matrix
    - all the remaining lines are filled with numbers representing the
      matrix elements.

   .. Note:: the matrix elements are printed in Fortran (column major)
             format. This means that if the matrix has elements H(i,j) for i and
             j in [1,N], then the matrix elements are given in the following
             order ::
                H(1,1)  H(2,1)  H(3,1)  H(4,1) ...
                H(j,1)  H(j+1,1) ...

   :Input(s):   a string corresponding to the name of the file
                containing the Hamiltonian matrix

   :Output(s):  a rank-2 numpy array containing the Hamiltonian matrix

  """

  import numpy

  # making sure that the matrix file is a file and is in the current
  # working directory
  if not path.isfile(matrix_file):
    print ""
    print " Error in read_htC :"
    print " %s is either not in the current working directory" %(file)
    print " or is not a file"
    print " Exiting the program..."
    print ""
    sys.exit(1)  

  # opening the matrix file and loading the matrix
  ham_file = open(matrix_file,'r')
  lines = ham_file.readlines()
  ham_file.close()

  # creating a progress bar
  load_widgets = [' Loading matrix : ', Percentage(), ' ',
                  Bar(marker="#",left="[",right="]"),' ',
                  ETA()]
  pgb = ProgressBar(widgets=load_widgets,maxval=len(lines),term_width=100).start()

  # now extracting the data
  try:
    matrix_size = int(lines[1].split()[0])
  except:
    print ""
    print " Error in read_htC :"
    print " Could not extract the matrix size from file %s" %(matrix_file)
    print " Was the matrix file obtained from Wannier90 ?"
    print " Exiting the program..."
    print ""
    sys.exit(1)
  pgb.update(2)

  hamiltonian = numpy.zeros((matrix_size,matrix_size),dtype='float')
  iterator = 0
  line_number = 2
  try:
    for i in xrange(2,len(lines)):
      line = lines[i].split()
      for item in line:
        col = iterator/matrix_size
        row = iterator%matrix_size
        hamiltonian[row,col] = float(item)
        iterator += 1
      line_number += 1
      pgb.update(i+1)
  except:
    print ""
    print " Error in read_htC :"
    print " Could not properly extract data on line : %i" %(line_number)
    print " in file %s" %(matrix_file)
    print " Exiting the program..."
    print ""
    sys.exit(1)

  pgb.finish()

  return hamiltonian

def write_htC(file_name,matrix):
  """
   This function will take a Hamiltonian matrix and write it into
   a file in a format usable by Wannier90.

   In particular, the matrix will be written to file in column-major
   format (see `read_htC` for more detail).

   :Input(s):   a string corresponding to the name of the file that
                will contain the Hamiltonian matrix.

                a rank-2 numpy array containing the matrix

   :Output(s):  a file containing the matrix

  """

  import time
  import numpy

  # creating the output file, making sure no other file with this
  # name exists
  answer = "dummy"
  if path.isfile(file_name):
    print ""
    print " Warning : %s already exists" %(file_name)
    print " Would you like to overwrite it ?"
    print ""
    while answer not in ["yes","no"]:
      answer = str(raw_input(" -> type 'yes' or 'no' : "))
    if answer=="no":
      return
    else:
      os.remove(file_name)
  output = open(file_name,'w')

  # creating a progress bar
  write_widgets = [' Writing matrix to file : ', Percentage(), ' ',
                   Bar(marker="#",left="[",right="]"),' ',
                   ETA()]
  pgb = ProgressBar(widgets=write_widgets,maxval=matrix.shape[0],term_width=100).start()

  output.write(" # Written by 'write_htC' on %s at %s \n" %(time.strftime("%A, %d %B %Y",time.localtime()),
                                                            time.strftime("%H:%M:%S",time.localtime())))
  output.write(" %6i \n" %(matrix.shape[0]))
  for col in xrange(matrix.shape[1]):
    for row in xrange(matrix.shape[0]):
      output.write(" %10.6f" %(matrix[row,col]))
    output.write(" \n")
    pgb.update(col+1)
  output.close()

  pgb.finish()

  return

### The sort_positions function ###

def sort_positions(atomic_coordinates,cond_dir):
  """
   This function takes a list of atom objects and sorts
   the atomic positions only in the direction of conduction

   :Input(s):  the list of `atom` objects

               a string which represents the direction of conduction
               (which serves as the main direction of sorting)

   :Output(s): a list of `atom` objects in sorted order

  """

  import numpy

  # extracting the coordinates in the direction of conduction
  first_coordinate = []
  if cond_dir=='x':
    for atome in atomic_coordinates:
      first_coordinate.append(atome.x)
  elif cond_dir=='y':
    for atome in atomic_coordinates:
      first_coordinate.append(atome.y)
  else:
    for atome in atomic_coordinates:
      first_coordinate.append(atome.z)
  first_coordinate = numpy.array(first_coordinate)

  # sorting in the direction of conduction
  first_sort = first_coordinate.argsort()

  # returning the list of atoms in proper sorted order
  final_sorted_positions = []
  for index in first_sort:
    final_sorted_positions.append(atomic_coordinates[index])

  return final_sorted_positions

### The "translate_to_zero" function ###

def translate_to_zero(atomic_coordinates,cond_dir):
  """
   This function takes a list of `atom` objects and
   translate the coordinates in the direction of conduction
   in such a way that the left-most atom is set to 0

   :Input(s):  the list of `atom` objects

               a string which represents the direction of conduction

   :Output(s): a list of `atom` objects with the properly
               translated coordinates

  """

  import numpy

  # getting all the coordinates in the direction of conduction
  first_coordinate = []
  if cond_dir=='x':
    for atome in atomic_coordinates:
      first_coordinate.append(atome.x)
  elif cond_dir=='y':
    for atome in atomic_coordinates:
      first_coordinate.append(atome.y)
  else:
    for atome in atomic_coordinates:
      first_coordinate.append(atome.z)
  first_coordinate = numpy.array(first_coordinate)

  # now finding the index of the atom with the lowest coordinate
  min_index = first_coordinate.argmin()

  # getting the lowest coordinate
  min_coordinate = first_coordinate[min_index]

  # translating all the coordinates
  if cond_dir=='x':
    for i in xrange(len(atomic_coordinates)):
      new_coord = atomic_coordinates[i].x - min_coordinate
      atomic_coordinates[i].setx(new_coord)
  elif cond_dir=='y':
    for i in xrange(len(atomic_coordinates)):
      new_coord = atomic_coordinates[i].y - min_coordinate
      atomic_coordinates[i].sety(new_coord)
  else:
    for i in xrange(len(atomic_coordinates)):
      new_coord = atomic_coordinates[i].z - min_coordinate
      atomic_coordinates[i].setz(new_coord)

  return atomic_coordinates

### The "read_lead" function ###

def read_lead(filename):
  """
   This function extracts the H00 and H01 matrices
   from a `*_htL.dat` or `*_htR.dat` Wannier90 formatted
   file.

   The H00 and H01 matrices correspond respectively
   to the **on-site** and **off-diagonal** matrix sub-blocks
   of a lead principal layer.

   :Input(s):  a string corresponding to the name (full
               path) of the file containing the principal
               layer matrices

   :Output(s): a tuple with two rank-2 numpy arrays in the
               following order (H00,H01)

  """

  import numpy

  temp_file = open(filename,'r')
  lines = temp_file.readlines()
  temp_file.close()

  starting_line_for_H00 = 1
  starting_line_for_H01 = int( (len(lines)+1)/2 )

  # extracting the size of H00
  try:
    H00_size = int(lines[starting_line_for_H00].split()[0])
  except:
    print ""
    print " Error in read_lead :"
    print " could not extract the size of H00 on line 1"
    print " in file %s" %(filename)
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # extracting the size of H01
  try:
    H01_size = int(lines[starting_line_for_H01].split()[0])
  except:
    print ""
    print " Error in read_lead :"
    print " could not extract the size of H01 on line %i" %(starting_line_for_H01)
    print " in file %s" %(filename)
    print " Exiting the program..."
    print ""
    sys.exit(1)

  H00 = numpy.zeros((H00_size,H00_size),dtype='float')
  H01 = numpy.zeros((H01_size,H01_size),dtype='float')

  # extracting the H00 matrix
  iterator = 0
  line_number = starting_line_for_H00+1
  try:
    for i in xrange(starting_line_for_H00+1,starting_line_for_H01):
      line = lines[i].split()
      for item in line:
        col = iterator/H00_size
        row = iterator%H00_size
        H00[row,col] = float(item)
        iterator += 1
      line_number += 1
  except:
    print ""
    print " Error in read_lead :"
    print " could not extract matrix H00 on line %i" %(line_number)
    print " in file %s" %(filename)
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # extracting the H01 matrix
  iterator = 0
  line_number = starting_line_for_H01+1
  try:
    for i in xrange(starting_line_for_H01+1,len(lines)):
      line = lines[i].split()
      for item in line:
        col = iterator/H01_size
        row = iterator%H01_size
        H01[row,col] = float(item)
        iterator += 1
      line_number += 1
  except:
    print ""
    print " Error in read_lead :"
    print " could not extract matrix H01 on line %i" %(line_number)
    print " in file %s" %(filename)
    print " Exiting the program..."
    print ""
    sys.exit(1)

  return (H00,H01)

### The "read_xyz" function ###

def read_xyz(filename):
  """
   This function reads an XYZ formatted file and returns
   a list of `atom` objects. If it can find it, the size
   of the system is also returned.

   :Input(s):  a string for the name of the xyz file

   :Output(s): a list of `atom` objects
               a float which represents the size of the system

  """

  # checking that the file exists
  if not path.isfile(filename):
    print ""
    print " Error in read_xyz :"
    print " %s is not a file or does not exist" %(filename)
    print " Exiting the program..."
    print ""
    sys.exit(1)

  file_unit = open(filename,'r')
  lines = file_unit.readlines()
  file_unit.close()

  # extracting the size of the system
  try:
    size = float(lines[1].split()[0])
  except:
    size = 0.0

  # now extracting the atomic coordinates
  atomic_coordinates = []
  line_number = 2
  try:
    for i in xrange(2,len(lines)):
      items = lines[i].split()
      if len(items)>=4:
        atom_to_add = atom()
        atom_to_add.setsymbol(str(items[0]))
        atom_to_add.setx(float(items[1]))
        atom_to_add.sety(float(items[2]))
        atom_to_add.setz(float(items[3]))
        atomic_coordinates.append(atom_to_add)
      line_number += 1
  except:
    print ""
    print " Error in read_xyz :"
    print " could not extract informations on line %i" %(line_number)
    print " in file %s" %(filename)
    print " Exiting the program..."
    print ""
    sys.exit(1)    

  return (atomic_coordinates,size)

### The "write_xyz" function ###

def write_xyz(filename,atomic_coordinates,size):
  """
   This function takes a list of `atom` objects and prints
   their cartesian coordinates in a file with xyz format

   :Input(s):  a string for the name of the xyz file

               the list of `atom` objects

               a float, the size of the system in Angstrom(s)

   :Output(s): a file with name `filename` containing the
               atomic positions in xyz format

  """

  import time

  # creating the output file, making sure no other file with this
  # name exists
  answer = "dummy"
  if path.isfile(filename):
    print ""
    print " Warning : %s already exists" %(file_name)
    print " Would you like to overwrite it ?"
    print ""
    while answer not in ["yes","no"]:
      answer = str(raw_input(" -> type 'yes' or 'no' : "))
    if answer=="no":
      return
    else:
      os.remove(filename)
  output = open(filename,'w')
  output.write('%6i \n' %(len(atomic_coordinates)))
  output.write('%10.5f ' %(size)) # the size will be used in 'system_builder.py'
  output.write(" # Written by 'write_xyz' on %s at %s \n" %(time.strftime("%A, %d %B %Y",time.localtime()),
                                                            time.strftime("%H:%M:%S",time.localtime())))
  for i in xrange(len(atomic_coordinates)):
    output.write(" %s  %10.6f  %10.6f  %10.6f \n" %(atomic_coordinates[i].symbol,
                                                    atomic_coordinates[i].x,
                                                    atomic_coordinates[i].y,
                                                    atomic_coordinates[i].z))
  output.close()

  return

### The "barycenter" function ###

def barycenter(wannier_functions):
  """
   This function takes a list of Wannier Functions and computes
   its barycenter.

   :Input(s):  a list of Wannier Function objects

   :Output(s): a rank-1 numpy array containing the x, y, z
               coordinates of the barycenter

  """

  import numpy

  x_coordinates = []
  y_coordinates = []
  z_coordinates = []
  for wf in wannier_functions:
    x_coordinates.append(wf.x)
    y_coordinates.append(wf.y)
    z_coordinates.append(wf.z)
  x_coordinates = numpy.array(x_coordinates)
  y_coordinates = numpy.array(y_coordinates)
  z_coordinates = numpy.array(z_coordinates)

  barycenter_x = x_coordinates.mean()
  barycenter_y = y_coordinates.mean()
  barycenter_z = z_coordinates.mean()
  barycenter = numpy.array([barycenter_x,barycenter_y,barycenter_z]) 

  return barycenter

### The "print_wf" function for pretty print of a list of WF ###

def print_wf(wf_list):
  """
   A utility function that allows the user to print the information
   contained in each Wannier Function of a list

   :Input(s):  a list of Wannier Function objects

   :Output(s): a string printed on the screen

  """
  print ""
  print " Index        x       y       z      on-site"
  string_to_print = ""
  for wf in wf_list:
    string_to_print += "%6i   %6.3f  %6.3f  %6.3f      %6.3f \n" %(wf.index,
                                                                   wf.x,
                                                                   wf.y,
                                                                   wf.z,
                                                                   wf.onsite)
  print string_to_print

  return

### The "sort" function ###

def sort(wf_list,cond_dir,second_dir,delta):
  """
   This function takes a list of Wannier Functions and sorts the
   atomic positions according to a scheme similar to the one in
   Wannier90.

   The first direction of sorting is `cond_dir` and the
   second direction of sorting is `second_dir`.

   :Input(s):   a list of Wannier Function objects

                a string that gives the first direction of sorting (which
                corresponds to the direction of conduction)

                a string that gives the second direction of sorting

                a float which gives the distance threshold for grouping WF

   :Output(s):  a list of Wannier Function objects properly sorted

  """

  import numpy

  # we check first that cond_dir and second_dir are distinct
  if cond_dir==second_dir:
    print ""
    print " Error in sort :"
    print " The direction of conduction and the second"
    print " sorting direction are equal. Something is"
    print " wrong."
    print " Exiting the program..."
    print ""
    sys.exit(1)

  # first we extract the coordinates
  coord1 = []
  coord2 = []
  coord3 = []
  if cond_dir=='x':
    if second_dir=='y':
      for wf in wf_list:
        coord1.append(wf.x)
        coord2.append(wf.y)
        coord3.append(wf.z)
    else:
      for wf in wf_list:
        coord1.append(wf.x)
        coord2.append(wf.z)
        coord3.append(wf.y)
  elif cond_dir=='y':
    if second_dir=='z':
      for wf in wf_list:
        coord1.append(wf.y)
        coord2.append(wf.z)
        coord3.append(wf.x)
    else:
      for wf in wf_list:
        coord1.append(wf.y)
        coord2.append(wf.x)
        coord3.append(wf.z)
  else:
    if second_dir=='x':
      for wf in wf_list:
        coord1.append(wf.z)
        coord2.append(wf.x)
        coord3.append(wf.y)
    else:
      for wf in wf_list:
        coord1.append(wf.z)
        coord2.append(wf.y)
        coord3.append(wf.x)
  coord1 = numpy.array(coord1)
  coord2 = numpy.array(coord2)
  coord3 = numpy.array(coord3)

  # a function to find the group structure of a list of WF
  # in a given direction
  def find_gp(wf_index_list,coord):
    groups = []
    in_group = numpy.zeros(len(wf_index_list),dtype='bool')
    for i in xrange(len(wf_index_list)):
      if in_group[i]==False:
        in_group[i]=True
        gp_i = []
        gp_i.append(wf_index_list[i])
        for j in xrange(min(i+1,len(wf_index_list)),len(wf_index_list)):
          if in_group[j]==False:
            if abs(coord[wf_index_list[j]]-coord[wf_index_list[i]])<=delta:
              in_group[j]=True
              gp_i.append(wf_index_list[j])
        groups.append(gp_i)
    return groups

  # now we sort in the first direction
  initial_order = coord1.argsort()

  # now we group the WF in the first direction
  group = find_gp(initial_order,coord1)

  # we now look at each group and sorts in the second direction
  # group again and sort  each sub_group in the third direction
  final_order = []
  for gp in group:
    positions = numpy.array([coord2[i] for i in gp])
    sorted_gp = [gp[i] for i in positions.argsort()]
    sub_gp = find_gp(sorted_gp,coord2)
    for sgp in sub_gp:
      spositions = numpy.array([coord3[i] for i in sgp])
      sorted_sgp = [sgp[i] for i in spositions.argsort()]
      final_order.extend(sorted_sgp)

  # now we can construct the final list of sorted Wannier Functions
  final_wf_list = []
  for pos in final_order:
    new_wf = wannier_function()
    new_wf.setindex(wf_list[pos].index)
    new_wf.setx(wf_list[pos].x)
    new_wf.sety(wf_list[pos].y)
    new_wf.setz(wf_list[pos].z)
    new_wf.setonsite(wf_list[pos].onsite)
    final_wf_list.append(new_wf)

  return final_wf_list

### The "atom" class ###$

class atom:
  """
   This is the class for Atoms

   An atom object has the following attributes :
     - a symbol (a string) (access through `*.symbol`)
     - an x coordinate (a float) (access through `*.x`)
     - a  y coordinate (a float) (access through `*.y`)
     - a  z coordinate (a float) (access through `*.z`)
  """
  def setsymbol(self,atom_symbol):
    self.symbol = str(atom_symbol)
  def setx(self,x_coord):
    self.x = float(x_coord)
  def sety(self,y_coord):
    self.y = float(y_coord)
  def setz(self,z_coord):
    self.z = float(z_coord)

### The "wannier_function" class ###

class wannier_function:
  """
   This is the class for Wannier Functions

   A wannier_function object has the following attributes :
     - an index (an integer) (access through `*.index`)
     - an x coordinate (a float) (access through `*.x`)
     - a  y coordinate (a float) (access through `*.y`)
     - a  z coordinate (a float) (access through `*.z`)
     - an on-site matrix element (a float) (access through `*.onsite`)
  """
  def setindex(self,wf_index):
    self.index = int(wf_index)
  def setx(self,x_coord):
    self.x = float(x_coord)
  def sety(self,y_coord):
    self.y = float(y_coord)
  def setz(self,z_coord):
    self.z = float(z_coord)
  def setonsite(self,wf_onsite):
    self.onsite = float(wf_onsite)

