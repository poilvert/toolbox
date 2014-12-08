#!/usr/bin/env python
# constraints.py

"""
   **Module information** :

   :Program:    module containing the `constraints`
                function needed by some of the programs
                like visualize_defect_influence.py or
                compactify_conductor.py
   :Purpose:    separate Wannier Functions into **defect**
                Wannier Functions and **backbone** Wannier
                Functions
   :Input(s):   3 cartesian coordinates `"x"`, `"y"` and `"z"`
   :Output(s):  a boolean

"""

def constraints(x,y,z):
  """
   This function has to be built by the user.

   One needs to construct the function so as to fit the
   following requirements :

    - the function needs to return *False* if the Wannier
      Function belongs to the `defect`
    - the function needs to return *True* if the Wannier
      Function belongs to the `backbone`

   :Input(s):  3 floats corresponding to the Wannier Function
               coordinates

   :Output(s): a boolean

  """

  from math import sqrt
  if sqrt(x**2+y**2)<=2.7:
    boolean = True
  else:
    boolean = False
  # user code here

  return boolean

