'''
McGill University, Montreal
MECH(447|652) - Dynamics of Combustion
Project #0

Submitted to Prof, Jeffrey Bergthorson
by
Student Name (Student Number)

(Month day, year)
Montreal, Quebec, Canada
'''

# =============================================================================
# Table of content
# =============================================================================

'''
Introduction
Part 1 - General-Information
Part 2 - Software
   Jupyter Notebook
   Cantera
Part 3 - Questions
Conclusion
'''

# =============================================================================
# Introduction
# =============================================================================

'''
Project 0 is an introduction to the tools used in the other projects during the semester. 
The objective of Project 0 is to ensure that everyone is able to run every element necessary for future projects:
- Spyder
- Python
- Cantera

You will also have to solve 3 basic problems.  (See below.)

The Python kernel is initialized for Cantera, Numpy, and MatplotLib with the commands below.

'''

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Part 1 - General Information
# =============================================================================

'''
Here are some useful links for more information about the installation of 
Anaconda and Cantera, and basic introductions to Python and Spyder. 
Cantera can also be installed with the Navigator provided with Anaconda.
 
Anaconda: 
* https://www.anaconda.com/download/ (download distribution)
 
Cantera: 
* https://cantera.org/documentation/ (installation instruction)
* https://anaconda.org/Cantera/cantera (Additional command from Anaconda documentation)

Spyder:
* https://medium.com/coderbyte/spyder-python-ide-for-absolute-beginners-89e4ea1832af (Introduction to Spyder and Pyhton)
* http://161.246.38.75/download/sdte/62_01Spyder.pdf  (Tutorial)

'''

# =============================================================================
# Part 2 - Software
# =============================================================================


# ### 2.1 Python

# Python 3 should be used for all projects. Ensure that the proper version of 
# Python is used with Anaconda.

import sys
print (sys.version)


# The following code snippet produces a simple plot.  Identify the different 
# variables and functions related to data generation and figure creation. Note 
# how comments can also be added to your code using the hash symbol #.


x = np.arange(0,2*np.pi,0.01)
y = np.sin(x)

# Plot Definition
figId = plt.figure(figsize=(10,4))
plt.plot(x, y, lw=2) # plot species colution

# Text Formating
plt.title('Sine Wave', fontsize = 18)
plt.xlabel('x')
plt.ylabel('y')

# Figure Formating
plt.xlim([0,2*np.pi])
plt.show()


# ### 2.2 Cantera

# Cantera 2.5.1, the latest stable release, is used in the projects throughout 
# the semester. Ensure the proper version is used.

# Display Cantera version used
print("Running Cantera Version: " + str(ct.__version__))


# The following snippet validates that Cantera is working properly. You should 
# see the equilibrium temperature printed, once the calculation is complete, 
# just below the code cell.


# import gas object with GRI30
gas = ct.Solution('gri30.cti')

# Define reactants condition
gas.TPX = 300, ct.one_atm, 'CH4:1, O2:2, N2:7.52'

# Solve equilibrium problem for adiabatic conditions
gas.equilibrate('HP')

print('Equilibrium temperature: %.2f K' % gas.T)


# =============================================================================
# Part 3 - Questions
# =============================================================================

# ### Finding Help. 

# Use the python command to display the help message of the gas command.

# From reading the help output, what is the function you would call if you wanted to get the species list?




# ### Debugging

# Correct the error(s) in the code below.

x = 3.o
  Print('x = ', x)



# ### Building a loop

# Write a for loop that iterates 5 times and prints the iteration value to the 
# screen (that is, on the first pass it should print 0, the second pass should 
# print 1, etc.



# ### Adding a second curve to the plot below.

#   Add a cosin wave with double the frequency of the sin wave.  
#   
#   Add a label for each curve and show the legend on the chart.

x = np.arange(0,2*np.pi,0.01)
y = np.sin(x)


# Plot Definition
figId = plt.figure(figsize=(10,4))
plt.plot(x, y, lw=2) # plot species colution

# Text Formating
plt.title('Sin-Cos Waves', fontsize = 18)
plt.xlabel('x')
plt.ylabel('y')

# Figure Formating
plt.xlim([0,2*np.pi])
plt.show()

# =============================================================================
# Conclusion
# =============================================================================

# Upon completion of Project 0, you should be able to run spyder with a Python 
# kernel and Cantera package. 

# Don't forget to save your script, and to save your figures and results as you 
# need. 