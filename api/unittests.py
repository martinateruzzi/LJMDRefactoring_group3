import unittest
from ljmd import *

class Test(unittest.TestCase):
	md = Ljmd("argon_003.inp")

	#TOREAD: file in which you can find the output to compare with reference is "ergfile" 

	def test_force(self): #use ergfile as frcfile!!! Otherwise mess with input
        md.force()
        mk.ekin()
        """Run test for the forces.

        This method will call the run of the simulation and print the forces
        on the C side. It will evolve the system for nstep times.

        """
        print("\nStarting test force with {} atoms for {} steps.".format(
            	self.natoms, self.nsteps))
        try:
        	with open("results/"+self.ergfile, "a") as file:
        		self.sys.nfi = 1
        		while self.sys.nfi <= self.sys.nsteps:
        			if(self.sys.nfi % self.nprint == 0):	
        				for i in range(self.natoms):
               				_ergstr = "%.3f, %.3f, %.3f\n" %(
                       		self.sys.fx[i], self.sys.fy[i], self.sys.fz[i])
                       		print(_ergstr)
                       		file.write(_ergstr)
           	self.velverlet()
           	self.ekin()
           	self.sys.nfi += 1
           	print("Forces have been calculated!\n")
        except Exception as err:
        	print("Error writing to file: {}".format(str(err)))

	def test_integration(self):
        md.force()
        mk.ekin()
        """Run test for integration.

        This method will call the run of the updates of velocities and positions
        on the C side. It will evolve the system for one step.

        """
        print("\nStarting test integration with {} atoms for {} steps.".format(
            	self.natoms, 1))
        try:
        	with open("results/"+self.ergfile, "a") as file:
        		
        		#setting forces manually
        		for i in range(self.natoms):
        			self.sys.fx[i] = i*10e3
        			self.sys.fy[i] = -i*10e3
        			self.sys.fz[i] = i*10e3

        		self.update_velocities_positions()

        		#setting forces manually
        		for i in range(self.natoms):
        			self.sys.fx[i] = i*10e4
        			self.sys.fy[i] = -i*10e3
        			self.sys.fz[i] = i*10e4

        		self.update_velocities()

        		for i in range(self.natoms):
               		_ergstr = "Ar  %20.8f %20.8f %20.8f\n" %(
                	self.sys.rx[i], self.sys.ry[i], self.sys.rz[i])
                    file.write(_ergstr)

                file.write("\n")

                for i in range(self.natoms):
               		_ergstr = "Ar  %20.8f %20.8f %20.8f\n" %(
                	self.sys.vx[i], self.sys.vy[i], self.sys.vz[i])
                    file.write(_ergstr)
                       		
           	           	
           	print("Test integration Done.\n")
        except Exception as err:
        	print("Error writing to file: {}".format(str(err)))


	def test_kinetic(self):
        md.force()
        mk.ekin()
        """Run test for kinetic.

        This method will call the run of the updates of velocities and positions
        on the C side. It will evolve the system for one step.

        """
        print("\nStarting test kinetic with {} atoms for {} steps.".format(
            	self.natoms, 1))
        try:
        	with open("results/"+self.ergfile, "a") as file:
        		
        		self.ekin()

           		_ergstr = " % 20.8f \n" %(self.sys.ekin)
				file.write(_ergstr)
                       		
           	           	
           	print(" Ekin test Done.\n")
        except Exception as err:
        	print("Error writing to file: {}".format(str(err)))
        
		
	def test_io(self):
        #md.force()
        #mk.ekin()
        #"""Run test for I/O.

        #This method will .

        #"""
        #print("\nStarting test kinetic with {} atoms for {} steps.".format(
        #    	self.natoms, 1))
        #try:
        #	with open("results/"+self.ergfile, "a") as file:
        #		
        #		self.ekin()

        #   		_ergstr = " % 20.8f \n" %(self.sys.ekin)
	#			file.write(_ergstr)
        #               		
        #   	           	
        #   	print(" Ekin test Done.\n")
        #except Exception as err:
        #	print("Error writing to file: {}".format(str(err)))
        
		
