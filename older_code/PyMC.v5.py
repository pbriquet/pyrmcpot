#!/usr/bin/env python3

'''
Import modules
'''
from __future__ import print_function
import random as rd
import math as math
import shutil as sh
import time as tm
from mpi4py import MPI
from lammps import lammps

'''
Define a class where the attributes and methods used in an off-lattice Metropolis
Monte Carlo simulation are defined
'''
class MonteCarlo:
	comm=MPI.COMM_WORLD	# Create an MPI communicator
	me=comm.Get_rank()	# Processor rank in a parallel simulation
	lmp=lammps()		# Create the LAMMPS object that allows to drive LAMMPS for energy calculations

	inistep=None		# Initial step
	ecurr=None		# Energy of the current configuration
	nsites=None		# Number of atoms
	interstitial=None	# Species of the interstitial atom, if any
	
	t=[]			# Store the type of the species occupying a given site
	boxlo=[]		# Store the lowest boundaries of the simulation box
	boxhi=[]		# Store the highest boundaries of the simulation box

	'''
	Methods:
	'''
	def __init__(self):	# Class initialization method... a lot of things in it
		'''
		---------------------------------------------------------------------------------
		List of the MC simulation parameters:

		init => Initial configuration generated at random ('from_scratch') or taken from a 
			checkpoint file ('restart')
		temp => Temperature in Kelvin (float)
		press => Pressure in bars (float)
		nspecsubs => Number of substitutional impurity species (integer)
		nintersatoms => Number of insterstitial atoms to be introduced in the simulation box (integer)
		concentration[*] => The fractional concentrations of the substitutional impurities (float)
		nsteps => Number of MC steps (integer)
		nrunequil => Number of MD steps for equilibrating the volume of the simulation box (integer)
		latfile => File containing the lattice points (string)
		production => Whether this is a production (True) or equilibration (False) run
		inptemplate => Template of the LAMMPS input file without run commands (string)
		maxdisp => Maximum displacement of a lattice atom (float)
		xcprob => Probability of performing position exchanges between different atom types on the lattice (float)
		dispprob => Probability of displacing atoms slightly  (float)
		enesave => Write energies at every this interval (integer)
		ncheckpoint => Write checkpoint data at every this interval (integer)
		ndump => Write LAMMPS dump files in production runs at every this interval (integer)
		xlo => Lowest system boundary in X (float)
		ylo => Lowest system boundary in Y (float)
		zlo => Lowest system boundary in Z (float)
		xhi => Highest system boundary in X (float)
		yhi => Highest system boundary in Y (float)
		zhi => Highest system boundary in Z (float)
		---------------------------------------------------------------------------------

		Below the default values of the MC control parameters are assigned as strings (or None) 
		to the entries of a Python dictionary
		'''
		self.params=dict()	# Python dictionary containing the MC simulation parameters
		self.params['init']='from_scratch'
		self.params['temp']='300.0'
		self.params['press']='0.0'
		self.params['nspecsubs']='1'
		self.params['nintersatoms']='0'
		self.params['nsteps']='1000000'
		self.params['nrunequil']='0'
		self.params['production']='False'
		self.params['latfile']='lattice.ini'
		self.params['inptemplate']='input.template'
		self.params['maxdisp']='0.02'
		self.params['xcprob']='0.1'
		self.params['dispprob']='0.9'
		self.params['enesave']='1000'
		self.params['ncheckpoint']='1000'
		self.params['ndump']='1000'
		self.params['xlo']=None
		self.params['ylo']=None
		self.params['zlo']=None
		self.params['xhi']=None
		self.params['yhi']=None
		self.params['zhi']=None

		MAXIMP=4	# Define the maximum number of substitutional impurity species

		self.params['concentration[0]']=1.0

		for i in range(MAXIMP):	# Initialize the concentrations
			self.params['concentration['+str(i+1)+']']='0.1'

		self.lmp.command("print '-----------------------------------------------' file mc.log")
		self.lmp.command("print '              Python Monte Carlo               ' append mc.log")
		self.lmp.command("print '    developed by Prof. Roberto G. A. Veiga     ' append mc.log")
		self.lmp.command("print '   at the Federal University of ABC, Brazil    ' append mc.log")
		self.lmp.command("print '      Funded by FAPESP grant 2014/10294-4      ' append mc.log")
		self.lmp.command("print '-----------------------------------------------' append mc.log")
		self.lmp.command("print '' append mc.log")
		self.lmp.command("print '*** MC simulation started! ***' append mc.log")
		self.lmp.command("print '' append mc.log")

		if self.me==0:	# I/O tasks are performed by the first processor only
			try:	# Try to open the file containing the simulation parameters
				f=open('mc.inp','r')
			except:
				print("File 'mc.inp' not found!")

			for line in f:	# Read the simulation parameters from the file
				if len(line.strip()) == 0: break

				l=line.split()

				if len(l)>=2:
					l[0]=l[0].lower()
					self.params[l[0]]=l[1]
				else:
					raise IOError("Entries in 'mc.inp' must be pairs 'keyword value'!")
					quit()

			f.close()
 
			'''
			Check the consistency of the parameters:
			'''
			self.params['init']=self.params['init'].lower()

			if not (self.params['init']=='from_scratch' or self.params['init']=='restart'):
				raise ValueError("'init' must either be 'from_scratch' or 'restart'!")
				quit()

			try:	# Is the temperature OK?
				self.params['temp']=float(self.params['temp'])
			except:
				print("Temperature must be a floating point number!")
			else:
				if self.params['temp']<=0.0:
					raise ValueError("Temperature must be greater than zero!")
					quit()

			kb=8.6173303e-5		# Boltzmann constant in eV/K
			self.beta=1.0/(kb*self.params['temp'])

			try:	# Is the pressure OK?
				self.params['press']=float(self.params['press'])
			except:
				print("Pressure must be a floating point number!")
			else:
				if self.params['press']<0.0:
					raise ValueError("Pressure must be greater than zero!")
					quit()
			
			try:	# Is the number of substitutional impurity species OK?
				self.params['nspecsubs']=int(self.params['nspecsubs'])
			except:
				print("The number of substitutional impurity species must be an integer!")
			else:
				if self.params['nspecsubs']<0 or self.params['nspecsubs']>MAXIMP:
					raise ValueError("The number of substitutional impurity species must be a positive number lower than "+str(MAXIMP)+"!")
					quit()			

			if self.params['nspecsubs']>0:			
				for i in range(self.params['nspecsubs']):	# Are the concentrations OK?
					try:
						self.params['concentration['+str(i+1)+']']=float(self.params['concentration['+str(i+1)+']'])
					except:
						print("Concentration must be a floating point number!")
					else:
						if self.params['concentration['+str(i+1)+']']<=0.0 or self.params['concentration['+str(i+1)+']']>0.5: 
							raise ValueError("Concentration of impurity "+str(i+1)+" is out of bounds!")
							quit()

					self.params['concentration[0]']=self.params['concentration[0]']-self.params['concentration['+str(i+1)+']']

			try:	# Is the number of interstitial impurities OK?
				self.params['nintersatoms']=int(self.params['nintersatoms'])
			except:
				print("The number of interstitial impurities must be an integer!")
			else:
				if self.params['nintersatoms']<0:
					raise ValueError("The number of interstitial impurities cannot be negative!")
					quit()
		
			try:	# Is the number of MC steps OK?
				self.params['nsteps']=int(self.params['nsteps'])
			except:
				print("Number of MC steps must be an integer greater than zero!")
			else:
				if self.params['nsteps']<=0: 
					raise ValueError("Number of MC steps must be an integer greater than zero!")
					quit()

			try:	# Is the number of MD steps for volume equilibration OK?
				self.params['nrunequil']=int(self.params['nrunequil'])
			except:
				print("Number of MD steps for volume equilibration must be a positive integer!")
			else:
				if self.params['nrunequil']<0: 
					raise ValueError("Number of MD steps for volume equilibration must be a positive integer!")
					quit()

			try:	# Is the maximum displacement of lattice atoms OK?
				self.params['maxdisp']=float(self.params['maxdisp'])
			except:
				print("Maximum atomic displacement must be a floating point number greater than zero!")
			else:
				if self.params['maxdisp']<=0.0:
					raise ValueError("Maximum atomic displacement must be a floating point number greater than zero!")
					quit()

			if self.params['nspecsubs']>0:
				try:	# Is the exchange probability of atom types OK?
					self.params['xcprob']=float(self.params['xcprob'])
				except:
					print("Exchange probability must be positive and less than one!")
				else:
					if self.params['xcprob']<0.0 or self.params['xcprob']>=1.0:
						raise ValueError("Exchange probability must be positive and less than one!")
						quit()

			if self.params['nintersatoms']>0:
				try:	# Is the probability of slightly displacing an atom OK?
					self.params['dispprob']=float(self.params['dispprob'])
				except:
					print("Displacement probability of lattice atoms must be greater than zero and less than one!")
				else:
					if self.params['dispprob']<=0.0 or self.params['dispprob']>=1.0:
						raise ValueError("Displacement probability of lattice atoms must be greater than zero and less than one!")
						quit()

				self.interstitial=self.params['nspecsubs']+2

			try:	# Is the energy save interval OK?
				self.params['enesave']=int(self.params['enesave'])
			except:
				print("Energy save interval must be specified as an integer greater than zero!")
			else:
				if self.params['enesave']<=0:
					raise ValueError("Energy save interval must be an integer greater than zero!")
					quit()

			try:	# Is the checkpoint interval OK?
				self.params['ncheckpoint']=int(self.params['ncheckpoint'])
			except:
				print("Checkpoint interval must be specified as an integer greater than zero!")
			else:
				if self.params['ncheckpoint']<=0:
					raise ValueError("Checkpoint interval must be an integer greater than zero!")
					quit()

			try:	# Is the dump interval OK?
				self.params['ndump']=int(self.params['ndump'])
			except:
				print("Dump interval must be specified as an integer greater than or equal to zero!")
			else:
				if self.params['ndump']<0:
					raise ValueError("Dump interval must be an integer greater than or equal to zero!")
					quit()

			self.params['production']=self.params['production'].lower()

			if self.params['production']=='true':	# Is the production keywork a Boolean?
				self.params['production']=True
			elif self.params['production']=='false': 
				self.params['production']=False
			else: 
				raise ValueError("production must be either True or False!")
				quit()

			try:	# Have the system boundaries been correctly set?
				self.params['xlo']=float(self.params['xlo'])
				self.params['ylo']=float(self.params['ylo'])
				self.params['zlo']=float(self.params['zlo'])
				self.params['xhi']=float(self.params['xhi'])
				self.params['yhi']=float(self.params['yhi'])
				self.params['zhi']=float(self.params['zhi'])
			except:
				print("The system boundaries were ill defined!")
			else:
				if(self.params['xhi']<=self.params['xlo']): 
					raise ValueError("XHi must be greater than XLo!")
					quit()
				elif(self.params['yhi']<=self.params['ylo']): 
					raise ValueError("YHi must be greater than YLo!")
					quit()
				elif(self.params['zhi']<=self.params['zlo']): 
					raise ValueError("ZHi must be greater than ZLo!")
					quit()

		self.comm.barrier()	# Force synchronization
		
		'''
		Broadcast dictionary values to all processors
		'''
		self.params=self.comm.bcast(self.params,root=0)

		self.lmp.command("print '   Temperature: %f K.' append mc.log" % self.params['temp'])

		if self.me==0:	# I/O tasks are performed by the first processor only				
			'''
			Next, the initialization function of the class will create 
			lists and tuples that will store information related to the 
			atomistic systems of interest
			'''
			x=[]			# Cartesian positions (x1,x2,x3) of the sites on the lattice
			t=[]			# Type of the atoms occupying the sites on the lattice
			spec=[[],[],[],[],[]]	# List of sites occupied by an atomic species on the lattice

			'''
			The lists containing the data of the system are populated
			with the data extracted from the lattice file
			'''
			try:	# Open (or try to open) the file containing the lattice
				f=open(self.params['latfile'],'r')
			except:
				print("File '%s' not found!" % self.params['latfile'])
		
			n=0	# Just a counter...

			for line in f:	# Loop over the lines of the lattice file, gathering the lattice data
				if len(line.strip()) == 0: break	# Use a blank line to indicate the end of file

				l=line.split()
				n=n+1
			
				if len(l)>=3:
					try:
						x1=float(l[0])
						x2=float(l[1])
						x3=float(l[2])
						x.append(tuple([x1,x2,x3]))
						t.append(1)
					except:
						print("Wrong data types in '%s'!" % self.params['latfile'])
				else: 
					raise IOError("Wrong number of columns at position %d in '%s'!" % (n,self.params['latfile']))
					quit()

			f.close()

			self.nsites=n

			'''
			In the following lines, the atoms of type 1 and atoms of other species are
			separated into distinct lists containing the indexes of the sites that
			they occupy
			'''
			if self.params['init']=='from_scratch':	# MC simulation starts from scratch
				self.inistep=0
				self.naccepted=0

				'''
				In a production run, the atomic species occupying each site are
				read from a file saved in the last iteration of an equilibration
				run
				'''
				if self.params['production']:
					try:
						f=open('last_equil.dat','r')
					except:
						print("File 'last_equil.dat' not found!")

					l=f.readline().split()	# Read the file header with information from the equilibration run

					if len(l)>=3:
						try:
							self.ecurr=float(l[1])	# Energy of the last configuration in the equilibration run
						except:
							print("Wrong data types at 'last_equil.dat' header!")
					else:
						raise IOError("Invalid number of columns at 'last_equil.dat' header!")
						quit()

					'''
					Read the atomic type of the sites from last_equil.dat and assign
					the site index to the appropriate list
					'''
					l=f.readline().split()

					if len(l)>=self.nsites:
						for i in range(self.nsites):
							try:
								v=int(l[i])
							except:
								print("Site types in 'last_equil.dat' must be integers!")
							else:
								if v>self.params['nspecsubs']+1: 
									raise ValueError("Site types in 'last_equil.dat' must be between 1 and "+str(self.params['nspecsubs']+1)+"!")
									quit()
							spec[v-1].append(i+1)
							t[i]=v				
					else:
						raise IOError("Invalid number of sites represented in 'last_equil.dat'!")
						quit()

					f.close()

					'''
					Read the information of the minimum energy configuration
					so far
					'''
					try:
						f=open('minconf.dat','r')
					except:
						print("File 'minconf.dat' not found!")

					l=f.readline().split()

					if len(l)>=3:
						try:
							self.emin=float(l[1])
						except:
							print("Wrong data types at 'minconf.dat' header!")
					else:
						raise IOError("Invalid number of columns at 'minconf.dat' header!")
						quit()
				else:	# In an equilibration run, select randomly the indices to be occupied
					for i in range(self.nsites):
						spec[0].append(i+1)
					
					for k in range(self.params['nspecsubs']):
						for i in range(round(self.nsites*self.params['concentration['+str(k+1)+']'])):
							j=rd.randrange(0,len(spec[0]))
							v=spec[0].pop(j)
							spec[k+1].append(v)
							t[spec[k+1][i]-1]=k+2
			elif self.params['init']=='restart':		# Restart a previous MC simulation
				try:
					f=open('checkpoint.dat','r')
				except:
					print("File 'checkpoint.dat' not found!")

				'''
				First, read the first line of 'checkpoint.dat', which
				contains the last saved step of the previous run as well
				as its current total energy and number of accepted trial
				moves 
				'''		
				l=f.readline().split()

				if len(l)>=3:
					try:
						self.inistep=int(l[0])+1
						self.ecurr=float(l[1])
						self.naccepted=int(l[2])
					except:
						print("Wrong data types at 'checkpoint.dat' header!")
				else:
					raise IOError("Invalid number of columns at 'checkpoint.dat' header!")
					quit()

				'''
				Read the atomic type of the sites from 'checkpoint.dat' and assign
				the site index to the appropriate list
				'''
				l=f.readline().split()

				if len(l)>=self.nsites:
					for i in range(self.nsites):
						try:
							v=int(l[i])
						except:
							print("Site types in 'checkpoint.dat' must be integers!")
						else:
							if v>self.params['nspecsubs']+1: 
								raise ValueError("Site types in 'checkpoint.dat' must be between 1 and "+str(self.params['nspecsubs']+1)+"!")
								quit()

						spec[v-1].append(i+1)
						t[i]=v		
				else:
					raise IOError("Invalid number of sites represented in 'checkpoint.dat'!")
					quit()

				f.close()

				'''
				Read the information of the minimum energy configuration
				so far
				'''
				try:
					f=open('minconf.dat','r')
				except:
					print("File 'minconf.dat' not found!")

				l=f.readline().split()

				if len(l)>=3:
					try:
						self.emin=float(l[1])
					except:
						print("Wrong data types at 'minconf.dat' header!")
				else:
					raise IOError("Invalid number of columns at 'minconf.dat' header!")
					quit()	
			
				f.close()
				
			'''
			Create the initial simulation box
			IMPORTANT: The template script MUST contain the command
			'read_data simbox.read_data'
			in order to open the coordinates file that will be created 
			by the Python driver script
			'''
			if self.params['production']==False and self.params['init']=='from_scratch':
				f=open('simbox.read_data','w')

				f.write("# Initial coordinate file provided to LAMMPS\n")

				f.write("%d atoms\n" % self.nsites)
				
				if self.params['nintersatoms']==0:
					f.write("%d atom types\n" % (self.params['nspecsubs']+1))
				else:
					f.write("%d atom types\n" % (self.params['nspecsubs']+2))

				f.write("%f %f xlo xhi\n" % (self.params['xlo'],self.params['xhi']))
				f.write("%f %f ylo yhi\n" % (self.params['ylo'],self.params['yhi']))
				f.write("%f %f zlo zhi\n" % (self.params['zlo'],self.params['zhi']))
				f.write("\n")
				f.write("Atoms\n")
				f.write("\n")

				for i in range(self.nsites):
					f.write("%d %d %f %f %f\n" % (i+1,t[i],x[i][0],x[i][1],x[i][2]))

				f.close()		

		self.comm.barrier()	# Force synchronization

		'''
		Broadcast variable value to all processors in the MPI pool
		'''
		self.inistep=self.comm.bcast(self.inistep,root=0)
		self.ecurr=self.comm.bcast(self.ecurr,root=0)
		self.nsites=self.comm.bcast(self.nsites,root=0)
		self.interstitial=self.comm.bcast(self.interstitial,root=0)

		'''
		Run line-by-line the first part of the template 
		script, which should contain simulation settings only, 
		passed to the LAMMPS library
		'''
		self.lmp.command("atom_modify map hash")	# This must be set in order to use the scatter_atoms() function later

		cmds=[]

		if self.me==0:	# I/O tasks are performed by the first processor only	
			try:
				f=open(self.params['inptemplate'],'r')
			except:
				print("File '%s' not found!" % self.params['inptemplate'])

			cmds=f.readlines()

			f.close()

		self.comm.barrier()	# Force synchronization

		cmds=self.comm.bcast(cmds,root=0)

		for cmd in cmds:
			self.lmp.command(cmd)

		'''
		Assign values to class' attributes related to 
		the atomistic system of interest
		'''
		if self.me==0:
			self.t=t
			self.spec=spec			

			if self.params['nintersatoms']>0:
				for i in range(self.params['nintersatoms']):
					self.t.append(self.interstitial)

		if self.params['init']=='restart': 
			self.lmp.command("print '   Restarting the MC simulation from step %d.' append mc.log" % self.inistep)

		'''
		Assign the values of box boundaries
		'''
		self.boxlo.append(float(self.lmp.extract_global("boxxlo",1)))
		self.boxlo.append(float(self.lmp.extract_global("boxylo",1)))
		self.boxlo.append(float(self.lmp.extract_global("boxzlo",1)))
		self.boxhi.append(float(self.lmp.extract_global("boxxhi",1)))
		self.boxhi.append(float(self.lmp.extract_global("boxyhi",1)))
		self.boxhi.append(float(self.lmp.extract_global("boxzhi",1)))

		self.comm.barrier()     # Force synchronization

		self.t=self.comm.bcast(self.t,root=0)
		
	'''
	Perform an MC trial move
	'''
	def move(self,step):
		r=None		# Initialize the variables involved in the 
		gamma=None	# acceptance condition to be broadcasted to the
		de=None		# MPI processors when running this function

		minconf=False	# Flag that defines whether this iteration is the minimum energy one	

		tfirst=None	# Type of the first site in an exchange trial move
		tsecond=None	# Type of the second site in an exchange trial move
		first=None	# Store the LAMMPS atomic indices that will the
		second=None	# types of which will be exchanged during an exchange trial move

		sel=None	# Store the index and the displacement of the atom
		dx=None		# that will be displaced during a displacement trial move
		dy=None
		dz=None

		trialprob=None	# Probability of choosing either a displacement or an exchange trial move
		dispprob=None	# Probability of choosing a slight atomic displacement or moving an interstitial atom to a new position

		xinternew=[]	# New position of the selected interstitial atom

		self.lmp.command("print '=> Step %d is running...' append mc.log" % step)
	
		if step>0 or self.params['production']:
			if self.me==0:	# Trial move selection is performed by the first processor only
				trialprob=rd.random()	

				if self.params['nspecsubs']>0 and trialprob<self.params['xcprob']:	# Perform an exchange trial move
					sellength=[]
					chosen=False
					totlength=0.0
					r=rd.random()

					for i in range(self.params['nspecsubs']+1):	# Select the type of the first atom
						length=totlength+self.params['concentration['+str(i)+']']

						if not chosen and r<length:
							sellength.append(0.0)
							tfirst=i
							chosen=True
						else:
							sellength.append(length)
							totlength=length		

					old=rd.randrange(0,len(self.spec[tfirst]))				
					first=self.spec[tfirst][old]

					r=rd.uniform(0.0,totlength)
					
					for i in range(self.params['nspecsubs']+1):	# Select the type of the other atom
						if r<sellength[i]:
							tsecond=i

							break

					new=rd.randrange(0,len(self.spec[tsecond]))
					second=self.spec[tsecond][new]
				else:	# Perform a displacement trial move
					dispprob=rd.random()

					if self.params['nintersatoms']==0 or dispprob<self.params['dispprob']:	# Slightly displace an atom from its current position
						sel=rd.randrange(0,self.nsites+self.params['nintersatoms'])
						dx=rd.uniform(-self.params['maxdisp'],self.params['maxdisp'])
						dy=rd.uniform(-self.params['maxdisp'],self.params['maxdisp'])
						dz=rd.uniform(-self.params['maxdisp'],self.params['maxdisp'])
					else:	# Move an interstitial atom to a random new position in the simulation box
						sel=rd.randrange(self.nsites,self.params['nintersatoms']+self.nsites)
						xinternew.append(rd.uniform(self.boxlo[0],self.boxhi[0]))
						xinternew.append(rd.uniform(self.boxlo[1],self.boxhi[1]))
						xinternew.append(rd.uniform(self.boxlo[2],self.boxhi[2]))

			self.comm.barrier()	# Force synchronization

			trialprob=self.comm.bcast(trialprob,root=0)	# Broadcast the trial probability to all processors
			dispprob=self.comm.bcast(dispprob,root=0)	# Broadcast the displacement probability to all processors

			if self.params['nspecsubs']>0 and trialprob<self.params['xcprob']:	# Perform an exchange trial move
				first=self.comm.bcast(first,root=0)
				second=self.comm.bcast(second,root=0)
				tfirst=self.comm.bcast(tfirst,root=0)
				tsecond=self.comm.bcast(tsecond,root=0)

				self.lmp.command("print '   Trial move: site %d -> %d, site %d -> %d' append mc.log" % (first,tsecond+1,second,tfirst+1))
				self.lmp.command("set atom %d type %d" % (first,tsecond+1))
				self.lmp.command("set atom %d type %d" % (second,tfirst+1))
			else:	# Perform a displacement trial move
				xold=self.lmp.gather_atoms("x",1,3)     # Get the atomic coordinates and store it into the old configuration
				x=self.lmp.gather_atoms("x",1,3)	# Get the atomic coordinates and store it into the current, transient configuration

				self.comm.barrier()	# Force synchronization

				sel=self.comm.bcast(sel,root=0)

				if self.params['nintersatoms']==0 or dispprob<self.params['dispprob']:	# Slightly displace an atom from its current position
					dx=self.comm.bcast(dx,root=0)
					dy=self.comm.bcast(dy,root=0)
					dz=self.comm.bcast(dz,root=0)

					self.comm.barrier()	# Force synchronization

					self.t=self.comm.bcast(self.t,root=0)

					x[sel*3]+=dx
					x[sel*3+1]+=dy
					x[sel*3+2]+=dz			

					self.lmp.command("print '   Trial move: atom %d, type %d, displaced by (%f,%f,%f)' append mc.log" % (sel+1,self.t[sel],dx,dy,dz))
				else:	# Move an interstitial site to a new position inside the simulation box
					self.comm.barrier()	# Force synchronization

					xinternew=self.comm.bcast(xinternew,root=0)

					x[sel*3]=xinternew[0]
					x[sel*3+1]=xinternew[1]
					x[sel*3+2]=xinternew[2]

					self.lmp.command("print '   Trial move: atom %d, interstitial, moved to (%f,%f,%f)' append mc.log" % (sel+1,xinternew[0],xinternew[1],xinternew[2]))

				self.lmp.scatter_atoms("x",1,3,x)								
		else:	# First MC step
			if self.params['nrunequil']>0:	# Relax the simulation box volume and positions with an initial MD run
				self.lmp.command("minimize 0 1e-4 5000 10000")
				self.lmp.command("reset_timestep 0")
				self.lmp.command("velocity all create %f %d dist gaussian" % (self.params['temp']*2,int(tm.time())))
				self.lmp.command("fix 1 all npt temp %f %f 100 aniso %f %f 1000 drag 2.0" % (self.params['temp'],self.params['temp'],self.params['press'],self.params['press']))
				self.lmp.command("run %d" % self.params['nrunequil'])
				self.lmp.command("unfix 1")
				self.lmp.command("reset_timestep 0")

				'''
				Assign the new values of box boundaries
				'''
				self.boxlo[0]=self.lmp.extract_global("boxxlo",1)
				self.boxlo[1]=self.lmp.extract_global("boxylo",1)
				self.boxlo[2]=self.lmp.extract_global("boxzlo",1)
				self.boxhi[0]=self.lmp.extract_global("boxxhi",1)
				self.boxhi[1]=self.lmp.extract_global("boxyhi",1)
				self.boxhi[2]=self.lmp.extract_global("boxzhi",1)

			if self.params['nintersatoms']>0:   # Create interstitial atoms at random positions
				self.lmp.command("create_atoms %d random %d %d NULL" % (self.interstitial,self.params['nintersatoms'],int(tm.time())))

			x=self.lmp.gather_atoms("x",1,3)        # Get the atomic coordinates and store it into the current, transient configuration

			if self.params['nrunequil']==0:	# Perturb atomic positions randomly
				for i in range(self.nsites+self.params['nintersatoms']):
					if self.me==0:
						dx=rd.uniform(-self.params['maxdisp'],self.params['maxdisp'])
						dy=rd.uniform(-self.params['maxdisp'],self.params['maxdisp'])
						dz=rd.uniform(-self.params['maxdisp'],self.params['maxdisp'])

					self.comm.barrier()	# Force synchronization

					dx=self.comm.bcast(dx,root=0)
					dy=self.comm.bcast(dy,root=0)
					dz=self.comm.bcast(dz,root=0)

					x[i*3]+=dx
					x[i*3+1]+=dy
					x[i*3+2]+=dz

				self.lmp.scatter_atoms("x",1,3,x) 

		self.lmp.command("variable e equal etotal")
		self.lmp.command("thermo_modify lost ignore")
		self.lmp.command("run 0")

		totatoms=int(self.lmp.get_natoms())

		if totatoms<(self.nsites+self.params['nintersatoms']):	
			e=0.0
		else:
			e=float(self.lmp.extract_variable("e","all",0))

			if math.isnan(e):
				e=0.0

		'''
		Apply the acceptance condition to decide whether the trial move
		will be accepted or not
		'''
		if step>0 or self.params['production']:
			if self.me==0:	# Serial computations are performed by the first processor only
				de=e-self.ecurr
				r=rd.random()

				try:
					gamma=math.exp(-self.beta*de)
				except:
					if de<=0.0:
						gamma=1.0
					else:
						gamma=0.0

			self.comm.barrier()	# Force synchronization

			de=self.comm.bcast(de,root=0)
			r=self.comm.bcast(r,root=0)
			gamma=self.comm.bcast(gamma,root=0)

			if r<gamma:
				self.lmp.command("print '   Energy variation: %f eV; trial move accepted.' append mc.log" % de)

				if self.me==0:	# I/O tasks are performed by the first processor only
					if self.params['nspecsubs']>0 and trialprob<self.params['xcprob']:
						'''
						Remove the entries from the corresponding lists, and then
						save the site indices in buffer variables
						'''
						buffer1=self.spec[tsecond].pop(new)
						buffer2=self.spec[tfirst].pop(old) 

						'''
						Update permanently the class attributes changed during the trial move
						'''
						self.t[buffer1-1]=tfirst+1
						self.t[buffer2-1]=tsecond+1
						self.spec[tsecond].append(buffer2)
						self.spec[tfirst].append(buffer1)

					self.ecurr=e
					self.naccepted+=1

					if e<self.emin:
						self.emin=e
						minconf=True

				self.comm.barrier()	# Force synchronization

				minconf=self.comm.bcast(minconf,root=0)
				self.ecurr=self.comm.bcast(self.ecurr,root=0)
			else:	# Undo the trial move
				self.lmp.command("print '   Energy variation: %f eV; trial move rejected.' append mc.log" % de)
				self.lmp.command("print '   Undoing the trial move...' append mc.log")

				if self.params['nspecsubs']>0 and trialprob<self.params['xcprob']:
					self.lmp.command("set atom %d type %d" % (first,tfirst+1))
					self.lmp.command("set atom %d type %d" % (second,tsecond+1))
				else:		
					self.lmp.scatter_atoms("x",1,3,xold)
		else:	# What if this is the first MC step of the equilibration run?
			if self.me==0:
				self.ecurr=e
				self.emin=e
				minconf=True

			minconf=self.comm.bcast(minconf,root=0)
			self.ecurr=self.comm.bcast(self.ecurr,root=0)

			if totatoms<(self.nsites+self.params['nintersatoms']):
				self.lmp.scatter_atoms("x",1,3,x)

		self.lmp.command("reset_timestep 0")

		'''
		Write information about the minimum energy configuration
		so far
		'''
		if minconf:
			if self.me==0:	# I/O tasks are performed by the first processor only
				try:
					sh.copy("minconf.dat","minconf.tmp")
				except:
					print("File 'minconf.dat' not found!")

				f=open('minconf.dat','w')
				f.write("%d %f %d\n" % (step,self.emin,self.naccepted))

				for i in range(self.nsites):
					f.write("%d " % self.t[i])

				f.write("\n")
				f.close()

			self.lmp.command("write_data minconf.read_data")

		'''
		Save the current energy to energies.dat and the dump files
		'''
		if self.params['production']:
			if (step+1)%self.params['enesave']==0:	# Save the energy into energies.dat
				self.lmp.command("print '%d %f' append energies.dat screen no" % (step,self.ecurr))

			if self.params['ndump']>0:	# Save dump files during production runs
				if (step+1)%self.params['ndump']==0:	
					self.lmp.command("write_dump all custom dump.%d.lammpstrj id type x y z" % step)

			if step==self.params['nsteps']-1:	# Write the last configuration in LAMMPS Data format
				self.lmp.command("write_data lastconf.read_data")
		else: # Save the current energy to energ_equil.dat and the initial configuration to iniconf.read_data
			if (step+1)%self.params['enesave']==0:
				self.lmp.command("print '%d %f' append energ_equil.dat screen no" % (step,self.ecurr))

			if step==0: 
				self.lmp.command("write_data iniconf.read_data")

		'''
		Write checkpoint at prescribed iterations
		'''
		if (step+1)%self.params['ncheckpoint']==0:
			'''
			Save the checkpoint coordinates
			'''
			self.lmp.command("write_data simbox.read_data")

			if self.me==0:	# I/O tasks are performed by the first processor only
				if (step+1)/self.params['ncheckpoint']>=2: 
					try:
						sh.copy("checkpoint.dat","checkpoint.tmp")
					except:
						print("File 'checkpoint.dat' not found!")

				f=open('checkpoint.dat','w')
				f.write("%d %f %d\n" % (step,self.ecurr,self.naccepted))

				for i in range(self.nsites):
					f.write("%d " % self.t[i])

				f.write("\n")
				f.close()

		'''
		Save the last accepted configuration in the equilibration run to be 
		later used in the production run
		'''
		if (not self.params['production']) and step==self.params['nsteps']-1:
			if self.me==0:	# I/O tasks are performed by the fist processor only
				f=open('last_equil.dat','w')
				f.write("%d %f %d\n" % (step,self.ecurr,self.naccepted))

				for i in range(self.nsites):
					f.write("%d " % self.t[i])

				f.write("\n")
				f.close()

			self.lmp.command("write_data simbox.read_data")

		self.comm.barrier()

	'''
	Below we have the class destructor. The only thing it does
	is to destroy the LAMMPS object
	'''
	def __del__(self):
		self.lmp.command("print '*** MC simulation finished! ***' append mc.log")
		self.lmp.close()

		self.comm.barrier()

		MPI.Finalize()

'''
MAIN code starts below
'''
mc=MonteCarlo()	# Instantiate a new Monte Carlo object

'''
Loop over the MC steps
'''
for i in range(mc.inistep,mc.params['nsteps']):
	mc.move(i)		# Perform the i-th trial move
	mc.comm.barrier()	# Force synchronization

del mc	# Destroy the Monte Carlo object

quit()	# Quit Python
