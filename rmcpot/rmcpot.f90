! Main program: RMCPOT
program rmcpot
	use globals
	use eamalloy
	use meam
	use tersoff
	use watmodel
	use buckingham
	use potfunc
	use datafunc
	use localmin

	implicit none

	! Constants
	integer,parameter	:: MAXSIGMADECAY=100
	integer,parameter	:: MINSIGMADECAYFOREXIT=4
	real,parameter		:: MINSTEPSIZE=0.1
	real,parameter		:: PERCDECAY=0.9

	! Work variables
	integer			:: i,j,k,shellstatus,mcsteps,shift,naccepted,parcycle,parindex
	integer			:: savebest,genpotstatus,getpotdatastatus,mcstatus
	real			:: accperpar,sigma
	double precision	:: chi,currchi,dchi,prevminchi,zeta0
	double precision	:: rnumber,defsign,perturb,tmpval1,tmpval2
	logical			:: tmpconst
	logical,allocatable	:: maxmovefixed(:) 

	! Input variables
	logical			:: performsa=.true.,performmin=.false.
	integer			:: nstepsperpar=100,maxitermin=100,resetevery=0
	real			:: sigma0=1.0,chitol=10.0,maxstepsize=100.0
	real			:: convthr=EPSILON,alphamin=0.1,maxperturb=1.0,inistepsize=1.0
	character*10		:: minalgo='stepmin'
	character*20		:: potmodel='eam/alloy'
	character*40		:: potfile='temp.eam'
	character*80		:: runpot='lmp < lammps.in >& lammps.out &'
	logical			:: randiniparams=.false.,verbose=.false.
	logical			:: runonce=.false.,autosigma=.false.,recalcsigma=.false.

	! Namelist definitions
	namelist /inpmain/ performsa,performmin,seed,runpot,potmodel,potfile,initialparams,iniparamsfile, &
		 & normweight,penalty,randiniparams,maxperturb,verbose
	namelist /inpsa/ sigma0,chitol,nstepsperpar,maxstepsize,inistepsize,runonce,resetevery,autosigma,recalcsigma
	namelist /inpmin/ minalgo,maxitermin,convthr,alphamin

	write(6,*) '---------------------------------------------------------'
	write(6,*)
	write(6,*) '                         RMCPOT:                         ' 
	write(6,*) '            a Reverse Monte Carlo code for the           ' 
	write(6,*) '      construction of interatomic POTential models,      '
	write(6,*) '                      version 0.27,                      '
	write(6,*) '             developed by Dr. Roberto Veiga              '
	write(6,*) '        at Universidade Federal do ABC, Brazil,          '
	write(6,*) '              as part of a project funded                '
	write(6,*) '           by São Paulo Research Foundation              '
	write(6,*)
	write(6,*) '---------------------------------------------------------'
	write(6,*)

	! Initialize a few global variables
	seed=123456789
	initialparams='random'
	iniparamsfile='paramset.dat.0'
	normweight=.true.
	penalty=1.0
	chosen=1

	! Open the input file
	open(66,file='pot.inp',status='old',err=519)

	goto 520

519	stop 'RMCPOT: File pot.inp not found!'

	! Read the general input variables
520	read(66,inpmain)

	if(performsa) read(66,inpsa)	! Read the input variables of simulated annealing
	if(performmin) read(66,inpmin)	! Read the input variables for the local minimization algorithm

	close(66)

	call inputpot(potmodel)	! Initialize variables specific of a potential model
	call inirefdata		! Initialize the reference dataset arrays
	call iniparams		! Initialize the parameter set array

	! Adjust the value of the weight applied to the zeta function
	if(normweight.and.penalty>1.0)then
		penalty=1.0
	elseif(penalty<0.0)then
		penalty=0.0
	end if

	! Get the parameters that have to be kept fixed during the fitting procedure
	if(initialparams=='file')then
		open(42,file='paramlist1',status='old',err=556)

555		read(42,*,end=554) parindex 

		if(parindex>=1.and.parindex<=nparams) fixpar(parindex)=.true.

		goto 555

554		close(42)

		write(6,'(A)') '   Parameter list (1) was found!'
		write(6,'(A)') '   The following parameters will not change during the fitting procedure:'
		
		do i=1,nparams
			if(fixpar(i)) write(6,'(I6,F20.7)') i,params(i)
		end do

		write(6,'(A)')

556		continue
	end if

	! Verify if some parameters have predetermined boundary values
	open(52,file='paramlist2',status='old',err=724)

723	read(52,*,end=722) parindex,tmpval1,tmpval2,tmpconst

	if(parindex<=nparams)then
		if(parindex>=1)then
			hardconst(parindex)=tmpconst

			if(abs(tmpval1)>EPSILON) minparval(parindex)=tmpval1
			if(abs(tmpval2)>EPSILON) maxparval(parindex)=tmpval2
		else
			hardconst=tmpconst

			if(abs(tmpval1)>EPSILON) minparval=tmpval1
			if(abs(tmpval2)>EPSILON) maxparval=tmpval2
		end if
	end if

	goto 723

722	close(52)

	write(6,'(A)') '   Parameter list (2) was found!'
	write(6,'(A)') '   Minimum (absolute) values of the parameters:'

	do i=1,nparams
		if(abs(minparval(i))>EPSILON)then
			write(6,'(I6,F20.7)') i,minparval(i)
		else
			write(6,'(I6,A)') i,' (Unbounded)'
		end if
	end do

	write(6,'(A)') '   Maximum (absolute) values of the parameters:'

	do i=1,nparams
		if(abs(maxparval(i))>EPSILON)then
			write(6,'(I6,F20.7)') i,maxparval(i)
		else
			write(6,'(I6,A)') i,' (Unbounded)'
		end if
	end do

	write(6,'(A)')

724	continue
	
	! Print the reference data with associated weight factors
	write(6,'(A)') '   Reference data and weight factors:'

	do i=1,ndata
		write(6,'(I6,2F20.7)') i,refdata(i),weight(i)
	end do
	
	write(6,'(A)')

	! Apply random variations to the parameters read from the file
	if(initialparams=='file'.and.randiniparams)then
		write(6,'(A)') '   Perturbed parameter set from file (and corresponding perturbation):'

		do i=1,nparams
			if(.not.fixpar(i))then
				rnumber=rand()

				if(rnumber>=0.5)then
					defsign=1.0
				else
					defsign=-1.0
				end if

				rnumber=rand()
				perturb=rnumber*maxperturb*params(i)
				params(i)=params(i)+defsign*perturb

				write(6,'(I6,2F20.7)') i,params(i),defsign*perturb
			else
				write(6,'(I6,2F20.7,A)') i,params(i),0.0e0,' (fixed parameter)'
			end if
		end do

		bestparams=params

		write(6,'(A)')
	end if

	! Obtain the initial potential data
	zeta=0.0				! Zero the value of zeta function
	genpotstatus=genpot(potmodel,potfile)	! Generate a potential with the current parameter set

	if(genpotstatus==-1) stop 'RMCPOT: Bad potential parameters!'	! Some NaN value(s) were generated, who can guess why?

	shellstatus=system(runpot)	! Run trial simulations with the potential

	if(shellstatus==0) then		! The trial simulations finished normally
		getpotdatastatus=getpotdata()

		if(getpotdatastatus==-1) stop 'RMCPOT: Bad data generated by the potential!'	! Again, for some reason, NaN values...

		call checkdata		! Check for qualitative disagreement (different signs) between reference and potential data 

		savebest=savebestpot(potmodel)		! Save the initial potential version as the best one (initializing things...)
		bestpotdata=potdata	! Initialize the array of best potential data with the data of the initial parameter set
		currchi=fchi()		! Obtain the first current chi value from the initial parameter set
		currchi=currchi+zeta	! Add the zeta contribution from the initial parameter set
		minchi=currchi		! Initialize the minimum value of the objective function
		zeta0=zeta		! Initialize the current zeta
		minzeta=zeta0		! Initialize the value of the zeta function at the minimum of the objective function

		call saveparamset(minchi,.true.)	! Save the initial parameter set as the best one
		call savebestpotdata			! Save the initial data generated by the potential
	else
		stop 'RMCPOT: Error running the external program!'	! For some reason, the call to the external code failed
	end if

	! What if simulated annealing (and, eventually, minimization) is disabled?
	if((.not.performsa).and.(.not.performmin))then
		write(6,'(A)') 'RMCPOT: Nothing else to do!'

		goto 79
	elseif(.not.performsa.and.performmin)then
		write(6,'(A)') 'RMCPOT: Simulation annealing disabled!'		

		goto 78	! Jump directly to minimization
	endif

	! File containing the convergence data of simulation annealing
	open(77,file='mc_chiconv.dat')

	if(verbose)then
		write(77,'(A)') '# Step   Chi   Zeta   Parameter   Displacement'
		write(77,'(A)') '# -------------------------------------------------'
		write(77,*) 0,minchi-minzeta,minzeta,0,0.0	
	else
		write(77,'(A)') '# Step   Chi   Zeta'
		write(77,'(A)') '# ----------------------------------------'
		write(77,*) 0,minchi-minzeta,minzeta
	end if

	! File containing the values of the parameters at each iteration of the Monte Carlo simulation
	if(verbose)then 
		open(222,file='mc_parevol.dat')

		write(222,*) 0,params
	end if

	! Verifiy if some parameters have customized step sizes
	allocate(maxmovefixed(nparams))	

	maxmove=inistepsize		! Initialize the array with the maximum step length of a parameter
	maxmovefixed=.true.

	open(43,file='paramlist3',status='old',err=643)

644	read(43,*,end=645) parindex,tmpval1
	
	if(parindex<=nparams)then
		if(parindex>=1)then
			if(tmpval1>0.0)then
				maxmove(parindex)=tmpval1
				maxmovefixed(parindex)=.false.
			end if
		else
			if(tmpval1>0.0)then
				maxmove=tmpval1
				maxmovefixed=.false.
			end if
		end if
	end if

	goto 644

645	close(43)

	write(6,'(A)') '   Parameter list (3) was found!'
	write(6,'(A)') '   Initial step sizes per parameter:'

	do i=1,nparams
		write(6,'(I6,F20.7)') i,maxmove(i)
	end do

	write(6,*)

643	continue

	! Initialize sigma
	if(autosigma)then
		sigma=fsigma(minchi)
	else
		sigma=sigma0
	end if
				
	! Initialize the remaining MC-specific variables
	parcycle=0				! Intialize the number of cycles through the parameter set
	mcsteps=nstepsperpar*nparams		! Number of MC steps per sigma

	! Loop over sigma
	do i=1,MAXSIGMADECAY
		shift=(i-1)*mcsteps
		naccepted=0

		write(6,'(A)') '---------------'
		Write(6,'(A,I5)') '=> ITERATION:',i
		write(6,'(A,F20.7)') '=> Performing MC for sigma=',sigma
		write(6,*)


		write(77,'(A,F21.7)') '# SIGMA=',sigma

		! MC loop for the current value of sigma	
		do j=1,mcsteps
			! Assign the best parameter set so far to the current parameter set if the MC iteration 
			! is multiple of 'resetevery'
			if(resetevery>0.and.mod(j,resetevery)==0)then
				params=bestparams
				currchi=minchi
			end if

			if(fixpar(chosen)) goto 77	! What if this parameter is fixed? Nothing changes...
				
			mcstatus=mcmove()		! Change the value of one parameter

			! In case the chosen parameter became NaN (why???) or is out of the allowed boundaries
			if(mcstatus==-1)then
				params(chosen)=prevparam
				
				goto 77	
			end if
			
			call testtoosmall		! Test for (absolute) too small numbers in the parameter set

			if(verbose) tmpval1=params(chosen)-prevparam	! Displacement of the parameter in this trial move 

			! Obtain a new value of chi to be compared to the previous one
			zeta=0.0				! Reset the value of zeta function
			genpotstatus=genpot(potmodel,potfile)	! Generate the potential

			! In case the chosen parameter is NaN (WTF! Why?) or is out of the allowed boundaries
			if(genpotstatus==-1)then
				params(chosen)=prevparam

				goto 77		! Anyway, if this is the case, the trial parameter set is obviously rejected...
			end if		

			shellstatus=system(runpot)		! Run trial simulations with the potential

			if(shellstatus==0)then
				getpotdatastatus=getpotdata() 	! Get the data generated by the potential

				! In case there are NaN data generated by the potential (this is unbelievable...)
				if(getpotdatastatus==-1)then
					params(chosen)=prevparam
					
					goto 77
				end if

				call checkdata

				chi=fchi()		! Obtain chi of the current parameter set
				chi=chi+zeta		! Add the contribution of the zeta function
				dchi=chi-currchi

				if(dchi<0.0)then
					currchi=chi		! Chi becomes the new chi for the next iteration
					zeta0=zeta
					naccepted=naccepted+1
					accepted(chosen)=accepted(chosen)+1

					if(chi<minchi)then
						minchi=chi				! Update the minimum chi
						minzeta=zeta				! Value of the zeta function at minimum chi
						savebest=savebestpot(potmodel)		! Save the best potential so far
						bestpotdata=potdata			! Best data from the potential so far
						bestparams=params			! Best parameter set so far

						call saveparamset(minchi,.true.)	! Save the best parameter set so far
						call savebestpotdata			! Save the data generated by the best parameter set so far
					end if
					
					call saveparamset(chi,.false.)			! Save the current parameter set
				else
					rnumber=rand()
			
					if(rnumber>exp(-dchi/sigma))then 
						params(chosen)=prevparam	! Keep unaltered the parameter set
					else
						currchi=chi			! Chi becomes the current chi for the next iteration
						zeta0=zeta
						naccepted=naccepted+1
						accepted(chosen)=accepted(chosen)+1
	
						call saveparamset(chi,.false.)	! Save the current parameter set
					endif
				endif
			else
				params(chosen)=prevparam			! Keep unaltered the parameter set
			end if
			
			! Output convergence data
77			if(verbose)then
				write(77,*) j+shift,currchi-zeta,zeta,chosen,tmpval1
				write(222,*) j+shift,params	! Output the values of the parameters at this MC iteration
			else			
				write(77,*) j+shift,currchi-zeta0,zeta0
			end if

			! Increment the parameter index for the current MC run
			chosen=chosen+1

			! What if the parameter index is out of bounds?
			if(chosen>nparams)then
				chosen=1		! Restart the cycle through the parameter set
				parcycle=parcycle+1

				! Finally, update the step size of each parameter according to Corana's formula
				do k=1,nparams
					if(.not.fixpar(k))then
						accperpar=real(accepted(k))/real(parcycle)

						if(accperpar>0.6)then
							maxmove(k)=maxmove(k)*(1.0+2.0*(accperpar-0.6)/0.4)

							if(maxmovefixed(k).and.maxmove(k)>maxstepsize) maxmove(k)=maxstepsize
						elseif(accperpar<0.4)then
							maxmove(k)=maxmove(k)/(1.0+2.0*(0.4-accperpar)/0.4)

							if(maxmove(k)<MINSTEPSIZE) maxmove(k)=MINSTEPSIZE
						end if
					end if
				end do				
			end if	
		end do

		params=bestparams	! Assign the best parameter set to the current parameter set
		currchi=minchi		! Assign the minimum chi so far to the current chi value

		write(6,'(A)') '   MC finished for this sigma value!'
		write(6,'(A,F10.6)') '=> Ratio of accepted moves:',real(naccepted)/real(mcsteps)
		write(6,*)
		write(6,'(A,F20.6)') '=> Minimum of the objective function:',minchi
		write(6,*)

		if(i>1)then
			write(6,'(A,F15.7)') '=> Improvement of the objective function:',prevminchi-minchi
			write(6,*)
		end if

		! Test some other criteria for finishing the run
		if(minchi<=chitol)then 					! The allowed level of convergence was attained
			write(6,'(A)') '   Level of convergence attained.'
			
			exit
		elseif(i>MINSIGMADECAYFOREXIT.and.(prevminchi-minchi)<EPSILON)then	! No significant improvement in chi
			write(6,'(A)') '   The objective function was not significantly improved.'

			exit
		elseif(runonce)then
			exit
		else
			if(recalcsigma)then
				sigma=fsigma(minchi)
			else
				sigma=sigma*PERCDECAY
			end if

			prevminchi=minchi

			write(6,'(A)') '   Ready for the next MC run...'
			write(6,*)
		end if
	end do

	close(77)

	potdata=bestpotdata	! Potential-generated data that best fitted the reference data
	
	write(6,*) 'Minimum value of chi:'
	write(6,*) '>>>>>>> ',minchi
	write(6,*)
	write(6,*) 'RMCPOT: Simulated annealing finished!'
	write(6,*)

	! Enter the local minimization step
78	if(performmin)then
		write(6,'(A)') '---------------'
		write(6,'(A)') '=> Initializing minimization towards a local minimum...'

		call minimize(minalgo,potmodel,potfile,runpot,maxitermin,convthr,alphamin)
	end if

79	call savebestpotdata

	! Deallocate global arrays
	deallocate(refdata)
	deallocate(weight)
	deallocate(tolerance)
	deallocate(potdata)
	deallocate(bestpotdata)
	deallocate(params)
	deallocate(bestparams)
	deallocate(maxparval)
	deallocate(minparval)
        deallocate(maxmove)
        deallocate(accepted)
	deallocate(qualiweight)
	deallocate(fixpar)

	write(6,*) 'RMCPOT: Simulation finished!'
end program rmcpot

