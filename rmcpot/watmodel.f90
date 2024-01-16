module watmodel
	use globals

	public		:: watmodel_readinput,watmodel_writepot

! Constants
	double precision,parameter,private	:: OHEXPBOND=0.942

! Work variables
	integer,private				:: watmodel_return

! Input variable
	character*10,public			:: wat_potential
	namelist /inpwatmodel/ wat_potential

contains
	! Read the input file for the water potential
	subroutine watmodel_readinput
		implicit none

		open(103,file='pot.inp',status='old',err=390)

		goto 391

390		stop 'READWATMODELINPUT: File pot.inp not found!'

391		read(103,inpwatmodel)

		close(103)		

		if(wat_potential=='tip4p')then
			nparams=11
		elseif(wat_potential=='spc')then
			nparams=10
		else
			stop 'READWATMODELINPUT: Potential type not supported!'
		end if
	end subroutine

	! Write into a file the water potential parameters that can be later used in H2O simulations
	function watmodel_writepot(potfile)
		character*40,intent(in)		:: potfile

		watmodel_return=0

		! Check if, for some weird reason, a parameter became NaN
		do i=1,nparams
			if(isnan(params(i)))then
				watmodel_return=-1

				goto 880
			end if
		end do

		open(454,file=potfile)
		rewind(454)

		! O-H bond length
		if(params(1)<0.9*OHEXPBOND)then
			zeta=zeta+fzeta(params(1),0.95d0*OHEXPBOND)
		elseif(params(1)>1.1*OHEXPBOND)then
			zeta=zeta+fzeta(params(1),1.05d0*OHEXPBOND)
		end if

		write(454,*) 'variable r0 equal ',params(1),' # r0 -- Equilibrium O-H distance'

		! H-O-H angle
		if(params(2)>109.5)then
			zeta=zeta+fzeta(params(2),109.5d0)
		elseif(params(2)<99.5)then
			zeta=zeta+fzeta(params(2),99.5d0)
		end if

		write(454,*) 'variable theta0 equal ',params(2),' # theta0 -- Equilibrium H-O-H angle'

		! LJ epsilon for O-O
		if(params(3)<0.0)then
			watmodel_return=-1

			goto 880
		end if

		write(454,*) 'variable eps_oo equal ',params(3),' # LJ eps_OO'

		! LJ epsilon for H-H
		if(params(4)<0.0)then
			watmodel_return=-1

			goto 880
		end if

		write(454,*) 'variable eps_hh equal ',params(4),' # LJ eps_HH'

		! LJ epsilon for O-H
		if(params(5)<0.0)then
			watmodel_return=-1

			goto 880
		end if

		write(454,*) 'variable eps_oh equal ',params(5),' # LJ eps_OH'

		! LJ sigma for O-O
		if(params(6)<0.0)then
			watmodel_return=-1

			goto 880
		end if

		write(454,*) 'variable sig_oo equal ',params(6),' # LJ sigma_OO'

		! LJ sigma for HH
		if(params(7)<0.0)then
			watmodel_return=-1

			goto 880
		end if

		write(454,*) 'variable sig_hh equal ',params(7),' # LJ sigma_HH'

		! LJ sigma for OH
		if(params(8)<0.0)then
			watmodel_return=-1

			goto 880
		end if

		write(454,*) 'variable sig_oh equal ',params(8),' # LJ sigma_OH'

		! Oxygen charge
		if(params(9)>0.0)then
			watmodel_return=-1

			goto 880
		end if

		write(454,*) 'variable qo equal ',params(9),' # O charge'
		write(454,*) 'variable qh equal ',-1.0*params(9)/2.0,' # H charge'

		! Potential cutoff
		if(params(10)<0.0)then
			watmodel_return=-1

			goto 880
		end if

		write(454,*) 'variable cutoff equal ',params(10),' # Potential cutoff'

		! Distance from O to the additional charge point in TIP4P
		if(wat_potential=='tip4p')then
			if(params(10)<0.0)then
				watmodel_return=-1

				goto 880
			end if

			write(454,*) 'variable r_om ',params(10),' # O-M distance'

                end if

		close(454)

880		watmodel_writepot=watmodel_return
	
		return
	end function
end module
