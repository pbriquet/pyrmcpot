module buckingham
	use globals

	public		:: buck_readinput,buck_writepot

! Constants
	integer,parameter,private	:: MAXELEBUCK=2
	integer,parameter,private	:: NBUCKPARAMS=3

! Work variables
	integer,private			:: buck_return

! Input variables
	integer,public			:: buck_nel		! Number of elements in the potential
	character*4,public		:: buck_el(MAXELEBUCK)	! Name of the elements
	real,public			:: buck_elratio		! Number of atoms of element 2 per atom of element 1
	namelist /inpbuck/ buck_nel,buck_el,buck_elratio

contains
	! Read input information for the Buckingham potential
	subroutine buck_readinput
		implicit none

		buck_nel=2
		buck_el(1)='Mg'
		buck_el(2)='O'
		buck_elratio=1.0
	
		open(32,file='pot.inp',status='old',err=440)

		goto 441

440		stop 'BUCK_READINPUT: File pot.inp not found!'

441		read(32,inpbuck)

		close(32)		

		if(buck_nel==1)then
			nparams=NBUCKPARAMS
		elseif(buck_nel==2)then
			nparams=NBUCKPARAMS*3+1
		else
			stop 'BUCK_READINPUT: Number of elements not allowed!'
		end if	
	end subroutine

	! Write a potential file with the parameters of the Buckingham potential
	function buck_writepot(potfile)
		implicit none

		character*40,intent(in)	:: potfile
		integer			:: i
		integer			:: buck_writepot

		buck_return=0

		! Perform a check in the parameter set, looking for unappropriate parameter values
		do i=1,nparams
			if(isnan(params(i)))then
				buck_return=-1

				goto 902
			end if

			if(i<=(nparams-1).and.params(i)<0.0)then
				buck_return=-1

				goto 902
			end if
		end do

		open(36,file=potfile)
		rewind(36)

		! Write down the parameters to a file
		write(36,'(F21.7,A14,2A5)') params(1),' # A for: ',buck_el(1),buck_el(1)
		write(36,'(F21.7,A14,2A5)') params(2),' # rho for: ',buck_el(1),buck_el(1)
		write(36,'(F21.7,A14,2A5)') params(3),' # B for: ',buck_el(1),buck_el(1)

		if(buck_nel==2)then
			write(36,'(F21.7,A14,2A5)') params(4),' # A for: ',buck_el(2),buck_el(2)
			write(36,'(F21.7,A14,2A5)') params(5),' # rho for: ',buck_el(2),buck_el(2)
			write(36,'(F21.7,A14,2A5)') params(6),' # B for: ',buck_el(2),buck_el(2)
			write(36,'(F21.7,A14,2A5)') params(7),' # A for: ',buck_el(1),buck_el(2)
			write(36,'(F21.7,A14,2A5)') params(8),' # rho for: ',buck_el(1),buck_el(2)
			write(36,'(F21.7,A14,2A5)') params(9),' # B for: ',buck_el(1),buck_el(2)
			write(36,'(F21.7,A17,A5)') params(10), ' # Charge for: ',buck_el(1)
			write(36,'(F21.7,A17,A5)') -params(10)/buck_elratio, ' # Charge for: ',buck_el(2)
		end if

		close(36)

902		buck_writepot=buck_return

		return
	end function
end module
