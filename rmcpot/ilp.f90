module ilp
	use globals

	public		:: ilp_readinput,ilp_writepot

! Constants
	integer,private,parameter		:: MAXELEMILP=4

! Work variables
	integer,private				:: ilp_return

! Input variable
	character*20,public			:: ilp_potential
	integer,public				:: ilp_nelements
	character*4,allocatable			:: ilp_elements(MAXELEMILP)
	namelist /inpilp/ ilp_potential

contains
	! Read the input file for the interlayer potential
	subroutine ilp_readinput
		implicit none

		open(117,file='pot.inp',status='old',err=342)

		goto 343

342		stop 'READILPINPUT: File pot.inp not found!'

343		read(117,inpilp)

		close(117)		

		if(ilp_potential=='lebedeva')then
			nparams=10*(ilp_nelements*(ilp_nelements+1)/2)
		elseif(ilp_potential=='kolmogorov')then
			nparams=9*(ilp_nelements*(ilp_nelements+1)/2)
		else
			stop 'READILPINPUT: Potential type not supported!'
		end if
	end subroutine

	! Write into a file the interlayer potential parameters
	function ilp_writepot(potfile)
		character*40,intent(in)		:: potfile

		ilp_return=0

		! Check if, for some weird reason, a parameter became NaN
		do i=1,nparams
			if(isnan(params(i)))then
				ilp_return=-1

				goto 802
			end if
		end do

		open(540,file=potfile)
		rewind(540)

		close(540)

802		ilp_writepot=ilp_return
	
		return
	end function
end module
