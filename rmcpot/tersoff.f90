! Obtain a Tersoff potential in LAMMPS format
module tersoff
	use globals

	private		:: check_params
	public		:: ter_readinput,ter_writepot

! Constants
	integer,parameter,private	:: MAXELETER=3	! A maximum of 3 elements can be considered in the fitting procedure
	integer,parameter,private	:: NTERPAR=13	! Number of Tersoff parameters per interaction

! Work variables
	integer,private			:: ter_return
	double precision		:: par_tmp(NTERPAR)

! Input variables
	integer,public			:: ter_nel=1			! Number of elements in the fitting procedure
	real,public			:: ter_m=3.0			! Value of parameter "m" in the Tersoff potential
	character*2,public		:: ter_el(MAXELETER)		! Elements in the fitting procedure
	logical,public			:: ter_ignore(MAXELETER)	! Elemental interactions will be not checked for consistency
	
	namelist /inptersoff/ ter_nel,ter_el,ter_m,ter_ignore

contains
	! Read the input file for the Tersoff potential
	subroutine ter_readinput
		implicit none

		ter_el(1)='C'
		ter_el(2)='Si'
		ter_el(3)='O'
		ter_ignore(1)=.false.
		ter_ignore(2)=.false.
		ter_ignore(3)=.false.

		open(112,file='pot.inp',status='old',err=332)

		goto 333

332		stop 'READTERSOFFINPUT: File pot.inp not found!'

333		read(112,inptersoff)

		close(112)

		if(ter_nel<1.or.ter_nel>MAXELETER)then
			 stop 'READTERSOFFINPUT: ter_nel is out of bounds!'
		end if

		if(ter_nel==1)then
			nparams=NTERPAR
		elseif(ter_nel==2)then
			nparams=74
		elseif(ter_nel==3)then
			nparams=NTERPAR*27
		end if
	end subroutine
	
	! Write a potential file in the format expected by LAMMPS
	function ter_writepot(potfile)
		implicit none

		integer			:: i
		character*40,intent(in)	:: potfile
		character*20		:: fmt1='(3A3,14E14.5)'
		integer			:: ter_writepot

		ter_return=0

		! Check if, for some weird reason, a parameter became NaN or if it has a value out of the prescribed boundaries
		do i=1,nparams
			if(isnan(params(i)))then
				ter_return=-1

				goto 334
			end if
		end do

		open(46,file=potfile)
		rewind(46)

		! Write the parameters for the interactions involving the first element
		if(.not.(ter_ignore(1)))then
			call check_params(1,13,.false.)
		endif

		par_tmp=params(1:13)		

		write(46,fmt1) ter_el(1),ter_el(1),ter_el(1),ter_m,par_tmp

		if(ter_nel>1)then	! There is, at least, another element in the alloy or compound
			if(.not.(ter_ignore(2)))then		! 2-2-2
				call check_params(14,13,.false.)
			endif

			par_tmp=params(14:26)

			write(46,fmt1) ter_el(2),ter_el(2),ter_el(2),ter_m,par_tmp

			call check_params(27,7,.false.)		! 1-1-2

			par_tmp(1:5)=params(27:31)
			par_tmp(6:9)=0.0d0
			par_tmp(10:11)=params(32:33)
			par_tmp(12:13)=0.0d0

			write(46,fmt1) ter_el(1),ter_el(1),ter_el(2),ter_m,par_tmp

			call check_params(34,7,.false.)		! 1-2-1

			par_tmp(1:5)=params(34:38)
			par_tmp(6:9)=0.0d0
			par_tmp(10:11)=params(39:40)
			par_tmp(12:13)=0.0d0

			write(46,fmt1) ter_el(1),ter_el(2),ter_el(1),ter_m,par_tmp

			call check_params(41,13,.false.) 	! 1-2-2

			par_tmp=params(41:53)

			write(46,fmt1) ter_el(1),ter_el(2),ter_el(2),ter_m,par_tmp

			call check_params(54,7,.true.)		! 2-1-1

			par_tmp(1:7)=params(54:60)

			write(46,fmt1) ter_el(2),ter_el(1),ter_el(1),ter_m,par_tmp

			call check_params(61,7,.false.)		! 2-2-1

			par_tmp(1:5)=params(61:65)
			par_tmp(6:9)=0.0d0
			par_tmp(10:11)=params(66:67)
			par_tmp(12:13)=0.0d0

			write(46,fmt1) ter_el(2),ter_el(2),ter_el(1),ter_m,par_tmp

			call check_params(68,7,.false.)		! 2-1-2

			par_tmp(1:5)=params(68:72)
			par_tmp(6:9)=0.0d0
			par_tmp(10:11)=params(73:74)
			par_tmp(12:13)=0.0d0

			write(46,fmt1) ter_el(2),ter_el(1),ter_el(2),ter_m,par_tmp

			if(ter_nel==3)then	! In case there is a third element and all cross interactions have to be defined
				call check_params(75,13,.false.)	! 3-3-3

				par_tmp=params(75:87)

				write(46,fmt1) ter_el(3),ter_el(3),ter_el(3),ter_m,par_tmp
			end if
		end if

		close(46)

334		ter_writepot=ter_return

		return
	end function

	! Keep the parameters inside acceptable boundaries 
	subroutine check_params(x,y,repeated)
		implicit none

		integer,intent(in)	:: x,y
		integer			:: i
		logical			:: repeated

		do i=x,x+y-1
			if(i==x+4)then	! This is cos(theta)
				if(params(i)<-1.0)then
					params(i)=-1.0
				elseif(params(i)>1.0)then
					params(i)=1.0
				end if
			else
				if(params(i)<0.0d0) params(i)=0.0d0
			end if
		end do

		! Treating D and R (R must be greater than D)
		if(.not.repeated)then	! This conditions is necessary because A, B, D, R, LAMBDA1 and LAMBDA2 for XYY and YXX are the same
			if(y==13)then
				if(params(x+10)>=params(x+9)) params(x+10)=params(x+9)-EPSILON
			elseif(y==7)then
				if(params(x+6)>=params(x+5)) params(x+6)=params(x+5)-EPSILON
			end if
		end if

		! Treating A, B, gamma1 and gamma2 (A > B and gamma1 > gamma2, for instance)
		if(y==13)then
			if(params(x+7)<=10.0*EPSILON) params(x+7)=10.0*EPSILON
			if(params(x+8)<=10.0*EPSILON) params(x+8)=10.0*EPSILON
			if(params(x+11)<=10.0*EPSILON) params(x+11)=10.0*EPSILON
			if(params(x+12)<=10.0*EPSILON) params(x+12)=10.0*EPSILON
			if(params(x+8)>=params(x+12)) params(x+8)=params(x+12)-EPSILON
			if(params(x+7)>=params(x+11)) params(x+7)=params(x+11)-EPSILON
		end if		
	end subroutine
end module
