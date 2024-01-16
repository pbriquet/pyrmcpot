! Define subroutines that perform potential-related tasks
module potfunc
	use globals
	use eamalloy
	use meam
	use tersoff
	use watmodel
	use buckingham

	public		:: inputpot,genpot,savebestpot

contains
	! Initialize variables that are specific for a potential model
	subroutine inputpot(potmodel)
		implicit none

		character*20,intent(in)		:: potmodel

		select case(potmodel)
			case('eam/alloy')		! Input variables for an EAM/alloy (or FS) potential
				call eam_readinput
			case('meam')			! Input variables for an MEAM potential
				call meam_readinput	
			case('tersoff')			! Input variables for a Tersoff potential
				call ter_readinput
			case('watmodel')		! Input variables for a water potential
				call watmodel_readinput	
			case('buckingham')		! Input variables for a Buckingham potential
				call buck_readinput
			case default
				stop 'GENPOT: Potential model not allowed!'
		endselect
	end subroutine

	! Generate the specified potential model
	function genpot(potmodel,potfile)
		implicit none

		character*20,intent(in)		:: potmodel
		character*40,intent(in)		:: potfile
		integer				:: genpot

		genpot=0

		select case(potmodel)
			case('eam/alloy')	! Generate a potential of EAM/alloy type
				genpot=eam_writepot(potfile,.false.)
			case('meam')		! Generate a potential of MEAM type (the library file or the parameter file or both) 
				genpot=meam_writepot(potfile,.false.)
			case('tersoff')		! Generate a potential of Tersoff type
				genpot=ter_writepot(potfile)
			case('watmodel')		! Generate a potential for water
				genpot=watmodel_writepot(potfile)
			case('buckingham')	! Generate a potential of Buckingham type
				genpot=buck_writepot(potfile)
			case default
				stop 'GENPOT: Potential model not allowed!'
		end select

		return
	end function

	! Save separately the best potential generated so far
	function savebestpot(potmodel)
		implicit none

		character*20,intent(in)		:: potmodel
		character*40			:: bestpotfile
		integer				:: savebestpot

		savebestpot=0

		! Write down a file containing the best interatomic potential
		select case(potmodel)
			case('eam/alloy')	! Save the best EAM/alloy potential so far
				bestpotfile='bestpot.eam'			

				savebestpot=eam_writepot(bestpotfile,.true.)
			case('meam')		! Save the best parameter file for alloys and, optionally, the best library file
				bestpotfile='bestparamfile.meam'

				savebestpot=meam_writepot(bestpotfile,.true.)
			case('tersoff')		! Save the best Tersoff potential so far
				bestpotfile='bestpot.tersoff'

				savebestpot=ter_writepot(bestpotfile)
			case('watmodel')		! Save the best water model potential so far
				bestpotfile='bestpot.watmodel'

				savebestpot=watmodel_writepot(bestpotfile)
			case('buckingham')
				bestpotfile='bestpot.buck'

				savebestpot=buck_writepot(bestpotfile)
			case default
				stop 'SAVEBESTPOT: Potential model not implemented yet!'
		endselect
	end function
end module
