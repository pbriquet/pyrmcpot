module meam
	use globals

	public		:: meam_readinput,meam_writepot
	private		:: check_meam

! Constants
	integer,parameter,private		:: NLIB=10
	integer,parameter,private		:: MAXELEMEAM=4
	real,parameter,private			:: MEAMRHOZERO=1.0

! Work variables
	integer,private			:: meam_return

! Input variables
	integer,public			:: meam_nel			! Number of element types
	logical,public			:: meam_existinglibrary		! Use an existing MEAM library (e.g., the LAMMPS one)
	character*6,public		:: meam_elt(MAXELEMEAM)		! Element names
	character*4,public		:: meam_lat(MAXELEMEAM)		! Lattice types
	integer,public			:: meam_z(MAXELEMEAM)		! Number of nearest neighbors in the atomic structure
	integer,public			:: meam_anumber(MAXELEMEAM)	! Atomic number of the elements
	real,public			:: meam_aweight(MAXELEMEAM)	! Atomic weight of the elements
	real,public			:: meam_alat(MAXELEMEAM)	! Lattice parameters
	real,public			:: meam_esub(MAXELEMEAM)	! Energies per atom in the reference structures
	integer,public			:: meam_ibar(MAXELEMEAM)	! Define the Gamma function (used to compute the electron density)
	real,public			:: meam_rc			! Cutoff radius
	real,public			:: meam_delr			! Smoothing distance
	real,public			:: meam_ec(MAXELEMEAM-1,MAXELEMEAM-1)	! Cohesive energy odf the reference structure 
	real,public			:: meam_delta(MAXELEMEAM-1,MAXELEMEAM-1)	! Heat of formation of the alloy
	real,public			:: meam_re(MAXELEMEAM-1,MAXELEMEAM-1)	! Distance between atoms in the reference structure
	character*4,public		:: meam_lattce(MAXELEMEAM-1,MAXELEMEAM-1)	! Lattice structure of the reference structure
	real,public			:: meam_nn(MAXELEMEAM,MAXELEMEAM)	! Turn on/off second nearest neighbor interactions
	integer				:: meam_augt1			! Augment or not t1 by 3/5*t3

	namelist /inpmeam/ meam_nel,meam_existinglibrary,meam_elt,meam_lat,meam_z,meam_anumber,meam_aweight,meam_alat, &
	& meam_esub,meam_ibar,meam_rc,meam_delr,meam_ec,meam_delta,meam_re,meam_lattce,meam_nn,meam_augt1

contains
	! Read the input file for the MEAM potential
	subroutine meam_readinput
		implicit none

		! Initialize input variables
		meam_existinglibrary=.false.
		meam_nel=1
		meam_elt='Si'
		meam_lat='dia'
		meam_z=4
		meam_anumber=14
		meam_aweight=28.086
		meam_alat=5.431
		meam_esub=4.63
		meam_ibar=0
		meam_rc=4.0
		meam_delr=0.1
		meam_ec=0.0
		meam_delta=0.0
		meam_re=2.0
		meam_lattce='b1'
		meam_nn=0
		meam_augt1=0

		! Reopen the input file
		open(25,file='pot.inp',status='old',err=537)

		goto 538

537		stop 'MEAM_READINPUT: File pot.inp not found!'

538		read(25,inpmeam)

		close(25)

		! If only one element is considered, a new library will be mandatorily created
		if(meam_nel==1) meam_existinglibrary=.false.

		! Obtain the total number of parameters
!		if(meam_existinglibrary)then
!			nparams=
!		else
!			nparams=
!		end if		
	end subroutine

	! Write the MEAM file(s) to be used with LAMMPS
	function meam_writepot(potfile,isthebest)
		implicit none

		character*40,intent(in)	:: potfile
		logical,intent(in)	:: isthebest
		integer			:: i,shift
		integer			:: meam_writepot

		meam_return=0

		call check_meam

		! Open the custom library file and write the parameters for the single elements
		if(.not.meam_existinglibrary)then 
			if(.not.isthebest)then
				open(38,file='mylibrary.meam')
				rewind(38)
			else
				open(38,file='mybestlibrary.eam')
				rewind(38)
			end if

			do i=1,meam_nel
				shift=(i-1)*NLIB

				write(38,'(2A8,2I5,F10.4)') meam_elt(i),meam_lat(i),meam_z(i),meam_anumber(i),meam_aweight(i)
				write(38,'(8F10.4)') params(1+shift:5+shift),meam_alat(i),meam_esub(i),params(6+shift)
				write(38,'(5F10.4,I3)') params(7+shift:10+shift),MEAMRHOZERO,meam_ibar(i)
			end do

			close(38)
		end if

		open(39,file=potfile)	! Open the parameter file
		rewind(39)

		meam_writepot=meam_return
	end function

	! Chek the parameter values
	subroutine check_meam
		implicit none

				
	end subroutine	
end module meam
