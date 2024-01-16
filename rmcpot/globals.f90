! Global variables, functions and subroutines
module globals
	public		:: rand,mcmove,testtoosmall,fchi,fzeta,iniparams,saveparamset
	public		:: fsinscr,fquartscr

! Global constants
	real,parameter,public			:: EPSILON=1.0e-6
	real,parameter,public			:: PI=3.14159265359

! Global variables
	integer,public				:: chosen
        integer*8,public	               	:: seed
	integer,public				:: nparams,ndata
	integer,allocatable,public		:: accepted(:)
	real,allocatable,public			:: tolerance(:)
	double precision,public			:: prevparam
	double precision,public			:: minchi
	double precision,public			:: penalty,zeta,minzeta
	double precision,allocatable,public	:: refdata(:),weight(:),potdata(:),bestpotdata(:)
	double precision,allocatable,public	:: params(:),bestparams(:),maxmove(:),maxparval(:),minparval(:)
	logical,allocatable,public		:: qualiweight(:),fixpar(:),hardconst(:)
	character*6,public			:: initialparams
	character*40,public			:: iniparamsfile
	logical,public				:: normweight

! Global functions and subroutines
contains
	! Compute a pseudo-random number
	function rand()
	        implicit none

		double precision,parameter    	:: m=714025.0d0
		double precision,parameter    	:: a=150889.0d0
		double precision,parameter    	:: c=1366.0d0
		double precision,parameter     	:: d=30629.0d0
		double precision,parameter	:: MINRAND=1.0d-10
		double precision               	:: tmp
		double precision		:: rand

		tmp=dble(seed)/d
		rand=tmp-int(tmp)

		if(rand<MINRAND) rand=MINRAND

		seed=mod((a*seed+c),m)
		return
	end function

	! Change a randomly chosen parameter
	function mcmove()
	        implicit none
	
		double precision		:: defsign,rnumber
		integer				:: mcmove

		mcmove=0
		rnumber=rand()
	
		if(rnumber>=0.5)then
			defsign=1.0
		else
			defsign=-1.0
		end if
	
		rnumber=rand()
		prevparam=params(chosen)
		params(chosen)=params(chosen)+defsign*rnumber*maxmove(chosen)

		! I cannot figure out how this can happen, but it seems that sometimes... it happens
		if(isnan(params(chosen)))then 
			mcmove=-1

			return
		elseif(abs(maxparval(chosen))>EPSILON.or.abs(minparval(chosen))>EPSILON)then
			mcmove=checkparambound(chosen)
		end if

		return
	end function

	! Test for (absolute) too small numbers in the parameter set
	subroutine testtoosmall
		implicit none

		integer				:: i
		double precision,parameter	:: MINVALUE=1.0d-10
	
		do i=1,nparams
			if((params(i)<MINVALUE).and.(params(i)>0.0))then 
				params(i)=MINVALUE
			elseif((params(i)>-MINVALUE).and.(params(i)<0.0))then
				params(i)=-MINVALUE
			end if
		end do
	end subroutine

	! Return chi (indeed, chi squared)
	function fchi()
		implicit none

		integer			:: i
		double precision	:: dchi=0.0
		double precision	:: fchi

		fchi=0.0

		do i=1,ndata
			if(tolerance(i)<-EPSILON)then		! Normalize the difference considering the reference value
				dchi=(potdata(i)-refdata(i))/refdata(i)
			elseif(tolerance(i)>EPSILON)then	! Normalize the difference considering the value of tolerance
				dchi=(potdata(i)-refdata(i))/tolerance(i)
			else					! Use the difference between the calculated and reference value as is
				dchi=potdata(i)-refdata(i)
			end if

			fchi=fchi+weight(i)*dchi*dchi
		end do

		return
	end function

	! Compute the squared deviation from boundary values that work as soft constraints during the fitting procedure
	function fzeta(x,y)
		implicit none

		double precision,intent(in)	:: x,y
		double precision		:: dzeta,fzeta

		fzeta=0.0

		if(abs(y)<EPSILON)then
			dzeta=x-y
		else
			dzeta=(x-y)/y
		end if

		fzeta=penalty*dzeta*dzeta

		return
	end function

	! Determine the value of sigma based on the value of chi
	function fsigma(x)
		double precision,intent(in)	:: x
		real,parameter			:: MINSIGMA=0.01,MAXSIGMA=1000000.0
		real				:: upperval
		real				:: fsigma

		fsigma=MINSIGMA

		if(x<=MAXSIGMA)then
			upperval=x
		else
			upperval=MAXSIGMA
		end if

		goto 661

662		fsigma=fsigma*10.0

661		if(fsigma<upperval) goto 662	

		return
	end function

	! Assign initial values to parameters
	subroutine iniparams
                implicit none

		integer				:: i
		double precision		:: defsign
		double precision,parameter	:: MAXVAL=10.0d0
		double precision		:: rnumber

		allocate(params(nparams))
		allocate(bestparams(nparams))
		allocate(maxmove(nparams))
		allocate(maxparval(nparams))
		allocate(minparval(nparams))
		allocate(hardconst(nparams))
		allocate(accepted(nparams))
		allocate(fixpar(nparams))

		select case(initialparams)
			case('random')
				do i=1,nparams	
					rnumber=rand()

					if(rnumber>=0.5)then
						defsign=1.0
					else
						defsign=-1.0
					end if
		
					params(i)=defsign*rnumber*MAXVAL

					! If for some weird reason the parameter is NaN
					if(isnan(params(i))) params(i)=1.0d0
				end do
			case('ones')
				params=1.0d0
			case('file')
				! Try to open the file containing the best parameter set so far, if it exists
				open(10,file='bestparamset.dat',status='old',err=586)

				goto 588

				! Try to open the file containing the initial parameter set, if it exists
586				open(10,file=iniparamsfile,status='old',err=587)

				goto 588
				
587				stop 'INIPARAMS: File not found!'

588				read(10,*)

				do i=1,nparams
					read(10,*) params(i)
				end do

				close(10)
			case default
				stop 'INIPARAMS: This value is not allowed!'
		end select

		bestparams=params	! Initialize the array containing the best parameter set
		accepted=0		! Initialize the array with the number of accepted trial moves per parameter
		maxparval=0.0		! Initialize the array with the maximum value that a parameter can have
		minparval=0.0		! Initialize the array with the minimum value that a parameter can have
		fixpar=.false.		! Initialize the array with the fixed parameters
		hardconst=.false.	! Initialize the array defining if either hard or soft constraints are applied
	end subroutine

	! Save a file containing the parameter set during the fitting procedure
	subroutine saveparamset(chi,isthebest)
		implicit none

		double precision,intent(in)	:: chi
		logical,intent(in)		:: isthebest
		character*16			:: paramfile
		integer				:: i

		if(isthebest)then
			paramfile='bestparamset.dat'
		else
			paramfile='paramset.dat'
		end if

		open(99,file=paramfile)
		rewind(99)

		write(99,*) chi,nparams

		do i=1,nparams
			write(99,*) params(i)
		end do

		close(99)
	end subroutine

	! Apply soft or hard constraints to a parameter, if it is listed in the paramlist2 file
	function checkparambound(x)
		implicit none

		integer				:: checkparambound
		integer,intent(in)		:: x

		checkparambound=0

		if(params(x)>maxparval(x))then
			if(hardconst(x))then
				checkparambound=-1
			else
				zeta=zeta+fzeta(params(x),maxparval(x))
			end if
		elseif(params(x)<minparval(x))then
			if(hardconst(x))then
				checkparambound=-1
			else
				zeta=zeta+fzeta(params(x),minparval(x))
			end if
		end if

		return
	end function

	! Sinus screening function
	function fsinscr(r,r0,d)
		implicit none

		double precision,intent(in)	:: r,r0,d
		double precision		:: fsinscr

		fsinscr=1.0d0

		if(r>=r0)then
			fsinscr=1.0d0-sin(PI/2.0d0*(r-r0)/d)
		end if

		return
	end function

	! Quartic screening function
	function fquartscr(d)
		implicit none

		double precision,intent(in)	:: d
		double precision		:: fquartscr	

		fquartscr=(d*d*d*d)/(1+d*d*d*d)

		return
	end function	
end module
