program mendel
	use globals,only:seed,rand

	implicit none

	integer			:: i
	character*40		:: file1,file2
	double precision	:: chi1,chi2,rnumber,alpha,tmp1,tmp2,sigma
	integer			:: nparams1,nparams2
	character*10		:: mixtype

	write(6,'(A)') '=> File names:'
	read(*,*) file1,file2

	open(10,file=file1,status='old',err=1000)
	open(11,file=file2,status='old',err=1000)
	open(12,file='newparamset.dat')
	rewind(12)

	read(10,*) chi1,nparams1
	read(11,*) chi2,nparams2

	if(.not.(nparams1==nparams2)) stop 'ERROR: Number of parameters is not the same!'

	write(6,'(A)') '=> How do you want do specify the mixing of the potential?'
	write(6,'(A)') '   Choose either "chi" or "user":'
87	read(*,*) mixtype

	if(mixtype=='chi')then
		write(6,'(A)') '=> Seed and sigma for the random number generator:'
		read(*,*) seed,sigma

		alpha=exp(-chi1/sigma)/(exp(-chi1/sigma)+exp(-chi2/sigma))
	elseif(mixtype=='user'
		write(6,'(A)') '=> % of the potential 1 in the new potential:'
		read(*,*) alpha
	else
		write(6,*) '   Option not allowed!'
		
		goto 87
	end if

	write(12,*) 1000000,nparams1

	do i=1,nparams1
		read(10,*) tmp1
		read(11,*) tmp2

		rnumber=rand()

		if(rnumber<=alpha)then
			write(12,*) tmp1
		else
			write(12,*) tmp2
		end if
	end do

	close(10)
	close(11)
	close(12)

	stop 'Done!'

1000	stop 'File not found!'	
end program mendel
