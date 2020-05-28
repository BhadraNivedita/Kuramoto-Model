	! computes the full lyapunov spectrum of the classical kuramoto model by the method of christiansen and rugh.

	implicit none

	integer, parameter :: nn=10, nm=2*nn+nn*nn
	integer, parameter :: nh=1, tran=9900, iter=100
	double precision   :: xx(nm),yy(nm),dy(nm),le(nm)
	double precision   :: pi, pi2, tau, dtau, beta, beta2, grnd
	integer            :: i, j, k, l, m, n, n1,n2, tran1,itlen
        double precision   :: w(nn), kk

	common w, kk, beta, beta2

	open(unit=1,file='kura_classical.dat')
	open(unit=5,file='kura_classical_le.dat')

	!call sgrnd(2701)
	
	
	!value of parameters

	kk=.50; beta=2.0d0	! maximum LE ~ 0.03
	beta2=2.0d0*beta
	tran1 = tran + 1; itlen = iter - tran
	
	

	pi   = 4.0d0*datan(1.0d0)
	pi2 = 2.0d0*pi
	tau  = 0.0d0
	dtau = 0.1d0
	n1=nn+nn*nn+1
        n2=2*nn+nn*nn

	do k=1,nh		! history

	yy = 0.0d0
	le = 0.0d0

	!Initialization
	do l = 1, nn
	yy(l) = pi*(-1 + 2*(l-1)/(nn-1)) !pi2*(grnd()-0.5d0)
	write(*,*) pi*(-1+2*(l-1)/float(nn-1))
	
	end do
        
	do i = 1, nn
        w(i) = -1 + 2*(i-1)/float(nn-1)
        end do

	do m = 1, nn
	yy(nn*m+m) = 1.0d0
	end do

	do i=1,tran	! transient
	call derivs(tau, yy, dy)
        call rk4(yy, dy, nm, tau, dtau, xx, derivs)
 	yy=xx
 	write(1,'(1024(xf12.6))') (xx(j),j=1,nn)
	write(5,'(1024(xf12.6))') (xx(j)/(dtau*float(i)),j=n1,n2)
	enddo	! transient

	do i=tran1,iter	! iteration
	call derivs(tau, yy, dy)
        call rk4(yy, dy, nm, tau, dtau, xx, derivs)
 	yy=xx
 	write(1,'(1024(xf12.6))') (xx(j),j=1,nn)
	write(5,'(1024(xf12.6))') (xx(j)/(dtau*float(i)),j=n1,n2)
	do j = n1, n2
	le(j) = le(j) + xx(j)/(dtau*float(i))
	end do
	enddo	! iteration

	!write(5,'(1024(xf12.6))') (le(j)/itlen,j=n1,n2)

        enddo   ! history

	end
	
	! calculates the time derivative 

        subroutine derivs(tau, yy, dy)

        integer, parameter :: nn=10, nm=2*nn+nn*nn
	double precision   :: yy(nm), dy(nm)
	double precision   :: jj(nn,nn), e(nn,nn), je(nn,nn)
	double precision   :: eje(nn,nn), ee(nn,nn), ll(nn,nn)
	double precision   :: tau, dtau, beta, beta2, tmp
	integer            :: i, j, k, l, m, n
        double precision   :: w(nn), kk, kknn

	common w, kk, beta, beta2

	kknn = kk/nn

        !Equations to be integrated
	do i=1,nn
	dy(i) = w(i)
	do j=1,nn
	dy(i) = dy(i) + kknn*dsin(yy(j)-yy(i))
	 
	end do
	end do

	jj = 0.0d0

	do i=1,nn
	do j=1,nn
	jj(i,j) = dcos(yy(j)-yy(i))
	end do
	end do

	do i=1,nn
	do j=1,nn
	jj(i,i) = jj(i,i) - dcos(yy(j)-yy(i))
	end do
	end do

	jj = kknn*jj

	!.......................................... do not edit the portion below....................................................

	k=nn+1
        do m=1,nn
        do n=1,nn
        e(n,m) = yy(k)
        k=k+1
        end do
        end do

	do m=1,nn
	do n=1,nn
	je(n,m) = 0.0d0
	do k=1,nn
	je(n,m) = je(n,m) + jj(n,k)*e(k,m)
	end do
	end do
	end do

	do m=1,nn
	do n=1,nn
	eje(m,n) = 0.0d0
	do k=1,nn
	eje(m,n) = eje(m,n) + e(k,m)*je(k,n)
	end do
	end do
	end do

	do m=1,nn
	do n=1,nn
	ee(m,n) = 0.0d0
	do k=1,nn
	ee(m,n) = ee(m,n) + e(k,m)*e(k,n)
	end do
	end do
	end do

	do m=1,nn
        do n=1,nn
        ll(m,n) = eje(m,n) + eje(n,m) + beta2*ee(m,n)
        end do
        ll(m,m) = eje(m,m) + beta*(ee(m,m)-1.0d0)
        end do

	i=nn+1
        do m=1,nn
        do n=1,nn
        tmp=0.0d0
        do k=1,m
        tmp = tmp + e(n,k)*ll(k,m)
        end do
        dy(i) = je(n,m) - tmp
        i=i+1
        end do
        end do
	
	do i = 1, nn
	dy(nn+nn*nn+i) = eje(i,i)
	end do

        return
        end
	
	! integrates by fourth order runge-kutta method
	
	SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)

	INTEGER n,NMAX
	REAL*8 h,x,dydx(n),y(n),yout(n)
	EXTERNAL derivs
	PARAMETER (NMAX=440)
	INTEGER i
	REAL*8 h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)

	hh=h*0.5
	h6=h/6.
	xh=x+hh

	do i=1,n
        yt(i)=y(i)+hh*dydx(i)
	end do

	call derivs(xh,yt,dyt)
	do i=1,n
        yt(i)=y(i)+hh*dyt(i)
	end do

	call derivs(xh,yt,dym)
	do i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
	end do

	call derivs(x+h,yt,dyt)
	do i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
	end do

	return
	END

	!include 'mt.f'
