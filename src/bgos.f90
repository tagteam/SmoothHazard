!=========================== bgos ======================
      subroutine bgos(sx,id,x1,x2,ro)
!      use ifport

! id=1:u(0,sx); id diff de 1 :n(0,sx)
      implicit none 
      double precision ro,sx
      integer id
      double precision f,v1,v2,s,dls,ro2
      double precision x1,x2,uniran,rand
!      write(*,*)'dans bgos'

5     continue

      x1=dble(rand(0)) 
      x2=dble(rand(0))   
      if(id.ne.1) go to 10
      f=2.d0*dsqrt(3.d0)
      x1=(x1-0.5d0)*f
      x2=(x2-0.5d0)*f
      go to 20
10    continue
      v1=2.d0*x1-1.d0
      v2=2.d0*x2-1.d0
      s=v1*v1+v2*v2
      if(s.ge.1.) go to 5
      dls=dsqrt(-2.d0*log(s)/s)
      x1=v1*dls
      x2=v2*dls
20    continue
      ro2=ro*ro
      if(dabs(ro).gt.1.d-10)x2=(x1+x2*dsqrt(1.d0/ro2-1.d0))*ro
      x1=x1*sx
      x2=x2*sx



      return
      end subroutine bgos



