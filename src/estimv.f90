
        
                                
        double precision function estimv(b2,b,omeg01,omeg12,omeg02,tra)
        
        use optim
        use tailles
        use commun,only:nz01,nz12,nz02,opt2,k0,rec,nn
        
        implicit none
        
        integer::i,k,j,ier,ni,l   
        integer,dimension(np)::indx
        double precision,dimension(np,np)::omeg,hes2,h,hessh,y
        double precision,dimension(np*(np+3)/2)::v,v12
        double precision,dimension(np)::b12,b 
        double precision::rl,idmPlLikelihood0,res1,tra,d,res
        double precision,dimension(nz01+2,nz01+2)::omeg01
        double precision,dimension(nz12+2,nz12+2)::omeg12
        double precision,dimension(nz02+2,nz02+2)::omeg02
        double precision,dimension(nn)::b2
        double precision::ca,cb,dd
        double precision,external::idmPlLikelihood,idmPlLikelihood2
        integer::n,istop

        estimv = 0.d0   
        tra = 0.d0 
        indx = 0
        i = 0
        ier = 0
        ni = 0
        l = 0
        k = 0
        ca=0.d0
        cb=0.d0
        dd=0.d0 
        res1 =0.d0
        d=0.d0
        res=0.d0
        omeg = 0.d0
        hes2 = 0.d0
        h = 0.d0
        hessh = 0.d0
        y = 0.d0
    

        n  = nz01+2+nz12+2+nz02+2


       if(rec.eq.1)then 

                k0(1) = dexp(b2(1))
                k0(2) = dexp(b2(2))
                k0(3) = dexp(b2(3)) 
       end if
  
        istop =0        
        call marq98(b,n,ni,v,res,ier,istop,ca,cb,dd,idmPlLikelihood)

        if (istop.ne.1) goto 100

        res1 =idmPlLikelihood0(b,n)

!----------------------------------
        do i=1,n
                b12(i) = b(i)*b(i)
                opt2(i) = 5.d-5
        end do

        call deriva(b12,n,v12,rl,idmPlLikelihood2)

        do i=1,n
                do j=i,n
                        k = (((j-1)*j)/2) +i
                        hes2(i,j)= v12(k)
                end do
        end do
!
        do i = 2,n
                do  j = 1,i-1
                        hes2(i,j)=hes2(j,i)
                end do
        end do

        do i = 1,n
                do j = 1,n
                        h(i,j) = -hes2(i,j)
                        omeg(i,j) = 0.d0
                end do
        end do

        do i=1,nz01+2
                do j=1,nz01+2
                        omeg(i,j)=k0(1)*omeg01(i,j)
                end do
        end do   
        k = 0
        do i=nz01+3,nz01+2+nz02+2
                k = k + 1 
                l = 0
                do j=nz01+3,nz01+2+nz02+2
                        l = l + 1
                        omeg(i,j)=k0(2)*omeg02(k,l)
                end do
        end do   

        k = 0
        do i=nz01+2+nz02+3,nz01+2+nz12+2+nz02+2
                k = k + 1
                l = 0
                do j=nz01+2+nz02+3,nz01+2+nz12+2+nz02+2
                        l = l + 1
                        omeg(i,j)=k0(3)*omeg12(k,l)
                end do
        end do   

        do i = 1,n
                do j = 1,n
                        hessh(i,j)= (h(i,j)-(2.d0*omeg(i,j)))
                end do   
        end do

        np = n
        do i=1,np
                do j=1,n
                        y(i,j)=0.d0
                end do
                y(i,i)=1.d0
        end do

        call ludcmp(hessh,n,indx,d)

        do j=1,n
                call lubksb(hessh,n,indx,y(1,j))
        end do

        tra =0.d0
        do k = 1,n
                do j = 1,n
                        tra = tra + y(k,j)*h(j,k)
                end do
        end do
        estimv = res1 - tra
 !          estimv = res+pe - tra
  !          write(*,*)k0,estimv,res1,tra
100    continue

        return

        end function estimv     
        
        
!========================          ESTIMV         ===================


        double precision function estimvSurv(k00,b,aux,niter,res)

        use optim
        use parameters
        use commun,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m,nz,k0Surv
        
        implicit none

        double precision,dimension((nz+2)*(nz+2+3)/2)::v
        double precision::res,k00,pe,aux,ca,cb,dd
        double precision,external::survPlLikelihood
        double precision,dimension(nz+2)::b 
        double precision,dimension(-2:(nz))::the
        integer::ier,istop,niter,np,i

        np = nz+2 
        estimvSurv = 0.d0
        k0Surv = k00*k00

        call marq98(b,np,niter,v,res,ier,istop,ca,cb,dd,survPlLikelihood)

        if(k0Surv.gt.0.d0)then
                pe = 0.d0
                do i=1,nz+2
                        the(i-3)=(b(i))*(b(i))
!       the01(i-3)=dexp(bh(i))
                end do

                do i=1,nz-1
                        pe = pe+(the(i-3)*the(i-3)*m3m3(i))+(the(i-2) &
                        *the(i-2)*m2m2(i))+(the(i-1)*the(i-1)*m1m1(i))+( &
                        the(i)*the(i)*mmm(i))+(2.d0*the(i-3)*the(i-2)* &
                        m3m2(i))+(2.d0*the(i-3)*the(i-1)*m3m1(i))+(2.d0* &
                        the(i-3)*the(i)*m3m(i))+(2.d0*the(i-2)*the(i-1)* &
                        m2m1(i))+(2.d0*the(i-2)*the(i)*m2m(i))+(2.d0*the(i-1) &
                        *the(i)*m1m(i))
                end do

                pe = -k0Surv*pe

                call test(the,k0Surv,np,aux)

                estimvSurv = - ((res-pe)) - aux
        else

                aux = -np
        endif

        return

        end function estimvSurv


