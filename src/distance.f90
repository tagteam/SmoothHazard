!================ Distance pour idmPl

        subroutine distance(nz01,nz02,nz12,b,t,a01,a01_l,a01_u,a02,a02_l,a02_u,a12,a12_l,a12_u)
         
        use tailles

        use commun,only:zi01,zi12,zi02,hessienne

        implicit none
         
        integer::i,j,n,k,l,nz01,nz12,nz02
        integer,dimension(nz01+2+nz12+2+nz02+2)::indx
        double precision::x1,x2,x3,h1,h2,h3,su,bsup,binf,lam,lbinf,lbsup,d
        double precision,dimension(np)::b
        double precision,dimension(np,np)::h,y
         double precision,dimension(-2:(nz01-1))::the01
        double precision,dimension(-2:(nz12-1))::the12
        double precision,dimension(-2:(nz02-1))::the02

        double precision,dimension(nz01+2,nz01+2)::hes01
        double precision,dimension(nz02+2,nz02+2)::hes02 
        double precision,dimension(nz12+2,nz12+2)::hes12

        double precision,dimension(99,3)::t
        double precision,dimension(99)::a01,a01_l,a01_u,a02,a02_l,a02_u,a12,a12_l,a12_u
!----------------------------------------------------------------------------------------

        n = nz01+2+nz12+2+nz02+2
        
        do i = 1,n
                do j = 1,n
                        h(i,j) = - hessienne(i,j)
                end do   
        end do

        np = n

        do i=1,n
                do j=1,n
                        y(i,j)=0.d0
                end do
                y(i,i)=1.d0
        end do
        
        call ludcmp(h,n,indx,d)

        do j=1,n
                call lubksb(h,n,indx,y(1,j))
        end do

        do i=1,nz01+2
                do j=1,nz01+2
                        hes01(i,j)=y(i,j)
                end do
        end do   

        k = 0
        do i=nz01+3,nz01+2+nz02+2
                k = k + 1 
                l = 0
                do j=nz01+3,nz01+2+nz02+2
                        l = l + 1
                        hes02(k,l)=y(i,j)
                end do
        end do   


        k = 0
        do i=nz01+2+nz02+3,nz01+2+nz12+2+nz02+2
                k = k + 1
                l = 0
                do j=nz01+2+nz02+3,nz01+2+nz12+2+nz02+2
                        l = l + 1
                        hes12(k,l)=y(i,j)
                end do
        end do   


        do i=1,nz01+2
                the01(i-3)=(b(i))*(b(i))
!               the01(i-3)=dexp(b(i))
        end do
         
        do i=1,nz02+2
                j = nz01+2+i
                the02(i-3)=(b(j))*(b(j))
!               the02(i-3)=dexp(b(j))
        end do
        
        do i=1,nz12+2
                j = nz01+2+nz02+2+i
                the12(i-3)=(b(j))*(b(j))
!               the12(i-3)=dexp(b(j))
        end do

        h1 = (zi01(nz01)-zi01(1))*0.01d0
        x1 = zi01(1)

        h2 = (zi12(nz12)-zi12(1))*0.01d0
        x2 = zi12(1)

        h3 = (zi02(nz02)-zi02(1))*0.01d0
        x3 = zi02(1)

        do i=1,99
                x1 = x1 + h1
                t(i,1)=x1
                call cosp(x1,the01,nz01+2,hes01,zi01,binf,su,bsup,lbinf,lam,lbsup)
! 0------->1
                if(binf.lt.0.d0)then
                        binf = 0.d0
                endif
                if(bsup.gt.1.d0)then
                        bsup = 1.d0
                endif
                if(lbinf.lt.0.d0)then
                        lbinf = 0.d0
                endif
                
                a01(i)=lam
                a01_l(i)=lbinf
                a01_u(i)=lbsup
! 1------->2
                x2 = x2 + h2
                t(i,2)=x2       
                call cosp(x2,the12,nz12+2,hes12,zi12,binf,su,bsup,lbinf,lam,lbsup)
                
                if(binf.lt.0.d0)then
                        binf = 0.d0
                endif
                if(bsup.gt.1.d0)then
                        bsup = 1.d0
                endif
                if(lbinf.lt.0.d0)then
                        lbinf = 0.d0
                endif   
                
                a12(i)=lam
                a12_l(i)=lbinf
                a12_u(i)=lbsup   
                     
! 0------->2            
                x3 = x3 + h3
                t(i,3) = x3
                call cosp(x3,the02,nz02+2,hes02,zi02,binf,su,bsup,lbinf,lam,lbsup)

                if(binf.lt.0.d0)then
                        binf = 0.d0
                endif
                if(bsup.gt.1.d0)then
                        bsup = 1.d0
                endif
                if(lbinf.lt.0.d0)then
                        lbinf = 0.d0
                endif

                a02(i) = lam
                a02_l(i) = lbinf
                a02_u(i) = lbsup         
        end do

        end subroutine distance


!==========================  DISTANCE survie   =================================
        subroutine distanceSurv(v,b,t,S,S_l,S_u,h,h_l,h_u)
         
        use tailles
        use optim
        use commun,only:nz,zi,k0Surv

         implicit none
         
        integer::i,j,n,k,kkk
        integer,dimension(nz+2)::indx
        double precision::x,pas,su,bsup,binf,lam,lbinf,lbsup,d
        double precision,dimension(np)::b
        double precision,dimension(nz+2,nz+2)::hes,hh,y,omeg
        double precision,dimension(-2:(nz))::the
        double precision,dimension(99)::t
        double precision,dimension(99)::S,S_l,S_u,h,h_l,h_u
        double precision,dimension(np*(np+3)/2)::v
        double precision,external::survPlLikelihood
        double precision::rl


        n  = nz+2

         do i=1,nz+2
            the(i-3)=(b(i))*(b(i))
!        the(i-3)=dexp(b(i))
         end do

!        write(*,*)'DISTANCE',n,np      
!------------------------------  Version 1 -----------------------------------

        call deriva(b,np,v,rl,survPlLikelihood)

 !        write(*,*)'DISTANCE 0',v(1),rl

        kkk=np*(np+1)/2
        
        do i=1,n
                kkk = kkk+1
                do j=i,n
                        k = (((j-1)*j)/2) +i
                        hes(i,j)= v(k)/(4.d0*b(i)*b(j))
                end do
                hes(i,i)= hes(i,i)+(v(kkk)/(4.d0*b(i)*b(i)*b(i)))
        end do
        
        do i = 2,n
                do j = 1,i-1
                        hes(i,j)=hes(j,i)
                end do
        end do

        do i = 1,n
                do j = 1,n
                        hh(i,j) = -hes(i,j)
                end do   
        end do
!         write(*,*)'DISTANCE hh',hh(1,1)

        do i=1,n
                do j=1,n
                        y(i,j)=0.d0
                end do
                y(i,i)=1.d0
        end do

        call ludcmp(hh,n,indx,d)

        do j=1,n
                call lubksb(hh,n,indx,y(1,j))
        end do

!------------------------------  Version 2 -----------------------------------

        do i = 1,n
                do j = 1,n
                        hes(i,j) = 0.d0 
                end do
        end do
   
        do i = 1,n
                do j = i,n
                        call mat(hes(i,j),the,i,j,n)
                end do
        end do

        do i = 2,n
                do j = 1,i-1
                        hes(i,j)=hes(j,i)
                end do
        end do

        call calcomegSurv(n,omeg)

        do i = 1,n
                do j = 1,n
                        hes(i,j) = hes(i,j) - (2.d0*k0Surv*omeg(i,j))  
                end do   
        end do

        do i=1,n
                do j=1,n
                        y(i,j)=0.d0
                end do
                y(i,i)=1.d0
        end do

        call ludcmp(hes,n,indx,d)

        do j=1,np
                call lubksb(hes,n,indx,y(1,j))
        end do

!------------------------------ fin  Version 2 -----------------------------------

        do i=1,nz+2
                do j=1,nz+2
                        hes(i,j)=y(i,j)
                end do
        end do   
        
        pas = (zi(nz)-zi(1))*0.01d0
        x = zi(1)
!       write(*,*)'DISTANCE 2',pas,x
        do i=1,99
                x = x + pas
                t(i)=x 
                call cosp(x,the,nz+2,hes,zi,binf,su,bsup,lbinf,lam,lbsup)
                if(binf.lt.0.d0)then
                        binf = 0.d0
                endif
                if(bsup.gt.1.d0)then
                        bsup = 1.d0
                endif
                if(lbinf.lt.0.d0)then
                        lbinf = 0.d0
                endif
                h(i)=lam
                h_l(i)=lbinf
                h_u(i)=lbsup          
                S(i)=su
                S_l(i)=binf 
                S_u(i)=bsup
        end do

        end subroutine distanceSurv


!==========================  COSP  ====================================
        subroutine cosp(x,the,n,y,zi,binf,su,bsup,lbinf,lam,lbsup)
       
        implicit none
      
        integer::k,j,n,i
        double precision::x,ht,ht2,h2,som,lam,su, &
        binf,bsup,lbinf,lbsup,pm,htm,h2t,h3,h2n,hn,im,im1,im2,mm1,mm3, & 
        ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2,h,gl,hh        
        double precision,dimension(-2:n-3)::the
        double precision,dimension(-2:(n+1))::zi
        double precision,dimension(n,n)::y
        
        ht3=0.d0
        hht=0.d0
        h4=0.d0
        h3m=0.d0
        hh3=0.d0
        hh2=0.d0
        mm=0.d0
        im3=0.d0
        mm2=0.d0
        h=0.d0
        gl=0.d0
        hh=0.d0
        som = 0.d0
        j=1
        
         do k = 2,n-2
            if ((x.ge.zi(k-1)).and.(x.lt.zi(k)))then
               j = k-1
               if (j.gt.1)then
                  do i=2,j
                     som = som+the(i-4)
                  end do  
               endif 
            ht = x-zi(j)
            htm= x-zi(j-1)
            h2t= x-zi(j+2)
            ht2 = zi(j+1)-x
            ht3 = zi(j+3)-x
            hht = x-zi(j-2)
            h = zi(j+1)-zi(j)
            hh= zi(j+1)-zi(j-1)
            h2= zi(j+2)-zi(j)
            h3= zi(j+3)-zi(j)
            h4= zi(j+4)-zi(j)
            h3m= zi(j+3)-zi(j-1)
            h2n=zi(j+2)-zi(j-1)
            hn= zi(j+1)-zi(j-2)
            hh3 = zi(j+1)-zi(j-3)
            hh2 = zi(j+2)-zi(j-2)
            mm3 = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
            mm2 = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
            *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
            mm1 = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
            h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
            mm  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
            im3 = (0.25d0*(x-zi(j-3))*mm3)+(0.25d0*hh2*mm2) &
            +(0.25d0*h3m*mm1)+(0.25d0*h4*mm)
            im2 = (0.25d0*hht*mm2)+(h3m*mm1*0.25d0)+(h4*mm*0.25d0)
            im1 = (htm*mm1*0.25d0)+(h4*mm*0.25d0)
            im  = ht*mm*0.25d0
            gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1) &
            +(the(j)*im)
            lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1) &
            +(the(j)*mm)
            end if
            end do
            
            if(x.ge.zi(n))then
                        som = 0.d0
                        do i=1,n
                        som = som+the(i-3)
                        end do
                        j = n
                        gl = som
            endif

         
   
         call conf(x,j,n,y,pm,zi)

         binf = dexp(-gl + 1.96d0*pm)
         su  = dexp(-gl)
         bsup = dexp(-gl - 1.96d0*pm)

         call conf1(x,j,n,y,pm,zi)
         lbinf = lam - 1.96d0*pm
         lbsup = lam + 1.96d0*pm
	
         return

         end subroutine cosp
!=====================  CONF1  =============================
        subroutine conf1(x,ni,n,y,pm,zi)
        
        implicit none
 
         integer::ni,i,n,j
         double precision::mmsp,x,pm,res
         double precision,dimension(n,n)::y
         double precision,dimension(n)::aux,vecti
         double precision,dimension(-2:(n+1)):: zi
           
            do i=1,n
               vecti(i) = mmsp(x,ni,i,zi,n)
            end do   

            do i=1,n
               aux(i) = 0.d0
               do j=1,n
                  aux(i) = aux(i) - y(i,j)*vecti(j)
               end do
            end do   

            res = 0.d0
            do i=1,n
               res = res + aux(i)*vecti(i)
            end do

            pm = dsqrt(res)
               
         end subroutine conf1

!=====================  CONF  =============================

      subroutine conf(x,ni,n,y,pm,zi)

         implicit none
         
         integer::ni,i,n,j
         double precision::isp,x,pm,res
         double precision,dimension(n,n)::y
         double precision,dimension(n)::aux,vecti
         double precision,dimension(-2:(n+1)):: zi


            do i=1,n
               vecti(i) = isp(x,ni,i,zi,n)
            end do   

            do i=1,n
               aux(i) = 0.d0
               do j=1,n
                  aux(i) = aux(i) - y(i,j)*vecti(j)
               end do
            end do   

            res = 0.d0
            do i=1,n
               res = res + aux(i)*vecti(i)
            end do

            pm = dsqrt(res)

         end subroutine conf



