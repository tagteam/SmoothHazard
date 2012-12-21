


!========================          MNBRAK         ===================
        subroutine mnbrak(ax,bx,cx,fa,fb,fc,b)

        use commun,only:nz

        implicit none

        double precision::ax,bx,cx,fa,fb,fc,aux,res
        double precision,dimension(nz+2)::b 
        double precision,parameter::gold=1.618034d0,glimit=100.d0,tiny=1.d-20
        double precision::estimvSurv,dum,fu,q,r,u,ulim
        integer::ni

        fa = estimvSurv(ax,b,aux,ni,res)
        fb = estimvSurv(bx,b,aux,ni,res)

        if(fb.gt.fa)then
                dum = ax
                ax = bx
                bx = dum
                dum = fb
                fb = fa
                fa = dum
        endif
        cx = bx + gold*(bx-ax)
        fc = estimvSurv(cx,b,aux,ni,res)
 1      if(fb.ge.fc)then
                r = (bx-ax)*(fb-fc)
                q = (bx-cx)*(fb-fa)
                u = bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(dabs(q-r),tiny),q-r))
                ulim = bx + glimit*(cx-bx)
                if((bx-u)*(u-cx).gt.0.d0)then
                        fu = estimvSurv(u,b,aux,ni,res)
                        if(fu.lt.fc)then
                                ax = bx
                                fa = fb
                                bx = u
                                fb = fu
                                return
                        else
                                if(fu.gt.fb)then
                                        cx = u
                                        fc = fu
                                        return
                        endif   
                endif
                u = cx + gold*(cx-bx)
                fu = estimvSurv(u,b,aux,ni,res)
        else
                if((cx-u)*(u-ulim).gt.0.d0)then
                        fu = estimvSurv(u,b,aux,ni,res)
                        if(fu.lt.fc)then
                                bx = cx
                                cx = u
                                u = cx + gold*(cx-bx)
                                fb = fc
                                fc = fu
                                fu = estimvSurv(u,b,aux,ni,res)
                        endif  
                 else
                        if((u-ulim)*(ulim-cx).ge.0.d0)then
                                u = ulim
                                fu = estimvSurv(u,b,aux,ni,res)
                        else
                                u = cx + gold*(cx-bx)
                                fu = estimvSurv(u,b,aux,ni,res)
                        endif
                endif   
        endif
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu
        goto 1
        endif
        return 
        end subroutine mnbrak

!========================      GOLDEN   =========================
        double precision function golden(ax,bx,cx,tol,xmin,b,aux)

        use commun,only:nz

        implicit none

        double precision,dimension(nz+2)::b 
        double precision::ax,bx,cx,tol,xmin
        double precision::aux,res
        double precision,parameter::r=0.61803399d0,c=1.d0-r
        double precision::f1,f2,x0,x1,x2,x3,estimvSurv
        integer ni

        x0 = ax
        x3 = cx
        if(dabs(cx-bx).gt.dabs(bx-ax))then
                x1 = bx
                x2 = bx + c*(cx-bx)
        else
                x2 = bx
                x1 = bx - c*(bx-ax)
        endif
        f1 = estimvSurv(x1,b,aux,ni,res)
        f2 = estimvSurv(x2,b,aux,ni,res)
         
 1      if(dabs(x3-x0).gt.tol*(dabs(x1)+dabs(x2)))then
                if(f2.lt.f1)then
                        x0 = x1
                        x1 = x2
                        x2 = r*x1 + c*x3
                        f1 = f2
                        f2 = estimvSurv(x2,b,aux,ni,res)
                else
                        x3 = x2
                        x2 = x1
                        x1 = r*x2+c*x0
                        f2 = f1
                        f1 = estimvSurv(x1,b,aux,ni,res)
                endif
        go to 1
        endif
        if(f1.lt.f2)then
                golden = f1
                xmin = x1
        else
                golden = f2
                xmin = x2
        endif
        return
        end function golden


!=================  calcul de la hessienne  et de omega  ==============
        subroutine test(the,k0,np,aux)
        
        implicit none

        integer::i,j,np
        double precision,dimension(np,np)::omeg,hess,hessh,y
        integer,dimension(np)::indx
        double precision::k0,tra,aux,d
        double precision,dimension(-2:np-2)::the


        do i = 1,np
                do j = 1,np
                        hess(i,j) = 0.d0 
                end do
        end do


        do i = 1,np
                do j = i,np
                        call mat(hess(i,j),the,i,j,np)
                end do
        end do

        do i = 2,np
                do j = 1,i-1
                        hess(i,j)=hess(j,i)
                end do
        end do

        call calcomegSurv(np,omeg)

        do i = 1,np
                do j = 1,np
                        hessh(i,j)=-hess(i,j)
                        hess(i,j) = hess(i,j) - (2.d0*k0*omeg(i,j))  
                end do   
        end do


        do i=1,np
                do j=1,np
                        y(i,j)=0.d0
                end do
                y(i,i)=1.d0
        end do

        call ludcmp(hess,np,indx,d)

        do j=1,np
                call lubksb(hess,np,indx,y(1,j))
        end do

        tra = 0.d0
        do i=1,np
                do j=1,np
                        tra = tra + y(i,j)*hessh(j,i)
                end do
        end do

        aux = (tra)

        end subroutine test


!====================  MAT  ==================================
        subroutine mat(res,the,k,l,np)
         
        use commun,only:zi,c,no,nz,t1,t2

        implicit none

        integer::i,k,l,j,ni,np,ni1
        double precision::res,mmsp,isp,u2,res1,aux2,su1,su2,ri
        double precision,dimension(-2:(nz))::the



!---------- calcul de la hessienne ij ------------------
        res = 0.d0
        res1 = 0.d0
        do 10 i=1,no
                if(c(i).eq.1)then
                        call susp(t1(i),the,nz,su1,ri,zi) 
                        u2 = ri
                        do 6 j = 2,np-2
                                if((t1(i).ge.zi(j-1)).and.(t1(i).lt.zi(j)))then
                                        ni = j-1
                                endif
 6                      continue 
                        if(t2(i).eq.zi(np-2))then
                                ni = np-2
                        endif   
!-------  attention numero spline 
                        aux2 = mmsp(t1(i),ni,k,zi,np)*mmsp(t1(i),ni,l,zi,np)
                        if (u2.le.0.d0)then
                                res1 = 0.d0
                        else   
                                res1 = - aux2/(u2*u2)
                        endif 
                else   
                        if(c(i).eq.2)then
                                do 7 j = 2,np-2
                                        if((t1(i).ge.zi(j-1)).and.(t1(i).lt.zi(j)))then
                                                ni1 = j-1
                                        endif
 7                              continue 
                                do 8 j = 2,np-2
                                        if((t2(i).ge.zi(j-1)).and.(t2(i).lt.zi(j)))then
                                                ni = j-1
                                        endif
 8                              continue 
                                if(t2(i).eq.zi(np-2))then
                                        ni = np-2
                                endif   
                                call susp(t1(i),the,nz,su1,ri,zi) 
                                call susp(t2(i),the,nz,su2,ri,zi) 
                                u2 = su1-su2
                                aux2 = -isp(t1(i),ni1,k,zi,np)*isp(t1(i),ni1,l,zi,np)
                                aux2 = aux2 + isp(t1(i),ni1,k,zi,np)*isp(t2(i),ni,l,zi,np)
                                aux2 = aux2 - isp(t2(i),ni,l,zi,np)*isp(t2(i),ni,k,zi,np)
                                aux2 = aux2 + isp(t2(i),ni,k,zi,np)*isp(t1(i),ni1,l,zi,np)
                                aux2 = aux2*su1*su2 
                                res1 = aux2/(u2*u2)
                        else 
                                res1 = 0.d0
                        endif   
                endif
                res = res + res1
 10     continue   

        end subroutine mat



!=======================  CALOMEG  ===========================
        subroutine calcomegSurv(n,omeg)

        use commun,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m
        implicit none
        !        remplissage de la matrice omega n*n
        !          elle a 7 diagonales
        integer::n,i,j
        double precision::calc00,calc01,calc02
        double precision,dimension(n,n)::omeg
        

        
        do i=1,n
                do j=1,n
                        omeg(i,j)=0.d0
                end do
        end do
        
        omeg(1,1)=calc00(1,n,m3m3,m2m2,m1m1,mmm)
        omeg(1,2)=calc01(1,n,m3m2,m2m1,m1m)
        omeg(1,3)=calc02(1,n,m3m1,m2m)
        omeg(1,4)=m3m(1)
        omeg(2,1)=omeg(1,2)
        omeg(2,2)=calc00(2,n,m3m3,m2m2,m1m1,mmm)
        omeg(2,3)=calc01(2,n,m3m2,m2m1,m1m)
        omeg(2,4)=calc02(2,n,m3m1,m2m)
        omeg(2,5)=m3m(2)
        omeg(3,1)=omeg(1,3)
        omeg(3,2)=omeg(2,3)
        omeg(3,3)=calc00(3,n,m3m3,m2m2,m1m1,mmm)
        omeg(3,4)=calc01(3,n,m3m2,m2m1,m1m)
        omeg(3,5)=calc02(3,n,m3m1,m2m)
        omeg(3,6)=m3m(3)
        do i=4,n-3
                omeg(i,i-3)=omeg(i-3,i)
                omeg(i,i-2)=omeg(i-2,i)
                omeg(i,i-1)=omeg(i-1,i)
                omeg(i,i)=calc00(i,n,m3m3,m2m2,m1m1,mmm)
                omeg(i,i+1)=calc01(i,n,m3m2,m2m1,m1m)
                omeg(i,i+2)=calc02(i,n,m3m1,m2m)
                omeg(i,i+3)=m3m(i)
        end do 
        omeg(n-2,n-5)=omeg(n-5,n-2)
        omeg(n-2,n-4)=omeg(n-4,n-2)
        omeg(n-2,n-3)=omeg(n-3,n-2)
        omeg(n-2,n-2)=calc00(n-2,n,m3m3,m2m2,m1m1,mmm)
        omeg(n-2,n-1)=calc01(n-2,n,m3m2,m2m1,m1m)
        omeg(n-2,n)=calc02(n-2,n,m3m1,m2m)
        omeg(n-1,n-4)=omeg(n-4,n-1)
        omeg(n-1,n-3)=omeg(n-3,n-1)
        omeg(n-1,n-2)=omeg(n-2,n-1)
        omeg(n-1,n-1)=calc00(n-1,n,m3m3,m2m2,m1m1,mmm)
        omeg(n-1,n)=calc01(n-1,n,m3m2,m2m1,m1m)
        omeg(n,n-3)=omeg(n-3,n)
        omeg(n,n-2)=omeg(n-2,n)
        omeg(n,n-1)=omeg(n-1,n)
        omeg(n,n)=calc00(n,n,m3m3,m2m2,m1m1,mmm)
        
        end subroutine calcomegSurv


