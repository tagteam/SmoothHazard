!=============================================================================================
!            Illness-death model with non censored or interval-censored or right censored
!                               left truncated
!                                 Weibull
!                              19/04/11
!            last              20/04/11
!=============================================================================================
        subroutine idmWeib(entrytime,l,r,d,idm,idd,x01,x02,x12,N,P01,P02,P12,truncated,eps &
                      ,maxiter0,loglik,basepar,regpar,v,converged,cv,&
                      niter,t,a01,a01_l,a01_u,a02,a02_l,a02_u,a12,a12_l,a12_u,prt,hess_tot) 
      
        use idmCommun  
        use parameters
        use optim
        use commun,only:pl
                
        implicit none
!  variables entrants
        integer,intent(in)::n,p01,p02,p12,prt
        double precision,dimension(n),intent(in)::l,r,d,entrytime
        integer,dimension(n),intent(in)::idm,idd
        double precision,dimension(n,p01),intent(in)::x01
        double precision,dimension(n,p02),intent(in)::x02
        double precision,dimension(n,p12),intent(in)::x12
        integer,intent(in)::truncated,maxiter0
        integer,dimension(3),intent(in)::eps
!  variables sortant
        integer,intent(out)::converged,niter
! istop 
        double precision,dimension(p01+p02+p12),intent(out)::regpar
        double precision,dimension(99),intent(out)::t,a01,a01_l,a01_u,a02,a02_l,a02_u,a12,a12_l,a12_u
        double precision,dimension(2),intent(out)::loglik
        double precision,dimension(6),intent(out)::basepar
        double precision,dimension(3),intent(out)::cv
        double precision,dimension((p01+p02+p12)*(p01+p02+p12)),intent(out)::v   
!  variables locales           
        double precision:: res,min,max,x1,x2
        double precision,dimension((p01+p02+p12)+6)::xi,ut,bh
        double precision,dimension(2)::the01,the02,the12
        double precision,dimension(((p01+p02+p12)+6),(p01+p02+p12)+6)::vinf,vsup
        double precision::su,gl,ri,tx,pas,ep,moyenne
        double precision,dimension(2000,99)::mate_ri01,mate_ri02,mate_ri12
        double precision,dimension(2000)::vect     
        double precision,dimension((p01+p02+p12)+6)::b
        double precision,dimension(((p01+p02+p12)+6)*((p01+p02+p12)+6+3)/2)::v1
        double precision,dimension(((p01+p02+p12)+6),((p01+p02+p12)+6))::hes
        double precision,dimension(((p01+p02+p12)+6),((p01+p02+p12)+6)),intent(out)::hess_tot
        double precision::ca,cb,dd,dsqrt
        integer::ier,npw,np,j,k,jj,i,istop
        double precision,external::idmLikelihood     
        integer::som_idd,som_idm
        double precision::ts
        integer,dimension(3)::noVar
!---------------------------------      

        if(P01.gt.0)then
        noVar(1)=0
        else
        noVar(1)=1
        end if
        if(P02.gt.0)then
        noVar(2)=0
        else
        noVar(2)=1
        end if
        if(P12.gt.0)then
        noVar(3)=0
        else
        noVar(3)=1
        end if
        
        pl=0
        loglik=0.d0       
        print_iter = prt
        no=N
        
        nva01=0 
        nva02=0
        nva12=0
        
        if(noVar(1).ne.1)nva01=P01
        if(noVar(2).ne.1)nva02=P02
        if(noVar(3).ne.1)nva12=P12
        
        troncature = truncated
        epsa = 0.1d0**eps(1)
        epsb = 0.1d0**eps(2)
        epsd = 0.1d0**eps(3)
        maxiter = maxiter0
        min = l(1)
        max = l(1)
        som_idd = 0
        som_idm = 0     
        ts = 0.d0
        if(truncated.eq.1) then
                allocate(t0(no),t1(no),t2(no),t3(no),ve01(no,nva01),ve02(no,nva02),ve12(no,nva12),c(no))
                do i=1,no
                        som_idd = som_idd + idd(i)
                        som_idm = som_idm + idm(i)
                        t0(i)=entrytime(i)
                        t1(i)=l(i)
                        t2(i)=r(i)
                        t3(i)=d(i)
                        ts = ts + (t3(i) - t0(i))
                        if (max.lt.t3(i))max = t3(i)
                        if (min.gt.t0(i))min = t0(i)
                end do 
        else
                allocate(t0(1),t1(no),t2(no),t3(no),ve01(no,nva01),ve02(no,nva02),ve12(no,nva12),c(no))
                min = 0.d0
                do i=1,no
                        som_idd = som_idd + idd(i)
                        som_idm = som_idm + idm(i)
                        t1(i)=l(i)
                        t2(i)=r(i)
                        t3(i)=d(i)
                        ts = ts + t3(i)
                        if (max.lt.t3(i))max = t3(i)
                end do 
        endif 
  
        do i=1,no
                if((idm(i).eq.0).and.(idd(i).eq.0).and.(t1(i).eq.t3(i)))then
!           censure à droite 01 et 02
                        c(i) = 1
                endif
                if((idm(i).eq.1).and.(idd(i).eq.0).and.(t1(i).lt.t2(i)))then
!           censure par intervalle 01 et droite pour 12
                        c(i) = 2
                endif 
                if((idm(i).eq.1).and.(idd(i).eq.0).and.(t1(i).eq.t2(i)))then
!           observation 01 et censure à droite pour 12
                        c(i) = 3
                endif
                if((idd(i).eq.1).and.(idm(i).eq.1).and.(t1(i).lt.t2(i)))then
!           censure par intervalle 01 et observation pour 12
                        c(i) = 4
                endif 
                if((idd(i).eq.1).and.(idm(i).eq.1).and.(t1(i).eq.t2(i)))then
!           observation 01 et observation pour 12
                        c(i) = 5
                endif 
                if((idm(i).eq.0).and.(idd(i).eq.0).and.(t1(i).lt.t3(i)))then
!           vivant
                        c(i) = 6
!               write(*,*)i,t0(i),t1(i),t2(i),t3(i)
                endif
                if((idm(i).eq.0).and.(idd(i).eq.1))then
!           mort
                        c(i) = 7
                endif
                
                do j=1,nva01
                        ve01(i,j)=x01(i,j)
                end do   
                do j=1,nva02
                        ve02(i,j)=x02(i,j)
                end do  
                do j=1,nva12
                        ve12(i,j)=x12(i,j)
                end do  
        end do

        npw = 6

        b(1) = 1.d0
        b(2) = dsqrt(dble(som_idm)/ts)
        b(3) = 1.d0
        b(4) = dsqrt(dble(som_idd)/ts)  
        b(5) = 1.d0     
        b(6) = b(4)

        nva01=0
        nva02=0
        nva12=0
        
        np = npw
        
        call marq98(b,np,niter,v1,res,ier,istop,ca,cb,dd,idmLikelihood)
        if (istop.ne.1) goto 1000

        loglik(1) = res
!            write(*,*)'niter1',niter,res


        if((noVar(1)==0).or.(noVar(2)==0).or.(noVar(3)==0)) then

                nva01=P01
                nva02=P02
                nva12=P12
                
                np = npw + nva01+nva02+nva12
                do i=npw+1,np
                        b(i)=0.d0
                end do  

                call marq98(b,np,niter,v1,res,ier,istop,ca,cb,dd,idmLikelihood)
                loglik(2) = res
                if (istop.ne.1) goto 1000
                
                do i=npw+1,np
                        regpar(i-npw) = b(i)
                end do        
        end if

!------ var cov matrix
        do i=1,np
                do j=i,np
                        k = (((j-1)*j)/2) +i
                        hes(i,j)= v1(k)
                end do
        end do

        do i = 2,np
                do j = 1,i-1
                        hes(i,j)=hes(j,i)
                end do
        end do

        hess_tot=hes
        K = 0
        do i = 7,np
                do j = 7,np
                        K=K+1
                        v(k)=hes(j,i)
                end do
        end do
        
!        converged = istop

        do i=1,npw
                basepar(i) = b(i)*b(i)
        end do

        cv(1) = ca
        cv(2) = cb
        cv(3) = dd
        
!----------------- sorties graphiques ------------------------
        
                do i=1,2
                        the01(i)=b(i)*b(i)
                end do

                do i=1,2
                        j = 2+i
                        the02(i)=b(j)*b(j)
                end do

                do i=1,2
                        j = 4 + i
                        the12(i)=b(j)*b(j)
                end do
        
                pas = (max-min)*0.01d0
                tx = min

                do i=1,99
                        tx = tx + pas
                        call fonct(tx,the01,ri,gl,su) 
                        t(i) = tx
                        a01(i) = ri
                        call fonct(tx,the02,ri,gl,su) 
                        a02(i) = ri
                        call fonct(tx,the12,ri,gl,su) 
                        a12(i) = ri
                end do

!------------------------- calcul de U  avec V1=UTU ---------------
                ep=10.d-10
                call dmfsd(v1,np,ep,ier)
        
                vsup=0.d0
                vinf=0.d0


                do i=1,np
                        do j=i,np
                                k = (((j-1)*j)/2) +i
                                vsup(i,j)= v1(k)
                        end do
                end do
        
                do i = 1,np
                        do j = 1,i
                                vinf(i,j)=vsup(j,i)
                        end do
                end do

!---------------- bootstrap ------------------------------
                x1 = 0.d0
                x2 = 0.d0
                do jj=1,2000
        
                        do i=1,np
                                call BGOS(1.d0,0,x1,x2,0.d0)
                                xi(i)=x1
                        end do   
        
                        do i=1,np
                                ut(i)=0.d0
                                do j=1,np ! ad avant =0,np
                                        ut(i) = ut(i)+(vinf(i,j)*xi(j))
                                end do   
                        end do 
        
                        do i=1,np
                                bh(i)= ut(i)+b(i)
                        end do 
        
                        do i=1,2
                                the01(i)=bh(i)*bh(i)
                                j = 2+i
                                the02(i)=bh(j)*bh(j)
                                j = 4 + i
                                the12(i)=bh(j)*bh(j)
                        end do

                        pas = (max-min)*0.01d0
                        tx = min
        
                        do i=1,99
                                tx = tx + pas
                                call fonct(tx,the01,ri,gl,su) 
                                mate_ri01(jj,i) = ri
                        end do
        
                        pas = (max-min)*0.01d0
                        tx = min
        
                        do i=1,99
                                tx = tx + pas
                                call fonct(tx,the02,ri,gl,su) 
                                mate_ri02(jj,i) = ri
                        end do
        
                        pas = (max-min)*0.01d0
                        tx = min
        
                        do i=1,99
                                tx = tx + pas
                                call fonct(tx,the12,ri,gl,su) 
                                mate_ri12(jj,i) = ri
                        end do
        
                end do

                do k = 1,99
                        do jj = 1,2000
                                vect(jj) = mate_ri01(jj,k)
                        end do
                        
                        call tri(vect,2000,moyenne)
                        a01_u(k) = vect(1950)
                        a01_l(k) = vect(51)
        
                        do jj = 1,2000
                                vect(jj) = mate_ri02(jj,k)
                        end do
        
                        call tri(vect,2000,moyenne)
                        a02_u(k) = vect(1950)
                        a02_l(k) = vect(51)
        
                        do jj = 1,2000
                        vect(jj) = mate_ri12(jj,k)
                        end do
        
                        call tri(vect,2000,moyenne)
                        a12_u(k) = vect(1950)
                        a12_l(k) = vect(51)
                end do  


1000    continue
        converged = istop
        deallocate(t0,t1,t2,t3,ve01,ve02,ve12,c)

        end subroutine idmWeib








