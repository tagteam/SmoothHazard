!            Survival model with non censored or interval-censored or right censored
!                               left truncated
!                                 Weibull
!                              20/05/10
!            last              08/04/11
        subroutine survWeib(entrytime,l,r,status,x,n,p,truncated,interval,eps,maxiter0,&
        loglik,basepar,regpar,v,converged,cv,niter,t,S,S_l,S_u,h,h_l,h_u,conf_bands,prt,hess_tot)
!noVar
        use survCommun  
        use parameters
        use optim
        use commun,only:pl,iconf
        
        implicit none

!  variables entrants
        integer,intent(in)::n,p,prt,conf_bands
        double precision,dimension(n),intent(in)::l,r,entrytime
        integer,dimension(n),intent(in)::status
        double precision,dimension(n,p),intent(in)::x
        integer,intent(in)::truncated,interval,maxiter0
        integer,dimension(3),intent(in)::eps
!  variables sortant
        integer,dimension(2),intent(out)::converged
        integer,intent(out)::niter
! istop 
        double precision,dimension(p),intent(out)::regpar
        double precision,dimension(100),intent(out)::t,S,S_l,S_u,h,h_l,h_u       
        double precision,dimension(2),intent(out)::loglik,basepar
        double precision,dimension(3),intent(out)::cv
        double precision,dimension(p*p),intent(out)::v   
!  variables locales             
        double precision:: res,tx0,min,max,x1,x2
        double precision,dimension(p+2)::xi,ut,bh
        double precision,dimension(2)::the
        double precision,dimension(p+2,p+2)::vinf
        double precision::su,gl,ri,tx,pas,ep
        double precision,dimension(2000,100)::mate_su,mate_ri
        double precision,dimension(52)::tab_su_i,tab_su_s,tab_ri_i,tab_ri_s        
        double precision,dimension(p+2,p+2)::vsup
        double precision,dimension(p+2)::b
        double precision,dimension((p+2)*(p+2+3)/2)::v1
        double precision,dimension(p+2,p+2)::hes
        double precision::ca,cb,dd
        integer::ier,npw,np,j,k,jj,i,kk,kkk,istop,noVar
        double precision,external::survLikelihood   
        double precision::som,ts
        double precision,dimension(p+2,p+2),intent(out)::hess_tot
!!! CT 26sept2012
!        integer,intent(in)::noVar
!!! fin CT 26sept2012

        if(p.gt.0)then
        noVar=0
        else
        noVar=1
        end if

!---------------------------------
        pl=0
        print_iter = prt
	iconf = conf_bands
        no=N
        
        if(noVar.ne.1) then
                nva=p
        else
                nva=0
        end if

        troncature = truncated


        epsa = 0.1d0**eps(1)
        epsb = 0.1d0**eps(2)
        epsd = 0.1d0**eps(3)

        maxiter = maxiter0

        som = 0.d0
        ts = 0.d0

        min = l(1)
        max = l(1)

!       write(*,*)'truncated    ',truncated
!       write(*,*)'interval     ',interval

        if(truncated.eq.1) then

                if(interval.eq.1) then
                        allocate(t0(no),t1(no),t2(no),ve(no,nva),c(no))

                        do i=1,no
                                som = som + status(i)
                                t0(i)=entrytime(i)
                                t1(i)=l(i)
                                t2(i)=r(i)
                                c(i)=status(i)
                                ts = ts + (t2(i)-t0(i))
                                if(noVar.ne.1) then
                                        do j=1,nva
                                                k = (j-1)*no+i
                                                !ve(i,j)=x(k)
                                                ve(i,j)=x(i,j)
                                        end do  
                                else
                                        ve=0.d0 
                                end if
                                if (max.lt.t2(i))max = t2(i)
                                if (min.gt.t0(i))min = t0(i)
                        end do 
                else  
                        allocate(t0(no),t1(no),t2(1),ve(no,nva),c(no))
                        do i=1,no
                                som = som + status(i)
                                t0(i)=entrytime(i)
                                t1(i)=l(i)
                                c(i)=status(i)
                                ts = ts + (t1(i)-t0(i))
                                if(noVar.ne.1) then
                                        do j=1,nva
                                                k = (j-1)*no+i
                                                !ve(i,j)=x(k)
                                                ve(i,j)=x(i,j)
                                        end do
                                else
                                        ve=0.d0 
                                end if
                                if (max.lt.t1(i))max = t1(i)
                                if (min.gt.t0(i))min = t0(i)
                        end do 
                endif   
        else
                if(interval.eq.1) then
                        
                        allocate(t0(1),t1(no),t2(no),ve(no,nva),c(no))
                        min = 0.d0
                        do i=1,no
                                som = som + status(i)
                                t1(i)=l(i)
                                t2(i)=r(i)
                                c(i)=status(i)
                                ts = ts + t2(i)
                                if(noVar.ne.1) then
                                        do j=1,nva
                                                k = (j-1)*no+i
                                                !ve(i,j)=x(k)
                                                ve(i,j)=x(i,j)
                                        end do  
                                else
                                        ve=0.d0
                                end if 
                                if (max.lt.t2(i))max = t2(i)
                        end do 

                else  

                        allocate(t0(1),t1(no),t2(1),ve(no,nva),c(no))
                        do i=1,no
                                som = som + status(i)
                                t1(i)=l(i)
                                c(i)=status(i)
                                ts = ts + t1(i)
                                if(noVar.ne.1) then
                                        do j=1,nva
                                                k = (j-1)*no+i
                                                !ve(i,j)=x(k)
                                                ve(i,j)=x(i,j)
                                        end do  
                                else
                                        ve=0.d0 
                                end if
                                if (max.lt.t1(i))max = t1(i)
                        end do 
                endif
        endif   
  
        npw = 2
        
        b(1) = 1.d0
        b(2) = dsqrt(som/ts)
!       b(1) = 1d-5! 0.1d0
!       b(2) = 1d-5! 0.1d0

        np = npw
        ca=0.d0
        cb=0.d0
        dd=0.d0
        niter=0
        res=0.d0
        ier=0
        
!       write(*,*)'first call of marquardt '
!       write(*,*)'np',np
!       write(*,*)(b(i),i=1,np)
        nva=0
        call marq98(b,np,niter,v1,res,ier,istop,ca,cb,dd,survLikelihood)
        
        if (istop.ne.1)then ! goto 1000
           b(1) = 1.d0
           b(2) = dsqrt(som/ts)
        endif
        converged(1) = istop
        loglik(1) = res
        
        if(noVar.ne.1) then
                nva=p
                np = npw + nva
                do i=npw+1,np
                        b(i)=0.d0
                end do
                ca=0.d0
                cb=0.d0
                dd=0.d0 
                ier=0
                niter=0
!               write(*,*)'second call of marquardt '
                call marq98(b,np,niter,v1,res,ier,istop,ca,cb,dd,survLikelihood)
                loglik(2) = res
                converged(2)=istop
                if (istop.ne.1) goto 1000
!c           write(*,*)'niter2',niter
                do i=npw+1,np
                        regpar(i-npw) = b(i)
                end do 
                
        else
                regpar=0.d0
                loglik(2) = 0.d0
        end if 

        do i=1,npw
                basepar(i) = b(i)**2
        end do
        
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
        do i = 3,np
                do j = 3,np
                        K=K+1
                        v(k)=hes(j,i)
                end do
        end do

        cv(1) = ca
        cv(2) = cb
        cv(3) = dd



!----------------- sorties graphiques ------------------------

                do i=1,2
                        the(i)=b(i)**2
                end do
        
                pas = (max-min)*0.01d0
                tx = min

                do i=1,100
                        if (tx.eq.0.d0)then
                                tx0 = 0.1d0*pas
                                call fonct(tx0,the,ri,gl,su) 
                                t(i) = tx0
                                S(i) = su
                                h(i) = ri
                        else   
                                call fonct(tx,the,ri,gl,su) 
                                t(i) = tx
                                S(i) = su
                                h(i) = ri
                        endif   
                        tx = tx + pas
                end do

!------------------------- calcul de U  avec V1=UTU ---------------
                ep=10.d-10
                call dmfsd(v1,np,ep,ier)
        
                do i=1,np
                        do j=1,np
                                vsup(i,j)=0.d0
                                vinf(i,j)=0.d0
                        end do
                end do   
        
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

		if (iconf.eq.1)then
                do jj=1,2000
                        do i=1,np
                                call bgos(1.d0,0,x1,x2,0.d0)
                                xi(i)=x1
                        end do
        
                        do i=1,np
                                ut(i)=0.d0
                                do j=1,np
                                        ut(i) = ut(i)+(vinf(i,j)*xi(j))
                                end do   
                        end do 
                
                        do i=1,np
                                bh(i)= ut(i)+b(i)
                        end do 
        
                        do i=1,2
                                the(i)=bh(i)**2
                        end do
        
                        pas = (max-min)*0.01d0
                        tx = min
        
                        do i=1,100
                                if (tx.eq.0.d0)then
                                        tx0 = 0.1d0*pas
                                        call fonct(tx0,the,ri,gl,su) 
                                else   
                                        call fonct(tx,the,ri,gl,su)  
                                endif 
                
                                mate_su(jj,i) = su
                                mate_ri(jj,i) = ri
                                tx = tx + pas
                        end do
                end do
    
                do k = 1,100
        
                        do kk=1,51
                                tab_su_i(kk)  = 10.d0
                                tab_su_s(kk)  = 0.d0
                                tab_ri_i(kk)  = 10.d0
                                tab_ri_s(kk)  = 0.d0
                        end do
        
                        do i=1,2000
        
                                do kk=1,51
                                        if(mate_su(i,k).lt.tab_su_i(kk))then
                                                if(kk.lt.51)then
                                                        do kkk = 51,kk+1,-1
                                                                tab_su_i(kkk) = tab_su_i(kkk-1)
                                                        end do   
                                                        tab_su_i(kk) = mate_su(i,k)
                                                        go to 426
                                                else
                                                        tab_su_i(kk) = mate_su(i,k)
                                                        go to 426
                                                endif
                                        endif 
                                end do
 426    continue   
                                do kk=1,51
                                        if(mate_ri(i,k).lt.tab_ri_i(kk))then
                                                if(kk.lt.51)then
                                                        do kkk = 51,kk+1,-1
                                                                tab_ri_i(kkk) = tab_ri_i(kkk-1)
                                                        end do   
                                                        tab_ri_i(kk) = mate_ri(i,k)
                                                        go to 430
                                                else
                                                        tab_ri_i(kk) = mate_ri(i,k)
                                                        go to 430
                                                endif
                                        endif 
                                end do
 430    continue
                                do kk=51,1,-1
                                        if(mate_su(i,k).gt.tab_su_s(kk))then
                                                if(kk.gt.1)then
                                                        do kkk = 1,kk-1
                                                                tab_su_s(kkk) = tab_su_s(kkk+1)
                                                        end do   
                                                        tab_su_s(kk) = mate_su(i,k)
                                                        go to 436
                                                else
                                                        tab_su_s(kk) = mate_su(i,k)
                                                        go to 436
                                                endif
                                        endif 
                                end do 
 436     continue

                                do kk=51,1,-1
                                        if(mate_ri(i,k).gt.tab_ri_s(kk))then
                                                if(kk.gt.1)then
                                                        do kkk = 1,kk-1
                                                                tab_ri_s(kkk) = tab_ri_s(kkk+1)
                                                        end do   
                                                        tab_ri_s(kk) = mate_ri(i,k)
                                                        go to 439
                                                else
                                                        tab_ri_s(kk) = mate_ri(i,k)
                                                        go to 439
                                                endif
                                        endif 
                                end do 
 439     continue
        
                        end do
        
                        S_l(k) = tab_su_i(51)
                        S_u(k) = tab_su_s(1)
                        h_l(k) = tab_ri_i(51)
                        h_u(k) = tab_ri_s(1)

                end do  
	else 
		S_l=0
		S_u=0
		h_l=0
		h_u=0
	end if

1000    continue
!        converged = istop
        deallocate(t0,t1,t2,ve,c)

        end subroutine survWeib








