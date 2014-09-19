!            Modele de survie 
!                avec censure par intervalles et troncature.
!                              22/02/12
!
!         Derniere modification : le         22/02/12
!
!         Commande de compilation:
!              ifort -Vaxlib -o exe commun.f90 estimv.f90 recKappa.f90 optim.f90 mark3.f90 main.f90
!         commande execution:
!              exe ou ./exe





!==============================================================================================
!====================================== Programme spline ======================================
!==============================================================================================
!add prt,noVar
        subroutine survPl(entrytime,l,r,id,x,N,nva,truncated,interval,eps &
        ,maxiter0,loglik,regpar,v,converged,cv,niter,t,S,S_l,S_u,h,h_l,h_u,&
        nknots,irec,kappa0,kappa,conf_bands,CVcrit,mdf,ti,theta,prt,hess_tot)
!!! CT 26sept2012
! noVar    
!!! fin CT 26sept2012
        use tailles
        use propre
        use parameters
        use optim
        use commun
        implicit none
        
        integer::npar,nvat,j,k,ii,i,ic,ier,interval,irep,truncated,maxiter0,nva
        double precision::min,max,estimvSurv,ca,cb,dd,res
        double precision::mdf_int
!!      double precision,dimension(:),allocatable::b2
        integer,intent(in)::conf_bands,N,nknots,irec
        integer,dimension(N)::id
        double precision,dimension(N)::entrytime,l,r
        double precision,intent(in)::kappa0
        double precision,intent(out)::kappa
        integer,dimension(3),intent(in)::eps
        double precision,dimension(N,nva)::x
        double precision,dimension(3),intent(out)::cv   
        double precision,dimension(nva),intent(out)::regpar
        double precision,dimension(nva*nva),intent(out)::v
        double precision,dimension(2),intent(out)::loglik       
        integer,dimension(2),intent(out)::converged
        integer,intent(out)::niter
! istop
        double precision,dimension(99),intent(out)::t,S,S_l,S_u,h,h_l,h_u
        double precision,external::survPlLikelihood
        double precision::golden
        double precision,intent(out)::CVcrit,mdf
        double precision::ax,bx,cx,fa,fb,fc,xmin,ddl,auxi,crit,tol
        integer::ni,istop,noVar
        double precision,dimension(nknots+6),intent(out)::ti
        double precision,dimension(nknots+2),intent(out)::theta
        double precision,dimension(nva+nknots+2,nva+nknots+2),intent(out)::hess_tot
        integer,intent(in)::prt
        !------------  entre parametres -----
        if(nva.gt.0)then
        noVar=0
        else
        noVar=1
        end if  

! pour le passage en R initialisation locale    
!add for general optim
        pl=1
        ind_hess=0
        print_iter = prt
	iconf = conf_bands
!end add
        npar = 0
        verSurv = 0
        j = 0
        k = 0
        ii = 0
        i = 0
        ic = 0
        
        kappa=0.d0
        mdf_int=0.d0
        mdf=0.d0
        loglik=0.d0
        
        cv=0.d0
        ca=0.d0
        cb=0.d0
        maxiter=maxiter0
        dd=0.d0

        epsa = 0.1d0**eps(1)
        epsb = 0.1d0**eps(2)
        epsd = 0.1d0**eps(3)

        troncature = truncated
        no=N
        k0Surv=kappa0
        
        if(noVar .ne. 1) verSurv=nva

        min = l(1)
        max = l(1)

        if(truncated.eq.1) then
                if(interval.eq.1) then
                        allocate(t0(no),t1(no),t2(no),ve(no,nva),c(no))
                        do i=1,no
                                t0(i)=entrytime(i)
                                t1(i)=l(i)
                                t2(i)=r(i)
                                c(i)=id(i)
                                if ((c(i).eq. 1).and.(t1(i).ne.t2(i)))c(i)=2 
                                do j=1,nva
                                        ve(i,j)=x(i,j)
                                end do   
                                if (max.lt.t2(i))max = t2(i)
                                if (min.gt.t0(i))min = t0(i)
                        end do 
                else  
                        allocate(t0(no),t1(no),t2(1),ve(no,nva),c(no))
                        do i=1,no
                                t0(i)=entrytime(i)
                                t1(i)=l(i)
                                c(i)=id(i)
                                do j=1,nva
                                        ve(i,j)=x(i,j)
                                end do 
                                if (max.lt.t1(i))max = t1(i)
                                if (min.gt.t0(i))min = t0(i)
                        end do 
                endif   
        else
                min = 0.d0
                if(interval.eq.1) then
                        allocate(t0(1),t1(no),t2(no),ve(no,nva),c(no))
                        do i=1,no
                                t1(i)=l(i)
                                t2(i)=r(i)
                                c(i)=id(i)
                                if ((c(i).eq. 1).and.(t1(i).ne.t2(i)))c(i)=2 
                                do j=1,nva
                                        ve(i,j)=x(i,j)
                                end do   
                                if (max.lt.t2(i))max = t2(i)
                        end do 
                else  
                        allocate(t0(1),t1(no),t2(1),ve(no,nva),c(no))
                        do i=1,no
                                t1(i)=l(i)
                                c(i)=id(i)
                                do j=1,nva
                                        ve(i,j)=x(i,j)
                                end do   
                                if (max.lt.t1(i))max = t1(i)
                        end do 
                endif
        endif   

        !--------------- zi ----------------------------------  
        !      construire vecteur zi (des noeuds)
        
        if(noVar.ne.1) ve = x
        nz = nknots

!-------------------------------------------------------------------
        allocate(zi(-2:(nz+3)))
!       write(*,*)'taille de zi',size(zi)
!-------------------------------------------------------------------

        zi(-2) = min
        zi(-1) = min
        zi(0) = min
        zi(1) = min

        do i=2,nz-1
                zi(i) =zi(i-1)+(max-min)/dble(nz-1)
        end do   

        zi(nz) = max
        zi(nz+1)=max
        zi(nz+2)=max
        zi(nz+3)=max
        
        ti = zi
!        write(*,*)zi
!---------- affectation des vecteurs de splines -----------------
        
        npar  = nz+2
        
!--------------------------------------------------------------------------------------------- 
!!      allocate(omeg(nz+2,nz+2))
!---------------------------------------------------------------------------------------------- 
        allocate(m3m3(nz),m2m2(nz),m1m1(nz),mmm(nz),m3m2(nz),m3m1(nz),&
        m3m(nz),m2m1(nz),m2m(nz),m1m(nz))
!---------------------------------------------------------------------------------------------- 
        call vecpen(nz+2,zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m)
!----------------------------------------------------------------------------------------------         

                
        np = npar
        allocate(b(np),bt(np))
        allocate(hes(np,np))

        do i=1,np
           b(i)=5.d-1
        end do
        
        rec = irec
!---------- newton --------------------------------------------

        nvat = nva
        verSurv = 0

!       write(*,*)'rec',rec
        if(rec.eq.1)then
!               write(*,*)'recherche automatique du parametre de lissage'
                xmin=dsqrt(kappa0)
                auxi = estimvSurv(xmin,b,ddl,ni,res)
                if(ddl.gt.-2.5d0)then
                xmin = dsqrt(xmin)
                auxi = estimvSurv(xmin,b,ddl,ni,res)
                if(ddl.gt.-2.5d0)then
                        xmin = dsqrt(xmin)
                        auxi = estimvSurv(xmin,b,ddl,ni,res)
                        if(ddl.gt.-2.5d0)then
                        xmin = dsqrt(xmin)
                        auxi = estimvSurv(xmin,b,ddl,ni,res)
                        if(ddl.gt.-2.5d0)then
                        xmin = dsqrt(xmin)
                        auxi = estimvSurv(xmin,b,ddl,ni,res)
                        if(ddl.gt.-2.5d0)then
                                xmin = dsqrt(xmin)
                        endif   
                        endif   
                        endif   
                endif
                endif 

                if (ni.ge.maxiter) then
                        do i=1,nz+2
                                b(i)=1.d-1
                        end do     
                        xmin = sqrt(10.d0)*xmin
                        auxi = estimvSurv(xmin,b,ddl,ni,res)
                        if (ni.ge.maxiter) then
                                do i=1,nz+2
                                        b(i)=1.d-1
                                end do     
                                xmin = sqrt(10.d0)*xmin
                        endif
                endif  
                ax = xmin
                bx = xmin*dsqrt(1.5d0)
                call mnbrak(ax,bx,cx,fa,fb,fc,b)
                tol = 0.001d0
!               write(*,*)ax,bx,cx
                res = golden(ax,bx,cx,tol,xmin,b,ddl)
!               write(*,*)' '
!               write(*,*)'Parametre de lissage optimum',real(xmin*xmin)  &
!                       , '  DDL :',-ddl
!               write(*,*)'================== Fin ======================'
                k0Surv =  xmin*xmin    
        endif
        crit = estimvSurv(dsqrt(k0Surv),b,ddl,ni,res)
!        write(*,*)'test ddl',ddl,crit
        mdf=-ddl

        CVcrit = crit/dble(no)

        allocate(v1(np*(np+3)/2))
     
        call marq98(b,np,niter,v1,res,ier,istop,ca,cb,dd,survPlLikelihood)
        converged(1) = istop

        if (istop.ne.1)then 
           do i=1,np
              b(i)=5.d-1
           end do
        endif
        loglik(1)=res


        if(noVar .ne. 1)then
!               write(*,*)' '
!               write(*,*)'       Recherche des parametres de regression'
                verSurv = nvat
                bt=b

                deallocate(b,hes,v1)
                npar = npar + nva
                allocate(b(npar),hes(npar,npar),v1(npar*(npar+3)/2))

                b(1:np)=bt
        !       deallocate(bt)  
                do i=np+1,npar
                        b(i)=0.d0
                end do  

                call marq98(b,npar,niter,v1,res,ier,istop,ca,cb,dd,survPlLikelihood)
                converged(2) = istop
                loglik(2) = res
                if (istop.ne.1) goto 1000
               
!               write(*,*)'log-vrais',res,' nb iter ',niter

                do i=1,nva
                        ii = npar-nva+i
                        regpar(i)=b(ii)
                end do 

        
!               write(*,*)' '
!               write(*,*)nomfichregr,' fichier variables explicatives'
!               write(*,*)' '

        else
                loglik(2) = 0.d0
                regpar = 0.d0
        endif
        
        !------ paramètres de coefficients de splines

        theta = b(1:(nz+2))  
        
        do i=1,npar
                do j=i,npar
                        k = (((j-1)*j)/2) +i
                        hes(i,j)= v1(k)
                end do
        end do

        do i = 2,npar
                do j = 1,i-1
                        hes(i,j)=hes(j,i)
                end do
        end do
        hess_tot=hes
        K = 0
        do i = np+1,npar
                do j = np+1,npar
                        K=K+1
                        v(k)=hes(j,i)
                end do
        end do
                
!       write(*,*)'taille v',size(v)
!       write(*,*)
!       write(*,*)'taille v1',size(v1)

        cv(1) = ca
        cv(2) = cb
        cv(3) = dd

!        irep = igraph

!        if(irep.eq.1)then
                call distanceSurv(v1,b,t,S,S_l,S_u,h,h_l,h_u)
!        endif   
        kappa=k0Surv

1000    continue  
              
        converged = istop        
        deallocate(b,hes,v1,bt)
        deallocate(t0,t1,t2,c,zi,&
        m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m)
        deallocate(ve)

        end subroutine survPl

!=====================================================================



