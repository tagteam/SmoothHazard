!            Modele illness-death a trois etats Markov non-homogene
!                avec censure par intervalles et troncature.
!                              03/08/98
!           corrige le         12/03/09
!
!         Derniere modification : le         07/06/11
!
!         Commande de compilation:
!              ifort -Vaxlib -o exe commun.f90 estimv.f90 recKappa.f90 optim.f90 mark3.f90 main.f90
!         commande execution:
!              exe ou ./exe




        module propre           
                double precision,dimension(:,:),allocatable,save:: omeg01,omeg12,omeg02         
                double precision,dimension(:),allocatable,save::b,bt,v1
                double precision,dimension(:,:),allocatable::hes
        end module propre
        
        
        module fichiers 
                character(len=14)::nomfich,nomfich01,nomfich12,nomfich02,nomfichregr, &
                fichparam
                integer,save::cpt1,cpt2,cpt3
                double precision,save::crit     
        end module fichiers

!==============================================================================================
!====================================== Programme spline ======================================
!==============================================================================================

        subroutine idmPl(entrytime,l,r,d,idm,idd,x01,x02,x12,N,P01,P02,P12,truncated,eps &
        ,maxiter0,loglik,regpar,v,converged,cv,niter,t,a01,a01_l,a01_u,a02,a02_l,a02_u,a12, &
        a12_l,a12_u,nknots01,knots01,nknots02,knots02,nknots12,knots12,irec,kappa0,kappa,&
	conf_bands,CVcrit,mdf,theta01,theta02,theta12,prt,hess_tot)
        
        use tailles
        use propre
        use parameters
        use optim
        use commun
        use recKappa
        use fichiers

        implicit none
        
        integer::npar,ver,nvat12,nvat02,j,k,ii,i,ic,ier,nvat01,truncated,maxiter0,istop
        double precision::min,max,maxt1,estimv,maxt2,maxt3,mint1,mint2,mint3,ca,cb,dd,res
        double precision::mdf_int
        double precision,dimension(:),allocatable::b2
        integer,intent(in)::P01,P12,P02,N,nknots01,nknots02,nknots12,irec,prt,conf_bands
        integer,dimension(N)::idm,idd
        double precision,dimension(N)::entrytime,l,r,d
        double precision,dimension(3),intent(in)::kappa0
        double precision,dimension(3),intent(out)::kappa
        double precision,dimension(N,P01)::x01
        double precision,dimension(N,P12)::x12
        double precision,dimension(N,P02)::x02
        integer,dimension(3),intent(in)::eps
        double precision,dimension(3),intent(out)::cv   
        double precision,dimension(p01+p02+p12),intent(out)::regpar
        double precision,dimension((p01+p02+p12)*(p01+p02+p12)),intent(out)::v
        double precision,dimension(2),intent(out)::loglik       
        integer,intent(out)::converged,niter
        double precision,dimension(99,3),intent(out)::t
        double precision,dimension(99),intent(out)::a01,a01_l,a01_u,a02,a02_l,a02_u,a12,a12_l,a12_u
        double precision,external::idmPlLikelihood
        double precision,intent(out)::CVcrit,mdf
        double precision,dimension(nknots01+6),intent(in)::knots01
        double precision,dimension(nknots12+6),intent(in)::knots12
        double precision,dimension(nknots02+6),intent(in)::knots02
        double precision,dimension(nknots01+2),intent(out)::theta01
        double precision,dimension(nknots12+2),intent(out)::theta12
        double precision,dimension(nknots02+2),intent(out)::theta02
        double precision,dimension(p01+p02+p12+nknots01+nknots12+nknots02+6,&
        p01+p02+p12+nknots01+nknots12+nknots02+6),intent(out)::hess_tot
        integer,dimension(3)::noVar

        nz01 = nknots01
        nz12 = nknots12
        nz02 = nknots02 

        allocate(zi01(-2:(nz01+3)))        
        allocate(zi02(-2:(nz02+3)))
        allocate(zi12(-2:(nz12+3)))
        
        zi01 = knots01
        zi12 = knots12
        zi02 = knots02 

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
        
        !------------  entre parametres -----
        
        print_iter = prt
	iconf = conf_bands
        npar = 0
        ver = 0
        pl=1
        nvat01 = 0      
        nvat12 = 0
        nvat02 = 0
        
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
   
        no=N
        k0=kappa0


!--------------- tableau de tailles ver --------------------------


        if(noVar(1) .ne. 1) nvat01=P01
        if(noVar(2) .ne. 1) nvat02=P02
        if(noVar(3) .ne. 1) nvat12=P12
        
        ver=nvat01+nvat02+nvat12        

        troncature = truncated  

        allocate(t1(no),t2(no),t3(no),t4(no),c(no),ve01(no,nvat01),ve12(no,nvat12),ve02(no,nvat02))
!!!! a voir avec Pierre et celia

        if(truncated.eq.1) then
                allocate(t0(no))
                do i=1,no
                        t0(i)=entrytime(i)
                        t1(i)=l(i)
                        t2(i)=r(i)
                        t3(i)=d(i)
                end do 
        else
                allocate(t0(1))
                do i=1,no
                        t1(i)=l(i)
                        t2(i)=r(i)
                        t3(i)=d(i)
                end do 
        endif 

! Prise des variables explicatives
        ve01=0.d0
        ve02=0.d0
        ve12=0.d0
        
        if(noVar(1).ne.1) ve01=x01
        if(noVar(2).ne.1) ve02=x02      
        if(noVar(3).ne.1) ve12=x12

        
        maxt1 = 0.d0
        maxt2 = 0.d0
        maxt3 = 0.d0
        cpt1 = 0
        cpt2 = 0
        cpt3 = 0

        do i = 1,no
        
                if((idm(i).eq.0).and.(idd(i).eq.0).and.(t1(i).eq.t3(i)))then
                        c(i) = 1
                endif

                if((idm(i).eq.1).and.(idd(i).eq.0).and.(t1(i).lt.t2(i)))then
                        c(i) = 2
                endif

                if((idm(i).eq.1).and.(idd(i).eq.0).and.(t1(i).eq.t2(i)))then
                        c(i) = 3
                endif

                if((idm(i).eq.1).and.(idd(i).eq.1).and.(t1(i).lt.t2(i)))then
                        c(i) = 4
                endif
        
                if((idm(i).eq.1).and.(idd(i).eq.1).and.(t1(i).eq.t2(i)))then
                        c(i) = 5
                endif
        
                if((idm(i).eq.0).and.(idd(i).eq.0).and.(t1(i).lt.t3(i)))then
                        c(i) = 6
                endif
        
                if((idm(i).eq.0).and.(idd(i).eq.1))then
                        c(i) = 7
                endif

        ! -----------------  c=1 censure a droite 01
                if(c(i).eq.1)then
                        c(i) = 1
                        if(truncated.eq.1) then                 
                                t0(i) = t0(i)
                                t1(i) = t1(i)
                                t2(i) = t1(i)
                                t3(i) = t1(i)
                                t4(i) = t1(i)
                        else            
                                t1(i) = t1(i)
                                t2(i) = t1(i)
                                t3(i) = t1(i)
                                t4(i) = t1(i)
                        end if

                else
! ----------------- c=2 censure par intervalle de 01 et a droite 12
                        if(c(i).eq.2)then

                                if(truncated.eq.1) then 
                                        t0(i) = t0(i)
                                        t1(i) = t1(i)
                                        t2(i) = t2(i)
                                        t3(i) = t3(i)
                                        t4(i) = t1(i)
                                else
                                        t1(i) = t1(i)
                                        t2(i) = t2(i)
                                        t3(i) = t3(i)
                                        t4(i) = t1(i)
                                end if

                        else
        ! ----------------- c=3 observation 01  et censure a droite 12
                                if(c(i).eq.3)then
                                        if(truncated.eq.1) then 
                                                t0(i) = t0(i)
                                                t1(i) = t1(i)
                                                t2(i) = t1(i)
                                                t3(i) = t3(i)
                                                t4(i) = t1(i)
                                        else
                                                t1(i) = t1(i)
                                                t2(i) = t1(i)
                                                t3(i) = t3(i)
                                                t4(i) = t1(i)
                                        end if
                ! ----------------- c=4 censure par intervalle de 01 et observation 12
                                else   
                                        if(c(i).eq.4)then
                                                cpt1 = cpt1 + 1
                                                cpt2 = cpt2 + 1
                                                if(truncated.eq.1) then 
                                                        t0(i) =  t0(i)
                                                        t1(i) =  t1(i)
                                                        t2(i) =  t2(i)
                                                        t3(i) =  t3(i)
                                                        t4(i) =  t1(i)
                                                else
                                                        t1(i) =  t1(i)
                                                        t2(i) =  t2(i)
                                                        t3(i) =  t3(i)
                                                        t4(i) =  t1(i)
                                                end if
                                        else
                                                if(c(i).eq.5)then
! ----------------- c=5 observation de 01 et observation 12
                                                        cpt1 = cpt1 + 1
                                                        cpt2 = cpt2 + 1
                                                        if(truncated.eq.1) then 
                                                                t0(i) =  t0(i)
                                                                t1(i) =  t1(i)
                                                                t2(i) =  t1(i)
                                                                t3(i) =  t3(i)
                                                                t4(i) = t1(i)
                                                        else
                                                                t1(i) =  t1(i)
                                                                t2(i) =  t1(i)
                                                                t3(i) =  t3(i)
                                                                t4(i) = t1(i)
                                                        end if
                                                else
                                                        if(c(i).eq.6)then
! ----------------- c=6 vivant 
                                                                if(truncated.eq.1) then 
                                                                        t0(i) = t0(i)
                                                                        t1(i) =  t1(i)
                                                                        t2(i) =  t3(i)
                                                                        t3(i) = t3(i)
                                                                        t4(i) = t3(i)
                                                                else
                                                                        t1(i) =  t1(i)
                                                                        t2(i) =  t3(i)
                                                                        t3(i) = t3(i)
                                                                        t4(i) = t3(i)
                                                                end if

                                                        else
! ----------------- c=7 passage 0 2
                                                                cpt3 = cpt3 + 1
                                                                if(truncated.eq.1) then 
                                                                        t0(i) = t0(i)
                                                                        t1(i) =  t1(i)
                                                                        t2(i) =  t3(i)
                                                                        t3(i) =  t3(i)
                                                                        t4(i) = t3(i)
                                                                else
                                                                        t1(i) =  t1(i)
                                                                        t2(i) =  t3(i)
                                                                        t3(i) =  t3(i)
                                                                        t4(i) = t3(i)
                                                                end if

                                                        endif
                                                endif
                                        endif
                                endif
                        endif
                endif 

                if (maxt1.lt.t2(i))then
                        maxt1 = t2(i)
                endif
        
                if (maxt3.lt.t4(i))then
                        maxt3 = t4(i)
                endif
        
                if (maxt2.lt.t3(i))then
                        maxt2 = t3(i)
                endif
        
        end do   !avoir

        !--------------- zi ----------------------------------  
        !      construire vecteur zi (des noeuds)
        
        min = 0.d0
        max = maxt1
        
        do k = 1,no
                if(truncated.eq.1) then
                        if(t0(k).ge.min)then 
                                if(t0(k).lt.max)then
                                        max = t0(k)
                                endif   
                        endif
                end if

                if((t1(k).ge.min))then
                        if(t1(k).lt.max)then
                                max = t1(k)
                        endif
                endif

                if(t2(k).ge.min)then 
                        if(t2(k).lt.max)then
                                max = t2(k)
                        endif   
                endif
        end do   
! ajout
        mint1=max
! fin ajout

!---------- affectation des vecteurs de splines -----------------
        
        npar  = nz01+2+nz02+2+nz12+2

                
!--------------------------------------------------------------------------------------------- 
        allocate(omeg01(nz01+2,nz01+2),omeg12(nz12+2,nz12+2),omeg02(nz02+2,nz02+2))
!---------------------------------------------------------------------------------------------
        allocate(mm3a(no),mm2a(no),mm1a(no),mma(no),im3a(no),im2a(no),im1a(no),ima(no))
        allocate(mm3b(no),mm2b(no),mm1b(no),mmb(no),im3b(no),im2b(no),im1b(no),imb(no))
        allocate(mm3c(no),mm2c(no),mm1c(no),mmc(no),im3c(no),im2c(no),im1c(no),imc(no))

!---------------------------------------------------------------------------------------------- 
        call vecspli(nz01+2,no,zi01,t1,mm3a,mm2a,mm1a,mma,im3a,im2a,im1a,ima)
        call vecspli(nz12+2,no,zi12,t2,mm3b,mm2b,mm1b,mmb,im3b,im2b,im1b,imb)
        call vecspli(nz02+2,no,zi02,t3,mm3c,mm2c,mm1c,mmc,im3c,im2c,im1c,imc)

!---------------------------------------------------------------------------------------------- 
        allocate(m3m3a(nz01),m2m2a(nz01),m1m1a(nz01),mmma(nz01),m3m2a(nz01),m3m1a(nz01),&
        m3ma(nz01),m2m1a(nz01),m2ma(nz01),m1ma(nz01))
        allocate(m3m3b(nz12),m2m2b(nz12),m1m1b(nz12),mmmb(nz12),m3m2b(nz12),m3m1b(nz12),&
        m3mb(nz12),m2m1b(nz12),m2mb(nz12),m1mb(nz12))
        allocate(m3m3c(nz02),m2m2c(nz02),m1m1c(nz02),mmmc(nz02),m3m2c(nz02),m3m1c(nz02),&
        m3mc(nz02),m2m1c(nz02),m2mc(nz02),m1mc(nz02))   
!----------------------------------------------------------------------------------------------         

        call vecpen(nz01+2,zi01,m3m3a,m2m2a,m1m1a,mmma,m3m2a,m3m1a,m3ma,m2m1a,m2ma,m1ma)
        call vecpen(nz12+2,zi12,m3m3b,m2m2b,m1m1b,mmmb,m3m2b,m3m1b,m3mb,m2m1b,m2mb,m1mb)
        call vecpen(nz02+2,zi02,m3m3c,m2m2c,m1m1c,mmmc,m3m2c,m3m1c,m3mc,m2m1c,m2mc,m1mc)
!----------------------------------------------------------------------------------------------         
        call calcomeg(nz01+2,omeg01,m3m3a,m2m2a,m1m1a,mmma,m3m2a,m3m1a,m3ma,m2m1a,m2ma,m1ma)    
        call calcomeg(nz12+2,omeg12,m3m3b,m2m2b,m1m1b,mmmb,m3m2b,m3m1b,m3mb,m2m1b,m2mb,m1mb)    
        call calcomeg(nz02+2,omeg02,m3m3c,m2m2c,m1m1c,mmmc,m3m2c,m3m1c,m3mc,m2m1c,m2mc,m1mc)
                
        np = npar
! Modif CT 19juin2012
!       allocate(b(np),bt(np))
        allocate(b(np))
! Fin modif changement de place CT 19juin2012   
        allocate(opt2(np))
        allocate(hes(np,np))

        
        do i=1,np
                b(i)=5.d-1
        end do
        
        rec = irec
!---------- newton --------------------------------------------

        nn = 3
        allocate(b2(nn))
        
        nva01 = 0
        nva12 = 0
        nva02 = 0       

        
!       write(*,*)'rec',rec
        if(rec.eq.1)then
!               write(*,*)'recherche automatique du parametre de lissage'
                ind_hess = 0
                call sub_rec(b,omeg01,omeg12,omeg02)
!               write(*,*)'================== Fin ======================'       
        endif
        allocate(v1(np*(np+3)/2))
!       write(*,*)'avant marq98'
        ind_hess = 1
        allocate(hessienne(np,np))
!       write(*,*) '====> marquardt sans variable expli'
        call marq98(b,np,niter,v1,res,ier,istop,ca,cb,dd,idmPlLikelihood)
!       write(*,*) 'apres marquardt'
        

!       write(*,*) 'apres marquardt'
!       write(*,*) '====> '


        if (istop.ne.1) then
                deallocate(hessienne)
                goto 1000
        end if

        loglik(1)=res

        ind_hess = 0    
        
        b2(1) = dlog(k0(1))
        b2(2) = dlog(k0(2))
        b2(3) = dlog(k0(3))

!       write(*,*) '====> estimv'       
        crit = estimv(b2,b,omeg01,omeg12,omeg02,mdf_int)
!       write(*,*) '====> '
        mdf=mdf_int
        
        CVcrit = crit/dble(no)
!       write(*,*)' '
!       write(*,*)'Critere de validation croisee : ',CVcrit


        if((noVar(1)==0).or.(noVar(2)==0).or.(noVar(3)==0)) then
        
!               write(*,*)' '
!               write(*,*)'       Recherche des parametres de regression'
                nva01 = nvat01
                nva12 = nvat12
                nva02 = nvat02
                 
                
                allocate(bt(np))
                bt=b
                deallocate(b,hes,v1)    
                npar = npar + nva01 + nva12 + nva02
                allocate(b(npar),hes(npar,npar),v1(npar*(npar+3)/2))
                
                deallocate(hessienne)
                allocate(hessienne(npar,npar))
                ind_hess=1
                 
                b(1:np)=bt
                deallocate(bt)
                do i=np+1,npar
                        b(i)=0.d0
                end do  
                
!               write(*,*) '====> marquardt avec variable expli'

                call marq98(b,npar,niter,v1,res,ier,istop,ca,cb,dd,idmPlLikelihood)
                ind_hess=0
!               write(*,*) '====> '
! Modif CT 4juillet2012
                if (istop.ne.1) then
                        deallocate(hessienne)
                        goto 1000
                end if
! Fin modif CT 4juillet2012
                
!               write(*,*)'log-vrais',res,' nb iter ',niter
                loglik(2) = res
                
                do i=1,nva01
                        ii = npar-nva01-nva12-nva02+i
                        regpar(i)=b(ii)
                end do 

                do i=1,nva02
                        ii = npar-nva12-nva02+i
                        regpar(i+nva01)=b(ii)
                end do 

                do i=1,nva12
                        ii = npar-nva12+i
                        regpar(i+nva01+nva02)=b(ii)
                end do 

                
                
        else
                t=0.d0
                a01=0.d0
                a01_l=0.d0
                a01_u=0.d0
                a02=0.d0
                a02_l=0.d0
                a02_u=0.d0
                a12=0.d0
                a12_l=0.d0
                a12_u=0.d0      
        
        end if  


!------ paramètres de coefficients de splines
                theta01 = b(1:(nz01+2))
                theta02 = b((nz01+3):(nz01+nz02+4))
                theta12 = b((nz01+nz02+5):(nz01+nz12+nz02+6))

                
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

            
                
                
!               write(*,*)' '
!               write(*,*)nomfichregr,' fichier variables explicatives'
!               write(*,*)' '

!!        if(igraph.eq.1)then
                call distance(nz01,nz02,nz12,b,t,a01,a01_l,a01_u,a02,a02_l,a02_u,a12,a12_l,a12_u)
!!       end if

        deallocate(hessienne)
!        if(igraph.ne.1)then 
!                t=0.d0
!                a01=0.d0
!                a01_l=0.d0
!                a01_u=0.d0
!                a02=0.d0
!                a02_l=0.d0
!                a02_u=0.d0
!                a12=0.d0
!                a12_l=0.d0
!                a12_u=0.d0
!        end if  

!        converged = istop
        cv(1) = ca
        cv(2) = cb
        cv(3) = dd      
        kappa=k0

1000    continue        
        converged = istop
        deallocate(b,hes,v1)
        deallocate(t0,t1,t2,t3,t4,c,zi01,zi12,zi02,omeg01,omeg12,omeg02,&
        mm3a,mm2a,mm1a,mma,im3a,im2a,im1a,ima,mm3b,mm2b,mm1b,mmb,im3b,im2b,im1b,&
        imb,mm3c,mm2c,mm1c,mmc,im3c,im2c,im1c,imc,m3m3a,m2m2a,m1m1a,mmma,m3m2a,&
        m3m1a,m3ma,m2m1a,m2ma,m1ma,m3m3b,m2m2b,m1m1b,mmmb,m3m2b,m3m1b,m3mb,m2m1b,&
        m2mb,m1mb,m3m3c,m2m2c,m1m1c,mmmc,m3m2c,m3m1c,m3mc,m2m1c,m2mc,m1mc,b2)
        deallocate(ve01,ve12,ve02,opt2) 
        end subroutine idmPl

!=====================================================================

        
!========================== VECSPLI ==============================

        subroutine vecspli(n,no,zi,t,mm3,mm2,mm1,mm,im3,im2,im1,im)
        
        implicit none
        
        integer::no,n
        double precision,dimension(no)::t        
        double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2, & 
        h3,h4,h3m,h2n,hn,hh3,hh2
        integer::i,j,k
        double precision,dimension(-2:(n+1))::zi
        double precision,dimension(no)::mm3,mm2,mm1,mm,im3,im2,im1,im

!----------  calcul de u(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)

!       write(*,*)'taille de zi dans vecspli',size(zi)
        i=0
        k=0
        j=0     
        do i=1,no
                do k = 2,n-2
                        if ((t(i).ge.zi(k-1)).and.(t(i).lt.zi(k))) then
                                j = k-1
                        endif
                end do 
        
                ht = t(i)-zi(j)
                htm= t(i)-zi(j-1)
                h2t= t(i)-zi(j+2)       
                ht2 = zi(j+1)-t(i)
                ht3 = zi(j+3)-t(i)
                hht = t(i)-zi(j-2)
        
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
                mm3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
                mm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))+((-4.d0*h2t*htm &
                *ht2)/(hh2*h2n*hh*h))+((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
                mm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h))+((-4.d0*htm*ht* &
                h2t)/(h3m*h2*h*h2n))+((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
                mm(i)  = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
        
                im3(i) = (0.25d0*(t(i)-zi(j-3))*mm3(i))+(0.25d0*hh2*mm2(i))+(0.25d0*h3m*mm1(i))+(0.25d0*h4*mm(i))
        
                im2(i) = (0.25d0*hht*mm2(i))+(h3m*mm1(i)*0.25d0) &
                        +(h4*mm(i)*0.25d0)
                im1(i) = (htm*mm1(i)*0.25d0)+(h4*mm(i)*0.25d0)
                im(i)  = ht*mm(i)*0.25d0
        end do

        end subroutine vecspli
!========================== VECPEN ==============================
      subroutine vecpen(n,zi,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m)
         implicit none
         integer::n,i
         double precision::h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2, &
          a3,a2,b2,c2,a1,b1,c1,a0,x3,x2,x
         double precision,dimension(-2:(n+1))::zi
         double precision,dimension(n-2):: m3m3,m2m2,m1m1,mmm,m3m2, &
         m3m1,m3m,m2m1,m2m,m1m


         do i=1,n-3
            h = zi(i+1)-zi(i)
            hh= zi(i+1)-zi(i-1)
            h2= zi(i+2)-zi(i)
            h3= zi(i+3)-zi(i)
            h4= zi(i+4)-zi(i)
            h3m= zi(i+3)-zi(i-1)
            h2n=zi(i+2)-zi(i-1)
            hn= zi(i+1)-zi(i-2)
            hh3 = zi(i+1)-zi(i-3)
            hh2 = zi(i+2)-zi(i-2)
            a3 = h*hh*hn*hh3
            a2 = hh2*hh*h*hn
            b2 = hh2*h2n*hh*h
            c2 = hh2*h2*h*h2n
            a1 = h3m*h2n*hh*h
            b1 = h3m*h2*h*h2n
            c1 = h3m*h3*h2*h
            a0 = h4*h3*h2*h
            x3 = zi(i+1)*zi(i+1)*zi(i+1)-zi(i)*zi(i)*zi(i)
            x2 = zi(i+1)*zi(i+1)-zi(i)*zi(i)
            x  = zi(i+1)-zi(i)
            m3m3(i) = (192.d0*h/(hh*hn*hh3*hh*hn*hh3))
            m2m2(i) = 64.d0*(((3.d0*x3-(3.d0*x2*(2.d0*zi(i+1)+zi(i-2) &
           ))+x*(4.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i-2)+4.d0*zi(i+1) &
           *zi(i-2)))/(a2*a2)))
            m2m2(i) = m2m2(i) + 64.d0*(((3.d0*x3-(3.d0*x2*(zi(i+2) &
           +zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1) &
           +zi(i+1)*zi(i+1)+2.d0*zi(i+2)*zi(i-1)+2.d0*zi(i+2) &
           *zi(i+1)+2.d0*zi(i-1)*zi(i+1)))/(b2*b2)))
            m2m2(i) = m2m2(i) +64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i+2) &
          +zi(i)))+x*(4.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i)+4.d0*zi(i+2) &
           *zi(i)))/(c2*c2))
            m2m2(i) = m2m2(i) +128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+2) &
           +zi(i-1)+3.d0*zi(i+1)+zi(i-2)))+x*(2.d0*zi(i+1)*zi(i+2) &
           +2.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1)*zi(i+1)+zi(i-2)*zi(i+2) &
           +zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*b2))
            m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0* &
           x2*(2.d0*zi(i+2)+zi(i)+2.d0*zi(i+1)+zi(i-2)))+x* &
           (4.d0*zi(i+1)*zi(i+2)+2.d0*zi(i+1)*zi(i)+2.d0*zi(i-2) &
           *zi(i+2)+zi(i-2)*zi(i)))/(a2*c2))
            m2m2(i) = m2m2(i) + 128.d0*((3.d0*x3-(1.5d0*x2  &
          *(3.d0*zi(i+2)+zi(i)+zi(i-1)+zi(i+1)))+x*(zi(i+2)*zi(i)+ &
           2.d0*zi(i-1)*zi(i+2)+zi(i)*zi(i-1)+2.d0*zi(i+1)*zi(i+2) &
           +zi(i+1)*zi(i)+2.d0*zi(i+2)*zi(i+2)))/(b2*c2))
            m1m1(i) = 64.d0*((3.d0*x3-(3.d0*x2*(2.d0*zi(i-1)+zi(i+1))) &
           +x*(4.d0*zi(i-1)*zi(i-1)+zi(i+1)*zi(i+1)+4.d0*zi(i-1) &
           *zi(i+1)))/(a1*a1))
            m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i-1)+zi(i) &    
           +zi(i+2)))+x*(zi(i-1)*zi(i-1)+zi(i)*zi(i)+zi(i+2)* &
           zi(i+2)+2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+2.d0* &
           zi(i)*zi(i+2)))/(b1*b1))
            m1m1(i) = m1m1(i) + 64.d0*((3.d0*x3-(3.d0*x2*(zi(i+3) &
           +2.d0*zi(i)))+x*(zi(i+3)*zi(i+3)+4.d0*zi(i)*zi(i) &
           +4.d0*zi(i+3)*zi(i)))/(c1*c1)) 
            m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(3.d0 &
           *zi(i-1)+zi(i)+zi(i+2)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i-1) &
           +2.d0*zi(i-1)*zi(i)+2.d0*zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1) &
           +zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(a1*b1))
            m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+ &
           2.d0*zi(i)+2.d0*zi(i-1)+zi(i+1)))+x*(2.d0*zi(i-1)*zi(i+3) &
           +4.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i))) &
            /(a1*c1))    
            m1m1(i) = m1m1(i) + 128.d0*((3.d0*x3-(1.5d0*x2*(zi(i+3)+3.d0 &
           *zi(i)+zi(i-1)+zi(i+2)))+x*(zi(i-1)*zi(i+3)+2.d0*zi(i-1)   &  
           *zi(i)+zi(i+3)*zi(i)+2.d0*zi(i)*zi(i)+zi(i+2)*zi(i+3) &
           +2.d0*zi(i+2)*zi(i)))/(b1*c1))
            mmm(i) = (192.d0*h/(h4*h3*h2*h4*h3*h2))
            m3m2(i) = 192.d0*(((-x3+(0.5d0*x2*(5.d0*zi(i+1)+zi(i-2) &
           ))-x*(2.d0*zi(i+1)*zi(i+1)+zi(i+1)*zi(i-2)))/(a3*a2)) &
           +((-x3+(0.5d0*x2*(4.d0*zi(i+1)+zi(i-1)+zi(i+2)))-x* &
            (zi(i+1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*b2)) &
           +((-x3+(0.5d0*x2*(3.d0*zi(i+1)+2.d0*zi(i+2)+zi(i)))-x* &
           (2.d0*zi(i+1)*zi(i+2)+zi(i+1)*zi(i)))/(a3*c2)) )
            m3m1(i) = 192.d0*(((x3-(0.5d0*x2*(4.d0*zi(i+1)+2.d0*zi(i-1) &
           ))+x*(2.d0*zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(a3*a1)) &
           +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+2)+zi(i-1)+zi(i))) &
           +x*(zi(i+1)*zi(i-1)+zi(i+1)*zi(i)+zi(i+1)*zi(i+2)))/(b1*a3)) &
          +((x3-(0.5d0*x2*(3.d0*zi(i+1)+zi(i+3)+2.d0*zi(i)))+x*(zi(i+1) &
           *zi(i+3)+2.d0*zi(i+1)*zi(i)))/(c1*a3)) )
            m3m(i) = 576.d0*((-(x3/3.d0)+(0.5d0*x2*(zi(i+1)+zi(i))) &
            -x*zi(i+1)*zi(i))/(a3*a0))
            m2m1(i) = 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1)+3.d0* &
           zi(i+1)+zi(i-2)))-x*(4.d0*zi(i+1)*zi(i-1)+2.d0*zi(i+1) &
           *zi(i+1)+2.d0*zi(i-2)*zi(i-1)+zi(i-2)*zi(i+1)))/(a2*a1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1)+ &
           zi(i)+zi(i+2)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i-1) &
           +2.d0*zi(i+1)*zi(i)+2.d0*zi(i+1)*zi(i+2)+zi(i-2)*zi(i-1)+ &
           zi(i-2)*zi(i)+zi(i-2)*zi(i+2)))/(a2*b1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3)+2.d0 & 
           *zi(i)+2.d0*zi(i+1)+zi(i-2)))-x*(2.d0*zi(i+1)*zi(i+3)+4.d0 &
           *zi(i+1)*zi(i)+zi(i-2)*zi(i+3)+2.d0*zi(i-2)*zi(i)))/(a2*c1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2* &
           (3.d0*zi(i-1)+2.d0*zi(i+1)+zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1) & 
           +zi(i+2)*zi(i+1)+2.d0*zi(i-1)*zi(i-1)+3.d0 &
            *zi(i+1)*zi(i-1)+zi(i+1)*zi(i+1)))/(b2*a1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0 &
           *zi(i-1)+zi(i)+2.d0*zi(i+2)+zi(i+1)))-x*(zi(i+2)*zi(i-1) & 
           +zi(i+2)*zi(i)+zi(i+2)*zi(i+2)+zi(i-1)*zi(i-1)+zi(i-1) &
           *zi(i)+zi(i-1)*zi(i+2)+zi(i+1)*zi(i-1)+zi(i+1)*zi(i) &
           +zi(i+1)*zi(i+2)))/(b2*b1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3) &
           +2.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1)))-x*(zi(i+2)*zi(i+3) &
           +2.d0*zi(i+2)*zi(i)+zi(i-1)*zi(i+3)+2.d0*zi(i-1)*zi(i) &
           +zi(i+1)*zi(i+3)+2.d0*zi(i+1)*zi(i)))/(b2*c1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(2.d0*zi(i-1) &
           +zi(i+1)+2.d0*zi(i+2)+zi(i)))-x*(4.d0*zi(i+2)*zi(i-1)+2.d0* &
           zi(i+2)*zi(i+1)+2.d0*zi(i)*zi(i-1)+zi(i)*zi(i+1)))/(c2*a1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i-1) &
           +2.d0*zi(i)+3.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i-1)+2.d0 &
           *zi(i+2)*zi(i)+2.d0*zi(i+2)*zi(i+2)+zi(i)*zi(i-1)+zi(i) &
           *zi(i)+zi(i)*zi(i+2)))/(c2*b1))
            m2m1(i) = m2m1(i) + 64.d0*((-3.d0*x3+(1.5d0*x2*(zi(i+3) &
           +3.d0*zi(i)+2.d0*zi(i+2)))-x*(2.d0*zi(i+2)*zi(i+3)+4.d0 &
            *zi(i+2)*zi(i)+zi(i)*zi(i+3)+2.d0*zi(i)*zi(i)))/(c2*c1))
            m2m(i) = 192.d0*(((x3-(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i+1) &
           +zi(i-2)))+x*(2.d0*zi(i+1)*zi(i)+zi(i-2)*zi(i)))/(a2*a0)) &
           +((x3-(0.5d0*x2*(3.d0*zi(i)+zi(i+2)+zi(i-1)+zi(i+1))) &
            +x*(zi(i+2)*zi(i)+zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(b2*a0)) &
           +((x3-(0.5d0*x2*(4.d0*zi(i)+2.d0*zi(i+2)))+x*(2.d0*zi(i+2) &
           *zi(i)+zi(i)*zi(i)))/(c2*a0)) )
            m1m(i) = 192.d0*(((-x3+(0.5d0*x2*(3.d0*zi(i)+2.d0*zi(i-1) &
           +zi(i+1)))-x*(2.d0*zi(i-1)*zi(i)+zi(i+1)*zi(i)))/(a1*a0)) &
           +((-x3+(0.5d0*x2*(4.d0*zi(i)+zi(i-1)+zi(i+2))) &
            -x*(zi(i-1)*zi(i)+zi(i)*zi(i)+zi(i+2)*zi(i)))/(b1*a0)) &
           +((-x3+(0.5d0*x2*(5.d0*zi(i)+zi(i+3)))-x*(zi(i+3)*zi(i) &
           +2.d0*zi(i)*zi(i)))/(c1*a0)) )
        end do

          end subroutine vecpen


!================================  QGAUS1   ==========================
      subroutine qgaussPL(cas,a,b,the01,the12,the02,res,v1,v2,v3)

        use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

         implicit none
         
         integer::j,cas
         double precision::a,b,dx,xm,xr,res,v1,v2,v3
         double precision,dimension(-2:(nz01-1))::the01
        double precision,dimension(-2:(nz12-1))::the12
        double precision,dimension(-2:(nz02-1))::the02
         double precision,dimension(5)::w,x
         double precision::xx,f1,su01,ri01,ri12,f2,su12,su02,ri02
         save w,x
         data w/0.2955242247d0,0.2692667193d0,0.2190863625d0, &
               0.1494513491d0,0.0666713443d0/
         data x/0.1488743389d0,0.4333953941d0,0.6794095682d0, &
               0.8650633666d0,0.9739065285d0/
         
         
         if((cas.eq.2).or.(cas.eq.6))then
            xm = 0.5d0*(b+a)
            xr = 0.5d0*(b-a)
            res = 0.d0
            do j=1,5
               dx=xr*x(j)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01)
               call susp(xx,the02,nz02,su02,ri02,zi02)
               call susp(xx,the12,nz12,su12,ri12,zi12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01)
               call susp(xx,the02,nz02,su02,ri02,zi02)
               call susp(xx,the12,nz12,su12,ri12,zi12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               res = res + w(j)*(f1+f2)
            end do
            res = res*xr
         else
            xm = 0.5d0*(b+a)
            xr = 0.5d0*(b-a)
            res = 0.d0
            do j=1,5
               dx=xr*x(j)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01)
               call susp(xx,the02,nz02,su02,ri02,zi02)
               call susp(xx,the12,nz12,su12,ri12,zi12)
               f1 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01)
               call susp(xx,the02,nz02,su02,ri02,zi02)
               call susp(xx,the12,nz12,su12,ri12,zi12)
               f2 = (su01**v1)*(su02**v3)*ri01*v1/(su12**v2)
               res = res + w(j)*(f1+f2)
            end do
            res = res*xr
         endif   

         end subroutine qgaussPL


!======================  LUDCMP  ======================================
        subroutine ludcmp(a,n,indx,d)

        implicit none
        
        integer::n,i,imax,j,k
        integer,dimension(n)::indx
        double precision::d,aamax,dum,sum
        double precision,dimension(n,n)::a
        integer,parameter::nmax=500
        double precision,parameter::tiny=1.d-20
        double precision,dimension(nmax)::vv

        d = 1.d0
        imax=0
        j=0
        k=0
        
        do i=1,n
                aamax=0.d0
                do j=1,n
                        if (dabs(a(i,j)).gt.aamax)then
                                aamax=dabs(a(i,j))
                        endif
                end do
                
                if (aamax.eq.0.d0) then
!                       write(*,*)'matrice singuliere'
                end if
                vv(i) = 1.d0/aamax
        end do
        
        do j = 1,n
                do i=1,j-1
                        sum = a(i,j)
                        do k=1,i-1
                                sum = sum - a(i,k)*a(k,j)
                        end do
                        a(i,j) = sum
                end do
                aamax = 0.d0
                do i = j,n
                        sum = a(i,j)
                        
                        do k=1,j-1
                                sum = sum -a(i,k)*a(k,j)
                        end do
                        
                        a(i,j) = sum
                        dum = vv(i)*dabs(sum)
                        
                        if (dum.ge.aamax) then
                                imax = i
                                aamax = dum
                        endif
                end do
                
                if(j.ne.imax)then
                        do k=1,n
                                dum = a(imax,k)
                                a(imax,k)=a(j,k)
                                a(j,k) = dum
                        end do
                        d = -d
                        vv(imax)=vv(j)
                endif
                
                indx(j)=imax
                if(a(j,j).eq.0.d0)then
                        a(j,j)=tiny
                endif
                if(j.ne.n)then
                        dum = 1.d0/a(j,j)
                        do i = j+1,n
                                a(i,j) = a(i,j)*dum
                        end do
                endif
        end do

        return

        end subroutine ludcmp
!======================  LUBKSB  ======================================

        subroutine lubksb(a,n,indx,b)
      
        implicit none
        
        integer::n,i,ii,j,ll
        double precision::sum
        integer,dimension(n)::indx
        double precision,dimension(n)::b
        double precision,dimension(n,n)::a

        ii = 0
        
        do i=1,n
        
                ll = indx(i)
                sum = b(ll)
                b(ll) = b(i)
                
                if(ii.ne.0)then
                        do j=ii,i-1
                                sum = sum -a(i,j)*b(j)
                        end do
                else
                        if(sum.ne.0.d0)then
                                ii=i
                        endif
                endif
                
                b(i)=sum
                
        end do
        
        do i=n,1,-1
                sum = b(i)
                do j = i+1,n
                        sum = sum-a(i,j)*b(j)
                end do
                b(i)=sum/a(i,i)
        end do
        
        return
        
        end subroutine lubksb
!==========================   ISP   ==================================
          double precision function isp(x,ni,ns,zi,n)
             implicit none
             integer::ni,ns,n
             double precision::val,mmsp,x
             double precision,dimension(-2:(n+1)):: zi       


          if(x.eq.zi(ni))then
             if(ni.le.ns-3)then
                val = 0.d0
             else
                if(ni.le.ns-2)then
           val = ((zi(ni)-zi(ni-1))*mmsp(x,ni,ns,zi,n))*0.25d0
                else
                   if (ni.eq.ns-1)then
                  val = ((zi(ni)-zi(ni-2))*mmsp(x,ni,ns,zi,n)+ &
              (zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+1,zi,n))*0.25d0
                   else
                      if(ni.eq.ns)then
                  val = ((zi(ni)-zi(ni-3))*mmsp(x,ni,ns,zi,n)+ &
                      (zi(ni+2)-zi(ni-2))*mmsp(x,ni,ns+1,zi,n) &
             +(zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+2,zi,n))*0.25d0
                      else
                         val = 1.d0
                      endif
                   endif
                endif   
             endif
          else   
          if(ni.lt.ns-3)then
             val = 0.d0
          else
             if(ni.eq.ns-3)then
                 val = (x-zi(ni))*mmsp(x,ni,ns,zi,n)*0.25d0
             else  
             if(ni.eq.ns-2)then
                   val = ((x-zi(ni-1))*mmsp(x,ni,ns,zi,n)+ &
             (zi(ni+4)-zi(ni))*mmsp(x,ni,ns+1,zi,n))*0.25d0
             else   
                if (ni.eq.ns-1)then
                   val =((x-zi(ni-2))*mmsp(x,ni,ns,zi,n)+ &
                   (zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+1,zi,n) &
             +(zi(ni+4)-zi(ni))*mmsp(x,ni,ns+2,zi,n))*0.25d0
                else
                   if(ni.eq.ns)then
                      val =((x-zi(ni-3))*mmsp(x,ni,ns,zi,n)+ &
                   (zi(ni+2)-zi(ni-2))*mmsp(x,ni,ns+1,zi,n) &
                   +(zi(ni+3)-zi(ni-1))*mmsp(x,ni,ns+2,zi,n) &
              +(zi(ni+4)-zi(ni))*mmsp(x,ni,ns+3,zi,n))*0.25d0
                   else
                      val = 1.d0
                   endif
                endif
             endif
             endif
          endif 
          endif
             isp = val

             return

             end function isp

!==========================  MMSP   ==================================

          double precision function mmsp(x,ni,ns,zi,n)
             integer::ni,ns,n
             double precision::val,x
             double precision,dimension(-2:(n+1)):: zi       


          if(ni.lt.ns-3)then
             val = 0.d0
          else
             if(ns-3.eq.ni)then
                if(x.eq.zi(ni))then
                   val = 0.d0
                else  
                   val = (4.d0*(x-zi(ni))*(x-zi(ni)) &
                   *(x-zi(ni)))/((zi(ni+4)-zi(ni))*(zi(ni+3) &
                   -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
                endif
             else 
                if(ns-2.eq.ni)then
                   if(x.eq.zi(ni))then
                      val = (4.d0*(zi(ni)-zi(ni-1))*(zi(ni)-zi(ni-1))) &
                   /((zi(ni+3)-zi(ni-1))*(zi(ni+2)-zi(ni-1)) &
                   *(zi(ni+1)-zi(ni-1)))
                   else  
                      val = (4.d0*(x-zi(ni-1))*(x-zi(ni-1)) &
                   *(zi(ni+1)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
                   -zi(ni-1))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)-zi(ni))) &
                   +   (4.d0*(x-zi(ni-1))*(x-zi(ni)) &
                   *(zi(ni+2)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+2) &
                   -zi(ni))*(zi(ni+1)-zi(ni))*(zi(ni+2)-zi(ni-1))) &
                   +   (4.d0*(x-zi(ni))*(x-zi(ni)) &
                   *(zi(ni+3)-x))/((zi(ni+3)-zi(ni-1))*(zi(ni+3) &
                   -zi(ni))*(zi(ni+2)-zi(ni))*(zi(ni+1)-zi(ni)))
                   endif
                else  
                   if (ns-1.eq.ni)then
                      if(x.eq.zi(ni))then
                         val = (4.d0*((zi(ni)-zi(ni-2))*(zi(ni+1) &
                       -zi(ni)))/((zi(ni+2)-zi(ni-2))*(zi(ni+1) &
                       -zi(ni-1))*(zi(ni+1)-zi(ni-2)))) &
                      +((4.d0*((zi(ni)-zi(ni-1))*(zi(ni+2)-zi(ni))) &
                      /((zi(ni+2)-zi(ni-2))*(zi(ni+2)-zi(ni-1)) &
                      *(zi(ni+1)-zi(ni-1)))))
                      else
                        val = (4.d0*((x-zi(ni-2))*(zi(ni+1) &
                       -x)*(zi(ni+1)-x))/((zi(ni+2) &
                       -zi(ni-2))*(zi(ni+1)-zi(ni-1))*(zi(ni+1)- &
                        zi(ni))*(zi(ni+1)-zi(ni-2)))) &
                      +((4.d0*((x-zi(ni-1))*(zi(ni+2)-x)  & 
                      *(zi(ni+1)-x))/((zi(ni+2)-zi(ni-2)) & 
                      *(zi(ni+2)-zi(ni-1))*(zi(ni+1)-zi(ni-1))* &
                      (zi(ni+1)-zi(ni))))) &
                      +((4.d0*((zi(ni+2)-x)*(zi(ni+2)-x)  &
                      *(x-zi(ni)))/((zi(ni+2)-zi(ni-2)) &
                      *(zi(ni+2)-zi(ni))*(zi(ni+2)-zi(ni-1))* &
                      (zi(ni+1)-zi(ni)))))
                      endif 
                   else
                      if(ni.eq.ns)then
                         if(x.eq.zi(ni))then
                            val =(4.d0*(x-zi(ni+1))*(x &
                         -zi(ni+1))/((zi(ni+1)-zi(ni-1))*(zi(ni+1) &
                         -zi(ni-2))*(zi(ni+1)-zi(ni-3)))) 
                         else   
                           val =(4.d0*(x-zi(ni+1))*(x &
                           -zi(ni+1))*(zi(ni+1)-x)/((zi(ni+1) &
                           -zi(ni-1))*(zi(ni+1)-zi(ni-2))*(zi(ni+1) &
                           -zi(ni))*(zi(ni+1)-zi(ni-3))))
                         endif
                      else
                         val = 0.d0
                      endif
                   endif
                endif
             endif
          endif

             mmsp = val
             return
             end function mmsp
!==========================  SUSP  ====================================
        subroutine susp(x,the,n,su,lam,zi)

        implicit none
        
        integer::j,k,n,i
        double precision::x,ht,ht2,h2,som,lam,su,htm,h2t,h3,h2n,hn, &
        im,im1,im2,mm1,mm3,ht3,hht,h4,h3m,hh3,hh2,mm,im3,mm2,h,gl,hh
        double precision,dimension(-2:(n+3))::zi
        double precision,dimension(-2:n-1)::the 

        som = 0.d0
        gl = 0.d0 
        
        do k = 2,n
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
                        gl = som +(the(j-3)*im3)+(the(j-2)*im2)+(the(j-1)*im1)+(the(j)*im)
                        lam = (the(j-3)*mm3)+(the(j-2)*mm2)+(the(j-1)*mm1)+(the(j)*mm)
                endif
        end do
   
        if(x.ge.zi(n))then
                som = 0.d0
                do i=1,n+2
                        som = som+the(i-3)
                end do
                gl = som
                lam = 4.d0*the(n-1)/(zi(n)-zi(n-1))
        endif

        su  = dexp(-gl)

        return

        end subroutine susp


!=======================  CALOMEG  ===========================
        subroutine calcomeg(n,omeg,m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m)
        
        implicit none
        !        remplissage de la matrice omega n*n
        !          elle a 7 diagonales
        integer::n,i,j
        double precision::calc00,calc01,calc02
        double precision,dimension(n,n)::omeg
        
        double precision,dimension(n-2)::m3m3,m2m2,m1m1,mmm,m3m2,m3m1, &     
        m3m,m2m1,m2m,m1m
        
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
        
        end subroutine calcomeg

!=========================  CALC00  =========================

        double precision function calc00(j,n,m3m3,m2m2,m1m1,mmm)
        
        implicit none
        
        integer::j,n
        double precision::part             
        double precision,dimension(n-2)::m3m3,m2m2,m1m1,mmm

        if(j.eq.1)then
                part = m3m3(j)
        else
                if(j.eq.2)then
                        part = m3m3(j) + m2m2(j-1)
                else
                        if(j.eq.3)then
                                part = m3m3(j) + m2m2(j-1) + m1m1(j-2)
                        else
                                if(j.eq.n-2)then
                                        part = m2m2(j-1) + m1m1(j-2) + mmm(j-3)
                                else   
                                        if(j.eq.n-1)then
                                                part = mmm(j-3) + m1m1(j-2)
                                        else
                                                if(j.eq.n)then
                                                        part = mmm(j-3)
                                                else   
                                                        part=mmm(j-3)+m1m1(j-2)+m2m2(j-1)+m3m3(j)
                                                endif
                                        endif
                                endif   
                        endif   
                endif   
        endif 

        calc00 = part

        return

        end function calc00
!=========================  CALC01  =========================
        double precision function calc01(j,n,m3m2,m2m1,m1m)
        
        implicit none
        
        integer::j,n
        double precision::part        
        double precision,dimension(n-2)::m3m2,m2m1,m1m

!        entre j+1 et j+2 ---> m1*m et       i = j-2
!        entre j+2 et j+3 ---> m2*m1 et    i = j -1
!        entre j+3 et j+4 ---> m3*m2 et     i = j 

        if(j.eq.1)then
                part = m3m2(j)
        else   
                if(j.eq.2)then
                        part = m3m2(j) + m2m1(j-1) 
                else
                        if(j.eq.n-2)then
                                part = m1m(j-2) + m2m1(j-1) 
                        else
                                if(j.ne.n-1)then
                                        part = m3m2(j) + m2m1(j-1) + m1m(j-2)
                                else
                                        part = m1m(j-2)
                                endif
                        endif   
                endif
        endif   

        calc01 = part

        return

        end function calc01
!=========================  CALC02  =========================
        double precision function calc02(j,n,m3m1,m2m)
        
        implicit none
        
        double precision::part
        integer::j,n
        double precision,dimension(n-2)::m3m1,m2m

!        entre j+3 et j+4 ---> m2*m et     i = j-1
!        entre j+2 et j+3 ---> m3*m1 et       i = j

        if(j.eq.1)then
                part = m3m1(j)
        else   
                if(j.ne.n-2)then
                part = m3m1(j) + m2m(j-1) 
                else
                part = m2m(j-1)
                endif
        endif   

        calc02 = part
        
        return
                
        end function calc02


!================================  QGAUS12   =========================
      subroutine qgauss12(cas,a,b,the01,the12,the02,res)

        use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

         implicit none
         
         integer::j,cas
         double precision::a,b,dx,xm,xr,res,xx,f1,su01,ri01,ri12,f2,su12,su02,ri02
         double precision,dimension(-2:(nz01+2))::the01
         double precision,dimension(-2:(nz12+2))::the12
         double precision,dimension(-2:(nz02+2))::the02
         double precision,dimension(5)::w,x

         save w,x
         data w/0.2955242247d0,0.2692667193d0,0.2190863625d0, &
               0.1494513491d0,0.0666713443d0/
         data x/0.1488743389d0,0.4333953941d0,0.6794095682d0, &
               0.8650633666d0,0.9739065285d0/
         
         
         if((cas.eq.2).or.(cas.eq.6))then
            xm = 0.5d0*(b+a)
            xr = 0.5d0*(b-a)
            res = 0.d0
            do j=1,5
               dx=xr*x(j)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01)
               call susp(xx,the02,nz02,su02,ri02,zi02)
               call susp(xx,the12,nz12,su12,ri12,zi12)
               f1 = (su01*su02*ri01)/(su12)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01)
               call susp(xx,the02,nz02,su02,ri02,zi02)
               call susp(xx,the12,nz12,su12,ri12,zi12)
               f2 = (su01*su02*ri01)/(su12)
               res = res + w(j)*(f1+f2)
            end do
            res = res*xr
         else
            xm = 0.5d0*(b+a)
            xr = 0.5d0*(b-a)
            res = 0.d0
            do j=1,5
               dx=xr*x(j)
               xx = xm+dx
               call susp(xx,the01,nz01,su01,ri01,zi01)
               call susp(xx,the02,nz02,su02,ri02,zi02)
               call susp(xx,the12,nz12,su12,ri12,zi12)
               f1 = (su01*su02*ri01)/(su12)
               xx = xm-dx
               call susp(xx,the01,nz01,su01,ri01,zi01)
               call susp(xx,the02,nz02,su02,ri02,zi02)
               call susp(xx,the12,nz12,su12,ri12,zi12)
               f2 = (su01*su02*ri01)/(su12)
               res = res + w(j)*(f1+f2)
            end do
            res = res*xr
         endif   

         end subroutine qgauss12


!================================  QGAUS13   =========================
        subroutine qgauss13(cas,a,b,the01,the12,the02,res)

        use commun,only:zi01,zi12,zi02,nz01,nz12,nz02

        implicit none
        
        integer::j,cas
        double precision::a,b,res,h,xx,f1,su01,ri01,ri12,f2,su12,su02,ri02
        double precision,dimension(-2:(nz01+2))::the01
        double precision,dimension(-2:(nz12+2))::the12
        double precision,dimension(-2:(nz02+2))::the02

        xx = a
        res = 0.d0
        h = 0.01d0*(b-a)
        
        if((cas.eq.2).or.(cas.eq.6))then
                call susp(xx,the01,nz01,su01,ri01,zi01)
                call susp(xx,the02,nz02,su02,ri02,zi02)
                call susp(xx,the12,nz12,su12,ri12,zi12)
                f1 =  (su01*su02*ri01)/(su12)
                do j=1,100
                        xx = xx + h
                        call susp(xx,the01,nz01,su01,ri01,zi01)
                        call susp(xx,the02,nz02,su02,ri02,zi02)
                        call susp(xx,the12,nz12,su12,ri12,zi12)
                        f2 = (su01*su02*ri01)/(su12)
                        res = res + 0.5d0*h*(f1+f2)
                        f1 = f2
                end do
        else
                call susp(xx,the01,nz01,su01,ri01,zi01)
                call susp(xx,the02,nz02,su02,ri02,zi02)
                call susp(xx,the12,nz12,su12,ri12,zi12)
                f1 = (su01*su02*ri01)/(su12)
                do j=1,100
                        xx = xx + h
                        call susp(xx,the01,nz01,su01,ri01,zi01)
                        call susp(xx,the02,nz02,su02,ri02,zi02)
                        call susp(xx,the12,nz12,su12,ri12,zi12)
                        f2 = (su01*su02*ri01)/(su12)
                        res = res + 0.5d0*h*(f1+f2)
                        f1 = f2
                end do
        endif   

        end subroutine qgauss13

