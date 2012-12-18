!            Survival model with non censored or interval-censored or right censored
!                               left truncated
!                                 Weibull
!                              20/05/10
!            last              08/04/11
         subroutine survWeib(entrytime,l,r,status,x,N,P,truncated,interval,eps,&
                             maxiter0,loglik,basepar,regpar,v,converged,&
                             cv,niter,t,S,S_l,S_u,h,h_l,h_u)
     
         use survCommun  
         use parameters
	 use optim
         implicit none
!  variables entrants
	 integer,intent(in)::n,p
	 double precision,dimension(n),intent(in)::l,r,entrytime
	 integer,dimension(n),intent(in)::status
	 double precision,dimension(n*p),intent(in)::x
	 integer,intent(in)::truncated,interval,maxiter0
	 integer,dimension(3),intent(in)::eps
!  variables sortant
         integer,intent(out)::converged,niter 
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
         double precision,dimension(P+2,P+2)::hes
	 double precision::ca,cb,dd
         integer::ier,istop,npw,irep,np,j,k,jj,ii,ij,i,ic,id,ide, kk,kkk,rec,nn
     
   
!---------------------------------
                  
         
         model=1
         no=N
         nva=P
         troncature = truncated
         epsa = 0.1d0**eps(1)
         epsb = 0.1d0**eps(2)
         epsd = 0.1d0**eps(3)
         ! write(*,*)epsa,epsb,epsd,eps
         maxiter = maxiter0

         min = l(1)
         max = l(1)
         if(truncated.eq.1) then
            if(interval.eq.1) then
               allocate(t0(no),t1(no),t2(no),ve(no,nva),c(no))
               do 12 i=1,no
                  t0(i)=entrytime(i)
                  t1(i)=l(i)
                  t2(i)=r(i)
                  c(i)=status(i)
                  do 11 j=1,nva
                     k = (j-1)*no+i
                     ve(i,j)=x(k)
 11               continue   
               if (max.lt.t2(i))max = t2(i)
               if (min.gt.t0(i))min = t0(i)
 12            continue 
             else  
               allocate(t0(no),t1(no),t2(1),ve(no,nva),c(no))
               do 14 i=1,no
                  t0(i)=entrytime(i)
                  t1(i)=l(i)
                  c(i)=status(i)
                  do 13 j=1,nva
                     k = (j-1)*no+i
                     ve(i,j)=x(k)
 13               continue 
                 if (max.lt.t1(i))max = t1(i)
                 if (min.gt.t0(i))min = t0(i)
 14            continue 
            endif   
         else
            if(interval.eq.1) then
               allocate(t0(1),t1(no),t2(no),ve(no,nva),c(no))
               min = 0.d0
               do 16 i=1,no
                  t1(i)=l(i)
                  t2(i)=r(i)
                  c(i)=status(i)
                  do 15 j=1,nva
                     k = (j-1)*no+i
                     ve(i,j)=x(k)
 15               continue   
                  if (max.lt.t2(i))max = t2(i)
 16            continue 
             else  
               allocate(t0(1),t1(no),t2(1),ve(no,nva),c(no))
               do 18 i=1,no
                  t1(i)=l(i)
                  c(i)=status(i)
                  do 17 j=1,nva
                     k = (j-1)*no+i
                     ve(i,j)=x(k)
 17               continue   
                 if (max.lt.t1(i))max = t1(i)
 18            continue 
            endif
         endif   
   
 
         npw = 2

         do 75 i=1,npw
            b(i)=1.d-1
 75      continue

         np = npw
         call marq98(b,np,niter,v1,res,ier,istop,ca,cb,dd)
         loglik(1) = res
         if((nva).gt.0)then
            np = npw + nva
            do 122 i=npw+1,np
               b(i)=0.d0
 122        continue  
            call marq98(b,np,niter,v1,res,ier,istop,ca,cb,dd)
!c           write(*,*)'niter2',niter
            do 124 i=npw+1,np
               regpar(i-npw) = b(i)
 124        continue 
            do 305 i=1,np
                do 304 j=i,np
                   k = (((j-1)*j)/2) +i
                   hes(i,j)= v1(k)
 304           continue
 305        continue
            do  307 i = 2,np
            do 306 j = 1,i-1
               hes(i,j)=hes(j,i)
 306        continue
 307        continue
            K = 0
            do  309 i = 3,np
               do 308 j = 3,np
                  K=K+1
                  v(k)=hes(j,i)
 308           continue
 309        continue
         endif   



         loglik(2) = res
         converged = istop
         do 127 i=1,npw
           basepar(i) = dexp(b(i))
 127    continue

        cv(1) = ca
        cv(2) = cb
        cv(3) = dd

      deallocate(t0,t1,t2,ve,c)

!----------------- sorties graphiques ------------------------

          do 390 i=1,2
           the(i)=dexp(b(i))
 390      continue


         pas = (max-min)*0.01d0
         tx = min
         do 391 i=1,100
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
 391      continue

!------------------------- calcul de U  avec V1=UTU ---------------
         ep=10.d-10
         call dmfsd(v1,np,ep,ier)


         do 393 i=1,np
            do 392 j=1,np
               vsup(i,j)=0.d0
               vinf(i,j)=0.d0
 392        continue
 393     continue   



          do 395 i=1,np
             do 394 j=i,np
                k = (((j-1)*j)/2) +i
                vsup(i,j)= v1(k)
 394        continue
 395        continue
       do  397 i = 1,np
            do 396 j = 1,i
               vinf(i,j)=vsup(j,i)
 396       continue
 397    continue




!---------------- bootstrap ------------------------------

       do 410 jj=1,2000

         do 401 i=1,np
            call BGOS(1.d0,0,x1,x2,0.d0)
            xi(i)=x1
 401     continue   

 

         do 403 i=1,np
            ut(i)=0.d0
            do 402 j=0,np
              ut(i) = ut(i)+(vinf(i,j)*xi(j))
 402        continue   
 403     continue 
 
 
         do 404 i=1,np
            bh(i)= ut(i)+b(i)
 404     continue 


         do 405 i=1,2
           the(i)=dexp(bh(i))
 405      continue

         pas = (max-min)*0.01d0
         tx = min
         do 406 i=1,100
            if (tx.eq.0.d0)then
               tx0 = 0.1d0*pas
               call fonct(tx0,the,ri,gl,su) 
            else   
               call fonct(tx,the,ri,gl,su)  
            endif   
            mate_su(jj,i) = su
            mate_ri(jj,i) = ri
            tx = tx + pas
 406     continue


 410     continue
    

         do 460 k = 1,100
            do 421 kk=1,51
               tab_su_i(kk)  = 10.d0
               tab_su_s(kk)  = 0.d0
               tab_ri_i(kk)  = 10.d0
               tab_ri_s(kk)  = 0.d0
 421        continue
            do 440 i=1,2000
               do 425 kk=1,51
                  if(mate_su(i,k).lt.tab_su_i(kk))then
                     if(kk.lt.51)then
                        do 424 kkk = 51,kk+1,-1
                           tab_su_i(kkk) = tab_su_i(kkk-1)
 424                    continue   
                        tab_su_i(kk) = mate_su(i,k)
                        go to 426
                     else
                        tab_su_i(kk) = mate_su(i,k)
                        go to 426
                     endif
                  endif 
 425           continue
 426           continue   
 
               do 429 kk=1,51
                  if(mate_ri(i,k).lt.tab_ri_i(kk))then
                     if(kk.lt.51)then
                        do 428 kkk = 51,kk+1,-1
                           tab_ri_i(kkk) = tab_ri_i(kkk-1)
 428                       continue   
                        tab_ri_i(kk) = mate_ri(i,k)
                        go to 430
                     else
                        tab_ri_i(kk) = mate_ri(i,k)
                        go to 430
                     endif
                  endif 
 429           continue
 430           continue

            do 435 kk=51,1,-1
                  if(mate_su(i,k).gt.tab_su_s(kk))then
                     if(kk.gt.1)then
                        do 434 kkk = 1,kk-1
                           tab_su_s(kkk) = tab_su_s(kkk+1)
 434                    continue   
                        tab_su_s(kk) = mate_su(i,k)
                        go to 436
                     else
                        tab_su_s(kk) = mate_su(i,k)
                        go to 436
                     endif
                  endif 
 435           continue 
 436           continue

            do 438 kk=51,1,-1
                  if(mate_ri(i,k).gt.tab_ri_s(kk))then
                     if(kk.gt.1)then
                        do 437 kkk = 1,kk-1
                           tab_ri_s(kkk) = tab_ri_s(kkk+1)
 437                    continue   
                        tab_ri_s(kk) = mate_ri(i,k)
                        go to 439
                     else
                        tab_ri_s(kk) = mate_ri(i,k)
                        go to 439
                     endif
                  endif 
 438              continue 
 439           continue

 440        continue
            S_l(k) = tab_su_i(51)
            S_u(k) = tab_su_s(1)
            h_l(k) = tab_ri_i(51)
            h_u(k) = tab_ri_s(1)
 460     continue  

         end subroutine survWeib








