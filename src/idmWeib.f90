!=============================================================================================
!            Illness-death model with non censored or interval-censored or right censored
!                               left truncated
!                                 Weibull
!                              19/04/11
!            last              20/04/11
!=============================================================================================
       subroutine idmWeib(entrytime,l,r,d,idm,idd,x01,x02,x12,N,P01,P02,P12,truncated,eps &
                      ,maxiter0,loglik,basepar,regpar,v,converged,cv,&
                      niter,t,a01,a01_l,a01_u,a02,a02_l,a02_u,a12,a12_l,a12_u) 
     
         use idmCommun  
         use parameters
	 use optim
         implicit none
!  variables entrants
	 integer,intent(in)::n,p01,p02,p12
	 double precision,dimension(n),intent(in)::l,r,d,entrytime
	 integer,dimension(n),intent(in)::idm,idd
	 double precision,dimension(n*p01),intent(in)::x01
	 double precision,dimension(n*p02),intent(in)::x02
	 double precision,dimension(n*p12),intent(in)::x12
	 integer,intent(in)::truncated,maxiter0
	 integer,dimension(3),intent(in)::eps
!  variables sortant
         integer,intent(out)::converged,niter 
         double precision,dimension(p01+p02+p12),intent(out)::regpar

         double precision,dimension(100),intent(out)::t,a01,a01_l,a01_u,a02,a02_l,a02_u,a12,a12_l,a12_u
	 double precision,dimension(2),intent(out)::loglik
	 double precision,dimension(6),intent(out)::basepar
	 double precision,dimension(3),intent(out)::cv
	 double precision,dimension((p01+p02+p12)*(p01+p02+p12)),intent(out)::v
	 
!  variables locales	 
	 
         
         double precision:: res,tx0,min,max,x1,x2
	 double precision,dimension((p01+p02+p12)+6)::xi,ut,bh
	 double precision,dimension(2)::the01,the02,the12
	 double precision,dimension((p01+p02+p12)+6,(p01+p02+p12)+6)::vinf,vsup
         double precision::su,gl,ri,tx,pas,ep,moyenne
         double precision,dimension(2000,100)::mate_ri01,mate_ri02,mate_ri12
        double precision,dimension(2000)::vect
	 double precision,dimension(52)::tab_su_i,tab_su_s,tab_ri_i,tab_ri_s           
	 double precision,dimension((p01+p02+p12)+6)::b
	 double precision,dimension(((p01+p02+p12)+6)*((p01+p02+p12)+6+3)/2)::v1
         double precision,dimension((p01+p02+p12)+6,(p01+p02+p12)+6)::hes
	 double precision::ca,cb,dd
         integer::ier,istop,npw,irep,np,j,k,jj,ii,ij,i,ic,id,ide,kk,kkk,rec,nn
        
     
   
!---------------------------------
                  
         
         model=3
         no=N
         nva01=P01
         nva02=P02
         nva12=P12
         troncature = truncated
         epsa = 0.1d0**eps(1)
         epsb = 0.1d0**eps(2)
         epsd = 0.1d0**eps(3)
         maxiter = maxiter0
         min = l(1)
         max = l(1)
         if(truncated.eq.1) then
            allocate(t0(no),t1(no),t2(no),t3(no),ve01(no,nva01),ve02(no,nva02),ve12(no,nva12),c(no))
            do 12 i=1,no
               t0(i)=entrytime(i)
               t1(i)=l(i)
               t2(i)=r(i)
               t3(i)=d(i)
               if (max.lt.t3(i))max = t3(i)
               if (min.gt.t0(i))min = t0(i)
 12         continue 
         else
            allocate(t0(1),t1(no),t2(no),t3(no),ve01(no,nva01),ve02(no,nva02),ve12(no,nva12),c(no))
            min = 0.d0
            do 16 i=1,no
               t1(i)=l(i)
               t2(i)=r(i)
               t3(i)=d(i)
               if (max.lt.t3(i))max = t3(i)
 16         continue 
         endif 
  
         do 24 i=1,no
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
            do 21 j=1,nva01
                  k = (j-1)*no+i
                  ve01(i,j)=x01(k)
 21         continue   
            do 22 j=1,nva02
                  k = (j-1)*no+i
                  ve02(i,j)=x02(k)
 22         continue  
            do 23 j=1,nva12
                  k = (j-1)*no+i
                  ve12(i,j)=x12(k)
 23         continue  
 24      continue

 
         npw = 6

         do 75 i=1,npw
            b(i)=1.d-1
 75      continue
         nva01=0
         nva02=0
         nva12=0
         np = npw
         call marq98(b,np,niter,v1,res,ier,istop,ca,cb,dd)
         loglik(1) = res
!            write(*,*)'niter1',niter,res
         nva01=P01
         nva02=P02
         nva12=P12
         if((nva01+nva02+nva12).gt.0)then
            np = npw + nva01+nva02+nva12
            do 122 i=npw+1,np
               b(i)=0.d0
 122        continue  
               call marq98(b,np,niter,v1,res,ier,istop,ca,cb,dd)
!           write(*,*)'niter2',niter,res
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
            do  309 i = 7,np
               do 308 j = 7,np
                  K=K+1
                  v(k)=hes(j,i)
 308           continue
 309        continue
         endif   



         loglik(2) = res
         converged = istop
         do 127 i=1,npw
           basepar(i) = b(i)*b(i)
 127    continue

        cv(1) = ca
        cv(2) = cb
        cv(3) = dd

      deallocate(t0,t1,t2,t3,ve01,ve02,ve12,c)

!----------------- sorties graphiques ------------------------

          do 390 i=1,2
           the01(i)=b(i)*b(i)
            j = 2+i
            the02(i)=b(j)*b(j)
            j = 4 + i
            the12(i)=b(j)*b(j)
 390      continue
 

         pas = (max-min)*0.01d0
         tx = min
         do 391 i=1,100
             if (tx.eq.0.d0)then
                tx0 = 0.1d0*pas
                call fonct(tx0,the01,ri,gl,su) 
                t(i) = tx0
                a01(i) = ri
                call fonct(tx0,the02,ri,gl,su) 
                a02(i) = ri
                call fonct(tx0,the12,ri,gl,su) 
                a12(i) = ri
             else 
                call fonct(tx,the01,ri,gl,su) 
                t(i) = tx
                a01(i) = ri
                call fonct(tx,the02,ri,gl,su) 
                a02(i) = ri
                call fonct(tx,the12,ri,gl,su) 
                a12(i) = ri
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
           the01(i)=bh(i)*bh(i)
            j = 2+i
            the02(i)=bh(j)*bh(j)
            j = 4 + i
            the12(i)=bh(j)*bh(j)
 405      continue

         pas = (max-min)*0.01d0
         tx = min
         do 406 i=1,100
            if (tx.eq.0.d0)then
               tx0 = 0.1d0*pas
               call fonct(tx0,the01,ri,gl,su) 
            else   
               call fonct(tx,the01,ri,gl,su)  
            endif   
            mate_ri01(jj,i) = ri
            tx = tx + pas
 406     continue

         pas = (max-min)*0.01d0
         tx = min
         do 407 i=1,100
            if (tx.eq.0.d0)then
               tx0 = 0.1d0*pas
               call fonct(tx0,the02,ri,gl,su) 
            else   
               call fonct(tx,the02,ri,gl,su)  
            endif   
            mate_ri02(jj,i) = ri
            tx = tx + pas
 407     continue


         pas = (max-min)*0.01d0
         tx = min
         do 408 i=1,100
            if (tx.eq.0.d0)then
               tx0 = 0.1d0*pas
               call fonct(tx0,the12,ri,gl,su) 
            else   
               call fonct(tx,the12,ri,gl,su)  
            endif   
            mate_ri12(jj,i) = ri
            tx = tx + pas
 408     continue

 410     continue
    

         do 460 k = 1,100


            do 452 jj = 1,2000
               vect(jj) = mate_ri01(jj,k)
 452        continue
             call tri(vect,2000,moyenne)
             a01_u(k) = vect(1950)
             a01_l(k) = vect(51)
 
           do 453 jj = 1,2000
               vect(jj) = mate_ri02(jj,k)
 453        continue
             call tri(vect,2000,moyenne)
             a02_u(k) = vect(1950)
             a02_l(k) = vect(51)


           do 454 jj = 1,2000
               vect(jj) = mate_ri12(jj,k)
 454        continue
             call tri(vect,2000,moyenne)
             a12_u(k) = vect(1950)
             a12_l(k) = vect(51)

 460     continue  

end subroutine idmWeib








