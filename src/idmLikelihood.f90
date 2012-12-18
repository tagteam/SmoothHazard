!======================== idmLikelihood ====================

      module idmCommun

      implicit none
      integer,save ::no,nva01,nva02,nva12,troncature
      double precision,dimension(:),allocatable,save::t0,t1,t2,t3
      double precision,dimension(:,:),allocatable,save::ve01,ve02,ve12
      integer,dimension(:),allocatable,save :: c

      end module idmCommun

      double precision function idmLikelihood(b,np,id,thi,jd,thj)
         use idmCommun
         implicit none
         integer,intent(in)::np,id,jd
	 double precision,intent(in)::thi,thj
         double precision,dimension(np),intent(in)::b
	 double precision,dimension(np)::bh
         double precision,dimension(2)::the01,the02,the12
         integer::i,j,k
         double precision::res,res1,res2,vet01,vet02,vet12,tronc
         double precision::ri01,gl01,su01,ri02,gl02,su02,ri12,gl12,su12


         do 2 i=1,np
            bh(i)=b(i)
 2       continue

         if (id.ne.0) bh(id)=bh(id)+thi
         if (jd.ne.0) bh(jd)=bh(jd)+thj    


         do 3 i=1,2
!            the01(i)=dexp(bh(i))
            the01(i)=bh(i)*bh(i)
 3       continue
         do 4 i=1,2
            j = 2+i
!            the02(i)=dexp(bh(j))
            the02(i)=bh(j)*bh(j)
 4       continue
         do 5 i=1,2
            j = 4 + i
!            the12(i)=dexp(bh(j))
            the12(i)=bh(j)*bh(j)
 5       continue


!---------- calcul de la vraisemblance ------------------

          res = 0.d0
          do 15 i=1,no

             vet01 = 0.d0
             vet02 = 0.d0
             vet12 = 0.d0
             if(nva01.gt.0)then
                do 12 j=1,nva01
                   vet01 =vet01 + bh(np-nva01-nva02-nva12+j)*dble(ve01(i,j))
 12             continue
             endif   
             if(nva02.gt.0)then
                do 13 j=1,nva02
                   vet02 =vet02 + bh(np-nva02-nva12+j)*dble(ve02(i,j))
 13             continue
             endif
             if(nva12.gt.0)then
                do 14 j=1,nva12
                   vet12 =vet12 + bh(np-nva12+j)*dble(ve12(i,j))
 14             continue
             endif
             vet01 = dexp(vet01)
             vet02 = dexp(vet02)
             vet12 = dexp(vet12)
             if(troncature.eq.1)then
                if(t0(i).eq.0.d0)then
                   tronc = 0.d0
                else   
                   call fonct(t0(i),the01,ri01,gl01,su01) 
                   call fonct(t0(i),the02,ri02,gl02,su02) 
                   tronc = (gl01*vet01)+(gl02*vet02)
                endif
             else
                tronc = 0.d0
             endif   
             res1 = 0.d0
!           censure à droite 01 et 02
             if(c(i).eq.1)then   
                call fonct(t1(i),the01,ri01,gl01,su01)
                call fonct(t3(i),the02,ri02,gl02,su02)
                res1 = -(gl01*vet01)-(gl02*vet02)
             else
!           censure par intervalle 01 et droite pour 12
                if(c(i).eq.2)then  
                    call fonct(t3(i),the12,ri12,gl12,su12)
                    call  qgauss1(t1(i),t2(i),the01,the02,the12,res2,vet01,vet02,vet12)
                    res1=dlog(res2*(su12**vet12))
                else
!           observation 01 et censure à droite pour 12
                   if(c(i).eq.3)then  
                      call fonct(t1(i),the01,ri01,gl01,su01)
                      call fonct(t1(i),the02,ri02,gl02,su02)
                      call fonct(t1(i),the12,ri12,gl12,su12)
                      res1 = -(gl01*vet01)-(gl02*vet02)+dlog(ri01*vet01)+(gl12*vet12)
                      call fonct(t3(i),the12,ri12,gl12,su12)
                      res1 = res1 -(gl12*vet12)
                   else   
!           censure par intervalle 01 et observation pour 12
                      if(c(i).eq.4)then  
                         call fonct(t3(i),the12,ri12,gl12,su12)
                         call  qgauss1(t1(i),t2(i),the01,the02,the12,res2,vet01,vet02,vet12)
                         res1=dlog(res2*(su12**vet12)*ri12*vet12)
                      else
!           observation 01 et observation pour 12
                         if(c(i).eq.5)then 
                            call fonct(t1(i),the01,ri01,gl01,su01)
                            call fonct(t1(i),the02,ri02,gl02,su02)
                            call fonct(t1(i),the12,ri12,gl12,su12)
                            res1 = -(gl01*vet01)-(gl02*vet02)+dlog(ri01*vet01)+(gl12*vet12)
                            call fonct(t3(i),the12,ri12,gl12,su12)
                            res1 = res1 -(gl12*vet12) + dlog(ri12*vet12)
                         else
!                alive   
                            if(c(i).eq.6)then
                               call fonct(t3(i),the01,ri01,gl01,su01)
                               call fonct(t3(i),the02,ri02,gl02,su02)
                               call fonct(t3(i),the12,ri12,gl12,su12)
                               call  qgauss1(t1(i),t3(i),the01,the02,the12,res2,vet01,vet02,vet12)
                               res1 = (res2*(su12**vet12))+((su01**vet01)*(su02**vet02))
                               res1 = dlog(res1)
                            else
!                passage 02 ?
                               call fonct(t3(i),the01,ri01,gl01,su01)
                               call fonct(t3(i),the02,ri02,gl02,su02)
                               call fonct(t3(i),the12,ri12,gl12,su12)
                               call  qgauss1(t1(i),t3(i),the01,the02,the12,res2,vet01,vet02,vet12)
                               res1 = (res2*(su12**vet12)*ri12*vet12)+&
                                ((su01**vet01)*(su02**vet02)*ri02*vet02)
                               res1 = dlog(res1)
                            endif
                         endif                        
                      endif
                   endif   
                endif   
             endif   
             res = res + res1 + tronc
 15       continue   

          idmLikelihood = res
          return
          end function idmLikelihood
