!======================== survLikelihood ====================

      module survCommun

      implicit none
      integer,save ::no,nva,troncature
      double precision,dimension(:),allocatable,save::t0,t1,t2
      double precision,dimension(:,:),allocatable,save::ve
      integer,dimension(:),allocatable,save :: c

      end module survCommun

      double precision function survLikelihood(b,np,id,thi,jd,thj)
         use survCommun
         implicit none
         integer,intent(in)::np,id,jd
	 double precision,intent(in)::thi,thj
         double precision,dimension(np),intent(in)::b
	 
	 double precision,dimension(np)::bh
         double precision,dimension(2)::the
         integer::i,j,k,vj,tronc
         double precision::pe,res,res1,gl,su,ri,vet,su1,su2


         do 2 i=1,np
            bh(i)=b(i)
 2       continue

         if (id.ne.0) bh(id)=bh(id)+thi
         if (jd.ne.0) bh(jd)=bh(jd)+thj    

         do 3 i=1,2
            the(i)=dexp(bh(i))
 3       continue




!---------- calcul de la vraisemblance ------------------

          res = 0.d0
          do 15 i=1,no

             vet = 0.d0
             if(nva.gt.0)then
                do 12 j=1,nva
                vet =vet + bh(np-nva+j)*dble(ve(i,j))
 12             continue
             endif   
             vet = dexp(vet)
             if(troncature.eq.1)then
                if(t0(i).eq.0.d0)then
                   tronc = 0.d0
                else   
                   call fonct(t0(i),the,ri,gl,su)
                   tronc = -(gl*vet)
                endif
             else
                tronc = 0.d0
             endif   
             res1 = 0.d0
             if(c(i).eq.0)then             
                call fonct(t1(i),the,ri,gl,su)
                res1 = -(gl*vet)
             else
                if(c(i).eq.1)then  
                   call fonct(t1(i),the,ri,gl,su)
                   res1 = -(gl*vet)+dlog(ri*vet)
                else
                   call fonct(t1(i),the,ri,gl,su1)
                   call fonct(t2(i),the,ri,gl,su2)
                   res1 = dlog((su1**vet)-(su2**vet)) 
                endif   
             endif   
             res = res + res1 + tronc
 15       continue   


          survLikelihood = res
          return
          end function survLikelihood
