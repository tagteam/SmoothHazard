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
        integer::i,j
        double precision::res,res1,gl,su,ri,vet,su1,su2,tronc


        bh=b
        the=0.d0
        gl=0.d0
        su=0.d0
        ri=0.d0
        su1=0.d0
        su2=0.d0

        if (id.ne.0) bh(id)=bh(id)+thi
        if (jd.ne.0) bh(jd)=bh(jd)+thj    

        do i=1,2
                the(i)=bh(i)*bh(i)
        end do

!---------- calcul de la vraisemblance ------------------

        res = 0.d0


        do i=1,no
                vet = 0.d0
                if(nva.gt.0)then
                        do j=1,nva
                                vet = vet + bh(np-nva+j)*dble(ve(i,j))
                        end do
                endif  
 
                vet = dexp(vet)

                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else   
                                call fonct(t0(i),the,ri,gl,su)
                                tronc = gl*vet
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
!               print*,i, ' res ',res,' res1',res1,' tronc',tronc
!               print*,' c(i) ',c(i),' t1(i)',t1(i)
!               print*,' gl ',gl, ' vet ',vet
                res = res + res1 + tronc
        

                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        survLikelihood=-1.d9
                        goto 123
                end if
        end do   

        survLikelihood = res
!       print*,' survLikelihood ',survLikelihood        
        return

123     continue

        end function survLikelihood
