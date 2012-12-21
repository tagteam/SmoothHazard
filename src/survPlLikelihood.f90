

!========================          survPlLikelihood         ====================
      double precision function survPlLikelihood(b,np,id,thi,jd,thj)

        use commun,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m, &
             zi,c,no,nz,ve,verSurv,pe,t0,t1,t2,k0surv,troncature

        implicit none
         
        double precision::thi,thj,res,res1,tronc,vet
        integer::np,id,jd,i,j,nva
        double precision,dimension(np)::b,bh
        double precision,dimension(-2:(nz))::the
        double precision::su,ri,su1,su2

        nva = verSurv 
        bh=b


         if (id.ne.0) bh(id)=bh(id)+thi
         if (jd.ne.0) bh(jd)=bh(jd)+thj    

         do i=1,nz+2
            the(i-3)=(bh(i))*(bh(i))
!       the01(i-3)=dexp(bh(i))
         end do

!---------- calcul de la vraisemblance ------------------

        res = 0.d0
        do i=1,no

           vet = 0.d0

           if(nva.gt.0)then
              do j=1,nva
                 vet =vet + bh(np-nva+j)*dble(ve(i,j))
              end do
           endif

           vet = dexp(vet)


                
             if(troncature.eq.1)then
                if(t0(i).eq.0.d0)then
                   tronc = 0.d0
                else   
                   call susp(t0(i),the,nz,su,ri,zi) 
                   tronc = - (dlog(su))*vet
                endif
             else
                tronc = 0.d0
             endif   
             res1 = 0.d0
             if(c(i).eq.0)then             
                call susp(t1(i),the,nz,su,ri,zi) 
                res1 = (dlog(su))*vet
             else
                if(c(i).eq.1)then  
                   call susp(t1(i),the,nz,su,ri,zi) 
                   res1 = (dlog(su)*vet)+dlog(ri*vet)
                else
                   call susp(t1(i),the,nz,su1,ri,zi) 
                   call susp(t2(i),the,nz,su2,ri,zi) 
                   res1 = dlog((su1**vet)-(su2**vet)) 
                endif   
             endif   
             res = res + res1 + tronc
             
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        survPlLikelihood=-1.d9
                        goto 123
                end if
                
          end do   

!---------- calcul de la penalisation -------------------

        pe = 0.d0

        do i=1,nz-1
                pe = pe+(the(i-3)*the(i-3)*m3m3(i))+(the(i-2) &
                *the(i-2)*m2m2(i))+(the(i-1)*the(i-1)*m1m1(i))+( &
                the(i)*the(i)*mmm(i))+(2.d0*the(i-3)*the(i-2)* &
                m3m2(i))+(2.d0*the(i-3)*the(i-1)*m3m1(i))+(2.d0* &
                the(i-3)*the(i)*m3m(i))+(2.d0*the(i-2)*the(i-1)* &
                m2m1(i))+(2.d0*the(i-2)*the(i)*m2m(i))+(2.d0*the(i-1) &
                *the(i)*m1m(i))
        end do
         

        pe = k0surv*pe
!!                write(*,*)'res pe',res,pe
        res = res - pe

        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                survPlLikelihood=-1.d9
                goto 123
        end if

        survPlLikelihood = res

123     continue

        return

        end function survPlLikelihood
