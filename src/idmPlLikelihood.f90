

!========================          idmPlLikelihood         ====================
      double precision function idmPlLikelihood(b,np,id,thi,jd,thj)

        use commun,only:m3m3a,m2m2a,m1m1a,mmma,m3m1a,m3ma,m2m1a,m2ma,m3m3b,m2m2b,m1m1b,mmmb, &
        m3m1b,m3mb,m2m1b,m2mb,m3m2a,m1ma,m3m2b,m1mb,m3m3c,m2m2c,m1m1c,mmmc,m3m1c,m3mc,m2m1c,m2mc, &
        m3m2c,m1mc,&
        zi01,zi12,zi02,c,no,nz01,nz12,nz02,ve01,ve12,ve02,&
        nva01,nva12,nva02,pe,t0,t1,t2,t3,k0,troncature

        implicit none
         
        double precision::thi,thj,pe01,pe12,pe02,res,res1,res2,tronc, &
        vet01,vet12,vet02
        integer::np,id,jd,i,j
        double precision,dimension(np)::b,bh
        double precision,dimension(-2:(nz01-1))::the01
        double precision,dimension(-2:(nz12-1))::the12
        double precision,dimension(-2:(nz02-1))::the02
        double precision::su01,ri01,su12,ri12,su02,ri02

         


        bh=b


         if (id.ne.0) bh(id)=bh(id)+thi
         if (jd.ne.0) bh(jd)=bh(jd)+thj    

         do i=1,nz01+2
            the01(i-3)=(bh(i))*(bh(i))
!       the01(i-3)=dexp(bh(i))
         end do
         do i=1,nz02+2
            j = nz01+2+i
            the02(i-3)=(bh(j))*(bh(j))
!       the12(i-3)=dexp(bh(j))
         end do
         do i=1,nz12+2
            j = nz02+2+nz01+2+i
            the12(i-3)=(bh(j))*(bh(j))
!       the02(i-3)=dexp(bh(j))
         end do


!---------- calcul de la vraisemblance ------------------

        res = 0.d0
        do i=1,no
                vet01 = 0.d0
                vet12 = 0.d0
                vet02 = 0.d0

                if(nva01.gt.0)then
                        do j=1,nva01
                                vet01 =vet01 + bh(np-nva01-nva12-nva02+j)*dble(ve01(i,j))
                        end do
                endif  
 
                if(nva02.gt.0)then
                        do j=1,nva02
                                vet02 =vet02 + bh(np-nva02-nva12+j)*dble(ve02(i,j))
                        end do
                endif

                if(nva12.gt.0)then
                        do j=1,nva12
                                vet12 =vet12 + bh(np-nva12+j)*dble(ve12(i,j))
                        end do
                endif

                vet01 = dexp(vet01)
                vet12 = dexp(vet12)
                vet02 = dexp(vet02)

                res1 = 0.d0
                
                if(troncature.eq.1)then
                        if(t0(i).eq.0.d0)then
                                tronc = 0.d0
                        else 
                                call susp(t0(i),the01,nz01,su01,ri01,zi01) 
                                call susp(t0(i),the02,nz02,su02,ri02,zi02) 
                                tronc=(-dlog(su01)*vet01)-(dlog(su02)*vet02)
                        end if
                else
                        tronc = 0.d0
                end if

                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
                        call susp(t1(i),the01,nz01,su01,ri01,zi01)
                        call susp(t3(i),the02,nz02,su02,ri02,zi02)
                        res1 = (dlog(su01)*vet01)+(dlog(su02)*vet02)
                else
                if(c(i).eq.2)then ! cpi 0-->1
                        call qgaussPL(c(i),t1(i),t2(i),the01,the12,the02,res2, &
                            vet01,vet12,vet02)
                        call susp(t3(i),the12,nz12,su12,ri12,zi12)
!               res1=dlog(res2)+dlog(su12)*vet12 (autre ecriture)
                        res1=dlog(res2*(su12**vet12))
                else  
                    if(c(i).eq.3)then ! obs 0-->1
                    call susp(t2(i),the01,nz01,su01,ri01,zi01)
                    call susp(t3(i),the12,nz12,su12,ri12,zi12)
                    call susp(t1(i),the02,nz02,su02,ri02,zi02)
                    res2=dlog(su01)*vet01+dlog(ri01*vet01) + dlog(su12)*vet12+dlog(su02)*vet02
                    call susp(t2(i),the12,nz12,su12,ri12,zi12)
                    res1=res2 - dlog(su12)*vet12
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2

                       call susp(t3(i),the12,nz12,su12,ri12,zi12)
                       call qgaussPL(c(i),t1(i),t2(i),the01,the12,the02,res2,&
                       vet01,vet12,vet02)
!                      res1=dlog(res2*ri12*vet12)+dlog(su12)*vet12 (autre ecriture)
                       res1=dlog(res2*(su12**vet12)*ri12*vet12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
                         call susp(t2(i),the01,nz01,su01,ri01,zi01)
                         call susp(t3(i),the12,nz12,su12,ri12,zi12)
                         call susp(t1(i),the02,nz02,su02,ri02,zi02)
                         res2=dlog(su01)*vet01 + dlog(ri01*vet01) + dlog(su12)*vet12 + dlog(ri12*vet12) + dlog(su02)*vet02
                         call susp(t2(i),the12,nz12,su12,ri12,zi12)
                         res1=res2-dlog(su12)*vet12
                         else
                            if(c(i).eq.6)then ! vivant ???
                            call qgaussPL(c(i),t1(i),t2(i),the01,the12,the02,&
                            res2,vet01,vet12,vet02)
                            call susp(t3(i),the12,nz12,su12,ri12,zi12)
                            call susp(t3(i),the01,nz01,su01,ri01,zi01)
                            call susp(t3(i),the02,nz02,su02,ri02,zi02)
                            res1=dlog(res2*(su12**vet12) + (su01**vet01)*(su02**vet02))  
                            else ! passage 0-->2  
                                call susp(t3(i),the12,nz12,su12,ri12,zi12)
                                call susp(t3(i),the02,nz02,su02,ri02,zi02)
                                call susp(t3(i),the01,nz01,su01,ri01,zi01)
                                call qgaussPL(c(i),t1(i),t3(i),the01,the12,the02,res2,vet01,vet12,vet02)
                                res1=dlog(res2*ri12*vet12*(su12**vet12) + ri02*vet02*(su02**vet02)*(su01**vet01))
                            endif
                          endif  
                        endif    
                    endif
                endif   
             endif   

                res = res + res1  + tronc
             
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        idmPlLikelihood=-1.d9
                        goto 123
                end if
                
          end do   

!---------- calcul de la penalisation -------------------

        pe01 = 0.d0
        pe12 = 0.d0
        pe02 = 0.d0
        do i=1,nz01-1
                pe01 = pe01+(the01(i-3)*the01(i-3)*m3m3a(i))+(the01(i-2) &
                *the01(i-2)*m2m2a(i))+(the01(i-1)*the01(i-1)*m1m1a(i))+( &
                the01(i)*the01(i)*mmma(i))+(2.d0*the01(i-3)*the01(i-2)* &
                m3m2a(i))+(2.d0*the01(i-3)*the01(i-1)*m3m1a(i))+(2.d0* &
                the01(i-3)*the01(i)*m3ma(i))+(2.d0*the01(i-2)*the01(i-1)* &
                m2m1a(i))+(2.d0*the01(i-2)*the01(i)*m2ma(i))+(2.d0*the01(i-1) &
                *the01(i)*m1ma(i))
        end do
         
        do i=1,nz12-1
                pe12 = pe12+(the12(i-3)*the12(i-3)*m3m3b(i))+(the12(i-2) &
                *the12(i-2)*m2m2b(i))+(the12(i-1)*the12(i-1)*m1m1b(i))+( &
                the12(i)*the12(i)*mmmb(i))+(2.d0*the12(i-3)*the12(i-2)* &
                m3m2b(i))+(2.d0*the12(i-3)*the12(i-1)*m3m1b(i))+(2.d0* &
                the12(i-3)*the12(i)*m3mb(i))+(2.d0*the12(i-2)*the12(i-1)* &
                m2m1b(i))+(2.d0*the12(i-2)*the12(i)*m2mb(i))+(2.d0*the12(i-1) &
                *the12(i)*m1mb(i))
        end do
         
        do i=1,nz02-1
                pe02 = pe02+(the02(i-3)*the02(i-3)*m3m3c(i))+(the02(i-2) &
                *the02(i-2)*m2m2c(i))+(the02(i-1)*the02(i-1)*m1m1c(i))+( &
                the02(i)*the02(i)*mmmc(i))+(2.d0*the02(i-3)*the02(i-2)* &
                m3m2c(i))+(2.d0*the02(i-3)*the02(i-1)*m3m1c(i))+(2.d0* &
                the02(i-3)*the02(i)*m3mc(i))+(2.d0*the02(i-2)*the02(i-1)* &
                m2m1c(i))+(2.d0*the02(i-2)*the02(i)*m2mc(i))+(2.d0*the02(i-1) &
                *the02(i)*m1mc(i))
        end do

        pe = k0(1)*pe01 + k0(2)*pe02+ k0(3)*pe12 
!                write(*,*)'res pe',res,pe

        res = res - pe

        if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                idmPlLikelihood=-1.d9
                goto 123
        end if
        
        idmPlLikelihood = res

123     continue

        return

        end function idmPlLikelihood
