          
!========================          idmPlLikelihood0         ====================
        double precision function idmPlLikelihood0(b,np)
        use commun,only:zi01,zi12,zi02, &
        c,no,nz01,nz12,nz02,ve01,ve12,ve02,nva01,nva12,nva02,t0,t1,t2,t3,troncature

        implicit none
        
        integer::np,i,j
        double precision,dimension(np)::b,bh
        double precision::res,res1,res2,tronc,vet01,vet12,vet02
        double precision,dimension(-2:(nz01-1))::the01
        double precision,dimension(-2:(nz12-1))::the12
        double precision,dimension(-2:(nz02-1))::the02
        double precision::su01,ri01,su12,ri12,su02,ri02
         



        do i=1,np
                bh(i)=b(i)
        end do

        do i=1,nz01+2
                the01(i-3)=bh(i)*bh(i)
!               the01(i-3)=dexp(bh(i))
        end do
        
        do i=1,nz02+2
                j = nz01+2+i
                the02(i-3)=bh(j)*bh(j)
!               the12(i-3)=dexp(bh(j))
        end do
        
        do i=1,nz12+2
                j = nz02+2+nz01+2+i
                the12(i-3)=bh(j)*bh(j)
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
                        endif
                else
                        tronc = 0.d0
                endif

                if(c(i).eq.1)then ! cad 0-->1 et 0-->2
                        call susp(t1(i),the01,nz01,su01,ri01,zi01)
                        call susp(t3(i),the02,nz02,su02,ri02,zi02)
                        res1 = (dlog(su01)*vet01)+(dlog(su02)*vet02)

                else
                if(c(i).eq.2)then ! cpi 0-->1
                        call qgaussPL(c(i),t1(i),t2(i),the01,the12,the02,res2, &
                        vet01,vet12,vet02)
                        call susp(t3(i),the12,nz12,su12,ri12,zi12)      
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
                                        res1=dlog(res2*(su12**vet12)*ri12*vet12)
                                else
                                        if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
                                                call susp(t2(i),the01,nz01,su01,ri01,zi01)
                                                call susp(t3(i),the12,nz12,su12,ri12,zi12)
                                                call susp(t1(i),the02,nz02,su02,ri02,zi02)
                                                res2=dlog(su01)*vet01 + dlog(ri01*vet01) + dlog(su12)*vet12 & 
                                                + dlog(ri12*vet12) + dlog(su02)*vet02
                                                call susp(t2(i),the12,nz12,su12,ri12,zi12)
                                                res1=res2-dlog(su12)*vet12
                                        else 
                                                if(c(i).eq.6)then ! vivant 
                                                        call qgaussPL(c(i),t1(i),t2(i),the01,the12,the02,&
                                                        res2,vet01,vet12,vet02)
                                                        call susp(t3(i),the12,nz12,su12,ri12,zi12)
                                                        call susp(t3(i),the01,nz01,su01,ri01,zi01)
                                                        call susp(t3(i),the02,nz02,su02,ri02,zi02)
                                                        res1=dlog(res2*(su12**vet12) + (su01**vet01)*(su02**vet02))  
                                                        else  ! passage 0-->2   
                                                                call qgaussPL(c(i),t1(i),t3(i),the01,the12,the02,&
                                                                res2,vet01,vet12,vet02)
                                                                call susp(t3(i),the12,nz12,su12,ri12,zi12)
                                                                call susp(t3(i),the02,nz02,su02,ri02,zi02)
                                                                call susp(t3(i),the01,nz01,su01,ri01,zi01)
                                                                res1=dlog(res2*ri12*vet12*(su12**vet12) + &
                                                                ri02*vet02*(su02**vet02)*(su01**vet01))

                                                        endif
                                                endif  
                                        endif    
                                endif
                        endif   
                endif   

                res = res + res1  + tronc
        end do   


        idmPlLikelihood0 = res

        return

        end function idmPlLikelihood0

