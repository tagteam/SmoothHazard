
!========================          idmPlLikelihood2         ====================
      double precision function idmPlLikelihood2(b,np,id,thi,jd,thj)
        use commun,only:mm3a,mm2a,mm1a,mma,im3a,im2a,im1a,ima,mm3b,mm2b,mm1b,mmb, &
        im3b,im2b,im1b,imb,mm3c,mm2c,mm1c,mmc,im3c,im2c,im1c,imc,zi01,zi12,zi02, &
        c,no,nz01,nz12,nz02,t0,t1,t2,t3,t4
         implicit none
         
         integer::np,id,jd,i,j
         double precision::thi,thj,res,res1,res2
         double precision,dimension(np)::b,bh
         double precision,dimension(-2:(nz01-1))::the01
         double precision,dimension(-2:(nz12-1))::the12
         double precision,dimension(-2:(nz02-1))::the02
         double precision::su01,ri01,su12,ri12,su02,ri02



            bh=b


         if (id.ne.0) bh(id)=bh(id)+thi
         if (jd.ne.0) bh(jd)=bh(jd)+thj    

         do i=1,nz01+2
            the01(i-3)=bh(i)
         end do
         do i=1,nz02+2
            j = nz01+2+i
            the02(i-3)=bh(j)
         end do
         do i=1,nz12+2
            j = nz02+2+nz01+2+i
            the12(i-3)=bh(j)
         end do


!---------- calcul de la vraisemblance ------------------

          res = 0.d0
          do i=1,no
             if(c(i).eq.1)then ! cad 0-->1 et 0-->2
                res1 = 0.d0
             else
                if(c(i).eq.2)then ! cpi 0-->1
                call qgauss12(c(i),t1(i),t2(i),the01,the12,the02,res1)
                res1 = dlog(res1)
                 else  
                    if(c(i).eq.3)then ! obs 0-->1
                       res1 =0.d0
                    else   
                       if(c(i).eq.4)then ! cpi 0-->1 et obs 1-->2

                       call  qgauss12(c(i),t1(i),t2(i),the01,the12,the02,res1)
                       call susp(t3(i),the12,nz12,su12,ri12,zi12)
                       res1=dlog(res1*ri12)
                       else
                         if(c(i).eq.5)then ! obs 0-->1 et obs 1-->2
                       res1 = 0.d0
                         else
                            if(c(i).eq.6)then ! vivant 

                            call  qgauss12(c(i),t1(i),t2(i),the01,the12,the02,res2)
                            call susp(t3(i),the12,nz12,su12,ri12,zi12)
                            call susp(t3(i),the01,nz01,su01,ri01,zi01)
                            call susp(t3(i),the02,nz02,su02,ri02,zi02)
                            res1 = dlog(res2*su12 + su01*su02)
                            else   

                            call  qgauss12(c(i),t1(i),t2(i),the01,the12,the02,res2)
                            call susp(t3(i),the12,nz12,su12,ri12,zi12)
                            call susp(t3(i),the02,nz02,su02,ri02,zi02)
                            call susp(t3(i),the01,nz01,su01,ri01,zi01)
                            res1=dlog(res2*ri12*su12 + ri02*su02*su01)

                            endif
                          endif  
                        endif    
                    endif
                endif   
             endif   

             res = res + res1 
         end do   

          idmPlLikelihood2 = res

          return

          end function idmPlLikelihood2

