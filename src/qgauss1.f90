!================================  QGAUS1   ==========================
      subroutine qgauss1(a,b,the01,the02,the12,res,v01,v02,v12)
        implicit none
         double precision a,b,the01(2),the02(2),the12(2)
         double precision dx,xm,xr,w(5),x(5),res,v01,v02,v12
         double precision xx,f1,su01,ri01,ri12,f2,su12,su02,ri02
         double precision gl01,gl12,gl02
         integer j,cas
         save w,x
         data w/0.2955242247d0,0.2692667193d0,0.2190863625d0,0.1494513491d0,0.0666713443d0/
         data x/0.1488743389d0,0.4333953941d0,0.6794095682d0,0.8650633666d0,0.9739065285d0/

         
            xm = 0.5d0*(b+a)
            xr = 0.5d0*(b-a)
            res = 0.d0
            if(a.eq.b)then
               res = 0.d0
            else
               do 11 j=1,5
                  dx=xr*x(j)
                  xx = xm+dx
                  call fonct(xx,the01,ri01,gl01,su01)
                  call fonct(xx,the02,ri02,gl02,su02)
                  call fonct(xx,the12,ri12,gl12,su12)
                  f1 = (su01**v01)*(su02**v02)*ri01*v01/(su12**v12)
                  xx = xm-dx
                  call fonct(xx,the01,ri01,gl01,su01)
                  call fonct(xx,the02,ri02,gl02,su02)
                  call fonct(xx,the12,ri12,gl12,su12)
                  f2 = ((su01**v01)*(su02**v02)*ri01*v01)/(su12**v12)
                  res = res + w(j)*(f1+f2)
 11            continue
            endif
            res = res*xr

            
          end subroutine qgauss1
