!------------------------  FONCT  -------------------------
       subroutine fonct(x,p,risq,glam,surv)
       
            implicit none
            double precision,dimension(2)::p
            double precision::x,dens,surv,risq,glam

              surv = dexp(-(p(2)*x)**p(1))

              glam = (p(2)*x)**p(1)

              risq = p(1)*(p(2)**p(1))*(x**(p(1)-1.d0))

              if (x.le.0.d0) then
                 surv = 1.d0
                 glam = 0.d0
                 risq = 0.d0
              endif 
	        
              return
        end subroutine fonct
