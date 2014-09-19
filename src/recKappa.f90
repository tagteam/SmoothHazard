

        
        module recKappa
        
        CONTAINS
        
        subroutine sub_rec(b,omeg01,omeg12,omeg02)
        
        use commun,only:nz01,nz12,nz02,k0    
        use tailles
        
        implicit none
        
        integer::iii,ii,j
        double precision,dimension(3)::u,fu
        double precision,dimension(4)::b2
        double precision::uh
        double precision,dimension(nz01+2,nz01+2)::omeg01       
        double precision,dimension(nz12+2,nz12+2)::omeg12
        double precision,dimension(nz02+2,nz02+2)::omeg02
        double precision,dimension(np)::b
        double precision::estimv,mdf
        
        b2(1) = dlog(k0(1))
        b2(2) = dlog(k0(2))
        b2(3) = dlog(k0(3))
                
        do iii=1,2
                do ii=1,3
!                       write(*,*)' '
!                       write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!                       write(*,*)'Recherche parametre ',ii,' num ',iii
!                       write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!                       write(*,*)' '
                        uh = 1.d0
        
                        if(b2(ii).gt.2.d0)then
                                u(2) = b2(ii)
                                u(1) = u(2) - uh
                                u(3) = u(2) + uh
                        else
                                u(2) = 3.d0
                                u(1) = u(2) - uh
                                u(3) = u(2) + uh
                        endif   
!                  write(*,*)'ici'
                        b2(ii) = u(1)
                        
                        fu(1) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
        
                        b2(ii) = u(2)
                        fu(2) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
        
                        b2(ii) = u(3)
                        fu(3) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)

                        do j=1,5
                                
                                if((fu(1).lt.fu(2)).and.((fu(2).lt.fu(3))))then
!                                       write(*,*)'cas 1',u(3)
                                        u(1) = u(2)
                                        fu(1) = fu(2)
                                        u(2) = u(3)
                                        fu(2) = fu(3)
                                        u(3) = u(2) + uh
                                        b2(ii) = u(3)
                                        fu(3) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
                                else 
                                        if((fu(1).lt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                               write(*,*)'cas 2',u(1),u(3)
                                        else
                                                if((fu(1).gt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                       write(*,*)'cas 3',u(1)
                                                        u(3) = u(2)
                                                        fu(3) = fu(2)
                                                        u(2) = u(1)
                                                        fu(2) = fu(1)
                                                        u(1) = u(3) - 2.d0*uh
                                                        b2(ii) = u(1)
                                                        fu(1) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
                                                else
!                                                       write(*,*)'PB'
                                                endif
                                                
                                        endif
                                endif  
                        end do
        
!                       write(*,*)' '
!                       write(*,*)'Direction'
!                       write(*,*)' '

                        if((fu(1).lt.fu(2)).and.((fu(2).lt.fu(3))))then
!                               write(*,*)'cas 1, apres ',u(3)
                                b2(ii) = u(3)
                        else 
                                if((fu(1).lt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                       write(*,*)'cas 2, ZOOM entre',u(1),' et ',u(3)
!                                       write(*,*)' '
                                        uh = 0.5d0
                                        u(3) = u(2)
                                        fu(3) = fu(2)
                                        u(2) = u(1) + uh
                                        b2(ii) = u(2)
                                        fu(2) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)

                                        do j=1,2
                                                if((fu(1).lt.fu(2)).and.((fu(2).lt.fu(3))))then
!                                                       write(*,*)'cas 1bis',u(3)
                                                        u(1) = u(2)
                                                        fu(1) = fu(2)
                                                        u(2) = u(3)
                                                        fu(2) = fu(3)
                                                        u(3) = u(2) + uh
                                                        b2(ii) = u(3)
                                                        fu(3) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
                                                else 
                                                        if((fu(1).lt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                               write(*,*)'cas 2 bis',u(1),u(3)
                                                        else
                                                                if((fu(1).gt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                                       write(*,*)'cas 3 bis, ????????',u(1)
                                                                        u(3) = u(2)
                                                                        fu(3) = fu(2)
                                                                        u(2) = u(1)
                                                                        fu(2) = fu(1)
                                                                        u(1) = u(3) - 2.d0*uh
                                                                        b2(ii) = u(1)
                                                                        fu(1) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
                                                                else
!                                                                       write(*,*)'PB'
                                                                endif
                                                                
                                                        endif
                                                endif 
                                        end do
                                else
                                        if((fu(1).gt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                               write(*,*)'cas 3, avant ',u(1)
                                                b2(ii) = u(1)
                                        else
!                                               write(*,*)'PB'
                                        endif
                                endif
                        endif 
        
!                       write(*,*)' '
!                       write(*,*)'ZOOM ZOOM'
!                       write(*,*)' '
                        
                        if((fu(1).lt.fu(2)).and.((fu(2).lt.fu(3))))then
!                               write(*,*)'cas 1, apres ',u(3)
                                b2(ii) = u(3)
                        else 
                                if((fu(1).lt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                       write(*,*)'cas 2, ZOOM ZOOM entre',u(1),' et ',u(3)
!                                       write(*,*)' '
                                        uh = 0.25d0
                                        u(3) = u(2)
                                        fu(3) = fu(2)
                                        u(2) = u(1) + uh
                                        b2(ii) = u(2)
                                        fu(2) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
        
                                        do j=1,2
                                                if((fu(1).lt.fu(2)).and.((fu(2).lt.fu(3))))then
!                                                       write(*,*)'cas 1Ter',u(3)
                                                        u(1) = u(2)
                                                        fu(1) = fu(2)
                                                        u(2) = u(3)
                                                        fu(2) = fu(3)
                                                        u(3) = u(2) + uh
                                                        b2(ii) = u(3)
                                                        fu(3) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
                                                else 
                                                        if((fu(1).lt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                       write(*,*)'cas 2 Ter',u(2)
!                                                           write(*,*)u
!                                                           write(*,*)fu
                                                            b2(ii) = u(2)
                                                        else
                                                                if((fu(1).gt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                                       write(*,*)'cas 3 Ter, ????????',u(1)
                                                                        u(3) = u(2)
                                                                        fu(3) = fu(2)
                                                                        u(2) = u(1)
                                                                        fu(2) = fu(1)
                                                                        u(1) = u(3) - 2.d0*uh
                                                                        b2(ii) = u(1)
                                                                        fu(1) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
                                                                else
!                                                                       write(*,*)'PB'
                                                                endif
                                                        endif
                                                endif
                                        end do
                                else
                                        if((fu(1).gt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                               write(*,*)'cas 3, avant ',u(1)
                                                b2(ii) = u(1)
                                        else
!                                               write(*,*)'PB'
                                        endif
                                endif
                        endif
        
!                       write(*,*)' '
!                       write(*,*)'ZOOM ZOOM ZOOM'
!                       write(*,*)' '

                        if((fu(1).lt.fu(2)).and.((fu(2).lt.fu(3))))then
!                               write(*,*)'cas 1, apres ',u(3)
                                b2(ii) = u(3)
                        else 
                                if((fu(1).lt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                       write(*,*)'cas 2, ZOOM ZOOM ZOOM entre',u(1),' et ',u(3)
!                                       write(*,*)' '

                                        uh = 0.125d0
                                        u(3) = u(2)
                                        fu(3) = fu(2)
                                        u(2) = u(1) + uh
                                        b2(ii) = u(2)
                                        fu(2) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
                
                                        do j=1,2
                                                if((fu(1).lt.fu(2)).and.((fu(2).lt.fu(3))))then
!                                                       write(*,*)'cas 1 Qua',u(3)
                                                        u(1) = u(2)
                                                        fu(1) = fu(2)
                                                        u(2) = u(3)
                                                        fu(2) = fu(3)
                                                        u(3) = u(2) + uh
                                                        b2(ii) = u(3)
                                                        fu(3) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
                                                else 
                                                        if((fu(1).lt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                               write(*,*)'cas 2 Qua',u(2)
!                                                          write(*,*)u
!                                                          write(*,*)fu
                                                           b2(ii) = u(2)
                                                        else
                                                                if((fu(1).gt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                                       write(*,*)'cas 3 Qua, ????????',u(1)
                                                                        u(3) = u(2)
                                                                        fu(3) = fu(2)
                                                                        u(2) = u(1)
                                                                        fu(2) = fu(1)
                                                                        u(1) = u(3) - 2.d0*uh
                                                                        b2(ii) = u(1)
                                                                        fu(1) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
                                                                else
!                                                                       write(*,*)'PB'
                                                                endif
                                                        endif
                                                endif
                                        end do
                                else
                
                                        if((fu(1).gt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                       write(*,*)'cas 3, avant ',u(1)
                                                        b2(ii) = u(1)
                                        else
!                                                       write(*,*)'PB'
                                        endif
                                endif
                        endif 
        
!                       write(*,*)' '
!                       write(*,*)'ZOOM FINAL'
!                       write(*,*)' '
                        if((fu(1).lt.fu(2)).and.((fu(2).lt.fu(3))))then
!                               write(*,*)'cas 1, apres ',u(3)
                                b2(ii) = u(3)
                        else 
                                if((fu(1).lt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                       write(*,*)'cas 2, ZOOM FINAL entre',u(1),' et ',u(3)
!                                       write(*,*)' '

                                        uh = 0.125d0*0.5d0
                                        u(3) = u(2)
                                        fu(3) = fu(2)
                                        u(2) = u(1) + uh
                                        b2(ii) = u(2)
                                        fu(2) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
        
                                        do j=1,2
                                                if((fu(1).lt.fu(2)).and.((fu(2).lt.fu(3))))then
!                                                       write(*,*)'cas 1 FINAL',u(3)
                                                        u(1) = u(2)
                                                        fu(1) = fu(2)
                                                        u(2) = u(3)
                                                        fu(2) = fu(3)
                                                        u(3) = u(2) + uh
                                                        b2(ii) = u(3)
                                                        fu(3) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
                                                else 
                                                        if((fu(1).lt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                               write(*,*)'cas 2 FINAL',u(2)
                                                        else
                                                                if((fu(1).gt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                                       write(*,*)'cas 3 FINAL, ????????',u(1)
                                                                        u(3) = u(2)
                                                                        fu(3) = fu(2)
                                                                        u(2) = u(1)
                                                                        fu(2) = fu(1)
                                                                        u(1) = u(3) - 2.d0*uh
                                                                        b2(ii) = u(1)
                                                                        fu(1) = estimv(b2,b,omeg01,omeg12,omeg02,mdf)
                                                                else
!                                                                       write(*,*)'PB'
                                                                endif
                                                        endif
                                                endif
                                        end do
        
        
!                                       write(*,*)
!                                       write(*,*)'-------------------------------'
                                        if((fu(1).lt.fu(2)).and.((fu(2).lt.fu(3))))then
!                                               write(*,*)'Term 1',u(3)
!                                               write(*,*)u
!                                               write(*,*)fu
                                                b2(ii) = u(3)
                                        else 
                                                if((fu(1).lt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                       write(*,*)'Term 2',u(2)
!                                                       write(*,*)u
!                                                       write(*,*)fu
                                                        b2(ii) = u(2)
                                                else
                                                        if((fu(1).gt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                                               write(*,*)'Term 3',u(1)
!                                                               write(*,*)u
!                                                               write(*,*)fu
                                                                b2(ii) = u(1)
                                                        else
!                                                               write(*,*)'PB'
                                                        endif
                                                endif
                                        endif
                                else
                
!                                       write(*,*)'oo'
                                        if((fu(1).gt.fu(2)).and.((fu(2).gt.fu(3))))then
!                                               write(*,*)'cas 3, avant ',u(1)
                                                b2(ii) = u(1)
                                        else
!                                               write(*,*)'PB'
                                        endif
                                endif
                        endif 
                end do
        
        end do
        
        k0(1) = dexp(b2(1))
        
        k0(2) = dexp(b2(2))
        
        k0(3) = dexp(b2(3))
        
        end subroutine sub_rec
                
        end module recKappa
