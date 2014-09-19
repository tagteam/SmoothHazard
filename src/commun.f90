module commun
        integer,save::nn
!-------------------- pen1 - pen6 -----------------------------------
         double precision,dimension(:),allocatable,save::m3m3a,m2m2a,m1m1a,mmma, & 
         m3m1a,m3ma,m2m1a,m2ma,m3m3b,m2m2b,m1m1b,mmmb, &
         m3m1b,m3mb,m2m1b,m2mb,m3m2a,m1ma,m3m2b,m1mb, &
         m3m3c,m2m2c,m1m1c,mmmc,m3m1c,m3mc,m2m1c,m2mc, &
         m3m2c,m1mc
!-------------------mem1 & mem2----------------------------------------
         double precision,dimension(:),allocatable,save::mm3a,mm2a,mm1a,mma, &
         im3a,im2a,im1a,ima
!-------------------  mem3 & mem4 ----------------------------------
         double precision,dimension(:),allocatable,save::mm3b,mm2b,mm1b,mmb, &
         im3b,im2b,im1b,imb
!-----------------  mem5 & mem6-----------------------------------------
         double precision,dimension(:),allocatable,save::mm3c,mm2c,mm1c,mmc, &
         im3c,im2c,im1c,imc
!-------------------  dace1 -------------------------------------------
        double precision,dimension(:),allocatable,save::zi01,zi12,zi02
!-------------------  dace2 -------------------------------------------
        integer::no
!-------------------  dace3 -------------------------------------------
        integer::nz01,nz12,nz02
!-------------------  ve1 -------------------------------------------
        double precision,dimension(:,:),allocatable,save::ve01,ve12,ve02

!-------------------  ve2 -------------------------------------------
        integer::nva01,nva12,nva02
!-------------------  pe -------------------------------------------
        double precision,save::pe
!-------------------  dace1new -------------------------------------------
        double precision,dimension(:),allocatable,save::t0,t1,t2,t3,t4
        integer,dimension(:),allocatable,save::c
!-------------------  dace1new -------------------------------------------
        double precision,dimension(:),allocatable,save::opt2
        double precision,dimension(3),save::k0
        integer,save::rec,troncature,ind_hess,pl,iconf
        double precision,dimension(:,:),allocatable,save::hessienne     

!add for survPl
!-------------------- pen -----------------------------------
        double precision,dimension(:),allocatable,save::m3m3,m2m2,m1m1,mmm, &
        m3m2,m3m1,m3m,m2m1,m2m,m1m
!-------------------  dace1 -------------------------------------------
        double precision,dimension(:),allocatable,save::zi
!-------------------  dace -------------------------------------------
        integer::nz,verSurv
!-------------------  ve1 -------------------------------------------
        double precision,dimension(:,:),allocatable,save::ve
        double precision::k0surv

end module commun

module tailles
        integer,save::np
end module tailles

