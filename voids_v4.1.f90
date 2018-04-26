! Author: Elena Ricciardelli
! Last modified: 12/05/16
! Version history:
! v1.1 : version used for results in submitted paper
! v1.2: added profiles for the barionyc matter
! v1.3: added profiles for the stellar component - write gas density maps
! v1.5? temporary ? 
! v2.0: added void merging history 
! 09/13/13 MODIFIED MALLA
! 09/18/13:  leer adapted to the new format of clus
!            diver_fina_gas: nl temporary set to 1; error for nl>1 occurs (line 2757)
! 09/24/2013: Added routine Haloes, finding dm haloes in voids and void walls
!            and commented Eulerian
! V4: - ADDED ROUTINE LEER_NBODY --> it also works with particle distribution
!  TODO:
!  - add flag: MASCLET/NBODY
!  - take into account periodic boundaries
!  - fix percolation in a more robust way ? whatershed ?
!  - void hierarchy form part distrib ???
MODULE COMMONDATA
      IMPLICIT NONE
      INCLUDE 'voids_parameters.dat'
!variables defined in the main unit
       REAL*4 ACHE, T0, RE0
       REAL*4 PI4ROD, REI, CGR, PI
!variables defined (and some of them allocated) in LEER
!DEFINIRE U1G(NHYX,NHYY,NHYZ,ILEV) CON ILEV=-2,0 
       REAL*4 U1G(NHYX,NHYY,NHYZ) !gas density contrast at ir=0
       REAL*4 U2G(NHYX,NHYY,NHYZ) !gas velocity in x
       REAL*4 U3G(NHYX,NHYY,NHYZ) !gas velocity in x
       REAL*4 U4G(NHYX,NHYY,NHYZ) !gas velocity in x
       REAL*4 U1S(NHYX,NHYY,NHYZ) !stellar density contrast at ir=0
       REAL*4, ALLOCATABLE:: U11G(:,:,:,:)  !gas density contrast at ir>0
       REAL*4, ALLOCATABLE:: U11S(:,:,:,:)  !stellar density contrast at ir>0
       REAL*4 U1DM(NHYX,NHYY,NHYZ) !DM density contrast at ir=0
       REAL*4 UU2DM(NHYX,NHYY,NHYZ) !DM X-VELOCITY FIELD 
       REAL*4 UU3DM(NHYX,NHYY,NHYZ) !DM Y-VELOCITY FIELD 
       REAL*4 UU4DM(NHYX,NHYY,NHYZ) !DM Z-VELOCITY FIELD 
       REAL*4, ALLOCATABLE:: U11DM(:,:,:,:)  !DM density contrast at ir>0
       REAL*4 U1(NHYX,NHYY,NHYZ) !total density contrast at ir=0
       REAL*4, ALLOCATABLE:: U11(:,:,:,:)  !total density contrast at ir>0
       REAL*4, ALLOCATABLE:: U12G(:,:,:,:) ! gas  x-velocity (eulerian) for refined levels
       REAL*4, ALLOCATABLE:: U13G(:,:,:,:) ! gas  y-velocity (eulerian) for refined levels
       REAL*4, ALLOCATABLE:: U14G(:,:,:,:) ! gas  z-velocity (eulerian) for refined levels
       REAL*4 U2DM(PARTIRED)  !DM particle velocity in x 
       REAL*4 U3DM(PARTIRED)  !DM particle velocity in y  
       REAL*4 U4DM(PARTIRED)  !DM particle velocity in z
       REAL*4 MASAP(PARTIRED) !DM particle mass      
       REAL*4 RXPA(PARTIRED)  !DM particle position in x
       REAL*4 RYPA(PARTIRED)  !DM particle position in x
       REAL*4 RZPA(PARTIRED)  !DM particle position in x
       INTEGER ORIPA(PARTIRED)
       INTEGER NPARTT
!defined/allocated in smooth
       REAL*4, ALLOCATABLE:: U1CO(:,:,:) !total density contrast in the coarse grid
       REAL*4, ALLOCATABLE:: U1GCO(:,:,:) !gas density contrast in the coarse grid
       REAL*4, ALLOCATABLE:: U1SCO(:,:,:) !gas density contrast in the coarse grid
       REAL*4, ALLOCATABLE:: U2CO(:,:,:) !DM x-velocity (eulerian) for the coarse cells 
       REAL*4, ALLOCATABLE:: U3CO(:,:,:) !DM y-velocity (eulerian)  for th871e coarse cells 
       REAL*4, ALLOCATABLE:: U4CO(:,:,:) !DM z-velocity (eulerian)  for the coarse cells 
       REAL*4, ALLOCATABLE:: WTCO(:,:,:) 
       REAL*4, ALLOCATABLE:: U2GCO(:,:,:) !DM x-velocity (eulerian) for the coarse cells 
       REAL*4, ALLOCATABLE:: U3GCO(:,:,:) !DM y-velocity (eulerian)  for the coarse cells 
       REAL*4, ALLOCATABLE:: U4GCO(:,:,:) !DM z-velocity (eulerian)  for the coarse cells 
!defined/allocated in eulerian
       REAL*4, ALLOCATABLE:: U2(:,:,:) !DM x-velocity  (eulerian) for level 0
       REAL*4, ALLOCATABLE:: U3(:,:,:) !DM y-velocity  (eulerian) 
       REAL*4, ALLOCATABLE:: U4(:,:,:) !DM z-velocity  (eulerian) 
       REAL*4, ALLOCATABLE:: WT0(:,:,:)
       REAL*4, ALLOCATABLE:: U12(:,:,:,:) !DM x-velocity  (eulerian) for refined levels
       REAL*4, ALLOCATABLE:: U13(:,:,:,:) !DM y-velocity  (eulerian) 
       REAL*4, ALLOCATABLE:: U14(:,:,:,:) !DM z-velocity  (eulerian) 
!variables defined in MALLA
       REAL*4 DX,DY,DZ
       REAL*4, ALLOCATABLE:: RADX(:), RADY(:),RADZ(:)
       !REAL*4  RADX(0:NCOX+1), RADY(0:NCOY+1),RADZ(0:NCOZ+1)
       REAL*4 DX0,DY0,DZ0
       REAL*4  RADX0(0:NHYX+1), RADY0(0:NHYY+1),RADZ0(0:NHYZ+1)
       REAL*4 DXRR,DYRR,DZRR
       REAL*4, ALLOCATABLE, DIMENSION(:):: RADXR, RADYR, RADZR
!variables defined in MARK
       REAL*4, ALLOCATABLE:: DIVERCO(:,:,:)
       REAL*4, ALLOCATABLE:: DIVERGCO(:,:,:)
!variables defined in MARK and then MARK_SUB
       INTEGER, ALLOCATABLE:: FLAGV(:,:,:)
       INTEGER, ALLOCATABLE:: FLAG_SUB(:,:,:)
       INTEGER, ALLOCATABLE:: MARCAP(:,:,:), MARCACO(:,:,:)
       INTEGER, ALLOCATABLE:: MARCA_OLD(:,:,:) !V2.0
!variables defined in void_find: centers(icx, icy, icz) and extensions of voids
       INTEGER, DIMENSION(NVOID_MAX) :: INICIOX, FINALX, INICIOY, FINALY, &
            INICIOZ, FINALZ, ICX, ICY, ICZ
       REAL*4, DIMENSION(NVOID_MAX):: VOL
       REAL*4, DIMENSION(NVOID_MAX):: RINIXCO, RFINXCO, RINIX0, RFINX0, RINIXR, RFINXR, &
            RINIYCO, RFINYCO, RINIY0, RFINY0, RINIYR, RFINYR, &
            RINIZCO, RFINZCO, RINIZ0, RFINZ0, RINIZR, RFINZR

       INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: MARCA

!variables defined in DIVER
       REAL*4, ALLOCATABLE:: DIVER0(:,:,:)
       REAL*4, ALLOCATABLE:: DIVER(:,:,:,:)
!variables defined in VMESH
       INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: FLAGAMR
!variables defined and allocated in profiles
       INTEGER, ALLOCATABLE:: NPV(:) !# part in each void
       INTEGER, ALLOCATABLE:: VIP(:) !# particle id in the void list
!variables defined and allocated in vmesh
       REAL*4, ALLOCATABLE:: UR(:,:,:),UGR(:,:,:), USR(:,:,:), DIVR(:,:,:)
       INTEGER, ALLOCATABLE:: MARCAR(:,:,:), FLAGR(:,:,:)
       logical :: flag_map
!PARAMETERS: --> put in voids_paramters.dat
       !REAL*4, PARAMETER:: RMIN=3. !minimum radius for void analysis --> put in voids.dat
       REAL*4, PARAMETER:: FR=0.5 ! --> put in voids.dat
       INTEGER, PARAMETER:: NGAL_MAX=2000.
       INTEGER, PARAMETER:: NCELL_GAL=10 !number of cells used for wall definition: should be always grater than Rmax*DR/DX (Rmax: largest void)  
       INTEGER, PARAMETER:: FLAG_GAS=1 !1=use diverGAS; 0=use diverDM
       INTEGER, PARAMETER:: FLAG_DENS=1 !0=TOTAL, 1=DM, 2=GAS
       INTEGER:: FLAG_GAL=1 !1= when ddens is above the threshold check the gradient also in the next cell 
       INTEGER, PARAMETER:: FLAG_REF=1 !1=use refined levels; 0=use only coarse (128^3) grid
       !REAL*4, PARAMETER:: FLAG_DIV=1 !1=use condition on diverV to define the edges 
       INTEGER, PARAMETER:: FLAG_DIV1=1 !1=use condition on diverV to define the edges of voids
       INTEGER, PARAMETER:: FLAG_DIV2=1 !1=use condition on diverV to define the edges of subvoids
       INTEGER, PARAMETER:: FLAG_PROF=0 !1=compute profiles, 0 speed up the code in the case  voids are searched just in levels coarser than levp
       INTEGER, PARAMETER:: FLAG_VEL=1 !0=do not use velocity field, 1=use velocity field
       !LOGICAL, PARAMETER :: FLAG_MAP = .FALSE.
       LOGICAL, PARAMETER :: FLAG_DMO = .TRUE.
END MODULE COMMONDATA

!!*********************************************************************** 
!*********************************************************************** 
       PROGRAM VOID_FINDER                                              
!*********************************************************************** 

!usare densita' tot=dm+gas
!smussare tutto su scala di 1 Mpc --> costruire nuova malla=128^3
!passare a vel dm euleriane --> TSC
!calc divergenza --> DIVER_FINA? DIVER_COARSE?
!usare analogo della routine shocked per marcare le celle che possono stare nei vuoti
! - divV > 0
! - delta< thresh
!estendere lo pseudovuoto fino a che sono soddisfatte queste condizioni 
! - divV>0
! - grad(delta)<thresh_grad
! - delta< thresh2 ?
!poi eventualemte refinare a livello 1-2 per avere i bordi piu' precisi
!cercare sottostrutture usando i livelli + alti

       USE COMMONDATA  
       IMPLICIT NONE                                                     

       INTEGER FIRST,LAST, EVERY, IFI, NFILE, NDXYZ0, LEVV, LEVS, LEVP, LEVMIN, LEVMAX
       INTEGER NHYX2,NHYY2,NHYZ2, NX,NY,NZ,ITER, CONTA, NVOIDP, NOVER 
       REAL*4 OMEGA0, ZI, LADO0, LADO, RMIN
       REAL*4 MINF, MAXF, MINF0, MAXF0
       REAL*4 UNTERCIO,CGR2,RODO,ROI, RETE, ROTE, CONST
       REAL*4 T,TEI, ROCRIT, UM, UV, ZETA
       REAL*4 DENS_THRE, DENS_THRE2, GRAD_THRE
       REAL*4 TIEMPOI
       REAL*4 VOLT, VOLT2, VOLT3, VOLT4, VOLM, VCELL	
       REAL*4 VELL, INVPOR, VMIN0, RMIN0
       INTEGER I,J,K,IX,JY,KZ, NCELL
       INTEGER NHYXP, NHYYP, NHYZP, NMAXP, NMAYP, NMAZP
       INTEGER NCOXP, NCOYP, NCOZP, DXP, DYP, DZP
       INTEGER NPALEVP, PARTIREDP
       INTEGER:: NXX, NYY, NZZ, NX2, NY2, NZ2, NXR, NYR, NZR
       REAL*4:: DXX, DYY, DZZ, DXCO, DYCO, DZCO, RX1CO, RY1CO, RZ1CO 
       REAL*4:: DXR, DYR, DZR, RX1R, RY1R, RZ1R, RX0, RY0, RZ0, DXO, DXR2, DYR2, DZR2
       CHARACTER(LEN=600) :: NBODYFILE
       INTEGER NL
       INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
       INTEGER NPART(0:NLEVELS)
       INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHZ(NPALEV)
       REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
       REAL*4  TIME_1, TIME_2
       INTEGER DATE(3), TIME(3)
       INTEGER NVOID, IV, NROT, NN, NVOIDT
       INTEGER INIX, FINX, INIY, FINY, INIZ, FINZ
       INTEGER:: IPA, IPA0, LEV, PX, PY, PZ, IR, IR0, IRR, I0, J0, K0, LEV0, LEVM
       INTEGER:: NVOID2,NVOID0, FLAG, NVOIDL
       INTEGER IND, IND0, I1, I2, I3, FLAGX, FLAGY, FLAGZ, NOV, NLL, INDV, INDP, ILEV0
       INTEGER:: FLAGT, N, NR, IXR1, JYR1, KZR1, IXR, JYR, KZR, IXV, JYV, KZV, IX0, JY0, KZ0
       INTEGER:: BASX, BASY, BASZ,  BASX0, BASY0, BASZ0, IXCO, JYCO, KZCO
       REAL*4 RX, RY, RZ, BAS, RX1, RX2, RY1, RY2, RZ1, RZ2       
       REAL*4 A,B,C, AA, BB, CC, UMEANCO, UMEANR, UMEAN0, REQPP
       REAL*4:: RR, DELTA
       CHARACTER*200 DIR, FILEO, FILEO1, FILEO2, FILEO3, FILEO4, FILEO5, FILEO6, FILEO0, FILEO7
       CHARACTER*1 S, LIR

       INTEGER:: NVV, IVV2, IND_OLD, IND_PROG1, IND_PROG2, NVOID_OLD, INDO, NTHRE !V2.0
       INTEGER:: NCELL_MAX, NCELL_MAXO
       REAL*4:: FLAGMAX, FL, VOL1 !V2.0
       !INTEGER, DIMENSION(NVOID_MAX) :: INDICE2, INDICE, INDICE0, INDICEL
       INTEGER, ALLOCATABLE, DIMENSION(:):: INDICE2, INDICE, INDICEL
       INTEGER, ALLOCATABLE, DIMENSION(:)::INDICE_OLD, INDICE_PROG, NCELLV_OLD !V2.0
       REAL, ALLOCATABLE, DIMENSION(:):: SHARED_VOL !V2.0
       INTEGER, ALLOCATABLE:: FLAGVO(:)
       REAL*4, ALLOCATABLE:: INDB(:) !V2.0       
       REAL*4, DIMENSION(NVOID_MAX):: VOL2
       INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: MARCA2, MARCAD!, MARCA
       INTEGER, ALLOCATABLE, DIMENSION(:) :: UVOID
       REAL*4, ALLOCATABLE, DIMENSION(:) :: XC, YC, ZC, VOLNEW, VOLMAX , UMEAN, MTOT, VOLCO, REQ, REQP
       INTEGER, ALLOCATABLE, DIMENSION(:) ::NMERG
       INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MERG
       REAL*4, ALLOCATABLE, DIMENSION(:,:,:):: INERTIA
       REAL*4, ALLOCATABLE, DIMENSION(:,:,:):: UU, UUG, UUS
       REAL*4, ALLOCATABLE,DIMENSION(:):: EPS, IP
       INTEGER, ALLOCATABLE, DIMENSION(:):: ILEV, NCELLV
       REAL*4, DIMENSION(3):: RK, BASEIGENVAL, AXIS
       REAL*4, DIMENSION(NHYX, NHYY, NHYz):: MARCA0
       INTEGER:: NUM,OMP_GET_NUM_THREADS
       INTEGER:: FLAG_VEL_P
       CHARACTER(LEN=50):: SLEV, COMM
       LOGICAL:: LEXIST
       LOGICAL:: FLAG_MASCLET !(True/False)
!   temp !!!!!!!!!
       INTEGER:: ITT,  NHAL, NPARTTOT, II , IH
       REAL*4:: XHAL, YHAL, ZHAL, X, DR, D, DMIN
       REAL*4, DIMENSION(1000) :: MH, XCH, YCH, ZCH, DIST, AGEM, METM 
       INTEGER, DIMENSION(1000):: IDH, NPARTH, LEVH
       INTEGER::  NVOIDH,  INDW, &
            JJ, KK, INDW1, IN, ICONTA, IC, KMIN, KMAX, JMIN, JMAX, IMIN, IMAX
       INTEGER, DIMENSION(1000):: FLAGW, IDHALV
       INTEGER, DIMENSION(NVOID_MAX):: IND2, INDH
       INTEGER, ALLOCATABLE, DIMENSION(:):: NHALV
       INTEGER, PARAMETER:: NCELL_MIN=10 !min # of cells for keeping voids; for few cell voids inertia tensor not reliable

!   temp !!!!!!!!!!!1
       INTEGER, ALLOCATABLE, DIMENSION(:):: COVERAGE, DDX, DDY, DDZ
       REAL*4, ALLOCATABLE, DIMENSION(:):: U0, REF
       REAL*4:: DD
!*********************************************************************


       CALL IDATE(DATE)
       CALL ITIME(TIME)
       WRITE(*,*) 'DATE=',DATE(1),'/',DATE(2),'/',DATE(3)
       WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3)


!*     NUMBER OF PROCESSORS      
      NUM=1
!$OMP PARALLEL SHARED(NUM)
!$OMP SINGLE
!$      NUM=OMP_GET_NUM_THREADS()
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL
      WRITE(*,*)'PROCESSOR =',NUM


!*      OPENING FILES 
       OPEN(1,FILE='voids.dat',STATUS='UNKNOWN',ACTION='READ')                     
                                                                       
!*      READING INITIAL DATA                                       
!****************************************************          
!*      NX,NY,NZ < or = NMAX,NMAY,NMAZ               *       
!****************************************************
       READ(1,*)
       READ(1,*)
       READ(1,*) FLAG_MASCLET
       READ(1,*)
       READ(1, *) NBODYFILE
       READ(1,*)
       read(1,*) FLAG_MAP
       READ(1,*)
       READ(1,*) FIRST,LAST,EVERY 
       READ(1,*)
       !READ(1,*) NHYX2,NHYY2,NHYZ2 !masclet coarse grid; should be equal to   NHYX,NHYY,NHYZ 
       !READ(1,*)
       !READ(1,*) NX,NY,NZ !new coarser grid --> used in MALLA; should be equal to NCOX, NCOY, NCOZ
       !READ(1,*)
       READ(1,*) LEVMIN, LEVMAX !level to be used by voidfind
       READ(1,*)
       READ(1,*) ACHE,OMEGA0
       READ(1,*)
       READ(1,*) ZI,LADO0
       READ(1,*)
       READ(1,*) DENS_THRE !used for the centers
       READ(1,*)
       READ(1,*) DENS_THRE2 !used for the edges
       READ(1,*)
       READ(1,*) GRAD_THRE
       READ(1,*)
       READ(1,*) MINF
       READ(1,*)
       READ(1,*) MAXF
       READ(1,*)
       READ(1,*) RMIN
       CLOSE(1)

       NX = NCOX
       NY = NCOY
       NZ = NCOZ 

       LEVS=0 !put in parameters file?
       LEVP=1 !level for profiles       
       LEVV=LEVMAX !LEVEL OF PARENT VOIDS

       WRITE(*,*) '************************************************'
       WRITE(*,*) '         Input                      '
       WRITE(*,*) '************************************************'
       write(*,*) 'Masclet base mesh:', NHYX,NHYY,NHYZ
       write(*,*) 'Coarser grid',NX,NY,NZ
       write(*,*) 'Threshold values', DENS_THRE, DENS_THRE2, GRAD_THRE


       WRITE(*,*) '************************************************'
       WRITE(*,*) '         Output directory:                     '
       WRITE(*,*) '************************************************'

       WRITE(SLEV, '(I2)') ABS(LEVMIN)
       IF(LEVMIN .LT. 0) SLEV='m'//TRIM(ADJUSTL(SLEV))
       !DIR='output_files_'//TRIM(ADJUSTL(SLEV))
       DIR='output_files'
       ! check wether dir exists or not, if does not create the dir
       INQUIRE(FILE=TRIM(ADJUSTL(DIR))//'/.', EXIST=LEXIST)
       !INQUIRE(FILE='output_files_'//TRIM(ADJUSTL(SLEV))//'/.', EXIST=LEXIST)
       COMM='mkdir '//TRIM(ADJUSTL(DIR))
       IF(LEXIST .EQv. .FALSE. ) CALL SYSTEM(COMM)
       
       WRITE(*,*) 'Output files will be written in: ',TRIM(ADJUSTL(DIR))
       WRITE(*,*) '************************************************'


       NXX=NHYX
       NYY=NHYY
       NZZ=NHYZ
       LADO= LADO0-(LADO0/NXX)
       ALLOCATE(RADX(0:NXX+1))
       ALLOCATE(RADY(0:NYY+1))
       ALLOCATE(RADZ(0:NZZ+1))

!*      GRID BUILDER  --> coarser grid: 128^3: cell size=1.1 Mpc         
          CALL MALLA(NXX,NYY,NZZ,LADO) !USE ONLY RADX
!IF NXX=NHYX --> DX0, DY0, DZ0, RADX0, RADY0, RADZ0


          WRITE(*,*) '************************************************'
          WRITE(*,*) '             GRID FOR LEVEL 0                          '
          WRITE(*,*) '************************************************'
          WRITE(*,*) 'SIDE LENGTH=',LADO 
          WRITE(*,*) 'NX,DX,RADX(1),RADX(NX)=',NXX,DX0,RADX0(1),RADX0(NX)  

                    
!***********************************************************************
!*      COSMOLOGICAL BACKGROUND                                           
!***********************************************************************
       PI=DACOS(-1.D0)                                                   
       UNTERCIO=1.D0/3.D0                                                
       CGR=1.D0/(8.D0*PI)                                                
       CGR2=2.D0*CGR                                                     
       ACHE=ACHE*3.66D-3                                                 
                                                                        
       T0=2.D0*UNTERCIO/ACHE                                             
       RODO=OMEGA0*3.D0*ACHE**2                                                 
       RE0=1.0/10.98                                                     
       ROI=RODO*(1+ZI)**3                                                
       PI4ROD=4.D0*PI*ROI                                                
       REI=RE0/(1.0+ZI)                                                  
      
       ROCRIT=0.0
       ROCRIT=(0.73*0.73)*1.879E-29    !en gramos/cm^3
                                                                  
       TEI=T0*(1.0+ZI)**(-1.5)                                           

       UV=299792.458
       UM=9.1717E18      !en Msol                                           
!***********************************************************************

!* initialization

       NVOIDP = 0
!* special variables for paralelization
       NHYXP=NHYX
       NHYYP=NHYY
       NHYZP=NHYZ
       NMAXP=NMAX
       NMAYP=NMAY
       NMAZP=NMAZ
       NPALEVP=NPALEV
       PARTIREDP=PARTIRED

!$OMP PARALLEL DO SHARED(NHYXP,NHYYP,NHYZP,U1DM,U1G, U1S),PRIVATE(I,J,K)
      DO K=1,NHYZP
      DO J=1,NHYYP
      DO I=1,NHYXP
       U1DM(I,J,K)=-1.0        !valores minimos                                   
       U1G(I,J,K)=-1.0
       U1S(I,J,K)=-1.0
      END DO
      END DO
      END DO


!$OMP PARALLEL DO SHARED(PARTIREDP,U2DM,U3DM,U4DM,RXPA,RYPA,RZPA, &
!$OMP            MASAP,ORIPA), &
!$OMP            PRIVATE(I)
      DO I=1,PARTIREDP
       U2DM(I)=0.0
       U3DM(I)=0.0
       U4DM(I)=0.0
       RXPA(I)=0.0
       RYPA(I)=0.0
       RZPA(I)=0.0
       MASAP(I)=0.0
       ORIPA(I)=0
      END DO 

       NFILE=INT((LAST-FIRST)/EVERY) + 1

!*////////////////////////////////////
       DO IFI=1,NFILE
!*////////////////////////////////////          

!* variable initialization

       PATCHNX=0
       PATCHNY=0
       PATCHNZ=0
       PATCHX=0
       PATCHY=0
       PATCHZ=0
       PATCHRX=0.0
       PATCHRY=0.0
       PATCHRZ=0.0
       NPATCH=0
       PARE=0
       
!* end initialization

       ITER=FIRST+EVERY*(IFI-1)

       WRITE(*,*) '************************************************'
       !WRITE(*,*) 'Reading iteration', ITER
       CALL CPU_TIME(TIME_1)
       !WRITE(*,*) 'Time before reading:', TIME_1


!*      READING DATA
!       CALL LEER_SPH(ITER, NDXYZ0, &!input
!                NL, T,ZETA, NPATCH, PARE, PATCHNX, PATCHNY, & !output
!                PATCHNZ,PATCHX, PATCHY, PATCHZ, &  !output
!                PATCHRX, PATCHRY, PATCHRZ,NPART) !output

!*********************** NBODY ********************88

       !--> ridefinire NL come livello con NPATCH(IR) > 0
       CALL LEER_NBODY(NBODYFILE, NDMXYZ, LADO,RADX0, RADY0, RADZ0, NXX, NYY, NZZ, DX0, DY0, DZ0)
       ALLOCATE(U1CO(NXX, NYY, NZZ))
       ALLOCATE(U2CO(NXX, NYY, NZZ))
       ALLOCATE(U3CO(NXX, NYY, NZZ))
       ALLOCATE(U4CO(NXX, NYY, NZZ))
       FLAG_VEL_P=FLAG_VEL

!$OMP PARALLEL DO SHARED(NXX,NYY,NZZ,U1CO,U2CO,U3CO, U4CO, &
!$OMP   U1DM, UU2DM, UU3DM, UU4DM,FLAG_VEL_P),PRIVATE(IX,JY,KZ)
       DO KZ=1, NZZ
          DO JY=1, NYY
             DO IX=1, NXX
                U1CO(IX,JY,KZ)=U1DM(IX,JY,KZ)
                IF(FLAG_VEL_P==1) THEN
                   U2CO(IX,JY,KZ)=UU2DM(IX,JY,KZ)
                   U3CO(IX,JY,KZ)=UU3DM(IX,JY,KZ)
                   U4CO(IX,JY,KZ)=UU4DM(IX,JY,KZ)
                   U2G(IX,JY,KZ)=UU2DM(IX,JY,KZ) !NO GAS --> GAS = DM
                   U3G(IX,JY,KZ)=UU3DM(IX,JY,KZ)
                   U4G(IX,JY,KZ)=UU4DM(IX,JY,KZ)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

!*********************** NBODY ********************88
       CALL CPU_TIME(TIME_2)
       WRITE(*,*) 'Time spent during reading:', TIME_2-TIME_1

       IF (ZETA.LT.0.0) ZETA=0.0
       ROTE=RODO*(1.0+ZETA)**3
       RETE=RE0/(1.0+ZETA)

       IF (FLAG_MASCLET .EQV. .TRUE.) THEN
          CALL NOMFILE3(DIR, ITER,FILEO,FILEO1,FILEO2, FILEO3, FILEO4, &
               FILEO5, FILEO6, FILEO7) !CHANGE FILENAME !!! 
       ELSE
          CALL NOMFILE_NBODY(DIR, NBODYFILE,FILEO, FILEO1,FILEO2, FILEO3, FILEO4, &
               FILEO5, FILEO6, FILEO7) 
       ENDIF

       OPEN(UNIT=10, FILE=FILEO) !void list
       WRITE(10,*) '# IR, NVOID, FILLING_FRACTION'
       WRITE(10,*) "# ID(1), X(2), Y(3), Z(4), VOLUME(5), REFF(6), DENS_MEAN(7), EPS(8), IP(9), IDP(10), REFFP(11), MASS(12), LEVEL0(13), LEVEL(14), NCELL(15), ID_PROG(16), VOL_PROG(17)"

       IF (FLAG_MAP .EQV. .TRUE.)  OPEN(UNIT=11, FILE=FILEO1, FORM='UNFORMATTED') !BINARY --> map 

       DO IR=LEVMIN, LEVMAX

          LEVP=IR ! o LEVP=1
          LEVP=1
          RMIN0=RMIN

          NXX=REAL(NHYX)*(2.**(IR))
          NYY=REAL(NHYY)*(2.**(IR))
          NZZ=REAL(NHYZ)*(2.**(IR))
          IF(ALLOCATED(RADX) .EQV. .TRUE.) DEALLOCATE(RADX, RADY, RADZ)
          ALLOCATE(RADX(0:NXX+1))
          ALLOCATE(RADY(0:NYY+1))
          ALLOCATE(RADZ(0:NZZ+1))

          WRITE(*,*) '************************************************'
          WRITE(*,*) '             NEW  LEVEL BEING ANALYSED:                              '
          WRITE(*,*) '************************************************'
          WRITE(*,*) 'Level:', IR
          write(*,*) 'Resolution of the grid:', NXX, NYY, NZZ

          LADO= LADO0-(LADO0/NXX)

!*      GRID BUILDER  --> coarser grid: 128^3: cell size=1.1 Mpc         
          CALL MALLA2(NXX,NYY,NZZ,LADO) ! --> RADX, RADY, RADZ, DX, DY, DZ 

          DXX=DX
          DYY=DY
          DZZ=DZ
          RX1=RADX(1)
          RY1=RADY(1)
          RZ1=RADZ(1)
          CONST=DXX*DYY*DZZ*ROTE*RETE**3

          WRITE(*,*) 'ADDITIONAL GRID INFO:'
          WRITE(*,*) 'SIDE LENGTH=',LADO 
          WRITE(*,*) 'NX,DX,RADX(1),RADX(NX)=',NXX,DXX,RADX(1),RADX(NX)  

          IF(ALLOCATED(U1CO) .EQV. .FALSE.) ALLOCATE(U1CO(NXX, NYY, NZZ)) !!! remove allocate in smooth/diver
          IF(ALLOCATED(U1GCO) .EQV. .FALSE.) ALLOCATE(U1GCO(NXX, NYY, NZZ)) 
          IF(ALLOCATED(U1SCO) .EQV. .FALSE.) ALLOCATE(U1SCO(NXX, NYY, NZZ))
          ALLOCATE(DIVERCO(NXX,NYY,NZZ))
          ALLOCATE(DIVERGCO(NXX,NYY,NZZ))
          ALLOCATE(FLAGAMR(NXX,NYY,NZZ))

          FLAGAMR(:,:,:)=1 !FLAGAMR updated only if VMESH is called


          IF(IR .LT. 0) CALL SMOOTH(NPART, NL) !--> U1GCO, U1SCO, U1CO, U2CO, U3CO, U4CO, generalizzarla per lev<-1

          WRITE(*,*) '--------------------------------------------------------------'
          WRITE(*,*) 'SOME CHECKS:' 
          WRITE(*,*) '  min/max U1CO:',MINVAL(U1CO), MAXVAL(U1CO)

          IF(IR .EQ. LEVMIN) THEN
             ALLOCATE(MARCAP(NXX,NYY,NZZ))
             ALLOCATE(REQP(NVOID_MAX))
             MARCAP(:,:,:)=-1 !MARCA for the lower level of the hier.; 0 whne IR=LEVMIN
             REQP(:)=0.
         
             !IF(IR .EQ. LEVMIN) CALL EULERIAN(NL, NDXYZ0, NPATCH, PARE, PATCHNX, PATCHNY, & 
             !     PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
             !     PATCHRX, PATCHRY, PATCHRZ,NPART) !--> dm eulerian vel: U2, U3, U4, U12, U13, U14
             CALL CPU_TIME(TIME_1)
             WRITE(*,*) '  Time spent during Eulerian (sec):', TIME_1-TIME_2

          ENDIF
!*-------------------------------------------------------------------------------*
          !*      Compute divergence : pass level !! --> in output only one var: diverco
!*-------------------------------------------------------------------------------*
         IF(FLAG_VEL == 1) THEN
            CALL DIVER_FINA_GAS(NL, NPATCH, PARE, PATCHNX, PATCHNY, & 
                 PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
                 PATCHRX, PATCHRY, PATCHRZ,NPART, IR) !--> DIVER0, DIVER, DIVERCO, DIVERGCO
            WRITE(*,*) '  Min/Max divergence:', MINVAL(DIVERCO), MAXVAL(DIVERCO), MINVAL(DIVER0), MAXVAL(DIVER0)
            CALL CPU_TIME(TIME_2)
            WRITE(*,*) '  Time spent during Diver_fina (sec):', TIME_2-TIME_1
         ELSE
             ALLOCATE(DIVER0(NHYX, NHYY, NHYZ)) 
             !NB: DIVER IS NOT ALLOCATED HERE, NEEDED ONLY FOR AMR 
         ENDIF

 
          IF(IR .EQ. 0) THEN 
             DO KZ=1, NZZ
                DO JY=1, NYY
                   DO IX=1, NXX
                      U1GCO(IX,JY,KZ)=U1G(IX,JY,KZ)
                      U1SCO(IX,JY,KZ)=U1S(IX,JY,KZ)
                      !U1CO(IX,JY,KZ)=U1(IX,JY,KZ) !CHECK THIS LINE !!FOR LEER_sph SHOULD BE COMMENTED BUT NOT FOR LEER
                      DIVERGCO(IX,JY,KZ)=DIVER0(IX,JY,KZ)
                      DIVERCO(IX,JY,KZ)=DIVER0(IX,JY,KZ)
                   ENDDO
                ENDDO
             ENDDO
             !DEALLOCATE(U1, DIVER0) ?
          ELSE IF (IR .GT. 0) THEN
              WRITE(*,*) '  Check on divergence (DIVER0):', ALLOCATED(DIVER0), DIVER0(50,50,50),  DIVER0(100,100,100)
             CALL VMESH(IR, NXX, NYY, NZZ, DXX, DYY, DZZ, RX1, RY1, RZ1, &
                  NHYX, NHYY, NHYZ, &
                  NPATCH, PARE, PATCHNX, PATCHNY, & 
                  PATCHNZ,PATCHX, PATCHY, PATCHZ, &
                  PATCHRX, PATCHRY, PATCHRZ) !--> UR, UGR, USR, DIVR, FLAGAMR (allocated within the subr)
       ! --> modify VMESH for IR>1 : use IR=1 instead of IR=0 where cells are not refined
             DO KZ=1, NZZ
                DO JY=1,NYY
                   DO IX=1, NXX
                      U1CO(IX,JY,KZ)=UR(IX,JY,KZ)
                      U1GCO(IX,JY,KZ)=UGR(IX,JY,KZ)
                      U1SCO(IX,JY,KZ)=USR(IX,JY,KZ)
                      DIVERGCO(IX,JY,KZ)=DIVR(IX,JY,KZ)
                      DIVERCO(IX,JY,KZ)=DIVR(IX,JY,KZ)
                      !MARCAR(IX,JY,KZ)=MARCA(IX,JY,KZ) !MARCA NOT YET DEFINED
                   ENDDO
                ENDDO
             ENDDO
             DEALLOCATE(UR, UGR, USR, DIVR, FLAGR) 
          ENDIF !IF IR < 1 U1C), DIVERCO already defined in SMOOTH/EULERIAN/DIVER

          !WRITE(*,*) 'Count potential void cells:',COUNT(U1CO .NE. 0), COUNT(DIVERCO .NE. 0)
          IF(FLAG_VEL == 1) THEN! velocity filed defined --> diverV defined
             CALL MARK_ALL(NXX, NYY, NZZ, DENS_THRE) !--> FLAGV, FLAG_SUB (mark cells suitable for being void centers)
          ELSE IF(FLAG_VEL == 0 ) THEN
             CALL MARK_ALL_DENS(NXX, NYY, NZZ, DENS_THRE) 
          ENDIF

          !CHECK: at the first level check if FLAG_SUB=1 in all cells
          CONTA=0
          !IF(IR .EQ. LEVMIN) THEN
             !FLAG_SUB(:,:,:)=0.
             DO KZ=1, NZZ
                DO JY=1, NYY
                   DO IX=1, NXX
                      !IF(FLAG_SUB(IX,JY,KZ) .NE. 0) WRITE(*,*) 'WARNING!: FLAG_SUB>0', FLAG_SUB(IX,JY,KZ), IR
                      IF(FLAGV(IX,JY,KZ).EQ.1) CONTA=CONTA+1
                   ENDDO
                ENDDO
             ENDDO
          !ENDIF

          WRITE(*,*) 'Potential cells in voids / Filling Fraction:',CONTA, CONTA/(NXX**3.)

!*-------------------------------------------*
!*      Find voids
!*-------------------------------------------*
          WRITE(*,*) '************************************************'
          CALL CPU_TIME(TIME_1)
          WRITE(*,*) 'Finding voids...'

!          ALLOCATE(MARCA(NXX,NYY,NZZ))

          VMIN0=0. !no limit for void finder if coarse level is used
          WRITE(*,*) 'RMIN0:', RMIN0
          VMIN0=4.*3.14*(RMIN0**3.)/3.
          CALL VOIDFIND_PAR(NUM, LEVMIN, NXX, NYY, NZZ, DXX, DYY, DZZ, RX1, RY1, RZ1, NVOID_MAX, REQP, &
               DENS_THRE, DENS_THRE2, GRAD_THRE, FLAG_DIV1, VMIN0, RMIN0, &
               NVOID, NPATCH, PARE, PATCHNX, PATCHNY, &
               PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
               PATCHRX, PATCHRY, PATCHRZ) ! UPDATE USING ONLY CELLS WITH FLAGS>1. 
          ! FOR PRINCIPAL VOIDS(AS HERE) PUT FLAGS=1 FOR EVERY CELLS
 
          DEALLOCATE(FLAGV, FLAG_SUB)
          !IF(IR .EQ. LEVMIN) DEALLOCATE(REQP)

          WRITE(*,*) 'Found ',NVOID,' protovoids'
          CALL CPU_TIME(TIME_2)                        
          WRITE(*,*) 'Time spent during  voidfind (sec):',  TIME_2-TIME_1 

!*-------------------------------------------*
!*      Remove overlapping voids   
!*-------------------------------------------*  
          ALLOCATE(INDICE(NVOID))
          ALLOCATE(INDICEL(NVOID))
          ALLOCATE(INDICE2(NVOID))

          DO I=1, NVOID
             INDICE(I)=0
             INDICE2(I)=0
             INDICEL(I)=0
             VOL2(I)=0.D0
          ENDDO
       
          CALL INDEXX(NVOID,VOL(1:NVOID),INDICE2) 

          I1=NVOID     
          DO I=1, NVOID
             VOL2(I1)=VOL(INDICE2(I))
             INDICE(I1)=INDICE2(I) 
             I1=I1-1
          ENDDO

          WRITE(*,*) 'Void sorted:'
          WRITE(*,*) 'Largest void (vol/Re):', vol(indice(1)), ((3.*vol(indice(1)))/(4.*pi))**0.333
          WRITE(*,*) 'Smallest void (vol/Re):', vol(indice(nvoid)), ((3.*vol(indice(nvoid)))/(4.*pi))**0.333

          ALLOCATE(XC(NVOID))
          ALLOCATE(YC(NVOID))
          ALLOCATE(ZC(NVOID))
          ALLOCATE(UVOID(NVOID))
          ALLOCATE(VOLNEW(NVOID))
          ALLOCATE(UMEAN(NVOID))
          ALLOCATE(MTOT(NVOID))
          ALLOCATE(NCELLV(NVOID))
          ALLOCATE(EPS(NVOID))
          ALLOCATE(ILEV(NVOID))
          ALLOCATE(IP(NVOID))
          ALLOCATE(REQ(NVOID))
          ALLOCATE(MARCA(NXX,NYY,NZZ))

          MARCA(:,:,:)=0

          NOVER=NVOID/100
          NOVER=20000
          MINF0=MINF
          MAXF0=MAXF

          WRITE(*,*) '************************************************'
          WRITE(*,*) 'Removing overlapping voids, using MINF,MAXF:', MINF0, MAXF0
          WRITE(*,*) 'Maximum number of overlapping allowed:', NOVER

          CALL OVERLAPPING(MINF0, MAXF0, NOVER, NVOID, INDICE, NXX, NYY, NZZ, DXX, DYY, DZZ, &
               RX1, RY1, RZ1, VOLNEW, UVOID, XC, YC, ZC, &
               NVOID2)

          !WRITE(*,*) 'MARCA > 0', COUNT(MARCA .NE. 0),  REAL(COUNT(MARCA .NE. 0))/(NXX*NYY*NZZ)
          WRITE(*,*) 'Number of non-overlapping voids:', NVOID2


!*-------------------------------------------*
!*      GLOBAL QUANTITIES FOR VOIDS   
!*-------------------------------------------* 
          !clean MARCA for cells not belonging to a real void
          NCELLV(:)=0
          UMEAN(:)=0.
          MTOT(:)=0.
          VOLT3=0.
          VOLT2=0.
          NN=0
          NCELL_MAX=0 !max num of cells in a void
          DO IX=1, NXX
             DO JY=1,NYY
                DO KZ=1,NZZ
                   IF(MARCA(IX,JY,KZ) .LT. 0 ) WRITE(*,*) 'WARNING! MARCA<0'
                   IF(MARCA(IX,JY,KZ) .GT. 0.) THEN !filter void cells
                      IND0= MARCA(IX,JY,KZ)
                      IF(UVOID(MARCA(IX,JY,KZ)) .EQ. -1 ) THEN !keep only marca of the real voids
                         IND0= MARCA(IX,JY,KZ)
                         VOLT3=VOLT3+DXX*DYY*DZZ
                         NCELLV(IND0)=NCELLV(IND0)+1
                         UMEAN(IND0)=UMEAN(IND0)+U1CO(IX,JY,KZ)
                         MTOT(IND0)=MTOT(IND0)+U1CO(IX,JY,KZ)*CONST
                         NN=NN+1
                      ELSE
                         MARCA(IX,JY,KZ)=0
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

          DO IX=1, NXX
             DO JY=1, NYY
                DO KZ=1, NZZ
                   IF(MARCA(IX, JY, KZ) .GT. 0) THEN
                      IND0=MARCA(IX,JY,KZ)
                      IF(NCELLV(IND0) .LT. NCELL_MIN) MARCA(IX,JY,KZ)=0
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

          NCELL_MAX=MAXVAL(NCELLV)
          WRITE(*,*) 'Total volume in voids / FF :', VOLT3/(LADO0**3.D0), COUNT(MARCA .NE. 0)

          CALL SHAPE(NXX, NYY, NZZ, NVOID, INDICE, NCELLV, UVOID, XC, YC, ZC, EPS, IP)

          !write map file: MARCA and density       
          IF (FLAG_MAP .EQV. .True.) THEN
             WRITE(11) (((MARCA(IX,JY,KZ), KZ=1,NZZ), JY=1, NYY), IX=1,NXX)
             WRITE(11) (((U1CO(IX,JY,KZ), KZ=1,NZZ), JY=1, NYY), IX=1,NXX)
             !WRITE(11) (((U1GCO(IX,JY,KZ), KZ=1,NZZ), JY=1, NYY), IX=1,NXX)
             WRITE(11) (((DIVERCO(IX,JY,KZ), KZ=1,NZZ), JY=1, NYY), IX=1,NXX)
          ENDIF

          NVOIDT=COUNT(UVOID .EQ. -1 .AND. NCELLV .GE. NCELL_MIN)
          WRITE(*,*) 'Number of voids:', NVOIDT, count(uvoid .eq. -1), count(ncellv .ge. ncell_min), ncell_min

          WRITE(10,*) IR, NVOIDT, VOLT3/(LADO0**3.D0)
 
          VOLT2=0.
          VOLT3=0.
          VOLNEW(:)=0.


          NTHRE=((4.*PI/2.)*(3./0.678)**3.)/(DX0*DX0*DX0) !Numb of cells for voids with REQ=10mpc/h

          IF(ITER .GT. FIRST) THEN !only if a previous snapshot has been analysed
             !if already run read info from map file
             ALLOCATE(INDICE_PROG(NVOID))
             ALLOCATE(SHARED_VOL(NVOID))
             INDICE_PROG(:)=0
             SHARED_VOL(:)=0

             DO IV=1, NVOID
                IND=INDICE(IV)
                NVV=0
                ALLOCATE(INDB(2*NCELL_MAXO)) !should be sized to the maximum num of cell in a void (sse ncellv in prev snap)
                ALLOCATE(FLAGVO(NVOID_OLD)) !num of cells for each void  of the previous snap
                INDB(:)=0
                FLAGVO(:)=0

                IF(NCELLV(IND) .GE. NTHRE) THEN ! tree only for large voids
                   DO IX=1, NXX
                      DO JY=1, NYY
                         DO KZ=1, NZZ
                            IF(MARCA(IX,JY,KZ) == IND) THEN
                               IF(MARCA_OLD(IX,JY, KZ) .GT. 0) THEN  !cell occcupied by a void in the previous snap
                                  NVV=NVV+1 !number of cells of the actual void occupied by voids in the previous snapshot
                                  IND_OLD=MARCA_OLD(IX,JY,KZ)
                                  INDB(NVV)=IND_OLD
                                  
                                  !FLAGV(MARCA_OLD(IX,IY,KZ))=1
                                  DO IVV2=1, NVV-1
                                     IF(INDB(IVV2) .EQ. INDB(NVV)) THEN
                                        FLAGVO(IND_OLD)=FLAGVO(IND_OLD)+1
                                        EXIT
                                     ENDIF
                                  ENDDO
                                  IF(FLAGVO(IND_OLD) .EQ. 0) FLAGVO(IND_OLD)=1
                               ENDIF
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO

!GET IND FOR WHICH FLAGVO IS MAX
                   FLAGMAX=0.
                   IND_PROG1=0
                   IND_PROG2=0
                   DO IVV2=1, NVOID_OLD
                      INDO=INDICE_OLD(IVV2)
                      FL=REAL(FLAGVO(INDO))!/REAL(NCELLV(IND)*NCELLV_OLD(INDO))
                      IF(FL .GT. FLAGMAX) THEN
                      !IF(FLAGVO(INDO) .GT. FLAGMAX) THEN
                         !FLAGMAX=FLAGVO(INDO)
                         IF(NCELLV_OLD(INDO) == 0) WRITE(*,*) 'WARNING: NCELLV_OLD=0 ', IND, INDO
                         !FLAGMAX=REAL(FLAGVO(INDO))/REAL(NCELLV(IND)*NCELLV_OLD(INDO))
                         FLAGMAX=FL
                         IND_PROG1=INDO
                      ENDIF
                   ENDDO
                   IF(IND_PROG1==0) WRITE(*,*) 'WARNING!! NO PROGENITOR FOUND !!',IND, REQ(IND)
                   INDICE_PROG(IND)=IND_PROG1
                   SHARED_VOL(IND)=FLAGMAX/NCELLV(IND)
                   WRITE(*,*) 'prog:', IND_PROG1, flagmax, flagmax/ncellv(ind)
                ENDIF !filter on void radius

                DEALLOCATE(FLAGVO, INDB)                
                
             ENDDO
             DEALLOCATE(INDICE_OLD, MARCA_OLD, NCELLV_OLD)
             NVOID_OLD=0
          ENDIF !if iter > first


          !save void properties for next snapshot
          ALLOCATE(INDICE_OLD(NVOID))
          ALLOCATE(MARCA_OLD(NXX,NYY,NZZ))
          ALLOCATE(NCELLV_OLD(NVOID))
          NVOID_OLD=NVOID

          NCELL_MAXO=NCELL_MAX
          INDICE_OLD=INDICE
          DO IV=1, NVOID
             IND=INDICE(IV)
             INDICE_OLD(IV)=IND
             NCELLV_OLD(IND)=NCELLV(IND)
          ENDDO
          DO KZ=1, NZZ
             DO JY=1, NYY
                DO IX=1, NXX
                   MARCA_OLD(IX,JY,KZ)=MARCA(IX,JY,KZ)
                ENDDO
             ENDDO
          ENDDO

!end V2.0


          DO IV=1, NVOID !O NVOID2?

             IND=INDICE(IV)
             IF(ITER .GT. FIRST) THEN
                IND_PROG1=INDICE_PROG(IND)
                VOL1=SHARED_VOL(IND)
             ELSE
                IND_PROG1=0
                VOL1=0.
             ENDIF

             !VOLT2=VOLT2+VOLNEW(IND)

             IF(UVOID(IND) .NE. -1 .OR. NCELLV(IND) .LT. NCELL_MIN) CYCLE !use only true voids with more than one cell

             UMEAN(IND)=UMEAN(IND)/NCELLV(IND)
             VOLM=REAL(NCELLV(IND))*DXX*DYY*DZZ
             VOLNEW(IND)=VOLM !Mpc^3
             REQ(IND)=((3.*VOLM)/(4.*PI))**(1./3.) !Mpc

             !find parent void:
             INDP=0
             REQPP=0
             IF(IR .GT. LEVMIN ) THEN
                IXCO=INT(((XC(IND)-RX1CO)/DXCO)+0.49999)+1
                JYCO=INT(((YC(IND)-RY1CO)/DYCO)+0.49999)+1
                KZCO=INT(((ZC(IND)-RZ1CO)/DZCO)+0.49999)+1
                IXCO=INT(((XC(IND)-RX1)/DXX)+0.49999)+1 !this level
                JYCO=INT(((YC(IND)-RY1)/DYY)+0.49999)+1
                KZCO=INT(((ZC(IND)-RZ1)/DZZ)+0.49999)+1
                INDP=MARCAP(IXCO,JYCO,KZCO)
                IF(INDP .GT. 0) REQPP=REQP(INDP) !DEFINE

             ENDIF


             !VOLT3=VOLT3+VOLM
             !       !find patch of the void center
             LEV0=0
             IPA=0
             DO IRR=1, NLEVELS
                DXR2=DXX/(2.**IRR)
                DYR2=DYY/(2.**IRR)
                DZR2=DZZ/(2.**IRR)
                CALL PAFIND(XC(IND),YC(IND),ZC(IND), IRR-1, DXR2, DYR2, DZR2, NLEVELS, NPATCH, PARE, PATCHNX, PATCHNY, & 
                     PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
                     PATCHRX, PATCHRY, PATCHRZ, &
                     IPA0, I0, J0, K0)
                IF(IPA0 .GT. 0) THEN
                   IPA=IPA0
                   LEV0=IRR
                ENDIF
                !
             ENDDO

             LEVM=0 !MAX LEVEL WITHIN THE VOIDS
             IPA=0
             DO IRR=1, NLEVELS
                DXR2=DXX/(2.**IRR)
                DYR2=DYY/(2.**IRR)
                DZR2=DZZ/(2.**IRR)
                
                DO IX=1, NXX
                   DO JY=1, NYY
                      DO KZ=1, NZZ
                         IF(MARCA(IX,JY,KZ) .EQ. IND) THEN
                            RX0=RX1+(IX-1)*DXR2+0.5*DXR2
                            RY0=RY1+(JY-1)*DYR2+0.5*DYR2
                            RZ0=RZ1+(KZ-1)*DZR2+0.5*DZR2
                            CALL PAFIND(RX0,RY0,RZ0, IRR-1, DXR2, DYR2, DZR2, NLEVELS, NPATCH, PARE, PATCHNX, PATCHNY, & 
                                 PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
                                 PATCHRX, PATCHRY, PATCHRZ, &
                                 IPA0, I0, J0, K0)
                            IF(IPA0 .GT. 0) THEN
                               IPA=IPA0
                               LEVM=IRR
                            ENDIF
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          

             WRITE(10,*) IND, XC(IND), YC(IND), ZC(IND), VOLM, REQ(IND), UMEAN(IND), EPS(IND), &
             IP(IND), INDP, REQPP, MTOT(IND)*UM, LEV0, LEVM, NCELLV(IND),IND_PROG1, VOL1

          ENDDO
          WRITE(*,*) 'Min-max radii:', MINVAL(REQ), MAXVAL(REQ)
          WRITE(*,*) 'Min-max OVERDENSITY:', MINVAL(UMEAN), MAXVAL(umean)

!************  test on refinement *************


          NL = NLEVELS
          WRITE(*,*) 'TEST ON REFINEMENTS...', NL
          IND=INDICE(1)
          ALLOCATE(COVERAGE(NL))
          COVERAGE(:)=0
          ALLOCATE(U0(NCELLV(IND)))
          ALLOCATE(DDX(NCELLV(IND)))
          ALLOCATE(DDY(NCELLV(IND)))
          ALLOCATE(DDZ(NCELLV(IND)))
          ALLOCATE(REF(NCELLV(IND)))
          
          N=0
          DO IX=1, NXX
             DO JY=1, NYY
                DO KZ=1, NZZ
                   IF(MARCA(IX,JY,KZ) .EQ. IND) THEN                      

                      N=N+1
                      U0(N)=U1CO(IX,JY,KZ)
                      DDX(N)=IX
                      DDY(N)=JY
                      DDY(N)=KZ
                      REF(N)=0 !LEVEL OF THE CELL

                      RX0=RX1+(IX-1)*DXX+0.5*DXX
                      RY0=RY1+(JY-1)*DYY+0.5*DYY
                      RZ0=RZ1+(KZ-1)*DZZ+0.5*DZZ
                      DD=SQRT((RX0-XC(IND))**2.+(RY0-YC(IND))**2.+(RZ0-ZC(IND))**2.)/REQ(IND)
                      LEV0=0
                      DXO=DXX
                      
                      LEVELS: DO IR0=1, NL
                         DXR=0.5*DXO
                         DYR=0.5*DXO
                         DZR=0.5*DXO
                         RX0=RX0-0.5*DXR
                         RY0=RY0-0.5*DYR
                         RZ0=RZ0-0.5*DZR
                         IPA0=0
                         CALL PAFIND(RX0, RY0, RZ0, IR0-1, DXR, DYR, DZR, NL, NPATCH, PARE, PATCHNX, PATCHNY, & 
                              PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
                              PATCHRX, PATCHRY, PATCHRZ, &
                              IPA0, I0, J0, K0) 

                         IF(IPA0 .GT. 0) THEN
                            REF(N)=IR0
                            LEV0=IR0
                            COVERAGE(LEV0)=COVERAGE(LEV0)+1
                            !EXIT LEVELS
                         ENDIF
                         DXO=DXR
                      ENDDO LEVELS
                      write(98,*) rx0, ry0, rz0,DD, ref(n), u0(n)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

          
          IMAX=MAXLOC(U0,1)

          !WRITE(*,*) 'MAX DENSITY:', U0(IMAX), MAXVAL(U0)
          DEALLOCATE(U0, COVERAGE, DDX, DDY, DDZ, REF)
!***********************************************

          WRITE(*,*) '************************************************'
          WRITE(*,*) 'Void list written in:', FILEO
          !CLOSE(10)

          !find galaxies/haloes in voids and void walls:
          !IF(IR .EQ. LEVMIN) CALL GALAXIES(ITER, FILEO6,  NXX, NYY, NZZ, NVOID, XC, YC, ZC, REQ)
          !IF(IR .EQ. LEVMIN) CALL HALOES(ITER, FILEO7,  NXX, NYY, NZZ, NVOID, XC, YC, ZC, REQ)

          !select voids larger than rmin (NVOIDL) and find level of refinement in the void center
          NVOIDL=0
          INDICEL=0
          CONTA=0
          DO IV=1, NVOID
             IND=INDICE(IV)
             IF(UVOID(IND) .EQ. -1 .AND. 0.73*REQ(IND) .GE. RMIN0 &
                  .AND. NCELLV(IND) .GE. NCELL_MIN ) THEN !use rmin condition only for voids at levmin
                NVOIDL=NVOIDL+1
                INDICEL(NVOIDL)=IND
                !add pafind to find levels where center is located
             ENDIF
          ENDDO


          IF(LEVP .EQ. 1 .AND. IR .LT. 1 .AND. FLAG_PROF .EQ. 1) THEN
             
             NXR=NHYX*(2.**LEVP)
             NYR=NHYY*(2.**LEVP)
             NZR=NHYZ*(2.**LEVP)
             WRITE(*,*) 'PROFILES: LEVEL=', LEVP, NXR, NYR, NZR

             IF(ALLOCATED(RADX) .EQV. .TRUE.) DEALLOCATE(RADX, RADY, RADZ)
             ALLOCATE(RADX(0:NXR+1))
             ALLOCATE(RADY(0:NYR+1))
             ALLOCATE(RADZ(0:NZR+1))

             LADO= LADO0-(LADO0/NXR)
             !*      GRID BUILDER  --> coarser grid: 128^3: cell size=1.1 Mpc         
             CALL MALLA2(NXR,NYR,NZR,LADO) ! --> RADX, RADY, RADZ, DX, DY, DZ

             DXR=DX
             DYR=DY
             DZR=DZ
             RX1R=RADX(1)
             RY1R=RADY(1)
             RZ1R=RADZ(1)

             IF(ALLOCATED(U2).EQV. .TRUE.) DEALLOCATE(U2, U3, U4, U12, U13, U14)
             IF(ALLOCATED(FLAGAMR).EQV. .TRUE.) DEALLOCATE(FLAGAMR)
             ALLOCATE(FLAGAMR(NXR, NYR, NZR))

             !IF(IR .EQ. LEVMIN) 
             !CALL EULERIAN(NL, NDXYZ0, NPATCH, PARE, PATCHNX, PATCHNY, & 
             !     PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
             !     PATCHRX, PATCHRY, PATCHRZ,NPART) !--> dm eulerian vel: U2, U3, U4, U12, U13, U14
             CALL CPU_TIME(TIME_1)
             WRITE(*,*) 'Time spent during Eulerian (sec):', TIME_1-TIME_2

             CALL DIVER_FINA_GAS(NL, NPATCH, PARE, PATCHNX, PATCHNY, & 
                  PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
                  PATCHRX, PATCHRY, PATCHRZ,NPART, LEVP) !--> DIVER0, DIVER, DIVERCO, DIVERGCO
             CALL VMESH(LEVP, NXR, NYR, NZR, DXR, DYR, DZR, RX1R, RY1R, RZ1R, &
                  NHYX, NHYY, NHYZ, &
                  NPATCH, PARE, PATCHNX, PATCHNY, & 
                  PATCHNZ,PATCHX, PATCHY, PATCHZ, &
                  PATCHRX, PATCHRY, PATCHRZ) !--> UR, UGR, DIVR, FLAGAMR (allocated within the subr)
             ! --> modify VMESH for IR>1 : use IR=1 instead of IR=0 where cells are not refined
             
             ALLOCATE(UU(NXR, NYR, NZR))
             ALLOCATE(UUG(NXR, NYR, NZR))
             ALLOCATE(UUS(NXR, NYR, NZR))
             DO KZ=1, NZR
                DO JY=1,NYR
                   DO IX=1, NXR
                      UU(IX,JY,KZ)=UR(IX,JY,KZ) !LEV LEVP
                      UUG(IX,JY,KZ)=UGR(IX,JY,KZ) !LEV LEVP
                      UUS(IX,JY,KZ)=USR(IX,JY,KZ) !LEV LEVP
                   ENDDO
                ENDDO
             ENDDO
             DEALLOCATE(UR, UGR, USR, DIVR, FLAGR, DIVER, DIVER0) 
             !DEALLOCATE(U2, U3, U4, U12, U13, U14)
             !IF(ALLOCATED(U2).EQV. .TRUE.) DEALLOCATE(U2, U3, U4, U12, U13, U14)

          ELSE

             NXR=NXX
             NYR=NYY
             NZR=NZZ
             DXR=DXX
             DYR=DYY
             DZR=DZZ
             RX1R=RX1
             RY1R=RY1
             RZ1R=RZ1
             
             ALLOCATE(UU(NXR, NYR, NZR))
             ALLOCATE(UUG(NXR, NYR, NZR))
             ALLOCATE(UUS(NXR, NYR, NZR))
             DO KZ=1, NZR
                DO JY=1,NYR
                   DO IX=1, NXR
                      UU(IX,JY,KZ)=U1CO(IX,JY,KZ) !LEV IR 
                      UUG(IX,JY,KZ)=U1GCO(IX,JY,KZ) 
                      UUS(IX,JY,KZ)=U1SCO(IX,JY,KZ)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          
          IF (FLAG_PROF == 1) CALL PROFILES_SQ(NPART, NL, RODO, ZETA, NVOID,  NVOIDL, &
               XC, YC, ZC, INDICE, INDICEL, UVOID,VOLNEW, REQ, EPS, UMEAN,  IR, LEVP, &
               UU, UUG, UUS, NXR,NYR,NZR, DXR, DYR, DZR, RX1R, RY1R, RZ1R, &
               NPATCH, PARE, PATCHNX, PATCHNY, & 
               PATCHNZ,PATCHX, PATCHY, PATCHZ, &
               PATCHRX, PATCHRY, PATCHRZ, FILEO2 ) !level 1

          DEALLOCATE(UU, UUG, UUS)

          WRITE(*,*) '************************************************'

          IF(IR .EQ. LEVMAX) THEN
             !DEALLOCATE(U12, U13, U14)
             !DEALLOCATE(U2, U3, U4) !deallocate u2g, u3g, u4g ?
             IF (FLAG_DMO .EQV. .FALSE.) DEALLOCATE(U11G, U11S, U12G, U13G, U14G)
             IF(LEVMIN .LT. 0) THEN
                DEALLOCATE(U2CO, U3CO, U4CO)
                DEALLOCATE(U2GCO, U3GCO, U4GCO)
             ENDIF
          ENDIF

!-------------------------------------------------------
! Store variables used in the next level
!--------------------------------------------------------
         !IF(ALLOCATED(MARCAP) .EQV. .TRUE) 
          DEALLOCATE(MARCAP,REQP)
          NX2=NXX*2
          NY2=NYY*2
          NZ2=NZZ*2
          WRITE(*,*) NX2, NY2, NZ2
          ALLOCATE(MARCAP(NX2,NY2,NZ2)) 
          ALLOCATE(REQP(NVOID))

          DO IV=1, NVOID !NVOIDL --> save only large voids ?
             IND=INDICE(IV)
             !VOLP(IND)=VOLNEW(IND) !volume of the parent void
             REQP(IND)=REQ(IND)
          ENDDO
          NVOIDP=NVOID

          DO KZ=1, NZ2
             DO JY=1, NY2
                DO IX=1, NX2
                   IXCO=INT((IX-1)/2.)+1
                   JYCO=INT((JY-1)/2.)+1
                   KZCO=INT((KZ-1)/2.)+1                 
                   MARCAP(IX,JY,KZ)=MARCA(IXCO,JYCO,KZCO)
                ENDDO
             ENDDO
          ENDDO

          DXCO=DXX
          DYCO=DYY
          DZCO=DZZ
          RX1CO=RX1
          RY1CO=RY1
          RZ1CO=RZ1

          DEALLOCATE(U1CO, U1GCO, U1SCO, DIVERCO, DIVERGCO, FLAGAMR,MARCA)
          DEALLOCATE(XC,YC,ZC, UVOID, VOLNEW, UMEAN, MTOT, NCELLV, EPS, ILEV, IP, REQ)
          DEALLOCATE(INDICE, INDICEL, INDICE2)
          DEALLOCATE(RADX, RADY, RADZ)
          IF(ALLOCATED(INDICE_PROG) .EQV. .TRUE.) DEALLOCATE(INDICE_PROG, SHARED_VOL) !V2.0
       ENDDO ! loop on hierarchies

       CLOSE(10)
       CLOSE(11)

       DEALLOCATE(MARCAP, REQP)
!??????????????????
        IF(ALLOCATED(DIVER) .EQV. .TRUE.) DEALLOCATE(DIVER)   !??????????????????
        IF(ALLOCATED(DIVER0) .EQV. .TRUE.) DEALLOCATE(DIVER0)  
        IF(ALLOCATED(U11) .EQV. .TRUE.) DEALLOCATE(U11)
        

        
      ENDDO !LOOP ON ITERATIONS: IFI=1, NFILE


       CALL IDATE(DATE)
       CALL ITIME(TIME)
       WRITE(*,*) 'DATE=',DATE(1),'/',DATE(2),'/',DATE(3)
       WRITE(*,*) 'TIME=',TIME(1),':',TIME(2),':',TIME(3)

    END PROGRAM VOID_FINDER



!************************************************************************* 
      SUBROUTINE MALLA(NX,NY,NZ,LADO)                                    
!************************************************************************* 
      USE COMMONDATA, ONLY: RADX0,RADY0,RADZ0, DX0,DY0,DZ0, &
         RADX, RADY, RADZ, DX, DY, DZ, RADXR, RADYR, RADZR, DXRR, DYRR, DZRR, NCOX, NHYX
      IMPLICIT NONE                                                     

!input variables
      INTEGER NX,NY,NZ
      REAL*4 LADO
!local variables
      INTEGER I,J,K
      REAL*4 A,B,C
!ouptut
      REAL*4 RADXX(0:NX+1), RADYY(0:NY+1), RADZZ(0:NZ+1)
      REAL*4 DXX, DYY, DZZ
                                                                                                                        
       A=-LADO/2.0                                                          
       B=LADO/2.0                                                           

!*     GRID                                                              
                                                                        
!*     X-AXIS                                                            
      C=(B-A)/(NX-1)                                                    
      RADXX(1)=A                                                         
      DO I=2,NX                                                        
        RADXX(I)=RADXX(1)+(I-1)*C                                            
      END DO                                                           
                                                                        
!*     FICTICIUS CELLS                                                   
      RADXX(0)=RADXX(1)-C                                                 
      RADXX(NX+1)=RADXX(NX)+C                                             
       
!*     Y-AXIS                                                            
      C=(B-A)/(NY-1)                                                    
      RADYY(1)=A                                                         
      DO J=2,NY                                                        
        RADYY(J)=RADYY(1)+(J-1)*C                                            
      END DO                                                           

!*     FICTICIUS CELLS                                                   
      RADYY(0)=RADYY(1)-C                                                 
      RADYY(NY+1)=RADYY(NY)+C                                             

!*     Z-AXIS                                                            
      C=(B-A)/(NZ-1)                                                    
      RADZZ(1)=A                                                         
      DO K=2,NZ                                                        
        RADZZ(K)=RADZZ(1)+(K-1)*C                                            
      END DO                                                           
                                                                        
!*     FICTICIUS CELLS                                                   
      RADZZ(0)=RADZZ(1)-C                                                 
      RADZZ(NZ+1)=RADZZ(NZ)+C                                             
                                                                                                                          
                                                                        
      DXX=RADXX(2)-RADXX(1)                                                        
      DYY=RADYY(2)-RADYY(1)                                                        
      DZZ=RADZZ(2)-RADZZ(1) 

      !update common variables
      IF(NX==NCOX) THEN !coarse grid
         RADX(0:NX+1)=RADXX(0:NX+1) 
         RADY(0:NY+1)=RADYY(0:NY+1)
         RADZ(0:NZ+1)=RADZZ(0:NZ+1)
         DX=DXX
         DY=DYY
         DZ=DZZ
      ENDIF
      IF(NX==NHYX) THEN !base grid
         RADX0(0:NX+1)=RADXX(0:NX+1)
         RADY0(0:NY+1)=RADYY(0:NY+1)
         RADZ0(0:NZ+1)=RADZZ(0:NZ+1)
         DX0=DXX
         DY0=DYY
         DZ0=DZZ
      ENDIF
      IF(NX==2*NHYX) THEN !LEV 1
         IF(ALLOCATED(RADXR) .EQV. .TRUE.) DEALLOCATE(RADXR, RADYR, RADZR)
         ALLOCATE(RADXR(0:NX+1))
         ALLOCATE(RADYR(0:NY+1))
         ALLOCATE(RADZR(0:NZ+1))
         RADXR(0:NX+1)=RADXX(0:NX+1)
         RADYR(0:NY+1)=RADYY(0:NY+1)
         RADYR(0:NZ+1)=RADZZ(0:NZ+1)
         DXRR=DXX
         DYRR=DYY
         DZRR=DZZ
      ENDIF

      RETURN                                                            
    END SUBROUTINE MALLA


!************************************************************************* 
      SUBROUTINE MALLA2(NX,NY,NZ,LADO)                                    
!************************************************************************* 
      USE COMMONDATA, ONLY: RADX, RADY, RADZ, DX, DY, DZ
!, ONLY: RADX0,RADY0,RADZ0, DX0,DY0,DZ0, &
         !RADX, RADY, RADZ, DX, DY, DZ, RADXR, RADYR, RADZR, DXRR, DYRR, DZRR, NCOX, NHYX
      IMPLICIT NONE                                                     

!input variables
      INTEGER NX,NY,NZ
      REAL*4 LADO
!local variables
      INTEGER I,J,K
      REAL*4 A,B,C
!ouptut
      REAL*4 RADXX(0:NX+1), RADYY(0:NY+1), RADZZ(0:NZ+1)
      REAL*4 DXX, DYY, DZZ
                                                                                                                        
       A=-LADO/2.0                                                          
       B=LADO/2.0                                                           
                                                                                                                                           
!*     GRID                                                              
                                                                        
!*     X-AXIS                                                            
      C=(B-A)/(NX-1)                                                    
      RADXX(1)=A                                                         
      DO I=2,NX                                                        
        RADXX(I)=RADXX(1)+(I-1)*C                                            
      END DO                                                           
                                                                        
!*     FICTICIUS CELLS                                                   
      RADXX(0)=RADXX(1)-C                                                 
      RADXX(NX+1)=RADXX(NX)+C                                             
       
!*     Y-AXIS                                                            
      C=(B-A)/(NY-1)                                                    
      RADYY(1)=A                                                         
      DO J=2,NY                                                        
        RADYY(J)=RADYY(1)+(J-1)*C                                            
      END DO                                                           

!*     FICTICIUS CELLS                                                   
      RADYY(0)=RADYY(1)-C                                                 
      RADYY(NY+1)=RADYY(NY)+C                                             

!*     Z-AXIS                                                            
      C=(B-A)/(NZ-1)                                                    
      RADZZ(1)=A                                                         
      DO K=2,NZ                                                        
        RADZZ(K)=RADZZ(1)+(K-1)*C                                            
      END DO                                                           
                                                                        
!*     FICTICIUS CELLS                                                   
      RADZZ(0)=RADZZ(1)-C                                                 
      RADZZ(NZ+1)=RADZZ(NZ)+C                                             
                                                                                                                          
                                                                        
      DXX=RADXX(2)-RADXX(1)                                                        
      DYY=RADYY(2)-RADYY(1)                                                        
      DZZ=RADZZ(2)-RADZZ(1) 

 
         RADX(0:NX+1)=RADXX(0:NX+1) 
         RADY(0:NY+1)=RADYY(0:NY+1)
         RADZ(0:NZ+1)=RADZZ(0:NZ+1)
         DX=DXX
         DY=DYY
         DZ=DZZ
 
      RETURN                                                            
    END SUBROUTINE MALLA2


!************************************************************************
       SUBROUTINE NOMFILE(ITER,FILNOM1,FILNOM2,FILNOM3, FILNOM4)
!************************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER*9 FILNOM1,FILNOM2,FILNOM4
       CHARACTER*10 FILNOM3
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT
      
       CONTA=0
      
       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO
      
       FILNOM1='clus'//NOM
       FILNOM2='cldm'//NOM
       FILNOM3='grids'//NOM
       FILNOM4='clst'//NOM
      
       RETURN
     END SUBROUTINE NOMFILE

                             
!**********************************************************************
       SUBROUTINE NOMFILE3(DIR, ITER,FILNOMO, FILNOMO1, FILNOMO2, &
            FILNOMO3, FILNOMO4, FILNOMO5, FILNOMO6, FILNOMO7) 
!**********************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER(LEN=*) DIR, FILNOMO, FILNOMO2, FILNOMO1, FILNOMO3, FILNOMO4, &
            FILNOMO5, FILNOMO6, FILNOMO7
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT
      
       CONTA=0
      
       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO
      

       FILNOMO=TRIM(ADJUSTL(DIR))//'/voids'//NOM !list of voids
       FILNOMO1=TRIM(ADJUSTL(DIR))//'/map'//NOM  
       FILNOMO2=TRIM(ADJUSTL(DIR))//'/profiles'//NOM 
       FILNOMO3=TRIM(ADJUSTL(DIR))//'/profiles_ref'//NOM 
       FILNOMO4=TRIM(ADJUSTL(DIR))//'/subvoids'//NOM !list of subvoids
       FILNOMO5=TRIM(ADJUSTL(DIR))//'/map_sub'//NOM
       FILNOMO6=TRIM(ADJUSTL(DIR))//'/galaxies'//NOM !galaxies in walls
       FILNOMO7=TRIM(ADJUSTL(DIR))//'/haloes'//NOM

       RETURN
     END SUBROUTINE NOMFILE3
      
!**********************************************************************
       SUBROUTINE NOMFILE_NBODY(DIR,INPUT, FILNOMO, FILNOMO1, FILNOMO2, &
            FILNOMO3, FILNOMO4, FILNOMO5, FILNOMO6, FILNOMO7) 
!**********************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER(LEN=*) DIR, INPUT, FILNOMO, FILNOMO2, FILNOMO1, FILNOMO3, FILNOMO4, &
            FILNOMO5, FILNOMO6, FILNOMO7
       CHARACTER(LEN=100):: INPUT_SHORT
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT, startin

       STARTIN = SCAN(INPUT,'/', BACK=.TRUE.)
       INPUT_SHORT=INPUT(STARTIN+1:)

       FILNOMO=TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(INPUT_SHORT))//'_voids' !list of voids
       FILNOMO1=TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(INPUT_SHORT))//'_map'
       FILNOMO2=TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(INPUT_SHORT))//'_profiles'
       FILNOMO3=TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(INPUT_SHORT))//'_profiles_ref'
       FILNOMO4=TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(INPUT_SHORT))//'_subvoids' !list of subvoids
       FILNOMO5=TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(INPUT_SHORT))//'_map_sub'
       FILNOMO6=TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(INPUT_SHORT))//'_galaxies' !galaxies in walls
       FILNOMO7=TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(INPUT_SHORT))//'_haloes'

       RETURN
     END SUBROUTINE NOMFILE_NBODY

!**********************************************************************
       SUBROUTINE NOMFILE2(DIR, ITER,FILNOMO) 
!**********************************************************************

       IMPLICIT NONE
       INTEGER ITER
       CHARACTER(LEN=*) DIR, FILNOMO
       CHARACTER*5 NOM
       INTEGER CONTA,I,N10,IT
      
       CONTA=0
      
       DO I=4,0,-1
         N10=10**I
         IT=ITER/N10 - CONTA
         CONTA=(CONTA+IT)*10
         NOM(5-I:5-I)=CHAR(48+IT)
       END DO
             
       FILNOMO=TRIM(ADJUSTL(DIR))//'/families'//NOM !list of voids


       RETURN
     END SUBROUTINE NOMFILE2

!************************************************************************* 
       SUBROUTINE LEER(ITER, NDXYZ0,  &!input
                NL, T,ZETA, NPATCH, PARE, PATCHNX, PATCHNY,  &!output
                PATCHNZ,PATCHX, PATCHY, PATCHZ,   &!output
                PATCHRX, PATCHRY, PATCHRZ,NPART) !output
!************************************************************************* 
 
      USE COMMONDATA 
      IMPLICIT NONE                                                     
      
!C input variables
       INTEGER ITER!, NHYX, NHYY, NHYZ

!C local variables
       INTEGER IR,  IRR, LOW1,LOW2, N1, N2, N3
       INTEGER I, IX, J, K
       INTEGER CONTA, PARTIBAS
       INTEGER:: DIM1, DIM2, DIM3, DIM4, NHYXP, NHYYP, NHYZP, FLAG_DENS_P
       REAL*4 AAA,BBB,CCC, MAP
       REAL*4, ALLOCATABLE::UBAS(:)
       INTEGER,ALLOCATABLE::UBAS2(:)
       CHARACTER*9 FILNOM1,FILNOM2, FILNOM4
       CHARACTER*10 FILNOM3
       CHARACTER*25 FIL1,FIL2, FIL4
       CHARACTER*26 FIL3 

!C output variables
       INTEGER NL, NDXYZ0, NST0
       REAL*4 T,ZETA
       INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
       INTEGER NPART(0:NLEVELS),  NPARTST(0:NLEVELS)
       INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHZ(NPALEV)
       REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
!TEMPORARY VAR
       INTEGER:: IPA, IX0, JY0, KZ0, IXR, JYR, KZR
       INTEGER:: IX01, IX02, JY01, JY02, KZ01, KZ02
       REAL*4:: UMEAN, UMEANG, UMEANDM
       REAL*4, ALLOCATABLE:: UPA(:,:,:)
       CHARACTER(LEN=3):: IL

!*      READING DATA
       CALL NOMFILE(ITER,FILNOM1,FILNOM2,FILNOM3,FILNOM4)

       !WRITE(*,*) 'Reading iter',ITER,' ',FILNOM1,FILNOM2,FILNOM3

       FIL1='simu_masclet/'//FILNOM1 ! gas 
       FIL2='simu_masclet/'//FILNOM2 ! dm
       FIL3='simu_masclet/'//FILNOM3 ! grid
       FIL4='simu_masclet/'//FILNOM4 ! stars

       OPEN (33,FILE=FIL3,STATUS='UNKNOWN',ACTION='READ')
       OPEN (31,FILE=FIL1, &
            STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')
       OPEN (32,FILE=FIL2, &
            STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')
       OPEN (34,FILE=FIL4, &
            STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')

!*      GRID DATA
       READ(33,*) IRR,T,NL,MAP
       READ(33,*) ZETA
       READ(33,*) IR,NDXYZ0, NST0
       WRITE(*,*) 'IR,NL,NDXYZ,NST0, MAP', IR,NL,NDXYZ0,NST0, MAP
       
       DO IR=1,NL
       READ(33,*) IRR,NPATCH(IR), NPART(IR),NPARTST(IR)
       !WRITE(*,*), 'NPATCH(IR), NPART(IR)',NPATCH(IR), NPART(IR) 
       READ(33,*)
       
       IF (IR.NE.IRR) WRITE(*,*)'Warning: fail in restart'
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        READ(33,*) PATCHNX(I),PATCHNY(I),PATCHNZ(I)
        READ(33,*) PATCHX(I),PATCHY(I),PATCHZ(I)
        READ(33,*) AAA,BBB,CCC
        PATCHRX(I)=AAA
        PATCHRY(I)=BBB
        PATCHRZ(I)=CCC
        READ(33,*) PARE(I)
       END DO
       END DO
       CLOSE(33)
       NPART(0)=NDXYZ0
       NPARTST(0)=NST0

       WRITE(*,*) 'NPART (LEV0, LEV>0, ALL):',NDXYZ0, SUM(NPART(1:NL)), SUM(NPART)


!***************** ALLOCATE AND INITIALIZE ARRAYS *******************************************
       DIM1=MAXVAL(PATCHNX)
       DIM2=MAXVAL(PATCHNY)
       DIM3=MAXVAL(PATCHNZ)
       !DIM1=64
       !DIM2=64
       !DIM3=64
       DIM4=SUM(NPATCH(0:NL))
       ALLOCATE(U11G(DIM1, DIM2, DIM3, DIM4))
       ALLOCATE(U12G(DIM1, DIM2, DIM3, DIM4))
       ALLOCATE(U13G(DIM1, DIM2, DIM3, DIM4))
       ALLOCATE(U14G(DIM1, DIM2, DIM3, DIM4))
       ALLOCATE(U11DM(DIM1, DIM2, DIM3, DIM4))
       ALLOCATE(U11S(DIM1, DIM2, DIM3, DIM4))

!$OMP PARALLEL DO SHARED(DIM1, DIM2, DIM3, DIM4, U11DM, U11G, U11S),PRIVATE(IX, I,J,K)
       DO I=1,DIM4
        DO K=1,DIM3
        DO J=1,DIM2
        DO IX=1,DIM1
         U11DM(IX,J,K,I)=-1.0 
         U11G(IX,J,K,I)=-1.0  
         U11S(IX,J,K,I)=-1.0 
        END DO
        END DO
        END DO
       END DO
!************************************************************************

!*      GAS
       READ(31) 
       IR=0
        N1=NHYX
        N2=NHYY
        N3=NHYZ
        READ(31) (((U1G(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U2G(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U3G(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U4G(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(31) !(((PRES(I,J,K),K=1,N3),J=1,N2),I=1,N1)
        READ(31) !(((POT(I,J,K),K=1,N3),J=1,N2),I=1,N1)
        READ(31) 
        READ(31)   !CAUTION with this line!! depends on MASCLET version
        READ(31)  !!new: metalicity!! depends on MASCLET version!
        READ(31) !NEW!! since September 2013
       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        IF(N1 .GT. DIM1 .OR. N2 .GT. DIM2 .OR. N3 .GT. DIM3 .OR. I .GT. DIM4) THEN
           WRITE(*,*) 'Bad dimension for U11G', N1,N2,N3,I
           STOP
        ENDIF
        READ(31) (((U11G(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U12G(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U13G(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        READ(31) (((U14G(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
        READ(31) !(((PRES21(IX,J,K,I),K=1,N3),J=1,N2),IX=1,N1)
        READ(31) !(((POT1(IX,J,K,I),K=1,N3),J=1,N2),IX=1,N1)
        READ(31) 
        READ(31) !CAUTION with this line!! depends on MASCLET version
        READ(31)  !!new: metalicity!! depends on MASCLET version!
        READ(31) !NEW!! since September 2013
        READ(31) !NEW!! since September 2013
       END DO
       END DO

      CLOSE(31)


!**     DARK MATTER
       READ(32) 
       IR=0
        N1=NHYX
        N2=NHYY
        N3=NHYZ
        
        READ(32) (((U1DM(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(32) (RXPA(I),I=1,NDXYZ0)
        READ(32) (RYPA(I),I=1,NDXYZ0)
        READ(32) (RZPA(I),I=1,NDXYZ0)
        
        READ(32) (U2DM(I),I=1,NDXYZ0)
        READ(32) (U3DM(I),I=1,NDXYZ0)
        READ(32) (U4DM(I),I=1,NDXYZ0)
        READ(32) (ORIPA(I),I=1,NDXYZ0)
        CONTA=NDXYZ0
        MASAP(1:NDXYZ0)=MAP
      
        write(*,*) 'MAX-MIN MP:',maxval(masap), minval(masap)
        !WRITE(*,*) 'NPART(IR)=',IR, NPART(IR),CONTA
        
       PARTIBAS=MAXVAL(NPART)
       ALLOCATE(UBAS(PARTIBAS))
       ALLOCATE(UBAS2(PARTIBAS))

       DO IR=1,NL
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1,LOW2                
        N1=PATCHNX(I)
        N2=PATCHNY(I)
        N3=PATCHNZ(I)
        READ(32) (((U11DM(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
       END DO 
        
        UBAS=0.0
        UBAS2=0
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RXPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RYPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        RZPA(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U2DM(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U3DM(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        U4DM(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS(IX),IX=1,NPART(IR))
        MASAP(CONTA+1:CONTA+NPART(IR))=UBAS(1:NPART(IR))
        READ(32) (UBAS2(IX),IX=1,NPART(IR))
        IF (NPART(IR).GT.0) &
        ORIPA(CONTA+1:CONTA+NPART(IR))=UBAS2(1:NPART(IR))
        
        IF (NPART(IR).GT.0) THEN
        !WRITE(*,*) 'ORIPA=',MAXVAL(ORIPA(CONTA+1:CONTA+NPART(IR))), &
        !                    MINVAL(ORIPA(CONTA+1:CONTA+NPART(IR)))
        END IF
        
        CONTA=CONTA+NPART(IR)
        WRITE(*,*) 'NPART(IR)=',IR,NPART(IR),CONTA
        NPARTT=CONTA
                
        END DO
       CLOSE(32)

       WRITE(*,*) 'MIN-MAX MASAP:', 9.1717E18*MINVAL(MASAP(1:npartt)), 9.1717E18*MAXVAL(MASAP(1:npartt))

!*** STARS
        READ(34)
        IR=0
        N1=NHYX
        N2=NHYY
        N3=NHYZ


        READ(34) (((U1S(I,J,K),I=1,N1),J=1,N2),K=1,N3)
        READ(34) !(RXPA(I),I=1,NST0)
        READ(34) !(RYPA(I),I=1,NST0)
        READ(34) !(RZPA(I),I=1,NST0)        
        READ(34) !(U2DM(I),I=1,NST0)
        READ(34) !(U3DM(I),I=1,NST0)
        READ(34) !(U4DM(I),I=1,NST0)
        READ(34) !(MASAP(I),I=1,NST0)
        READ(34) !(TEST(I),I=1,NST0)
        READ(34) !(MET(I),I=1,NST0)


        DO IR=1, NL
           LOW1=SUM(NPATCH(0:IR-1))+1
           LOW2=SUM(NPATCH(0:IR))
           WRITE(*,*) LOW1,LOW2
           DO I=LOW1,LOW2                
              N1=PATCHNX(I)
              N2=PATCHNY(I)
              N3=PATCHNZ(I)
              READ(34) (((U11S(IX,J,K,I),IX=1,N1),J=1,N2),K=1,N3)
           END DO
           READ(34) !(RXPA(I),I=1,NPARTST(IR))
           READ(34) !(RYPA(I),I=1,NPARTST(IR))
           READ(34) !(RZPA(I),I=1,NPARTST(IR))        
           READ(34) !(U2DM(I),I=1,NPARTST(IR))
           READ(34) !(U3DM(I),I=1,NPARTST(IR))
           READ(34) !(U4DM(I),I=1,NPARTST(IR))
           READ(34) !(MASAP(I),I=1,NPARTST(IR))
           READ(34) !(TEST(I),I=1,NPARTST(IR))
           READ(34) !(MET(I),I=1,NPARTST(IR))
           READ(34) !(ORIPA1(I),I=1,NPARTST(IR))
        ENDDO
        CLOSE(34)

!****  TOTAL DENSITY CONTRAST: U1G+U1DM

        DEALLOCATE(UBAS, UBAS2)
        ALLOCATE(U11(DIM1, DIM2, DIM3, DIM4))    

        FLAG_DENS_P=FLAG_DENS
  
!$OMP PARALLEL DO SHARED(DIM1, DIM2, DIM3, DIM4, U11DM, U11G, U11S, U11, FLAG_DENS_P),PRIVATE(IX, I,J,K) 
        DO I=1,DIM4
        DO K=1,DIM3
        DO J=1,DIM2
        DO IX=1,DIM1
           IF( FLAG_DENS_P .EQ. 0) THEN
              U11(IX,J,K,I)=U11DM(IX,J,K,I)+1.+U11G(IX,J,K,I)+1.
           ELSE IF(FLAG_DENS_P .EQ. 1) THEN
              U11(IX,J,K,I)=U11DM(IX,J,K,I)+1.
           ELSE IF(FLAG_DENS_P .EQ. 2) THEN
              U11(IX,J,K,I)=U11G(IX,J,K,I)+1.
           ELSE
              WRITE(*,*) 'FLAG_DENS must be 0/1/2'
              STOP
           ENDIF
           U11G(IX,J,K,I)=U11G(IX,J,K,I)+1.
           U11S(IX,J,K,I)=U11S(IX,J,K,I)+1.
        END DO
        END DO
        END DO
        END DO        

      DEALLOCATE(U11DM)

 !* special variables for paralelization
        NHYZP=NHYZ
        NHYYP=NHYY
        NHYXP=NHYX


!$OMP PARALLEL DO SHARED(NHYXP, NHYYP, NHYZP, U1DM, U1G, U1, FLAG_DENS_P), PRIVATE(IX,J,K)
        DO K=1, NHYZP
        DO J=1, NHYYP
        DO IX=1, NHYXP
           IF( FLAG_DENS_P .EQ. 0) THEN ! Total density
              U1(IX,J,K)=U1DM(IX,J,K)+1.+U1G(IX,J,K)+1.
           ELSE IF(FLAG_DENS_P .EQ. 1) THEN ! DM only
              U1(IX,J,K)=U1DM(IX,J,K)+1.
           ELSE IF(FLAG_DENS_P .EQ. 2) THEN ! gas only
              U1(IX,J,K)=U1G(IX,J,K)+1.
           ELSE
              WRITE(*,*) 'FLAG_DENS must be 0/1/2'
              STOP
           ENDIF
           U1G(IX,J,K)=U1G(IX,J,K)+1.
           U1S(IX,J,K)=U1S(IX,J,K)+1.
        ENDDO
        ENDDO
        ENDDO


       !WRITE(*,*) 'U1DM:', MINVAL(U1DM)+1., MAXVAL(U1DM)+1.
       !WRITE(*,*) 'U1G:', MINVAL(U1G)+1., MAXVAL(U1G)+1.
       !WRITE(*,*) 'U1:', MINVAL(U1), MAXVAL(U1)

       CONTA=SUM(NPART(0:NL))
       WRITE(*,*) 'TOTAL PARTICLES IN ITER=',CONTA
       
       RETURN

     END SUBROUTINE LEER


!************************************************************************* 
       SUBROUTINE LEER_NBODY(NBODYFILE, NPART_TOT, LADO, RADX0, RADY0, RADZ0, NX0, NY0, NZ0, DX0, DY0, DZ0)  !input
!************************************************************************* 
!************************************************************************* 
!       SUBROUTINE LEER_NBODY(ITER, NDXYZ0,  &!input
!                NL, T,ZETA, NPATCH, PARE, PATCHNX, PATCHNY,  &!output
!                PATCHNZ,PATCHX, PATCHY, PATCHZ,   &!output
!                PATCHRX, PATCHRY, PATCHRZ,NPART) !output
!************************************************************************* 
        USE COMMONDATA, ONLY: U1DM, UU2DM, UU3DM, UU4DM, RXPA, RYPA, RZPA, U2DM, U3DM, U4DM, FLAG_VEL
!C input variables
       INTEGER ITER!, NHYX, NHYY, NHYZ
       INTEGER:: NPART_TOT
       INTEGER:: NX0, NY0, NZ0
       REAL*4:: LADO
       REAL*4 DX0, DY0, DZ0
       REAL*4 RADX0(0:NX0+1), RADY0(0:NY0+1), RADZ0(0:NZ0+1)
       CHARACTER(LEN=*) :: NBODYFILE
!C local variables
       INTEGER:: IP
       INTEGER IR,  IRR, LOW1,LOW2
       INTEGER*4 N1, N2, N3, N4, N5, N6
       INTEGER I, IX, J, K
       INTEGER CONTA, PARTIBAS
       INTEGER:: DIM1, DIM2, DIM3, DIM4, NHYXP, NHYYP, NHYZP, FLAG_DENS_P
       REAL*4 AAA,BBB,CCC, MAP
       REAL*4 MASS, DENS_MEAN, BASS
       REAL*4, ALLOCATABLE::U4SPH(:,:,:), U1SPH(:,:,:), U2SPH(:,:,:), U3SPH(:,:,:)
       REAL*4, ALLOCATABLE::U4GSPH(:,:,:), U1GSPH(:,:,:), U2GSPH(:,:,:), U3GSPH(:,:,:)
       INTEGER,ALLOCATABLE::UBAS2(:)
       INTEGER UBAS(NX0, NY0, NZ0)
       CHARACTER*9 FILNOM1,FILNOM2, FILNOM4
       CHARACTER*10 FILNOM3
       CHARACTER*80 FIL1,FIL2, FIL4
       CHARACTER*26 FIL3 
       CHARACTER*26 FILEIN
       REAL*4 VX(-1:1),VY(-1:1),VZ(-1:1)
       REAL*4 WT(NX0, NY0, NZ0)
       INTEGER FLAG_VEL_P
       REAL*4 LENG, RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ
!C output variables
!       INTEGER NL, NDXYZ0, NST0
!       REAL*4 T,ZETA
!       INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
!       INTEGER NPART(0:NLEVELS),  NPARTST(0:NLEVELS)
!       INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
!       INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHZ(NPALEV)
!       REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
!TEMPORARY VAR
       INTEGER:: IPA, IX0, JY0, KZ0, IXR, JYR, KZR
       INTEGER:: IX01, IX02, JY01, JY02, KZ01, KZ02
       REAL*4:: UMEAN, UMEANG, UMEANDM
       REAL*4, ALLOCATABLE:: UPA(:,:,:)
       CHARACTER(LEN=3):: IL


       !FILEIN='link2reconstructed128' !link to file with halo info: X (3xfloat) V (3xfloat)
       !FILEIN=NBODYFILE

       WRITE(*,*) 'READING FILE: ', TRIM(ADJUSTL(NBODYFILE))
       OPEN(UNIT=31, FILE=NBODYFILE, STATUS='UNKNOWN',ACTION='READ')

       READ(31,*)
       DO IP=1, NPART_TOT
          IF(FLAG_VEL ==1) THEN
             READ(31,*) RXPA(IP), RYPA(IP), RZPA(IP), U2DM(IP), U3DM(IP), U4DM(IP)
          ELSE
            READ(31,*) RXPA(IP), RYPA(IP), RZPA(IP)
         ENDIF
       ENDDO
       !MASAP(:)=1. !FIXED RESOLUTION

       WRITE(*,*) 'Min/Max galaxy coordinates: ', MINVAL(rxpa), MAXVAL(rxpa)

       !RXPA=RXPA*1.E-3-0.5*LADO0 ! move to Mpc/h units
       !RYPA=RYPA*1.E-3-0.5*LADO0 ! move to Mpc/h units
       !RZPA=RZPA*1.E-3-0.5*LADO0 ! move to Mpc/h units
       RXPA=RXPA-0.5*LADO0 ! move to Mpc/h units
       RYPA=RYPA-0.5*LADO0 ! move to Mpc/h units
       RZPA=RZPA-0.5*LADO0 ! move to Mpc/h units

       RADX1=RADX0(1)-0.5*DX0
       RADXNX=RADX0(NX0)+0.5*DX0
       RADY1=RADX0(1)-0.5*DX0
       RADYNY=RADY0(NY0)+0.5*DX0
       RADZ1=RADZ0(1)-0.5*DX0
       RADZNZ=RADZ0(NZ0)+0.5*DX0
       LENG=RADXNX-RADX1


       U1DM(:,:,:)=0.
       UU2DM(:,:,:)=0.
       UU3DM(:,:,:)=0.
       UU4DM(:,:,:)=0.
       WT(:,:,:)=0.
       UBAS(:,:,:)=0.

       MASS=1. !valid for equal mass particles 

       ! TSC method to get continous density field from Nbody

       FLAG_VEL_P=FLAG_VEL

! $OMP PARALLEL DO SHARED(NPART_TOT,IR,NL,RXPA,RYPA,RZPA, &
! $OMP            U2DM, U3DM, U4DM, LENG,DX0,DY0,DZ0,RADX0,RADY0,RADZ0,NX0,NY0,NZ0,MASS, FLAG_VEL_P),&
! $OMP   PRIVATE(IP,VX,VY,VZ,BAS,I,J,K,IX,JY,KZ,I3,J3,K3),  &
! $OMP   REDUCTION(+:U1DM,UU2DM, UU3DM, UU4DM, WT, UBAS)  
      DO IP=1,NPART_TOT                                                    
       
         IF(RXPA(IP).LT.RADX1) RXPA(IP)=RXPA(IP)+LENG
         IF(RXPA(IP).GE.RADXNX) RXPA(IP)=RXPA(IP)-LENG
         IF(RYPA(IP).LT.RADY1) RYPA(IP)=RYPA(IP)+LENG
         IF(RYPA(IP).GE.RADYNY) RYPA(IP)=RYPA(IP)-LENG
         IF(RZPA(IP).LT.RADZ1) RZPA(IP)=RZPA(IP)+LENG
         IF(RZPA(IP).GE.RADZNZ) RZPA(IP)=RZPA(IP)-LENG

         VX=0.0
         VY=0.0
         VZ=0.0

         BAS=RXPA(IP)                                                    
         I=INT(((BAS-RADX0(1))/DX0)+0.4999) + 1                               
         BAS=RYPA(IP)                                                    
         J=INT(((BAS-RADY0(1))/DY0)+0.4999) + 1                               
         BAS=RZPA(IP)                                                  
         K=INT(((BAS-RADZ0(1))/DZ0)+0.4999) + 1                               
         
         BAS=ABS(RADX0(I-1)-RXPA(IP))
         VX(-1)=0.5*(1.5-BAS/DX0)**2
         BAS=ABS(RADX0(I)-RXPA(IP))
         VX(0)=0.75-(BAS/DX0)**2   
         BAS=ABS(RADX0(I+1)-RXPA(IP))
         VX(1)=0.5*(1.5-BAS/DX0)**2
               
         BAS=ABS(RADY0(J-1)-RYPA(IP))
         VY(-1)=0.5*(1.5-BAS/DY0)**2
         BAS=ABS(RADY0(J)-RYPA(IP))
         VY(0)=0.75-(BAS/DY0)**2
         BAS=ABS(RADY0(J+1)-RYPA(IP))
         VY(1)=0.5*(1.5-BAS/DY0)**2      
       
         BAS=ABS(RADZ0(K-1)-RZPA(IP))
         VZ(-1)=0.5*(1.5-BAS/DZ0)**2
         BAS=ABS(RADZ0(K)-RZPA(IP))
         VZ(0)=0.75-(BAS/DZ0)**2
         BAS=ABS(RADZ0(K+1)-RZPA(IP))
         VZ(1)=0.5*(1.5-BAS/DZ0)**2

         UBAS(I,J,K)=UBAS(I,J,K)+1

         DO KZ=-1,1
            DO JY=-1,1
               DO IX=-1,1
                  I3=I+IX
                  J3=J+JY
                  K3=K+KZ
                  IF (I3.LT.1) I3=I3+NX0
                  IF (I3.GT.NX0) I3=I3-NX0
                  IF (J3.LT.1) J3=J3+NY0
                  IF (J3.GT.NY0) J3=J3-NY0
                  IF (K3.LT.1) K3=K3+NZ0
                  IF (K3.GT.NZ0) K3=K3-NZ0

                  U1DM(I3,J3,K3)=MASS*VX(IX)*VY(JY)*VZ(KZ)+U1DM(I3,J3,K3)     
                  IF(FLAG_VEL_P == 1) THEN
                     UU2DM(I3,J3,K3)=U2DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+UU2DM(I3,J3,K3)
                     UU3DM(I3,J3,K3)=U3DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+UU3DM(I3,J3,K3)
                     UU4DM(I3,J3,K3)=U4DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+UU4DM(I3,J3,K3)
                     WT(I3,J3,K3)=WT(I3,J3,K3)+VX(IX)*VY(JY)*VZ(KZ)
                  ENDIF
               END DO
            END DO
         END DO


      END DO !loop over particles


      DENS_MEAN=MASS*NPART_TOT/(LADO**3.)
      BASS=(1./DENS_MEAN)/(DX0*DY0*DZ0)
      WRITE(*,*) 'Mean Density of the simulation:', DENS_MEAN

!$OMP PARALLEL DO  SHARED(U1DM,UU2DM, UU3DM,UU4DM, WT,BASS,NX0,NY0,NZ0, FLAG_VEL_P),PRIVATE(I,J,K, WEI)
      DO K=1,NZ0
      DO J=1,NY0
      DO I=1,NX0
       U1DM(I,J,K)=U1DM(I,J,K)*BASS - 1.0
       IF(FLAG_VEL_P == 1) THEN
          WEI=WT(I,J,K)
          IF(WT(I,J,K)==0.) WEI=1.D0
          UU2DM(I,J,K)=UU2DM(I,J,K)/WEI
          UU3DM(I,J,K)=UU3DM(I,J,K)/WEI
          UU4DM(I,J,K)=UU4DM(I,J,K)/WEI
       ENDIF
      END DO
      END DO
      END DO
      

      WRITE(*,*) 'Min-max density contrast (rho/rho_b-1):', minval(U1DM), MAXVAL(U1DM)
      WRITE(*,*) 'Min-max x-velocity:', minval(UU2DM), MAXVAL(UU2DM)

    END SUBROUTINE LEER_NBODY


!************************************************************************* 
       SUBROUTINE LEER_SPH(ITER, NDXYZ0,  &!input
                NL, T,ZETA, NPATCH, PARE, PATCHNX, PATCHNY,  &!output
                PATCHNZ,PATCHX, PATCHY, PATCHZ,   &!output
                PATCHRX, PATCHRY, PATCHRZ,NPART) !output
!************************************************************************* 
 
      USE COMMONDATA 
      IMPLICIT NONE                                                     
      
!C input variables
       INTEGER ITER!, NHYX, NHYY, NHYZ

!C local variables
       INTEGER IR,  IRR, LOW1,LOW2
       INTEGER*4 N1, N2, N3, N4, N5, N6
       INTEGER I, IX, J, K
       INTEGER CONTA, PARTIBAS
       INTEGER:: DIM1, DIM2, DIM3, DIM4, NHYXP, NHYYP, NHYZP, FLAG_DENS_P
       REAL*4 AAA,BBB,CCC, MAP
       REAL*4, ALLOCATABLE::U4SPH(:,:,:), U1SPH(:,:,:), U2SPH(:,:,:), U3SPH(:,:,:)
       REAL*4, ALLOCATABLE::U4GSPH(:,:,:), U1GSPH(:,:,:), U2GSPH(:,:,:), U3GSPH(:,:,:)
       INTEGER,ALLOCATABLE::UBAS2(:), UBAS(:)
       CHARACTER*9 FILNOM1,FILNOM2, FILNOM4
       CHARACTER*10 FILNOM3
       CHARACTER*80 FIL1,FIL2, FIL4
       CHARACTER*26 FIL3 

!C output variables
       INTEGER NL, NDXYZ0, NST0
       REAL*4 T,ZETA
       INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
       INTEGER NPART(0:NLEVELS),  NPARTST(0:NLEVELS)
       INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHZ(NPALEV)
       REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
!TEMPORARY VAR
       INTEGER:: IPA, IX0, JY0, KZ0, IXR, JYR, KZR
       INTEGER:: IX01, IX02, JY01, JY02, KZ01, KZ02
       REAL*4:: UMEAN, UMEANG, UMEANDM
       REAL*4, ALLOCATABLE:: UPA(:,:,:)
       CHARACTER(LEN=3):: IL

!*      READING DATA
       CALL NOMFILE(ITER,FILNOM1,FILNOM2,FILNOM3,FILNOM4)

       !WRITE(*,*) 'Reading iter',ITER,' ',FILNOM1,FILNOM2,FILNOM3
       !FIL1='../unigrid_dmo/snapshot_0032_N128.bin'
       fil1='../DensityField/ChallengeHaloPos_LCDMDEMNUni_N256.bin'
       OPEN (31,FILE=FIL1, &
            STATUS='UNKNOWN',ACTION='READ',FORM='UNFORMATTED')

       !READ(31) N1, N2, N3
       !N1=N2
       N1=NHYX
       N2=NHYY
       N3=NHYZ

       ALLOCATE(U1SPH(N1,N2,N3))
       ALLOCATE(U2SPH(N1,N2,N3))
       ALLOCATE(U3SPH(N1,N2,N3))
       ALLOCATE(U4SPH(N1,N2,N3))
       ALLOCATE(U1GSPH(N1,N2,N3))
       ALLOCATE(U2GSPH(N1,N2,N3))
       ALLOCATE(U3GSPH(N1,N2,N3))
       ALLOCATE(U4GSPH(N1,N2,N3))
       ALLOCATE(U1CO(N1,N2,N3))
       ALLOCATE(U2CO(N1,N2,N3))
       ALLOCATE(U3CO(N1,N2,N3))
       ALLOCATE(U4CO(N1,N2,N3))

       DIM1=N1
       DIM2=N2
       DIM3=N3
       
!$OMP PARALLEL DO SHARED(DIM1, DIM2, DIM3,  U1SPH, U2SPH, U3SPH, U4SPH, U1GSPH, U2GSPH, U3GSPH, U4GSPH),PRIVATE(I,J,K)
        DO K=1,DIM3
        DO J=1,DIM2
        DO I=1,DIM1
         U1SPH(I,J,K)=0. 
         U2SPH(I,J,K)=0.  
         U3SPH(I,J,K)=0. 
         U4SPH(I,J,K)=0.
         U1GSPH(I,J,K)=0. 
         U2GSPH(I,J,K)=0.  
         U3GSPH(I,J,K)=0. 
         U4GSPH(I,J,K)=0.
        END DO
        END DO
        END DO
       
        !N1=1
        !N2=2
        !N3=1
       READ(31) (((U1SPH(I,J,K),K=1,N1),J=1,N2),I=1,N3)
       
       WRITE(*,*) 'MIN/MAX density:', MINVAL(U1SPH), MAXVAL(U1SPH), U1SPH(16,16,16)

       IF(FLAG_VEL == 1 ) THEN
          READ(31) (((U2SPH(I,J,K),K=1,N1),J=1,N2),I=1,N3)
          WRITE(*,*) 'U[1,1,.:]', U1SPH(1,1,:)
          WRITE(*,*) U2SPH(1,1,:)
          READ(31) (((U3SPH(I,J,K),K=1,N1),J=1,N2),I=1,N3)
          READ(31) (((U4SPH(I,J,K),K=1,N1),J=1,N2),I=1,N3)
          ! read particles
      
          WRITE(*,*) 'MIN/MAX Vx:', MINVAL(U2SPH), MAXVAL(U2SPH), U2SPH(11,11,11)
          WRITE(*,*) 'MIN/MAX Vy:', MINVAL(U3SPH), MAXVAL(U3SPH), U3SPH(11,11,11)
          WRITE(*,*) 'MIN/MAX Vz:', MINVAL(U4SPH), MAXVAL(U4SPH), U4SPH(11,11,11)
       ENDIF

       FLAG_DENS_P=FLAG_DENS

!$OMP PARALLEL DO SHARED(NHYXP, NHYYP, NHYZP, U1DM, U1G, U1, FLAG_DENS_P), PRIVATE(IX,J,K)
        DO K=1, N3
        DO J=1, N2
        DO IX=1, N1
           IF( FLAG_DENS_P .EQ. 0) THEN ! Total density
              U1CO(IX,J,K)=U1SPH(IX,J,K)+U1GSPH(IX,J,K)
           ELSE IF(FLAG_DENS_P .EQ. 1) THEN ! DM only
              U1CO(IX,J,K)=U1SPH(IX,J,K)
              IF(FLAG_VEL==1) THEN
                 U2CO(IX,J,K)=U2SPH(IX,J,K)
                 U3CO(IX,J,K)=U3SPH(IX,J,K)
                 U4CO(IX,J,K)=U4SPH(IX,J,K)
              ENDIF
           ELSE IF(FLAG_DENS_P .EQ. 2) THEN ! gas only
              U1CO(IX,J,K)=U1GSPH(IX,J,K)
           ELSE
              WRITE(*,*) 'FLAG_DENS must be 0/1/2'
              STOP
           ENDIF
           !U1G(IX,J,K)=U1GSPH(IX,J,K)
        ENDDO
        ENDDO
        ENDDO
        WRITE(*,*) 'MIN/MAX U1CO:', MINVAL(U1CO), MAXVAL(U1CO)
        U2G(:,:,:)=U2CO(:,:,:)
        U3G(:,:,:)=U3CO(:,:,:)
        U4G(:,:,:)=U4CO(:,:,:)
       WRITE(*,*) 'TOTAL PARTICLES IN ITER=',CONTA
       
       RETURN

     END SUBROUTINE LEER_SPH

!******************************************************************************** 
     SUBROUTINE SMOOTH(NPART, NL)
!Smooth quantities to a coarse grid: (NCOX*NCOY*NCOZ), with NCOX > NHYX
! U1CO: overdensity
! U2CO, U3CO, U3CO: Eulerian velocities of DM
! U2GCO, U3GCO, U3GCO : gas velocity 
!********************************************************************************
     USE COMMONDATA
          !, ONLY: U1C0, RXPA, RYPA, RZPA, U2DM, U3DM, U4DM, &
          !RADX, RADY, RADZ, DX, DY, DZ,RADX0, RADY0, RADZ0, DX0, DY0, DZ0, &
          !U1CO, U2C0, U3CO, U4CO, U2DM, U3DM, U4DM
     IMPLICIT NONE
!input variables
     INTEGER NL!, LEVV
     INTEGER NPART(0:NLEVELS)
!local variables
     INTEGER I, J, K, IP, LOW1, NCOXP,NCOYP,NCOZP
     INTEGER IX, JY, KZ, I3, J3, K3, II, JJ, KK
     REAL*4 RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ, LENG
     REAL*4 VX(-1:1), VY(-1:1), VZ(-1:1)
     REAL*4 BAS, WEI, WEIV, RX, RY, RZ
     REAL*4 DT, DRX, DRY, DRZ
     !REAL*4, ALLOCATABLE:: WT(:,:,:)

     write(*,*) 'Starting smoothing...'

     !NCOX=NHYX*(2**LLEV)
     !NCOY=NHYY*(2**LLEV)
     !NCOZ=NHYZ*(2**LLEV)

     !ALLOCATE(U1CO(NCOX, NCOY, NCOZ)) !128^3, density contrast
     ALLOCATE(U2CO(NCOX, NCOY, NCOZ)) !128^3, velx_dm
     ALLOCATE(U3CO(NCOX, NCOY, NCOZ)) !128^3, vely_dm
     ALLOCATE(U4CO(NCOX, NCOY, NCOZ)) !128^3, velz_dm
     ALLOCATE(WTCO(NCOX, NCOY, NCOZ)) 
     ALLOCATE(U2GCO(NCOX, NCOY, NCOZ)) !128^3, velx_gas
     ALLOCATE(U3GCO(NCOX, NCOY, NCOZ)) !128^3, vely_gas
     ALLOCATE(U4GCO(NCOX, NCOY, NCOZ)) !128^3, velz_gas

!* special variables for paralelization
     NCOXP=NCOX
     NCOYP=NCOY
     NCOZP=NCOZ

!$OMP PARALLEL DO SHARED(NCOXP,NCOYP,NCOZP,U1CO,U2CO,U3CO, U4CO, &
!$OMP   WTCO, U1GCO, U1SCO, U2GCO, U3GCO, U4GCO),PRIVATE(I,J,K)
     DO K=1,NCOZP
     DO J=1,NCOYP
     DO I=1,NCOXP
       U1CO(I,J,K)=0.0
       U1GCO(I,J,K)=0.0
       U1SCO(I,J,K)=0.0
       U2CO(I,J,K)=0.0
       U3CO(I,J,K)=0.0
       U4CO(I,J,K)=0.0
       WTCO(I,J,K)=0.0
       U2GCO(I,J,K)=0.0
       U3GCO(I,J,K)=0.0
       U4GCO(I,J,K)=0.0
     END DO
     END DO
     END DO 

     !smallest/largest coordinates of the coarse grid
     RADX1=RADX(1)-0.5*DX 
     RADXNX=RADX(NCOX)+0.5*DX
     RADY1=RADY(1)-0.5*DY
     RADYNY=RADY(NCOY)+0.5*DY
     RADZ1=RADZ(1)-0.5*DZ
     RADZNZ=RADZ(NCOZ)+0.5*DZ
     LENG=RADXNX-RADX1


     IF(FLAG_VEL == 1) THEN

!****   from lagrangian velocities (U2DM, U3DM, U4DM)  to eulerian (U1CO, U2CO, U3CO)
        DO I=-1,1
           VX(I)=0.0
           VY(I)=0.0
           VZ(I)=0.0
        END DO

        !LOW1=SUM(NPART(0:NL))
        LOW1=NPART(0)

!$OMP PARALLEL DO SHARED(LOW1, RXPA, RYPA, RZPA, RADX1, RADY1, RADZ1, &
!$OMP     RADXNX, RADYNY, RADZNZ, LENG, DX, DY, DZ, NCOXP, NCOYP, NCOZP, &
!$OMP     U2DM, U3DM, U4DM), &
!$OMP PRIVATE(IP, BAS, I, J, K, VX, VY, VZ, IX, JY, KZ, I3, J3, K3), &
!$OMP REDUCTION(+: U2CO, U3CO, U4CO, WTCO)
        DO IP=1, LOW1
 
           !move by LENG particles outside the grid
           IF(RXPA(IP).LT.RADX1) RXPA(IP)=RXPA(IP)+LENG
           IF(RXPA(IP).GE.RADXNX) RXPA(IP)=RXPA(IP)-LENG
           IF(RYPA(IP).LT.RADY1) RYPA(IP)=RYPA(IP)+LENG
           IF(RYPA(IP).GE.RADYNY) RYPA(IP)=RYPA(IP)-LENG
           IF(RZPA(IP).LT.RADZ1) RZPA(IP)=RZPA(IP)+LENG
           IF(RZPA(IP).GE.RADZNZ) RZPA(IP)=RZPA(IP)-LENG

           !initial position
           !I,J,K: cell of the particle in the coarse grid
           BAS=RXPA(IP)                                                    
           I=INT(((BAS-RADX(1))/DX)+0.49999) + 1                               
           BAS=RYPA(IP)                                                    
           J=INT(((BAS-RADY(1))/DY)+0.49999) + 1                               
           BAS=RZPA(IP)                                                    
           K=INT(((BAS-RADZ(1))/DZ)+0.49999) + 1   


           !***     using a cell in cloud scheme (TSC)  to transform lagrangian velocities (DM part) into eulerian
           BAS=ABS(RADX(I-1)-RXPA(IP))
           VX(-1)=0.5*(1.5-BAS/DX)**2
           BAS=ABS(RADX(I)-RXPA(IP))
           VX(0)=0.75-(BAS/DX)**2
           BAS=ABS(RADX(I+1)-RXPA(IP))
           VX(1)=0.5*(1.5-BAS/DX)**2

           BAS=ABS(RADY(J-1)-RYPA(IP))
           VY(-1)=0.5*(1.5-BAS/DY)**2
           BAS=ABS(RADY(J)-RYPA(IP))
           VY(0)=0.75-(BAS/DY)**2
           BAS=ABS(RADY(J+1)-RYPA(IP))
           VY(1)=0.5*(1.5-BAS/DY)**2
        
           BAS=ABS(RADZ(K-1)-RZPA(IP))
           VZ(-1)=0.5*(1.5-BAS/DZ)**2
           BAS=ABS(RADZ(K)-RZPA(IP))
           VZ(0)=0.75-(BAS/DZ)**2
           BAS=ABS(RADZ(K+1)-RZPA(IP))
           VZ(1)=0.5*(1.5-BAS/DZ)**2

        
           DO KZ=-1,1
              DO JY=-1,1
                 DO IX=-1,1

                    I3=I+IX
                    J3=J+JY
                    K3=K+KZ
                    IF (I3.LT.1) I3=I3+NCOXP
                    IF (I3.GT.NCOX) I3=I3-NCOXP
                    IF (J3.LT.1) J3=J3+NCOYP
                    IF (J3.GT.NCOY) J3=J3-NCOYP
                    IF (K3.LT.1) K3=K3+NCOZP
                    IF (K3.GT.NCOZ) K3=K3-NCOZP

                    U2CO(I3,J3,K3)=U2DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U2CO(I3,J3,K3)
                    U3CO(I3,J3,K3)=U3DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U3CO(I3,J3,K3)  
                    U4CO(I3,J3,K3)=U4DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U4CO(I3,J3,K3)
                    WTCO(I3,J3,K3)=WTCO(I3,J3,K3)+VX(IX)*VY(JY)*VZ(KZ)

                 END DO
              END DO
           END DO

        ENDDO !loop on particles


        WRITE(*,*) 'Cell without particles:', COUNT(WTCO .EQ. 0.) 

        DO K=1,NCOZ
           DO J=1,NCOY
              DO I=1,NCOX
                 WEI=WTCO(I,J,K)
                 IF(WTCO(I,J,K)==0.) WEI=1.D0
                 U2CO(I,J,K)=U2CO(I,J,K)/WEI
                 U3CO(I,J,K)=U3CO(I,J,K)/WEI
                 U4CO(I,J,K)=U4CO(I,J,K)/WEI
              END DO
           END DO
        END DO
        !DEALLOCATE(WT)
     ENDIF

! now smooth the density field: U1 (DM+GAS)  --> U1CO

      WEI=(DX0*DY0*DZ0)/(DX*DY*DZ) !dilution factor
      WEIV=2.**3.

      DO KK=1, NHYZ !level 0 grid
      DO JJ=1, NHYY
      DO II=1, NHYX

         RX=RADX0(II)
         RY=RADY0(JJ)
         RZ=RADZ0(KK)
         
     !move by LENG particles outside the grid, check if happens!
        IF(RX.LT.RADX1) WRITE(*,*) 'WARNING: RX<RADX1 in SMOOTH'
        IF(RX.GE.RADXNX) WRITE(*,*) 'WARNING: RX>=RADXNX in SMOOTH'
        IF(RY.LT.RADY1) WRITE(*,*) 'WARNING: RY<RADY1 in SMOOTH'
        IF(RY.GE.RADYNY) WRITE(*,*) 'WARNING: RY>=RADYNY in SMOOTH'
        IF(RZ.LT.RADZ1) WRITE(*,*) 'WARNING: RZ<RADZ1 in SMOOTH'
        IF(RZ.GE.RADZNZ) WRITE(*,*)'WARNING: RZ>RADZNZ in SMOOTH'

        I=INT(((RX-RADX(1))/DX)+0.5) + 1 !coarse (NXCO^3) grid  
        J=INT(((RY-RADY(1))/DY)+0.5) + 1
        K=INT(((RZ-RADZ(1))/DZ)+0.5) + 1  
 
        U1CO(I,J,K)=U1CO(I,J,K)+U1(II,JJ,KK)*WEI
        U1GCO(I,J,K)=U1GCO(I,J,K)+U1G(II,JJ,KK)*WEI
        U1SCO(I,J,K)=U1SCO(I,J,K)+U1S(II,JJ,KK)*WEI
        IF(FLAG_VEL ==1 ) THEN
           U2GCO(I,J,K)=U2GCO(I,J,K)+U2G(II,JJ,KK) !smooth gas velocity field
           U3GCO(I,J,K)=U3GCO(I,J,K)+U3G(II,JJ,KK)
           U4GCO(I,J,K)=U4GCO(I,J,K)+U4G(II,JJ,KK)
        ENDIF
     ENDDO
     ENDDO
     ENDDO
        
     IF(FLAG_VEL ==1 ) THEN
        DO I=1, NCOX
           DO J=1,NCOY
              DO K=1, NCOZ
                 U2GCO(I,J,K)=U2GCO(I,J,K)/WEIV
                 U3GCO(I,J,K)=U3GCO(I,J,K)/WEIV
                 U4GCO(I,J,K)=U4GCO(I,J,K)/WEIV
              ENDDO
           ENDDO
        ENDDO
     ENDIF

     WRITE(*,*) 'MEAN DENSITY, LEVEL -1:', sum(U1CO)/(NCOX*NCOY*NCOZ)

     DEALLOCATE(WTCO)
     WRITE(*,*) 'End of smoothing'

     END SUBROUTINE SMOOTH

!******************************************************************************** 
     SUBROUTINE EULERIAN(NL, NDXYZ0, NPATCH, PARE, PATCHNX, PATCHNY, & 
                PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
                PATCHRX, PATCHRY, PATCHRZ,NPART)
!********************************************************************************
! Compute Eulerian velocities for the DM for levels >= 0 (i.e., levels not done in SMOOTH)
     USE COMMONDATA !ONLY: RXPA, RYPA, RZPA, U2D, U3DM, U4DM, U1, U2, U3, &
     !U12, U13, U14, RADX, RADY, RADZ, RADX0, RADY0, RADZ0, DX, DY, DZ, DX0, DY0, DZ0, &
     !NHYX, NHYY, NHYZ
     IMPLICIT NONE
!input variables
     INTEGER NL, NDXYZ0
     INTEGER NPART(0:NLEVELS)
     INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
     INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
     INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHz(NPALEV)
     REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
!local variables
     INTEGER I, J, K, IP, LOW1, LOW2, N1, N2, N3, NPNOPA
     INTEGER IX, JY, KZ, I3, J3, K3, II, JJ, KK
     REAL*4 RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ, LENG
     REAL*4 VX(-1:1), VY(-1:1), VZ(-1:1)
     REAL*4 BAS, WEI, RX, RY, RZ, PX, PY, PZ
     REAL*4 DT, DXR, DYR, DZR
     INTEGER DIM1, DIM2, DIM3, DIM4, IPA, IPA0, IPA1, IPA10, IPA00, IR
     !REAL*4, ALLOCATABLE:: WT0(:,:,:)
     REAL*4, ALLOCATABLE:: WT(:,:,:,:)
     REAL*4, ALLOCATABLE:: U120(:,:,:,:), U130(:,:,:,:), U140(:,:,:,:)

!*****   LEVEL 0 *****
     write(*,*) 'Eulerian: level 0'

     N1=NHYX
     N2=NHYY
     N3=NHYZ
     ALLOCATE(U2(N1, N2, N3)) 
     ALLOCATE(U3(N1, N2, N3)) 
     ALLOCATE(U4(N1, N2, N3)) 
     ALLOCATE(WT0(N1, N2, N3)) 


!$OMP PARALLEL DO SHARED(N1,N2,N3,U2,U3, U4, WT0),PRIVATE(I,J,K)
     DO K=1,N3
     DO J=1,N2
     DO I=1,N1
       U2(I,J,K)=0.0
       U3(I,J,K)=0.0
       U4(I,J,K)=0.0
       WT0(I,J,K)=0.0
     END DO
     END DO
     END DO 

     !smallest/largest coordinates of the coarse grid
     RADX1=RADX0(1)-0.5*DX0
     RADXNX=RADX0(N1)+0.5*DX0
     RADY1=RADY0(1)-0.5*DY0
     RADYNY=RADY0(N2)+0.5*DY0
     RADZ1=RADZ0(1)-0.5*DZ0
     RADZNZ=RADZ0(N3)+0.5*DZ0
     LENG=RADXNX-RADX1

!****   from lagrangian velocities (U2DM, U3DM, U4DM)  to eulerian (U1CO, U2CO, U3CO)
      DO I=-1,1
        VX(I)=0.0
        VY(I)=0.0
        VZ(I)=0.0
      END DO  

     LOW1=NDXYZ0
     !LOW1=SUM(NPART(0:NL))
!$OMP PARALLEL DO SHARED(LOW1, RXPA, RYPA, RZPA, RADX1, RADY1, RADZ1, &
!$OMP     RADXNX, RADYNY, RADZNZ, LENG, DX0, DY0, DZ0, N1, N2, N3, &
!$OMP     U2DM, U3DM, U4DM, RADX0, RADY0, RADZ0), &
!$OMP PRIVATE(IP, BAS, I, J, K, VX, VY, VZ, IX, JY, KZ, I3, J3, K3), &
!$OMP REDUCTION(+: U2, U3, U4, WT0)
     DO IP=1, LOW1
 
     !move by LENG particles outside the grid
        IF(RXPA(IP).LT.RADX1) RXPA(IP)=RXPA(IP)+LENG
        IF(RXPA(IP).GE.RADXNX) RXPA(IP)=RXPA(IP)-LENG
        IF(RYPA(IP).LT.RADY1) RYPA(IP)=RYPA(IP)+LENG
        IF(RYPA(IP).GE.RADYNY) RYPA(IP)=RYPA(IP)-LENG
        IF(RZPA(IP).LT.RADZ1) RZPA(IP)=RZPA(IP)+LENG
        IF(RZPA(IP).GE.RADZNZ) RZPA(IP)=RZPA(IP)-LENG

        !initial position
        !I,J,K: cell of the particle in the coarse grid
        BAS=RXPA(IP)                                                    
        I=INT(((BAS-RADX0(1))/DX0)+0.49999) + 1                               
        BAS=RYPA(IP)                                                    
        J=INT(((BAS-RADY0(1))/DY0)+0.49999) + 1                               
        BAS=RZPA(IP)                                                    
        K=INT(((BAS-RADZ0(1))/DZ0)+0.49999) + 1   

!***     using a cell in cloud scheme (TSC)  to transform lagrangian velocities (DM part) into eulerian
        BAS=ABS(RADX0(I-1)-RXPA(IP))
        VX(-1)=0.5*(1.5-BAS/DX0)**2
        BAS=ABS(RADX0(I)-RXPA(IP))
        VX(0)=0.75-(BAS/DX0)**2
        BAS=ABS(RADX0(I+1)-RXPA(IP))
        VX(1)=0.5*(1.5-BAS/DX0)**2

        BAS=ABS(RADY0(J-1)-RYPA(IP))
        VY(-1)=0.5*(1.5-BAS/DY0)**2
        BAS=ABS(RADY0(J)-RYPA(IP))
        VY(0)=0.75-(BAS/DY0)**2
        BAS=ABS(RADY0(J+1)-RYPA(IP))
        VY(1)=0.5*(1.5-BAS/DY0)**2
        
        BAS=ABS(RADZ0(K-1)-RZPA(IP))
        VZ(-1)=0.5*(1.5-BAS/DZ0)**2
        BAS=ABS(RADZ0(K)-RZPA(IP))
        VZ(0)=0.75-(BAS/DZ0)**2
        BAS=ABS(RADZ0(K+1)-RZPA(IP))
        VZ(1)=0.5*(1.5-BAS/DZ0)**2

       DO KZ=-1,1
        DO JY=-1,1
        DO IX=-1,1
         I3=I+IX
         J3=J+JY
         K3=K+KZ
         IF (I3.LT.1) I3=I3+N1
         IF (I3.GT.N1) I3=I3-N1
         IF (J3.LT.1) J3=J3+N2
         IF (J3.GT.N2) J3=J3-N2
         IF (K3.LT.1) K3=K3+N3
         IF (K3.GT.N3) K3=K3-N3

         U2(I3,J3,K3)=U2DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U2(I3,J3,K3)
         U3(I3,J3,K3)=U3DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U3(I3,J3,K3)  
         U4(I3,J3,K3)=U4DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U4(I3,J3,K3)
         WT0(I3,J3,K3)=WT0(I3,J3,K3)+VX(IX)*VY(JY)*VZ(KZ)

        END DO
        END DO
     END DO

  ENDDO !loop on particles
 
    DO K=1,N3
     DO J=1,N2
     DO I=1,N1
       WEI=WT0(I,J,K)
       IF(WT0(I,J,K)==0.) WEI=1.D0
       U2(I,J,K)=U2(I,J,K)/WEI
       U3(I,J,K)=U3(I,J,K)/WEI
       U4(I,J,K)=U4(I,J,K)/WEI
      END DO
     END DO
     END DO
     !DEALLOCATE(WT0)
       write(*,*) 'IR=0, MINMAX(U2)',  MINVAL(U2(:,:,:)), MAXVAL(U2(:,:,:))
       write(*,*) 'IR=0, MINMAX(U3)',  MINVAL(U3(:,:,:)), MAXVAL(U3(:,:,:))
       write(*,*) 'IR=0, MINMAX(U3)',  MINVAL(U4(:,:,:)), MAXVAL(U4(:,:,:))
   

     write(*,*) 'Eulerian: levels>0'



!LEVEL IR>0 (from eulerian_par)
     DIM1=MAXVAL(PATCHNX)
     DIM2=MAXVAL(PATCHNY)
     DIM3=MAXVAL(PATCHNZ)
     DIM4=SUM(NPATCH(0:NL))
     ALLOCATE(U12(DIM1, DIM2, DIM3, DIM4))
     ALLOCATE(U13(DIM1, DIM2, DIM3, DIM4))
     ALLOCATE(U14(DIM1, DIM2, DIM3, DIM4))

     DO IR=1, 0
     DIM4=NPATCH(IR)
     !write(*,*) 'IR, DIM4', IR, DIM4

     ALLOCATE(U120(DIM1, DIM2, DIM3, DIM4)) !temporary arrays, used to save memory in the parallel loop
     ALLOCATE(U130(DIM1, DIM2, DIM3, DIM4))
     ALLOCATE(U140(DIM1, DIM2, DIM3, DIM4))
     ALLOCATE(WT(DIM1, DIM2, DIM3, DIM4))

       LOW1=1
       LOW2=NPATCH(IR)
       DO IPA=LOW1, LOW2
          DO I=1,PATCHNX(IPA)
             DO J=1,PATCHNY(IPA)
                DO K=1, PATCHNZ(IPA)
                   U120(I,J,K,IPA)=0.0
                   U130(I,J,K,IPA)=0.0
                   U140(I,J,K,IPA)=0.0
                   WT(I,J,K,IPA)=0.0
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       DXR=DX0/(2.**IR)
       DYR=DY0/(2.**IR)
       DZR=DZ0/(2.**IR)

       RADX1=RADX0(1)-0.5*DX0 
       RADXNX=RADX0(NHYX)+0.5*DX0
       RADY1=RADY0(1)-0.5*DY0
       RADYNY=RADY0(NHYY)+0.5*DY0
       RADZ1=RADZ0(1)-0.5*DZ0
       RADZNZ=RADZ0(NHYZ)+0.5*DZ0
       LENG=RADXNX-RADX1

       LOW1=SUM(NPART(0:IR-1))+1
       LOW2=SUM(NPART(0:IR))

       NPNOPA=0
 
!$OMP PARALLEL DO SHARED(LOW1, LOW2, RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ, & 
!$OMP     RXPA, RYPA, RZPA, LENG, NPATCH, PATCHRX, PATCHRY, PATCHRZ, DXR, DYR, DZR, &
!$OMP     U2DM, U3DM, U4DM), &
!$OMP PRIVATE(IP, IPA0, I, J, K, IPA, BAS, PX, PY, PZ, VX, VY, VZ, IX, JY, KZ, & 
!$OMP     I3, J3, K3, RX, RY, RZ, II, JJ, KK, IPA1, IPA10, IPA00), &
!$OMP REDUCTION(+: NPNOPA, U120, U130, U140, WT), &
!$OMP IF(NPATCH(IR) .LT. 2000) !otherwise too much memory is required
      DO IP=LOW1, LOW2 !only part in level IR

        IF(RXPA(IP).LT.RADX1) RXPA(IP)=RXPA(IP)+LENG
        IF(RXPA(IP).GE.RADXNX) RXPA(IP)=RXPA(IP)-LENG
        IF(RYPA(IP).LT.RADY1) RYPA(IP)=RYPA(IP)+LENG
        IF(RYPA(IP).GE.RADYNY) RYPA(IP)=RYPA(IP)-LENG
        IF(RZPA(IP).LT.RADZ1) RZPA(IP)=RZPA(IP)+LENG
        IF(RZPA(IP).GE.RADZNZ) RZPA(IP)=RZPA(IP)-LENG

          IPA0=0
          I=0 ! cell in the patch with the particle
          J=0
          K=0

          OUT: DO IPA1=1, NPATCH(IR) !IPA1 USED ONLY IN U120, U130, U140
             IPA=IPA1+SUM(NPATCH(0:IR-1))

             BAS=(RXPA(IP)-(PATCHRX(IPA)-0.5*DXR))/DXR
             IF(BAS .GE. -0.5) &
                  I=INT(BAS+0.49999)+1 
             IF(BAS .LT. -0.5) I=0

             BAS=(RYPA(IP)-(PATCHRY(IPA)-0.5*DYR))/DYR
             IF(BAS .GE. -0.5) &
                  J=INT(BAS+0.49999)+1
             IF(BAS .LT. -0.5) J=0

             BAS=(RZPA(IP)-(PATCHRZ(IPA)-0.5*DZR))/DZR 
             IF(BAS .GE. -0.5) &
                  K=INT(BAS+0.49999)+1
             IF(BAS .LT. -0.5) K=0

             IF(I.GE.1 .AND. I.LE. PATCHNX(IPA) .AND.  &
                  J.GE.1 .AND. J.LE. PATCHNY(IPA) .AND. &
                  K.GE.1 .AND. K.LE. PATCHNZ(IPA)) THEN
                IPA0=IPA
 
                EXIT OUT
             ENDIF
          ENDDO OUT
  
          IPA10=IPA0-SUM(NPATCH(0:IR-1))

          IF(IPA0==0) THEN
             
             NPNOPA=NPNOPA+1

          ELSE

          PX=PATCHRX(IPA0)-0.5*DXR+(I-1)*DXR !center of the cell (i,j,k, ipa0)
          PY=PATCHRY(IPA0)-0.5*DYR+(J-1)*DYR
          PZ=PATCHRZ(IPA0)-0.5*DZR+(K-1)*DZR

          BAS=ABS(PX-DXR-RXPA(IP))         
          VX(-1)=0.5*(1.5-BAS/DXR)**2
          BAS=ABS(PX-RXPA(IP))
          VX(0)=0.75-(BAS/DXR)**2
          IF(VX(0)<0) WRITE(*,*) VX(0), BAS, PX, RXPA(IP), DXR
          BAS=ABS(PX+DXR-RXPA(IP))
          VX(1)=0.5*(1.5-BAS/DXR)**2

          BAS=ABS(PY-DYR-RYPA(IP))
          VY(-1)=0.5*(1.5-BAS/DYR)**2
          BAS=ABS(PY-RYPA(IP))
          VY(0)=0.75-(BAS/DYR)**2
          BAS=ABS(PY+DYR-RYPA(IP))
          VY(1)=0.5*(1.5-BAS/DYR)**2

          BAS=ABS(PZ-DZR-RZPA(IP))
          VZ(-1)=0.5*(1.5-BAS/DZR)**2
          BAS=ABS(PZ-RZPA(IP))
          VZ(0)=0.75-(BAS/DZR)**2
          BAS=ABS(PZ+DZR-RZPA(IP))
          VZ(1)=0.5*(1.5-BAS/DZR)**2

        DO KZ=-1,1
        DO JY=-1,1
        DO IX=-1,1
         I3=I+IX
         J3=J+JY
         K3=K+KZ

         IF(I3.GE.1 .AND. I3.LE. PATCHNX(IPA0) .AND.  & 
            J3.GE.1 .AND. J3.LE. PATCHNY(IPA0) .AND. &
            K3.GE.1 .AND. K3.LE. PATCHNZ(IPA0)) THEN
            !neighbour cell to whom i'm assigning the weights is still in the patch

            U120(I3,J3,K3, IPA10)=U2DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U120(I3,J3,K3, IPA10)
            U130(I3,J3,K3, IPA10)=U3DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U130(I3,J3,K3, IPA10)
            U140(I3,J3,K3, IPA10)=U4DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U140(I3,J3,K3, IPA10)
            WT(I3,J3,K3, IPA10)=WT(I3,J3,K3, IPA10)+VX(IX)*VY(JY)*VZ(KZ)

            ELSE !find patch with cell i3, j3, k3

            RX=PATCHRX(IPA0)-0.5*DXR+(I3-1)*DXR !center of the cell i3, j3, k3 in physical coord
            RY=PATCHRY(IPA0)-0.5*DYR+(J3-1)*DYR
            RZ=PATCHRZ(IPA0)-0.5*DZR+(K3-1)*DZR

            !OUT2: DO IPA=SUM(NPATCH(0:IR-1))+1, SUM(NPATCH(0:IR))
            OUT2: DO IPA1=1, NPATCH(IR) !IPA1 USED ONLY IN U120, U130, U140
             IPA=IPA1+SUM(NPATCH(0:IR-1))

             BAS=(RX-(PATCHRX(IPA)-0.5*DXR))/DXR
             IF(BAS .GE. 0) II=INT(BAS)+1
             IF(BAS .LT. 0) II=0
             BAS=(RY-(PATCHRY(IPA)-0.5*DYR))/DYR
             IF(BAS .GE. 0) JJ=INT(BAS)+1
             IF(BAS .LT. 0) JJ=0
             BAS=(RZ-(PATCHRZ(IPA)-0.5*DZR))/DZR
             IF(BAS .GE. 0) KK=INT(BAS)+1
             IF(BAS .LT. 0) K=0
             IF(II.GE.1 .AND. II.LE. PATCHNX(IPA) .AND.  &
                  JJ.GE.1 .AND. JJ.LE. PATCHNY(IPA) .AND. &
                  KK.GE.1 .AND. KK.LE. PATCHNZ(IPA)) THEN
                IPA00=IPA1
                U120(II,JJ,KK, IPA00)=U2DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U120(II,JJ,KK, IPA00)
                U130(II,JJ,KK, IPA00)=U3DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U130(II,JJ,KK, IPA00)
                U140(II,JJ,KK, IPA00)=U4DM(IP)*VX(IX)*VY(JY)*VZ(KZ)+U140(II,JJ,KK, IPA00)
                WT(II,JJ,KK, IPA00)=WT(II,JJ,KK, IPA00)+VX(IX)*VY(JY)*VZ(KZ)

                EXIT OUT2
             ENDIF

          ENDDO OUT2
       ENDIF

       ENDDO
       ENDDO
       ENDDO

    ENDIF!IF IPA0=0
    ENDDO !LOOP ON THE PARTICLES, IP - END PARALLEL LOOP


    !ENDDO !LOOP ON LEVELS

   !NORMALIZE VEL AND STORE GLOBAL ARRAYS: U12, U13, U14
   !DIM4=SUM(NPATCH(0:NL))
   !ALLOCATE(U12(DIM1, DIM2, DIM3, DIM4))
   !ALLOCATE(U13(DIM1, DIM2, DIM3, DIM4))
   !ALLOCATE(U14(DIM1, DIM2, DIM3, DIM4))
   !write(*,*) 'allocate ok', dim4

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO IPA=LOW1, LOW2
          DO I=1,PATCHNX(IPA)
             DO J=1,PATCHNY(IPA)
                DO K=1, PATCHNZ(IPA)
                   U12(I,J,K,IPA)=0.0
                   U13(I,J,K,IPA)=0.0
                   U14(I,J,K,IPA)=0.0
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       DO IPA1=1, NPATCH(IR) !
          IPA=IPA1+SUM(NPATCH(0:IR-1))
          DO I=1,PATCHNX(IPA)
             DO J=1,PATCHNY(IPA)
                DO K=1, PATCHNZ(IPA)
                   IF(WT(I,J,K,IPA1)==0.) WT(I,J,K,IPA1)=1.D0
                   U12(I,J,K,IPA)=U120(I,J,K,IPA1)/WT(I,J,K,IPA1)
                   U13(I,J,K,IPA)=U130(I,J,K,IPA1)/WT(I,J,K,IPA1)
                   U14(I,J,K,IPA)=U140(I,J,K,IPA1)/WT(I,J,K,IPA1)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !write(*,*) 'IR, MINMAX(U12)',   IR, MINVAL(U12(:,:,:,LOW1:LOW2)), MAXVAL(U12(:,:,:,LOW1:LOW2))
       !write(*,*) 'IR, MINMAX(U13)',   IR, MINVAL(U13(:,:,:,LOW1:LOW2)), MAXVAL(U13(:,:,:,LOW1:LOW2))
       !write(*,*) 'IR, MINMAX(U13)',   IR, MINVAL(U14(:,:,:,LOW1:LOW2)), MAXVAL(U14(:,:,:,LOW1:LOW2))
   
    DEALLOCATE(U120, U130, U140)
    DEALLOCATE(WT)

    ENDDO!LOOP ON LEVELS

    DEALLOCATE(WT0)
    WRITE(*,*)'END of levels'


   END SUBROUTINE EULERIAN


!******************************************************************************** 
     SUBROUTINE MARKR(IR,NX0, NY0, NZ0, DENS_THRE) !update for IR>0
!********************************************************************************
     USE COMMONDATA, ONLY: U1CO, U2CO, U3CO, U4CO, DX,DY,DZ, DIVERCO, FLAGV, &
          U2GCO, U3GCO, U4GCO, DIVERGCO
     IMPLICIT NONE
!input variables
     INTEGER NCOX,NCOY,NCOZ, NX0, NY0, NZ0, IR
     REAL*4 DENS_THRE
!local variables
     INTEGER NX, NY, NZ, IX, JY, KZ
     REAL*4 BAS21, BAS32, BAS43

     NX=NX0*(2.**IR) !IR=1 --> 2*NHYX=512; IR=-1 --> NHYX/2=128
     NY=NY0*(2.**IR)
     NZ=NZ0*(2.**IR)
     ALLOCATE(DIVERCO(NX,NY,NZ))
     ALLOCATE(DIVERGCO(NX,NY,NZ))
     ALLOCATE(FLAGV(NX,NY,NZ))
     DO KZ=1,NZ 
     DO JY=1,NY
     DO IX=1,NX
        DIVERCO(IX,JY,KZ)=0.0 !IF IR>=0 i need to call routine eulerian
        DIVERGCO(IX,JY,KZ)=0.0
        FLAGV(IX,JY,KZ)=0
     ENDDO
     ENDDO
     ENDDO


   END SUBROUTINE MARKR


!******************************************************************************** 
     SUBROUTINE MARK_SUB(NX,NY,NZ, DENS_THRE) 
!********************************************************************************
     USE COMMONDATA, ONLY: U1CO, U2CO, U3CO, U4CO, DX,DY,DZ, DIVERCO, FLAGV, &
          U2GCO, U3GCO, U4GCO, DIVERGCO, FLAG_SUB, FLAGAMR, MARCACO, FLAGV
     IMPLICIT NONE
!input variables
     INTEGER NX,NY,NZ
     REAL*4 DENS_THRE
!local variables
     INTEGER IX, JY, KZ

     ALLOCATE(FLAGV(NX,NY,NZ))
     ALLOCATE(FLAG_SUB(NX,NY,NZ))
     FLAGV(:,:,:)=0
     FLAG_SUB(:,:,:)=0

     DO KZ=1, NZ
        DO JY=1, NY
           DO IX=1, NX

              IF(MARCACO(IX,JY,KZ) .GT. 0 .AND. FLAGAMR(IX,JY,KZ) .EQ. 1) THEN !only refined cells !!
              !IF(MARCACO(IX,JY,KZ) .GT. 0 ) THEN !all cells in voids 
                   FLAG_SUB(IX,JY,KZ)=1 !can be part of a subvoid
                   IF(DIVERGCO(IX,JY,KZ) .GT. 0  .AND. U1CO(IX,JY,KZ) .LT. DENS_THRE+1.) FLAGV(IX,JY,KZ)=1
                   !can be center of a subvoid
                ENDIF
           ENDDO
        ENDDO
     ENDDO

   END SUBROUTINE MARK_SUB

!******************************************************************************** 
     SUBROUTINE MARK_ALL(NX,NY,NZ, DENS_THRE) !valid for all hierarchies
!********************************************************************************
     USE COMMONDATA, ONLY: U1CO, U2CO, U3CO, U4CO, DX,DY,DZ, DIVERCO, FLAGV, &
          U2GCO, U3GCO, U4GCO, DIVERGCO, FLAG_SUB, FLAGAMR, MARCAP, FLAGV
     IMPLICIT NONE
!input variables
     INTEGER NX,NY,NZ
     REAL*4 DENS_THRE
!local variables
     INTEGER IX, JY, KZ

     ALLOCATE(FLAGV(NX,NY,NZ)) !for all hier
     ALLOCATE(FLAG_SUB(NX,NY,NZ)) !for subvoids
     FLAGV(:,:,:)=0
     FLAG_SUB(:,:,:)=0

     DO KZ=1, NZ
        DO JY=1, NY
           DO IX=1, NX
! @lev=LEVMIN MARCAP=-1 in every cell
!             FLAGAMR=1 in every cell 
! @lev>LEVMIN MARCA>0 for void cells, =0 for non-void cells
!   && lev>0  FLAGAMR=1 only in refined cells
              !IF(MARCAP(IX,JY,KZ) .NE. 0 .AND. FLAGAMR(IX,JY,KZ) .EQ. 1) THEN !only refined cells !!
              IF(MARCAP(IX,JY,KZ) .NE. 0 ) THEN !all cells in voids 
                   FLAG_SUB(IX,JY,KZ)=1 !can be part of a subvoid
                   IF(DIVERGCO(IX,JY,KZ) .GT. 0  .AND. U1CO(IX,JY,KZ) .LT. DENS_THRE+1.) FLAGV(IX,JY,KZ)=1
                   !can be center of a subvoid
                ENDIF
           ENDDO
        ENDDO
     ENDDO

!--> if lev=LEVMIN: FLAG_SUB=1 in every cell

   END SUBROUTINE MARK_ALL

!******************************************************************************** 
     SUBROUTINE MARK_ALL_DENS(NX,NY,NZ, DENS_THRE) !valid for all hierarchies
!********************************************************************************
     USE COMMONDATA, ONLY: U1CO, U2CO, U3CO, U4CO, DX,DY,DZ, DIVERCO, FLAGV, &
          U2GCO, U3GCO, U4GCO, DIVERGCO, FLAG_SUB, FLAGAMR, MARCAP, FLAGV
     IMPLICIT NONE
!input variables
     INTEGER NX,NY,NZ
     REAL*4 DENS_THRE
!local variables
     INTEGER IX, JY, KZ

     ALLOCATE(FLAGV(NX,NY,NZ)) !for all hier
     ALLOCATE(FLAG_SUB(NX,NY,NZ)) !for subvoids
     FLAGV(:,:,:)=0
     FLAG_SUB(:,:,:)=0
     WRITE(*,*) 'In MARK_ALL_DENS:',minval(u1co), maxval(u1co),DENS_THRE+1.
     DO KZ=1, NZ
        DO JY=1, NY
           DO IX=1, NX
! @lev=LEVMIN MARCAP=-1 in every cell
!             FLAGAMR=1 in every cell 
! @lev>LEVMIN MARCA>0 for void cells, =0 for non-void cells
!   && lev>0  FLAGAMR=1 only in refined cells
              !IF(MARCAP(IX,JY,KZ) .NE. 0 .AND. FLAGAMR(IX,JY,KZ) .EQ. 1) THEN !only refined cells !!
              IF(MARCAP(IX,JY,KZ) .NE. 0 ) THEN !all cells in voids 
                   FLAG_SUB(IX,JY,KZ)=1 !can be part of a subvoid
                   IF(U1CO(IX,JY,KZ) .LT. DENS_THRE+1.) FLAGV(IX,JY,KZ)=1
                ENDIF
           ENDDO
        ENDDO
     ENDDO

!--> if lev=LEVMIN: FLAG_SUB=1 in every cell

   END SUBROUTINE MARK_ALL_DENS

!******************************************************************************** 
     SUBROUTINE MARK(NCOX,NCOY,NCOZ,DENS_THRE) !update for IR>0
!********************************************************************************
     USE COMMONDATA, ONLY: U1CO, U2CO, U3CO, U4CO, DX,DY,DZ, DIVERCO, FLAGV, &
          U2GCO, U3GCO, U4GCO, DIVERGCO, FLAG_SUB
     IMPLICIT NONE
!input variables
     INTEGER NCOX,NCOY,NCOZ
     REAL*4 DENS_THRE
!local variables
     INTEGER NX, NY, NZ, IX, JY, KZ
     REAL*4 BAS21, BAS32, BAS43

     NX=NCOX
     NY=NCOY
     NZ=NCOZ
 
     !ALLOCATE(DIVERCO(NX,NY,NZ))
     !ALLOCATE(DIVERGCO(NX,NY,NZ))
     ALLOCATE(FLAGV(NX,NY,NZ))
     ALLOCATE(FLAG_SUB(NX,NY,NZ))
     DO KZ=1,NZ 
     DO JY=1,NY
     DO IX=1,NX
        !DIVERCO(IX,JY,KZ)=0.0
        !DIVERGCO(IX,JY,KZ)=0.0
        FLAGV(IX,JY,KZ)=0
        FLAG_SUB(IX,JY,KZ)=1
     ENDDO
     ENDDO
     ENDDO

   !  !divergence only for cells 2,N-1, not  the border
   !  DO KZ=2,NZ-1 
   !  DO JY=2,NY-1
   !  DO IX=2,NX-1

   !  !velocity divergence for DM 
   !    BAS21=U2CO(IX+1,JY,KZ)-U2CO(IX-1,JY,KZ)
   !    BAS21=BAS21/(2.0*DX)       
       
   !    BAS32=U3CO(IX,JY+1,KZ)-U3CO(IX,JY-1,KZ)
   !    BAS32=BAS32/(2.0*DY)
       
   !    BAS43=U4CO(IX,JY,KZ+1)-U4CO(IX,JY,KZ-1)
   !    BAS43=BAS43/(2.0*DZ)

   !    DIVERCO(IX,JY,KZ)=BAS21+BAS32+BAS43

     !velocity divergence for gas
   !    BAS21=U2GCO(IX+1,JY,KZ)-U2GCO(IX-1,JY,KZ)
   !    BAS21=BAS21/(2.0*DX)       
       
   !    BAS32=U3GCO(IX,JY+1,KZ)-U3GCO(IX,JY-1,KZ)
   !    BAS32=BAS32/(2.0*DY)
       
   !    BAS43=U4GCO(IX,JY,KZ+1)-U4GCO(IX,JY,KZ-1)
   !    BAS43=BAS43/(2.0*DZ)

   !    DIVERGCO(IX,JY,KZ)=BAS21+BAS32+BAS43

   !  ENDDO
   !  ENDDO


   !  ENDDO

    
     DO KZ=2,NZ-1  !exclude edges for center, use NSEC ?
     DO JY=2,NY-1
     DO IX=2,NX-1
        !IF( DIVERCO(IX,JY,KZ) .GT. 0 .AND. U1CO(IX,JY,KZ) .LT. DENS_THRE+1.) FLAGV(IX,JY,KZ)=1 !USE DIVER DM
       IF( DIVERGCO(IX,JY,KZ) .GT. 0  .AND. U1CO(IX,JY,KZ) .LT. DENS_THRE+1.) FLAGV(IX,JY,KZ)=1
        !IF( DIVERGCO(IX,JY,KZ) .GT. 0 .AND. DIVERCO(IX,JY,KZ) .GT. 0 .AND. &
        !     U1CO(IX,JY,KZ) .LT. DENS_THRE+1.) FLAGV(IX,JY,KZ)=1 !USE DIVER DM
     ENDDO
     ENDDO
     ENDDO
            
     !diverV for cells at the borders equal to neighbouring cell
     FLAGV(1,:,:)=FLAGV(2,:,:)
     FLAGV(NX,:,:)=FLAGV(NX-1,:,:)
     FLAGV(:,1,:)=FLAGV(:,2,:)
     FLAGV(:,NY,:)=FLAGV(:,NY-1,:)
     FLAGV(:,:,1)=FLAGV(:,:,2)
     FLAGV(:,:,NZ)=FLAGV(:,:,NZ-1)
     !DIVERCO(1,:,:)=DIVERCO(2,:,:)
     !DIVERCO(NX,:,:)=DIVERCO(NX-1,:,:)
     !DIVERCO(:,1,:)=DIVERCO(:,2,:)
     !DIVERCO(:,NY,:)=DIVERCO(:,NY-1,:)
     !DIVERCO(:,:,1)=DIVERCO(:,:,2)
     !DIVERCO(:,:,NZ)=DIVERCO(:,:,NZ-1)
     !DIVERGCO(1,:,:)=DIVERGCO(2,:,:)
     !DIVERGCO(NX,:,:)=DIVERGCO(NX-1,:,:)
     !DIVERGCO(:,1,:)=DIVERGCO(:,2,:)
     !DIVERGCO(:,NY,:)=DIVERGCO(:,NY-1,:)
     !DIVERGCO(:,:,1)=DIVERGCO(:,:,2)
     !DIVERGCO(:,:,NZ)=DIVERGCO(:,:,NZ-1)
     
     RETURN

     END SUBROUTINE MARK


!***************************************************************************
      SUBROUTINE indexx(n,arr,indx)
!***************************************************************************     
      
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) then 
         write(*,*) 'NSTACK too small in indexx'
         stop
        endif 
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

!******************************************************************************** 

!*-------------------------------------
     SUBROUTINE DIVER_FINA_GAS(NL, NPATCH, PARE, PATCHNX, PATCHNY, & 
                PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
                PATCHRX, PATCHRY, PATCHRZ,NPART, IRR)
       USE COMMONDATA
       IMPLICIT NONE
!input variables
       INTEGER NL, IRR
       INTEGER NPART(0:NLEVELS)
       INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
       INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHz(NPALEV)
       REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
!local variables
       INTEGER IR, I, IX, JY, KZ, LOW1, LOW2, N1, N2, N3, NX, NY, NZ
       INTEGER L1, L2, L3, CR1, CR2, CR3, DIM1, DIM2, DIM3, DIM4
       INTEGER IXCO, JYCO, KZCO
       REAL*4 DXR, DYR, DZR
       REAL*4 BAS21, BAS32, BAS43
       REAL*4, ALLOCATABLE:: DIVERG(:,:,:,:), DIVER00(:,:,:)


       NX=REAL(NHYX)*(2.**IRR)
       NY=REAL(NHYX)*(2.**IRR)
       NZ=REAL(NHYX)*(2.**IRR)

       WRITE(*,*) '  MIN/MAX eulerian velocity:'
       WRITE(*,*) '  vx: ',MINVAL(U2CO), MAXVAL(U2CO)
       WRITE(*,*) '  vy: ',MINVAL(U3CO), MAXVAL(U3CO)
       WRITE(*,*) '  vz: ', MINVAL(U4CO), MAXVAL(U4CO)
!*-------------------------------------
!*      Divergencia coarse (DIVERCO)
!*-------------------------------------

       IF(IRR .LT. 0) THEN

          !     NX=NCOX
          !     NY=NCOY
          !     NZ=NCOZ
          !ALLOCATE(DIVERCO(NX,NY,NZ))
          !ALLOCATE(DIVERGCO(NX,NY,NZ))
          DO KZ=1,NZ 
             DO JY=1,NY
                DO IX=1,NX
                   !DIVERCO(IX,JY,KZ)=0.0
                   DIVERGCO(IX,JY,KZ)=0.0
                ENDDO
             ENDDO
          ENDDO

          !divergence only for cells 2,N-1, not  the border
          DO KZ=2,NZ-1 
             DO JY=2,NY-1
                DO IX=2,NX-1

                   !velocity divergence for DM 
                   BAS21=U2CO(IX+1,JY,KZ)-U2CO(IX-1,JY,KZ)
                   BAS21=BAS21/(2.0*DX)       
                   
                   BAS32=U3CO(IX,JY+1,KZ)-U3CO(IX,JY-1,KZ)
                   BAS32=BAS32/(2.0*DY)
                   
                   BAS43=U4CO(IX,JY,KZ+1)-U4CO(IX,JY,KZ-1)
                   BAS43=BAS43/(2.0*DZ)

                   DIVERCO(IX,JY,KZ)=BAS21+BAS32+BAS43
                   
                   !velocity divergence for gas
                   BAS21=U2GCO(IX+1,JY,KZ)-U2GCO(IX-1,JY,KZ)
                   BAS21=BAS21/(2.0*DX)       
       
                   BAS32=U3GCO(IX,JY+1,KZ)-U3GCO(IX,JY-1,KZ)
                   BAS32=BAS32/(2.0*DY)
       
                   BAS43=U4GCO(IX,JY,KZ+1)-U4GCO(IX,JY,KZ-1)
                   BAS43=BAS43/(2.0*DZ)

                   DIVERGCO(IX,JY,KZ)=BAS21+BAS32+BAS43
                   
                ENDDO
             ENDDO


          ENDDO

          DIVERCO(1,:,:)=DIVERCO(2,:,:)
          DIVERCO(NX,:,:)=DIVERCO(NX-1,:,:)
          DIVERCO(:,1,:)=DIVERCO(:,2,:)
          DIVERCO(:,NY,:)=DIVERCO(:,NY-1,:)
          DIVERCO(:,:,1)=DIVERCO(:,:,2)
          DIVERCO(:,:,NZ)=DIVERCO(:,:,NZ-1)
          DIVERGCO(1,:,:)=DIVERGCO(2,:,:)
          DIVERGCO(NX,:,:)=DIVERGCO(NX-1,:,:)
          DIVERGCO(:,1,:)=DIVERGCO(:,2,:)
          DIVERGCO(:,NY,:)=DIVERGCO(:,NY-1,:)
          DIVERGCO(:,:,1)=DIVERGCO(:,:,2)
          DIVERGCO(:,:,NZ)=DIVERGCO(:,:,NZ-1)

          WRITE(*,*) MINVAL(DIVERCO), MAXVAL(DIVERCO)

       ELSE IF (IRR .GE. 0) THEN !I need DIVER0 for IR>0 in roder to compute fix mesh
!*-------------------------------*       
!*      LEVEL 0 
!*-------------------------------*
          IF(ALLOCATED(DIVER0) .EQV. .TRUE.) DEALLOCATE(DIVER0)
          ALLOCATE(DIVER0(NHYX, NHYY, NHYZ)) !COMMON

          DO KZ=2,NHYZ-1
             DO JY=2,NHYY-1
                DO IX=2,NHYX-1
        
                   BAS21=0.0
                   BAS32=0.0
                   BAS43=0.0
        

                   BAS21=U2G(IX+1,JY,KZ)-U2G(IX-1,JY,KZ)
                   BAS21=BAS21/(2.0*DX0)       
                   
                   BAS32=U3G(IX,JY+1,KZ)-U3G(IX,JY-1,KZ)
                   BAS32=BAS32/(2.0*DY0)
       
                   BAS43=U4G(IX,JY,KZ+1)-U4G(IX,JY,KZ-1)
                   BAS43=BAS43/(2.0*DZ0)
   
                   DIVER0(IX,JY,KZ)=BAS21+BAS32+BAS43

                END DO
             END DO
          END DO


          DIVER0(1,:,:)=DIVER0(2,:,:)
          DIVER0(NHYX,:,:)=DIVER0(NHYX-1,:,:)
          DIVER0(:,1,:)=DIVER0(:,2,:)
          DIVER0(:,NHYY,:)=DIVER0(:,NHYY-1,:)
          DIVER0(:,:,1)=DIVER0(:,:,2)
          DIVER0(:,:,NHYZ)=DIVER0(:,:,NHYZ-1)


          IF (IRR .GT. 0) THEN 
!*-------------------------------------
!*      Divergencia fina  (DIVER)
!*-------------------------------------
             DIM1=MAXVAL(PATCHNX)
             DIM2=MAXVAL(PATCHNY)
             DIM3=MAXVAL(PATCHNZ)
             DIM4=SUM(NPATCH(0:NL))
             ALLOCATE(DIVER(DIM1, DIM2, DIM3, DIM4)) !COMMON 
             !ALLOCATE(DIVER0(NHYX, NHYY, NHYZ)) !COMMON
             ALLOCATE(DIVER00(NHYX, NHYY, NHYZ)) !LOCAL
             ALLOCATE(DIVERG(DIM1, DIM2, DIM3, DIM4)) !LOCAL
             
             DIVER=0.0
             DIVER00=0.0
             
             DO IR=1,NL
       
                DXR=0.0
                DYR=0.0
                DZR=0.0
       
                DXR=DX0/(2.0**IR)
                DYR=DY0/(2.0**IR)
                DZR=DZ0/(2.0**IR)
      
                LOW1=SUM(NPATCH(0:IR-1))+1
                LOW2=SUM(NPATCH(0:IR))
                DO I=LOW1, LOW2
       
                   N1=PATCHNX(I)
                   N2=PATCHNY(I)
                   N3=PATCHNZ(I)
       
                   DO KZ=2,N3-1 !exclude the borders of the patches
                      DO JY=2,N2-1
                         DO IX=2,N1-1
        
                            BAS21=0.0
                            BAS32=0.0
                            BAS43=0.0

                            BAS21=U12G(IX+1,JY,KZ,I)-U12G(IX-1,JY,KZ,I)
                            BAS21=BAS21/(2.0*DXR)       
                            
                            BAS32=U13G(IX,JY+1,KZ,I)-U13G(IX,JY-1,KZ,I)
                            BAS32=BAS32/(2.0*DYR)
                            
                            BAS43=U14G(IX,JY,KZ+1,I)-U14G(IX,JY,KZ-1,I)
                            BAS43=BAS43/(2.0*DZR)
                            
                            DIVER(IX,JY,KZ,I)=BAS21+BAS32+BAS43     
                            
                         END DO
                      END DO
                   END DO
       
                END DO
             END DO
             
             !*-------------------------------------
             !*      Divergencia grosera  (MACH)
             !*-------------------------------------
       
             DIVERG=0.0        !SOLO AHORA ES DIVERGENCIA
             
!*      NIVEL NL=1

             IR=1
      
             LOW1=SUM(NPATCH(0:IR-1))+1
             LOW2=SUM(NPATCH(0:IR))
             DO I=LOW1, LOW2

                N1=PATCHNX(I)
                N2=PATCHNY(I)
                N3=PATCHNZ(I)
                L1=PATCHX(I)
                L2=PATCHY(I)
                L3=PATCHZ(I)
                
                DO KZ=1,N3
                   DO JY=1,N2
                      DO IX=1,N1
        
                         BAS21=0.0
                         BAS32=0.0
                         BAS43=0.0

                         IF (IX.LT.3.OR.IX.GT.N1-2.OR. &
                              JY.LT.3.OR.JY.GT.N2-2.OR. &
                              KZ.LT.3.OR.KZ.GT.N3-2) THEN
         
                            CR1=INT((IX+1)/2)+L1-1
                            CR2=INT((JY+1)/2)+L2-1
                            CR3=INT((KZ+1)/2)+L3-1

                            BAS21=U2G(CR1+1,CR2,CR3)-U2G(CR1-1,CR2,CR3)
                            BAS21=BAS21/(2.0*DX0)       
       
                            BAS32=U3G(CR1,CR2+1,CR3)-U3G(CR1,CR2-1,CR3)
                            BAS32=BAS32/(2.0*DY0)
       
                            BAS43=U4G(CR1,CR2,CR3+1)-U4G(CR1,CR2,CR3-1)
                            BAS43=BAS43/(2.0*DZ0)

                            DIVERG(IX,JY,KZ,I)=BAS21+BAS32+BAS43
                            
                         END IF
                         
                      END DO
                   END DO
                END DO
             END DO
             
!*      NIVELES AMR>1
             DO IR=2,0 !NL
                
                DXR=0.0
                DYR=0.0
                DZR=0.0
                DXR=DX0/(2.0**IR)
                DYR=DY0/(2.0**IR)
                DZR=DZ0/(2.0**IR)
                
                LOW1=SUM(NPATCH(0:IR-1))+1
                LOW2=SUM(NPATCH(0:IR))
                DO I=LOW1, LOW2
                   
                   N1=PATCHNX(I)
                   N2=PATCHNY(I)
                   N3=PATCHNZ(I)
                   L1=PATCHX(I)
                   L2=PATCHY(I)
                   L3=PATCHZ(I)
                   
                   DO KZ=1,N3
                      DO JY=1,N2
                         DO IX=1,N1
                            
                            BAS21=0.0
                            BAS32=0.0
                            BAS43=0.0
                            
                            !*       shock detection 
                            IF (IX.LT.3.OR.IX.GT.N1-2.OR. &
                                 JY.LT.3.OR.JY.GT.N2-2.OR. &
                                 KZ.LT.3.OR.KZ.GT.N3-2) THEN

                               CR1=INT((IX+1)/2)+L1-1
                               CR2=INT((JY+1)/2)+L2-1
                               CR3=INT((KZ+1)/2)+L3-1
                               IF(CR2 .EQ. 1) WRITE(*,*) CR2,IR, I,  JY, L2, ix, jy, kz
                               BAS21=U12G(CR1+1,CR2,CR3,PARE(I))- &
                                    U12G(CR1-1,CR2,CR3,PARE(I)) 
                               BAS21=BAS21/(4.0*DXR)
                               
                               BAS32=U13G(CR1,CR2+1,CR3,PARE(I))- &
                                    U13G(CR1,CR2-1,CR3,PARE(I))
                               BAS32=BAS32/(4.0*DYR)
                               
                               BAS43=U14G(CR1,CR2,CR3+1,PARE(I))- &
                                    U14G(CR1,CR2,CR3-1,PARE(I))
                               BAS43=BAS43/(4.0*DZR)
                               
                               DIVERG(IX,JY,KZ,I)=BAS21+BAS32+BAS43
                               
                            END IF
                            
                         END DO
                      END DO
                   END DO
                END DO
             END DO !LOOP ON LEVELS
             
             !*      DERIVADA MAS SUAVE PARA LAS CELDAS 2 Y N-1
             DO IR=1, 1 !NL
                
                DXR=DX0/(2.0**IR)
                DYR=DY0/(2.0**IR)
                DZR=DZ0/(2.0**IR)
                
                LOW1=SUM(NPATCH(0:IR-1))+1
                LOW2=SUM(NPATCH(0:IR))
                DO I=LOW1, LOW2
                   
                   N1=PATCHNX(I)
                   N2=PATCHNY(I)
                   N3=PATCHNZ(I)
                   
                   DO KZ=1, N3
                      DO JY=1, N2
                         DO IX=1, N1
                            

                            IF(IX.EQ.2.OR.IX.EQ.N1-1.OR.JY.EQ.2.OR.JY.EQ.N2-1.OR. &
                                 KZ.EQ.2.OR.KZ.EQ.N3-1) THEN 
                               

                               DIVER(IX,JY,KZ,I)=0.5*(DIVERG(IX,JY,KZ,I)+ &
                                    DIVER(IX,JY,KZ,I)) 
                               
                            END IF
       
                            
                            IF(IX.EQ.1.OR.IX.EQ.N1.OR.JY.EQ.1.OR.JY.EQ.N2.OR. &
                                 KZ.EQ.1.OR.KZ.EQ.N3) THEN 
                               DIVER(IX,JY,KZ,I)=DIVERG(IX,JY,KZ,I)
                            END IF
                            
                         END DO
                      END DO
                   END DO
                END DO
             END DO !LOOP ON LEVELS
             
             DEALLOCATE(DIVERG)
          
          ENDIF !IRR > 0
       
          ENDIF !IR test
          

          RETURN
          !

          
     END SUBROUTINE DIVER_FINA_GAS
!*-------------------------------------

     SUBROUTINE DIVER_FINA(NL, NPATCH, PARE, PATCHNX, PATCHNY, & 
                PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
                PATCHRX, PATCHRY, PATCHRZ,NPART)
       USE COMMONDATA
       IMPLICIT NONE
!input variables
       INTEGER NL
       INTEGER NPART(0:NLEVELS)
       INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
       INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHz(NPALEV)
       REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
!local variables
       INTEGER IR, I, IX, JY, KZ, LOW1, LOW2, N1, N2, N3
       INTEGER L1, L2, L3, CR1, CR2, CR3, DIM1, DIM2, DIM3, DIM4
       REAL*4 DXR, DYR, DZR
       REAL*4 BAS21, BAS32, BAS43
       REAL*4, ALLOCATABLE:: DIVERG(:,:,:,:)

!*-------------------------------------
!*      Divergencia fina  (DIVER)
!*-------------------------------------
       DIM1=MAXVAL(PATCHNX)
       DIM2=MAXVAL(PATCHNY)
       DIM3=MAXVAL(PATCHNZ)
       DIM4=SUM(NPATCH(0:NL))
       ALLOCATE(DIVER(DIM1, DIM2, DIM3, DIM4)) !COMMON 
       ALLOCATE(DIVER0(NHYX, NHYY, NHYZ)) !COMMON
       ALLOCATE(DIVERG(DIM1, DIM2, DIM3, DIM4)) !LOCAL

       DIVER=0.0
       DIVER0=0.0

       DO IR=1,NL
       
       DXR=0.0
       DYR=0.0
       DZR=0.0
       
       DXR=DX0/(2.0**IR)
       DYR=DY0/(2.0**IR)
       DZR=DZ0/(2.0**IR)
      
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1, LOW2
       
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)
       
       DO KZ=2,N3-1 !exclude the borders of the patches
       DO JY=2,N2-1
       DO IX=2,N1-1
        
        BAS21=0.0
        BAS32=0.0
        BAS43=0.0

        BAS21=U12(IX+1,JY,KZ,I)-U12(IX-1,JY,KZ,I)
        BAS21=BAS21/(2.0*DXR)       
       
        BAS32=U13(IX,JY+1,KZ,I)-U13(IX,JY-1,KZ,I)
        BAS32=BAS32/(2.0*DYR)
       
        BAS43=U14(IX,JY,KZ+1,I)-U14(IX,JY,KZ-1,I)
        BAS43=BAS43/(2.0*DZR)
       
        DIVER(IX,JY,KZ,I)=BAS21+BAS32+BAS43     

       END DO
       END DO
       END DO
       
       END DO
       END DO
       
!*-------------------------------------
!*      Divergencia grosera  (MACH)
!*-------------------------------------
       
       DIVERG=0.0        !SOLO AHORA ES DIVERGENCIA
       
!*      NIVEL NL=1

       IR=1
      
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1, LOW2

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)
       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)
       
       DO KZ=1,N3
       DO JY=1,N2
       DO IX=1,N1
        
        BAS21=0.0
        BAS32=0.0
        BAS43=0.0

        IF (IX.LT.3.OR.IX.GT.N1-2.OR. &
           JY.LT.3.OR.JY.GT.N2-2.OR. &
           KZ.LT.3.OR.KZ.GT.N3-2) THEN
         
          CR1=INT((IX+1)/2)+L1-1
          CR2=INT((JY+1)/2)+L2-1
          CR3=INT((KZ+1)/2)+L3-1

        BAS21=U2(CR1+1,CR2,CR3)-U2(CR1-1,CR2,CR3)
        BAS21=BAS21/(2.0*DX0)       
       
        BAS32=U3(CR1,CR2+1,CR3)-U3(CR1,CR2-1,CR3)
        BAS32=BAS32/(2.0*DY0)
       
        BAS43=U4(CR1,CR2,CR3+1)-U4(CR1,CR2,CR3-1)
        BAS43=BAS43/(2.0*DZ0)

        DIVERG(IX,JY,KZ,I)=BAS21+BAS32+BAS43

      END IF

       END DO
       END DO
       END DO
       END DO

!*      NIVELES AMR>1
       DO IR=2,NL

       DXR=0.0
       DYR=0.0
       DZR=0.0
       DXR=DX0/(2.0**IR)
       DYR=DY0/(2.0**IR)
       DZR=DZ0/(2.0**IR)
       
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1, LOW2

       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)
       L1=PATCHX(I)
       L2=PATCHY(I)
       L3=PATCHZ(I)
       
       DO KZ=1,N3
       DO JY=1,N2
       DO IX=1,N1
        
        BAS21=0.0
        BAS32=0.0
        BAS43=0.0

!*       shock detection 
        IF (IX.LT.3.OR.IX.GT.N1-2.OR. &
           JY.LT.3.OR.JY.GT.N2-2.OR. &
           KZ.LT.3.OR.KZ.GT.N3-2) THEN
          
          CR1=INT((IX+1)/2)+L1-1
          CR2=INT((JY+1)/2)+L2-1
          CR3=INT((KZ+1)/2)+L3-1

        BAS21=U12(CR1+1,CR2,CR3,PARE(I))- &
             U12(CR1-1,CR2,CR3,PARE(I)) 
        BAS21=BAS21/(4.0*DXR)
              
        BAS32=U13(CR1,CR2+1,CR3,PARE(I))- &
             U13(CR1,CR2-1,CR3,PARE(I))
        BAS32=BAS32/(4.0*DYR)
       
        BAS43=U14(CR1,CR2,CR3+1,PARE(I))- &
             U14(CR1,CR2,CR3-1,PARE(I))
        BAS43=BAS43/(4.0*DZR)
   
        DIVERG(IX,JY,KZ,I)=BAS21+BAS32+BAS43
     
        END IF

        END DO
        END DO
        END DO
       END DO
    END DO !LOOP ON LEVELS

!*      DERIVADA MAS SUAVE PARA LAS CELDAS 2 Y N-1
       DO IR=1,NL

       DXR=DX0/(2.0**IR)
       DYR=DY0/(2.0**IR)
       DZR=DZ0/(2.0**IR)

       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))
       DO I=LOW1, LOW2
 
       N1=PATCHNX(I)
       N2=PATCHNY(I)
       N3=PATCHNZ(I)
       
       DO KZ=1, N3
       DO JY=1, N2
       DO IX=1, N1


       IF(IX.EQ.2.OR.IX.EQ.N1-1.OR.JY.EQ.2.OR.JY.EQ.N2-1.OR. &
         KZ.EQ.2.OR.KZ.EQ.N3-1) THEN 

       DIVER(IX,JY,KZ,I)=0.5*(DIVERG(IX,JY,KZ,I)+ &
                                DIVER(IX,JY,KZ,I)) 
     
       END IF


       IF(IX.EQ.1.OR.IX.EQ.N1.OR.JY.EQ.1.OR.JY.EQ.N2.OR. &
          KZ.EQ.1.OR.KZ.EQ.N3) THEN 
          DIVER(IX,JY,KZ,I)=DIVERG(IX,JY,KZ,I)
       END IF

       END DO
       END DO
       END DO
       END DO
    END DO !LOOP ON LEVELS

    DEALLOCATE(DIVERG)

!*-------------------------------*       
!*      LEVEL 0
!*-------------------------------*

       DO KZ=2,NHYZ-1
       DO JY=2,NHYY-1
       DO IX=2,NHYX-1
        
        BAS21=0.0
        BAS32=0.0
        BAS43=0.0
        

        BAS21=U2(IX+1,JY,KZ)-U2(IX-1,JY,KZ)
        BAS21=BAS21/(2.0*DX0)       
       
        BAS32=U3(IX,JY+1,KZ)-U3(IX,JY-1,KZ)
        BAS32=BAS32/(2.0*DY0)
       
        BAS43=U4(IX,JY,KZ+1)-U4(IX,JY,KZ-1)
        BAS43=BAS43/(2.0*DZ0)
   
        DIVER0(IX,JY,KZ)=BAS21+BAS32+BAS43

       END DO
       END DO
       END DO


       DIVER0(1,:,:)=DIVER0(2,:,:)
       DIVER0(NHYX,:,:)=DIVER0(NHYX-1,:,:)
       DIVER0(:,1,:)=DIVER0(:,2,:)
       DIVER0(:,NHYY,:)=DIVER0(:,NHYY-1,:)
       DIVER0(:,:,1)=DIVER0(:,:,2)
       DIVER0(:,:,NHYZ)=DIVER0(:,:,NHYZ-1)
       
!*-------------------------------*       
!*      COARSE LEVEL
!*-------------------------------*

       RETURN

     END SUBROUTINE DIVER_FINA



!*----------------------------------------------------------------------*  
     SUBROUTINE PAFIND(RX,RY,RZ, IR0, DXR, DYR, DZR, NL, NPATCH, PARE, PATCHNX, PATCHNY, & 
                PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
                PATCHRX, PATCHRY, PATCHRZ, &
                IPA0, I0, J0, K0)
       USE COMMONDATA, ONLY: NLEVELS, NPALEV, DX0, DY0, DZ0, RADX0, RADY0, RADZ0, U11 !then remove U11
       IMPLICIT NONE
       !input variables
       INTEGER NL
       INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
       INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
       INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHz(NPALEV)
       REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
       INTEGER:: IR0
       REAL*4 RX, RY, RZ, DXR, DYR, DZR
       !local variables
       INTEGER IR, LOW1, LOW2, IPA, II, JJ, KK, FLAG
       REAL*4 BAS
       REAL*4 RXP, RYP, RZP
       !output
       INTEGER IPA0, I0, J0, K0

       IR=IR0+1

       IPA0=0
       II=0
       JJ=0
       KK=0
       FLAG=0
       LOW1=SUM(NPATCH(0:IR-1))+1
       LOW2=SUM(NPATCH(0:IR))

       DO IPA=LOW1, LOW2
          IF(FLAG .EQ. 0 ) THEN 
          RXP=PATCHRX(IPA)-0.5*DXR
          RYP=PATCHRY(IPA)-0.5*DYR
          RZP=PATCHRZ(IPA)-0.5*DZR

          BAS=(RX-RXP)/DXR
          IF(BAS .GT. -0.5) II=INT(BAS+0.49999)+1
          IF(BAS .LT. -0.5) II=0

          BAS=(RY-RYP)/DYR
          IF(BAS .GT. -0.5) JJ=INT(BAS+0.49999)+1
          IF(BAS .LE. -0.5) JJ=0

          BAS=(RZ-RZP)/DZR
          IF(BAS .GT. -0.5) KK=INT(BAS+0.49999)+1
          IF(BAS .LE. -0.5) KK=0

          IF(II .GE. 1 .AND. II .LE. PATCHNX(IPA) .AND. &
              JJ .GE. 1 .AND. JJ .LE. PATCHNY(IPA) .AND. & 
              KK .GE. 1 .AND. KK .LE. PATCHNZ(IPA)) THEN

             IPA0=IPA            
             I0=II
             J0=JJ
             K0=KK
             FLAG=1
             !EXIT !use for serial run

          ENDIF
       ENDIF

       ENDDO


       RETURN


     END SUBROUTINE PAFIND


!*----------------------------------------------------------------------*  
     SUBROUTINE OVERLAPPING( MINF, MAXF, NOVER, NVOID, INDICE, NXX, NYY, NZZ, &
          DXX, DYY, DZZ, RXX1, RYY1, RZZ1, VOLNEW, &
          UVOID, XC, YC, ZC, NVOID2)
     USE COMMONDATA, ONLY: INICIOX, FINALX, INICIOY, FINALY, INICIOZ, FINALZ, ICX, ICY, ICZ, &
          RINIXCO, RFINXCO, RINIYCO, RFINYCO, RINIZCO, RFINZCO, RADX, RADY, RADZ, DX, DY, DZ, VOL, &
          MARCA
     IMPLICIT NONE
!input variables
     INTEGER:: NVOID, NOVER
     INTEGER:: NXX, NYY, NZZ
     REAL*4:: MINF, MAXF, DXX, DYY, DZZ, RXX1, RYY1, RZZ1
!local variables
     INTEGER I, III, J, IND0, IND2, IND1, IND3, INIX, FINX, INIY, FINY, INIZ, FINZ, IO, IOO, K, KK
     INTEGER IM, INDM, IND00, BASIND00, INDF, INDMIN, INDMAX, NU2
     INTEGER NVOID2, NOVT, NOVT2, IO2, IM1,IM2, INDM1, INDM2, IO3
     INTEGER:: IX, JY, KZ
     INTEGER:: FLAGX, FLAGY, FLAGZ
     INTEGER:: NOV(NVOID), OVERL(NVOID,NOVER), FLAGJ(NVOID, NOVER)
     INTEGER:: NMERG(NVOID), MERG(NVOID,NOVER)!, AA(NVOID), OVERL2(NVOID, NVOID), NOV2(NVOID)
     INTEGER, DIMENSION(NVOID) ::  INDICE
     REAL*4, DIMENSION(4)::  AX, AY, AZ
     INTEGER, DIMENSION(4)::  INDX, INDY,  INDZ
     REAL*4:: RX1, RX2, RY1, RY2, RZ1, RZ2, RX3, RX4, RY3, RY4, RZ3, RZ4
     REAL*4:: DV, DV2, VOL1, VOL2, VOL22, FRC, VOL11, SUMDX, SUMDY, SUMDZ, VT, FRAC
!output variables
     REAL*4, DIMENSION(NVOID)::XC,YC,ZC
     REAL*4, DIMENSION(NVOID):: VOLNEW
     INTEGER, DIMENSION(NVOID):: UVOID
     !INTEGER, DIMENSION(NXX, NYY, NZZ) :: MARCA

     !DX--> DDX
     ! NCOX --> NXX
!*----------------------------------------------------------*       
!*      cross-match all voids to find all the overlappings
!*----------------------------------------------------------*

     NOV(:)=0
     XC(:)=0
     YC(:)=0
     ZC(:)=0
     VOLNEW(:)=0.
     UVOID(:)=-1
     FLAGJ(:,:)=0
     NVOID2=NVOID
     MARCA(:,:,:)=0

     write(*,*) 'Routine overlapping:'
!NEW
     DO I=1, NVOID

        IND0=INDICE(I) !INDICE(1) --> largest void
        IF(IND0==0) WRITE(*,*) 'WARNING: IND0=0', I

        IF(UVOID(IND0) .EQ. 0) CYCLE !only uvoid=-1, -2

        RX1=RINIXCO(IND0)-0.5*DXX !EDGES OF THE CELL
        RX2=RFINXCO(IND0)+0.5*DXX
        RY1=RINIYCO(IND0)-0.5*DYY
        RY2=RFINYCO(IND0)+0.5*DYY
        RZ1=RINIZCO(IND0)-0.5*DZZ
        RZ2=RFINZCO(IND0)+0.5*DZZ
        XC(IND0)=RADX(1)+(ICX(IND0)-1)*DXX
        YC(IND0)=RADY(1)+(ICY(IND0)-1)*DYY
        ZC(IND0)=RADZ(1)+(ICZ(IND0)-1)*DZZ

        IF(VOLNEW(IND0) .EQ. 0 ) VOLNEW(IND0)=(RX2-RX1)*(RY2-RY1)*(RZ2-RZ1)



        DO J=I+1, NVOID
           IND2=INDICE(J) 
           IF(UVOID(IND2) .EQ. 0) CYCLE
           IF(IND2==0) WRITE(*,*) 'WARNING: IND2=0', j
           FLAGX=0
           FLAGY=0
           FLAGZ=0
           RX3=RINIXCO(IND2)-0.5*DXX !EDGES OF THE CELL
           RX4=RFINXCO(IND2)+0.5*DXX
           RY3=RINIYCO(IND2)-0.5*DYY
           RY4=RFINYCO(IND2)+0.5*DYY
           RZ3=RINIZCO(IND2)-0.5*DZZ
           RZ4=RFINZCO(IND2)+0.5*DZZ
 
          IF(VOLNEW(IND2) .EQ. 0 ) VOLNEW(IND2)=(RX4-RX3)*(RY4-RY3)*(RZ4-RZ3)
          !DO IX=INICIOX(IND2), FINALX(IND2) !metterlo prima
          !   DO JY=INICIOY(IND2), FINALY(IND2)
          !      DO KZ=INICIOZ(IND2), FINALZ(IND2)
          !         IF(MARCA(IX,JY, KZ) .EQ. 0)  MARCA(IX,JY, KZ)=IND2
          !      ENDDO
          !   ENDDO
          !ENDDO

           SUMDX=RX2-RX1+RX4-RX3
           SUMDY=RY2-RY1+RY4-RY3
           SUMDZ=RZ2-RZ1+RZ4-RZ3
           AX=(/RX1, RX2, RX3, RX4/)
           AY=(/RY1, RY2, RY3, RY4/)
           AZ=(/RZ1, RZ2, RZ3, RZ4/)
           CALL INDEXX(4, AX, INDX)
           CALL INDEXX(4, AY, INDY)
           CALL INDEXX(4, AZ, INDZ)
           IF( (AX(INDX(4))-AX(INDX(1))) .LT. SUMDX ) FLAGX=1
           IF( (AY(INDY(4))-AY(INDY(1))) .LT. SUMDY ) FLAGY=1
           IF( (AZ(INDZ(4))-AZ(INDZ(1))) .LT. SUMDZ ) FLAGZ=1
 
          IF(FLAGX==1 .AND. FLAGY==1 .AND. FLAGZ==1) THEN !overlapping
                CALL VOLUME(IND0, IND2, NXX, NYY, NZZ, DXX, DYY, DZZ, &
                     RXX1, RYY1, RZZ1,  &
                     VOL11, VOL22, VOL2,  DV, DV2, FRAC) !find overlapping volume and mark cells in DV
                DV2=MAX(DV2,0.)
                IF(FRAC .GT. MAXF) THEN   
                   UVOID(IND2)=0  
                   NVOID2=NVOID2-1
                ELSE IF( FRAC .GE. MINF .AND. FRAC .LE. MAXF) THEN   
                  !sept2012! VOLNEW(IND2)=VOLNEW(IND2)-DV  !solo se  DV/VOL2 .GE. MINF .AND. DV/VOL2 .LE. MAXF                   
                   IF(VOLNEW(IND2) .LE. 0.) THEN !define DVV?
                      VOLNEW(IND2)=0.
                      UVOID(IND2)=0
                      NVOID2=NVOID2-1
                   ELSE
                      NOV(IND0)=NOV(IND0)+1
                      NOV(IND2)=NOV(IND2)+1
                      OVERL(IND0,NOV(IND0))=IND2
                      OVERL(IND2,NOV(IND2))=IND0
                      FLAGJ(IND0, NOV(IND0))=1 ! merger IND0/IND2
                      FLAGJ(IND2, NOV(IND2))=1 ! merge  IND0/IND2
                   ENDIF

                ELSE !considered not overlapping, but if already merged into others update the volume
                   !sept2012! VOLNEW(IND2)=VOLNEW(IND2)-DV !? DEFINIRE VOLNEW2 ?
                   UVOID(IND2)=-2
                   IF(VOLNEW(IND2) .LE. 0) THEN
                      VOLNEW(IND2)=0.
                      UVOID(IND2)=0
                      NVOID2=NVOID2-1
                   ENDIF
                ENDIF
           ENDIF
        ENDDO
ENDDO

NU2=0
DO I=1, NVOID
IND0=INDICE(I)
IF(UVOID(I) .EQ. -2) THEN
NU2=NU2+1
UVOID(I)=-1
ENDIF
ENDDO

WRITE(*,*) 'Maximum number of overlappings per void:', MAXVAL(NOV), NOVER, NVOID

!*----------------------------------------------------------*       
!*     Redefine voids:
!case a) DV/V2: [0.4,0.8] --> join voids
!case b) DV/V2 < 0.4 --> keep 2 separate voids
!case c) DV/V2 > 0.8 --> remove the smallest 
!*----------------------------------------------------------*
     NOVT=SUM(NOV) !total number of overlappings
     !write(*,*) 'NOVT:',NOVT, COUNT(UVOID .EQ. -1), COUNT(UVOID .EQ. 0 ), sum(volnew)
     
     NMERG(:)=0
     MERG(:,:)=0
     MARCA(:,:,:)=0
     VT=0.
     DO I=1, NVOID

        IND0=INDICE(I)

        IF(UVOID(IND0) .NE. -1) CYCLE   !IF UVOID(IND0)>1 DEVO CONTROLLARE CHE I SUOI OVERLAPPING SIANO STATI ELIMINATI  

        DO IX=INICIOX(IND0), FINALX(IND0) !metterlo prima
           DO JY=INICIOY(IND0), FINALY(IND0)
              DO KZ=INICIOZ(IND0), FINALZ(IND0)
                 IF(MARCA(IX,JY, KZ) .EQ. 0)  MARCA(IX,JY, KZ)=IND0
              ENDDO
           ENDDO
        ENDDO
        VT=VT+VOLNEW(IND0)

        IND00=IND0
 
       DO IO=1, NOV(IND0) !--> nb IND0 PUO' MERGERE CON ALTRO MASTER VOID: UPDATE VALUES OF IND00

           J=OVERL(IND0,IO) 
           IND2=J

           IF(UVOID(IND2) .EQ. 0 ) CYCLE !ALREADY REMOVED: it can appear as overlapped if it has been removed later
             
           IF(FLAGJ(IND0,IO) .EQ. 0) WRITE(*,*) 'WARNING: FLAGJ=0 FOR OVERL>0 !', IND0, IO, IND2

              IF(UVOID(IND2) .EQ. -1) THEN !in this case IND2 is the smallest 
                 !if uvoid(ind2)=-1 means that it's the first time it appears and nmerg(ind2)=0, otherwise
                 ! the overlapping IND2-IND0 would already have been considered
                 NVOID2=NVOID2-1
                 UVOID(IND2)=IND00

                 VOL1=VOLNEW(IND00) !/= VOL1 IF NMERG(IND0)>0
                 VOL2=VOLNEW(IND2) !/= VOL2 IF NMERG(IND2)>0

                 VOLNEW(IND00)=VOL1+VOL2

                 XC(IND00)=(XC(IND00)*VOL1+XC(IND2)*VOL2)/(VOL1+VOL2)
                 YC(IND00)=(YC(IND00)*VOL1+YC(IND2)*VOL2)/(VOL1+VOL2)
                 ZC(IND00)=(ZC(IND00)*VOL1+ZC(IND2)*VOL2)/(VOL1+VOL2)

                 DO IX=INICIOX(IND2), FINALX(IND2)
                    DO JY=INICIOY(IND2), FINALY(IND2)
                       DO KZ=INICIOZ(IND2), FINALZ(IND2)
                          !IF(MARCA(IX,JY, KZ) .GT. 0) WRITE(*,*) 'WARNING: MARCA>0:',MARCA(IX,JY, KZ), IND2, IND0, IND00
                          IF(MARCA(IX,JY, KZ) .EQ. 0) MARCA(IX,JY, KZ)=IND00
                       ENDDO
                    ENDDO
                 ENDDO
                 
                 NMERG(IND00)=NMERG(IND00)+1 !forse lo posso aggiungere dopo
                 MERG(IND00,NMERG(IND00))=IND2
 
                 IF(NMERG(IND2) .GT. 0) WRITE(*,*) 'WARNING: NMERG( UVOID=-1 ) > 0!!', IND00, IND0, IND2, UVOID(IND2) 

                 
              ELSE IF(UVOID(IND2) .GT. 0 ) THEN !IND2 ALREADY MERGED INTO LARGE VOID --> MERGE EVERYTHING INTO MASTER VOID
                 
                 INDF=UVOID(IND2) !IND2'S FATHER VOID      

                 IF(INDF .EQ. IND00) CYCLE !TOGLIERE ? cosi' non sto considerando il primo overl indo-ind2 perche' uvoid=id0 per il ciclo iniziale

                 NVOID2=NVOID2-1 

                 IF (VOLNEW(IND00) .GE. VOLNEW(INDF) ) THEN !IND00 IS THE LARGEST

                    BASIND00=IND00 !DOESN'T CHANGE
                    
                    NMERG(IND00)=NMERG(IND00)+1
                    MERG(IND00,NMERG(IND00))=INDF
                    UVOID(INDF)=IND00                    
                    DO IM=1, NMERG(INDF) !--> loop over NOV(indm)
                       INDM=MERG(INDF,IM)
                       IF(NMERG(INDM) .GT. 0) WRITE(*,*) 'WARNING 2!! NMERG(INDM)>0', NMERG(INDM), &
                            MERG(INDM,1), INDM, INDF, IND00, IND0, IO
                       NMERG(INDM)=0 
                       UVOID(INDM)=IND00
                       NMERG(IND00)=NMERG(IND00)+1
                       MERG(IND00,NMERG(IND00))=INDM
                       DO IX=INICIOX(INDM), FINALX(INDM)
                          DO JY=INICIOY(INDM), FINALY(INDM)
                             DO KZ=INICIOZ(INDM), FINALZ(INDM)
                                IF(MARCA(IX,JY,KZ) .EQ. INDF) MARCA(IX,JY, KZ)=IND00 
                                !there might be cells in the INDM limits belonging to a unique void not overlapping with INDM (FRC<MINF)
                                !MARCA(IX,JY, KZ)=IND00
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                    NMERG(INDF)=0

                    DO IX=INICIOX(IND2), FINALX(IND2)
                       DO JY=INICIOY(IND2), FINALY(IND2)
                          DO KZ=INICIOZ(IND2), FINALZ(IND2)
                             IF(MARCA(IX,JY,KZ) .EQ. INDF) MARCA(IX,JY, KZ)=BASIND00 !non e' gia' stata definita se uvoid(ind2)>0 ?
                             !MARCA(IX,JY, KZ)=BASIND00
                          ENDDO
                       ENDDO
                    ENDDO
                    DO IX=INICIOX(INDF), FINALX(INDF) !NEW
                       DO JY=INICIOY(INDF), FINALY(INDF)
                          DO KZ=INICIOZ(INDF), FINALZ(INDF)
                             IF(MARCA(IX,JY,KZ) .EQ. INDF) MARCA(IX,JY, KZ)=BASIND00 !non e' gia' stata definita se uvoid(ind2)>0 ?
                             !MARCA(IX,JY, KZ)=BASIND00
                          ENDDO
                       ENDDO
                    ENDDO

                 ELSE !INDF IS THE LARGEST -- > BECOME MASTER VOID

                    BASIND00=INDF !CHANGE MASTER VOID
                    NMERG(INDF)=NMERG(INDF)+1
                    MERG(INDF,NMERG(INDF))=IND00
                    UVOID(IND00)=INDF

                    DO IM=1, NMERG(IND00) 
                       INDM=MERG(IND00,IM)
                     !  IF(NOV(INDM) .GT. 0 ) WRITE(*,*) 'WARNING!!!NOV(INDM)>0', IND00, INDM, NOV(INDM)
                       IF(NMERG(INDM) .GT. 0) WRITE(*,*) 'WARNING!! NMERG(INDM)>0', NMERG(INDM), MERG(INDM,1),&
                            INDM,  IND00, INDF, IND0, IO
                       NMERG(INDM)=0 
                       UVOID(INDM)=INDF
                       NMERG(INDF)=NMERG(INDF)+1
                       MERG(INDF,NMERG(INDF))=INDM
                       DO IX=INICIOX(INDM), FINALX(INDM)
                          DO JY=INICIOY(INDM), FINALY(INDM)
                             DO KZ=INICIOZ(INDM), FINALZ(INDM)
                                IF(MARCA(IX,JY,KZ) .EQ. IND00) MARCA(IX,JY, KZ)=INDF
                                !MARCA(IX,JY, KZ)=INDF
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO

                    DO IX=INICIOX(IND2), FINALX(IND2)!NEW
                       DO JY=INICIOY(IND2), FINALY(IND2)
                          DO KZ=INICIOZ(IND2), FINALZ(IND2)
                             IF(MARCA(IX,JY,KZ) .EQ. IND00) MARCA(IX,JY, KZ)=INDF !non e' gia' stata definita se uvoid(ind2)>0 ?
                             !MARCA(IX,JY, KZ)=BASIND00
                          ENDDO
                       ENDDO
                    ENDDO
                    DO IX=INICIOX(IND00), FINALX(IND00) !togliere quello iniziale?
                       DO JY=INICIOY(IND00), FINALY(IND00)
                          DO KZ=INICIOZ(IND00), FINALZ(IND00)
                             IF(MARCA(IX,JY,KZ) .EQ. IND00) MARCA(IX,JY, KZ)=INDF
                          ENDDO
                       ENDDO
                    ENDDO
                    
                    NMERG(IND00)=0

                 ENDIF !test INDF/IND00

                 VOLNEW(BASIND00)=VOLNEW(IND00)+VOLNEW(INDF)

                 XC(BASIND00)=(XC(IND00)*VOLNEW(IND00)+XC(INDF)*VOLNEW(INDF))/(VOLNEW(IND00)+VOLNEW(INDF))
                 YC(BASIND00)=(YC(IND00)*VOLNEW(IND00)+YC(INDF)*VOLNEW(INDF))/(VOLNEW(IND00)+VOLNEW(INDF))
                 ZC(BASIND00)=(ZC(IND00)*VOLNEW(IND00)+ZC(INDF)*VOLNEW(INDF))/(VOLNEW(IND00)+VOLNEW(INDF))
                 !DO IX=INICIOX(IND2), FINALX(IND2)
                 !   DO JY=INICIOY(IND2), FINALY(IND2)
                 !      DO KZ=INICIOZ(IND2), FINALZ(IND2)
                 !          IF(MARCA(IX,JY,KZ) .EQ. IND2) MARCA(IX,JY, KZ)=BASIND00 !non e' gia' stata definita se uvoid(ind2)>0 ?
                 !      ENDDO
                 !   ENDDO
                 !ENDDO
                 IND00=BASIND00

              ENDIF !test UVOID(IND2)

           
        ENDDO !LOOP ON OVERLAPPINGS


     ENDDO !LOOP ON ALL VOIDS       


     END SUBROUTINE OVERLAPPING

!*----------------------------------------------------------------------*  
     SUBROUTINE VOLUME(IND1, IND2, NXX, NYY, NZZ, DXX, DYY, DZZ, &
          RXX1, RYY1, RZZ1, &
          VOL1, VOL22, VOLMIN,  DV, DV2, FRAC)
     USE COMMONDATA, ONLY: RINIXCO, RFINXCO, RINIYCO, RFINYCO, RINIZCO, RFINZCO, &
          DX, DY, DZ, RADX, RADY, RADZ, MARCA
     IMPLICIT NONE
     !input variables
     INTEGER IND1, IND2
     INTEGER:: NXX, NYY, NZZ
     REAL*4:: DXX, DYY, DZZ, RXX1, RYY1, RZZ1
     !INTEGER, DIMENSION(NXX,NYY,NZZ):: MARCA!INOUT
     !local variables
     INTEGER:: IX,IX1,IX2, JY, JY1, JY2,KZ,KZ1,KZ2, INDMIN, INDMAX
     REAL*4 RX1, RX2, RY1, RY2, RZ1, RZ2,  RX3, RX4, RY3, RY4, RZ3, RZ4
     REAL*4 DDX, DDY, DDZ, VCELL, LENG
     REAL*4, DIMENSION(4)::  AX, AY, AZ
     INTEGER, DIMENSION(4)::  INDX, INDY,  INDZ
     !output variables
     REAL*4  DV, DV2, VOL1, VOLMIN, FRAC, VOL22

           RX1=RINIXCO(IND1)-0.5*DX
           RX2=RFINXCO(IND1)+0.5*DX
           RY1=RINIYCO(IND1)-0.5*DY
           RY2=RFINYCO(IND1)+0.5*DY
           RZ1=RINIZCO(IND1)-0.5*DZ
           RZ2=RFINZCO(IND1)+0.5*DZ
           VOL1=(RX2-RX1)*(RY2-RY1)*(RZ2-RZ1) 

           RX3=RINIXCO(IND2)-0.5*DX
           RX4=RFINXCO(IND2)+0.5*DX
           RY3=RINIYCO(IND2)-0.5*DY
           RY4=RFINYCO(IND2)+0.5*DY
           RZ3=RINIZCO(IND2)-0.5*DZ
           RZ4=RFINZCO(IND2)+0.5*DZ
           VOL22=(RX4-RX3)*(RY4-RY3)*(RZ4-RZ3)

           AX=(/RX1, RX2, RX3, RX4/)
           AY=(/RY1, RY2, RY3, RY4/)
           AZ=(/RZ1, RZ2, RZ3, RZ4/)
           CALL INDEXX(4, AX, INDX)
           CALL INDEXX(4, AY, INDY)
           CALL INDEXX(4, AZ, INDZ)

           DDX=AX(INDX(3))-AX(INDX(2))
           DDY=AY(INDY(3))-AY(INDY(2))
           DDZ=AZ(INDZ(3))-AZ(INDZ(2))
           DV=DDX*DDY*DDZ
           
           VOLMIN=VOL22
           INDMIN=IND2
           INDMAX=IND1
           IF(VOL22 .GT. VOL1) THEN
              VOLMIN=VOL1
              INDMIN=IND1
              INDMAX=IND2
           ENDIF
           
           FRAC=DV/VOLMIN !O DV/VOLMIN?
           FRAC=DV/VOL22

           !IF(VOL2 .GT. VOL1) WRITE(*,*) 'WARNING !! VOL2 > VOL1 IN VOLUME, UPDATE MARCA'

           !LENG=RADX(NCOX)-RADX(1)+DX
           !VCELL=(LENG/NCOX)**3. !UPDATE !!!
           VCELL=DX*DY*DZ

           !marking cell in DV
           IX1=INT(((AX(INDX(2))+0.5*DX)-RXX1)/DX+0.5)+1
           IX2=INT(((AX(INDX(3))-0.5*DX)-RXX1)/DX+0.5)+1
           JY1=INT(((AY(INDY(2))+0.5*DY)-RYY1)/DY+0.5)+1
           JY2=INT(((AY(INDY(3))-0.5*DY)-RYY1)/DY+0.5)+1
           KZ1=INT(((AZ(INDZ(2))+0.5*DZ)-RZZ1)/DZ+0.5)+1
           KZ2=INT(((AZ(INDZ(3))-0.5*DZ)-RZZ1)/DZ+0.5)+1
           

           DV2=DV
           IF(FRAC .LE. 0.6) THEN !USE MAXF
           DO IX=IX1, IX2              
              DO JY=JY1, JY2
                 DO KZ=KZ1, KZ2
                    IF(MARCA(IX,JY,KZ) .EQ. 0) THEN
                      MARCA(IX,JY,KZ)= IND1 !check if ind2 is always the smallest
                    ELSE
                      DV2=DV2-VCELL
                    ENDIF 
                 ENDDO
              ENDDO
           ENDDO
           ENDIF

         END SUBROUTINE VOLUME

!*----------------------------------------------------------------------* 


!************************************************************* 
!* Esta subrutina ordena en orden descendente los autovalores*
!* obtenidos por JACOBI  (ver diagonalizando2.f)             *
!************************************************************* 
       SUBROUTINE SORT(D,N,NP)

       IMPLICIT NONE

       integer N,NP,I,J,K
       real*4  D(NP)
       real*4 P

       DO I=1,N-1
       K=I
       P=D(I)
       DO J=I+1, N
               IF(D(J).GE.P) THEN
                  K=J
                  P=D(J)
               END IF
       END DO
       IF(K.NE.I) THEN
       D(K)=D(I)
       D(I)=P

       END IF

       END DO
       RETURN

     END SUBROUTINE SORT



!*************************************************************
!* This subroutine computes all eigenvalues and eigenvectors *
!* of a real symmetric square matrix A(N,N). On output, ele- *
!* ments of A above the diagonal are destroyed. D(N) returns *
!* the eigenvalues of matrix A. V(N,N) contains, on output,  *
!* the eigenvectors of A by columns. THe normalization to    *
!* unity is made by main program before printing results.    *
!* NROT returns the number of Jacobi matrix rotations which  *
!* were required.                                            *
!* --------------------------------------------------------- *
!* Ref.:"NUMERICAL RECIPES, Cambridge University Press, 1986,*
!*       chap. 11, pages 346-348" [BIBLI 08].                *
!*************************************************************
       Subroutine JACOBI(A,N,D,NROT)

       implicit none
       integer N,NROT,ip,iq,ialloc,i,j
       real*4  A(1:N,1:N),D(1:N)  
       real*4, pointer :: B(:), Z(:)
       real*4  c,g,h,s,sm,t,tau,theta,tresh

       allocate(B(1:100))   !,stat=ialloc)
       allocate(Z(1:100))   !,stat=ialloc)
       do ip=1, N
         B(ip)=A(ip,ip)
         D(ip)=B(ip)
         Z(ip)=0.d0
       end do
       NROT=0
       do i=1, 50
         sm=0.d0
         do ip=1, N-1           !sum off-diagonal elements
            do iq=ip+1, N
               sm=sm+ABS(A(ip,iq))
            end do
         end do
         if(sm==0.d0) return    !normal return
         if(i.lt.4) then
            tresh=0.2d0*sm**2
         else
            tresh=0.d0
         end if
         do ip=1, N-1
            do iq=ip+1, N
               g=100.d0*ABS(A(ip,iq))
!      after 4 sweeps,skip the rotation if the off-diag element is small      

       if((i.gt.4).and.(ABS(D(ip))+g.eq.ABS(D(ip))) &
          .and.(ABS(D(iq))+g.eq.ABS(D(iq)))) then

       A(ip,iq)=0.d0
       else if(ABS(A(ip,iq)).gt.tresh) then
                  h=D(iq)-D(ip)
       if(ABS(h)+g.eq.ABS(h)) then
          t=A(ip,iq)/h
       else
          theta=0.5d0*h/A(ip,iq)
          t=1.d0/(ABS(theta)+SQRT(1.d0+theta**2))
          if(theta.lt.0.d0) t=-t
       end if
       c=1.d0/SQRT(1.d0+t**2)
       s=t*c
          tau=s/(1.d0+c)
       h=t*A(ip,iq)
       Z(ip)=Z(ip)-h
       Z(iq)=Z(iq)+h
       D(ip)=D(ip)-h
       D(iq)=D(iq)+h
       A(ip,iq)=0.d0
       do j=1, ip-1
          g=A(j,ip)
          h=A(j,iq)
          A(j,ip)=g-s*(h+g*tau)
          A(j,iq)=h+s*(g-h*tau)
       end do
       do j=ip+1, iq-1
          g=A(ip,j)
          h=A(j,iq)
          A(ip,j)=g-s*(h+g*tau)
          A(j,iq)=h+s*(g-h*tau)
       end do
       do j=iq+1, N
          g=A(ip,j)
          h=A(iq,j)
          A(ip,j)=g-s*(h+g*tau)
          A(iq,j)=h+s*(g-h*tau)
       end do

       NROT=NROT+1
    end if                   !if ((i.gt.4)...
 end do                    !main iq loop
end do                    !main ip loop      
do ip=1, N
   B(ip)=B(ip)+Z(ip)
       D(ip)=B(ip)
       Z(ip)=0.d0
    end do
 end do                    !main i loop 
!c       pause ' 50 iterations !'
       return
    END 

!     end of file ujacobi.f90

!    SUBROUTINE PARTICLES(NPART, NL, MARCA, NPV)
!    USE COMMONDATA  
!    IMPLICIT NONE
!    INTEGER:: IP, LOW1, I, J, K, IVOID
!    REAL*4 BAS
!
!     LOW1=SUM(NPART(0:NL))
!!-------->  update parallel do !!!!!
!!$OMP PARALLEL DO SHARED(LOW1, RXPA, RYPA, RZPA, RADX1, RADY1, RADZ1, &
!!$OMP     RADXNX, RADYNY, RADZNZ, LENG, DX, DY, DZ, NCOXP, NCOYP, NCOZP, &
!!$OMP     U2DM, U3DM, U4DM), &
!!$OMP PRIVATE(IP, BAS, I, J, K, VX, VY, VZ, IX, JY, KZ, I3, J3, K3), &
!!$OMP REDUCTION(+: U2CO, U3CO, U4CO, WTCO)
!     DO IP=1, LOW1
! 
!     !move by LENG particles outside the grid
!     !   IF(RXPA(IP).LT.RADX1) RXPA(IP)=RXPA(IP)+LENG
!     !   IF(RXPA(IP).GE.RADXNX) RXPA(IP)=RXPA(IP)-LENG
!     !   IF(RYPA(IP).LT.RADY1) RYPA(IP)=RYPA(IP)+LENG
!     !   IF(RYPA(IP).GE.RADYNY) RYPA(IP)=RYPA(IP)-LENG
!     !   IF(RZPA(IP).LT.RADZ1) RZPA(IP)=RZPA(IP)+LENG
!     !   IF(RZPA(IP).GE.RADZNZ) RZPA(IP)=RZPA(IP)-LENG
!
!        !initial position
!        !I,J,K: cell of the particle in the coarse grid
!        BAS=RXPA(IP)                                                    
!        I=INT(((BAS-RADX(1))/DX)+0.49999) + 1                               
!        BAS=RYPA(IP)                                                    
!        J=INT(((BAS-RADY(1))/DY)+0.49999) + 1                               
!        BAS=RZPA(IP)                                                    
!        K=INT(((BAS-RADZ(1))/DZ)+0.49999) + 1   
!
!        IVOID=MARCA(I,J,K)
!        NPV(IVOID)=NPV(IVOID)+1
!        VIP(NPV(IVOID))=IP
!
!
!    END SUBROUTINE PARTICLES
!

!-------------------------------------------------------------------------------------------------


  SUBROUTINE PROFILES_SPH(NPART, NL, RODO, ZETA, NVOID, NVOIDL, XC, YC, ZC, &
       INDICE, INDICEL, UVOID, VOLNEW, FILEO) 
    USE COMMONDATA
    IMPLICIT NONE
    !input variables
    INTEGER:: NVOID, NVOIDL, NL
    INTEGER NPART(0:NL) 
    INTEGER, DIMENSION(NVOID):: UVOID
    INTEGER, DIMENSION(NVOID_MAX)::INDICE, INDICEL
    !INTEGER:: MARCA(NCOX, NCOY, NCOZ)
    REAL*4, DIMENSION(NVOID):: XC, YC, ZC, VOLNEW
    REAL*4:: RODO, ZETA, FRACR
    CHARACTER(LEN=*) :: FILEO
    !local variables
    INTEGER:: IVOID, NPARTBIN, LASTP, IP1, IP2, IP, IPP, IPAR, IPAR1, IPAR2
    INTEGER:: NMAXR, IND,  NPARTV, I, J, K, IRAD, IPV, CONTA
    INTEGER:: NCOXP, NCOYP, NCOZP, LOW1
    INTEGER:: BASX, BASY, BASZ
    INTEGER:: IX, JY, KZ, FLAG_LASTP
    INTEGER:: INIX, FINX, INIY, FINY, INIZ, FINZ
    REAL*4:: TIME_1, TIME_2
    REAL*4:: RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ, LENG
    REAL*4:: DENSC, BAS, XV, YV, ZV, UM, RETE, ROTE, FACT
    REAL*4:: MAXX, MAXY, MAXZ, DISTX, DISTY, DISTZ, REQ
    REAL*4:: VSHELL, VSHELL2, RSHELL, VCIRC, VCIRC2, D
    INTEGER, ALLOCATABLE:: INDP(:), NP0(:)
    REAL*4, ALLOCATABLE:: SDIST(:), MP0(:)
    !output variables
    REAL*4, ALLOCATABLE:: RAD(:), URAD(:), URADC(:)
    INTEGER, PARAMETER:: FLAG_C=1

    
    OPEN(UNIT=50, FILE=FILEO)

!get mass profiles
    NPARTBIN=10000 !--> put in parameters
    FRACR=0.1
    DENSC=2.7E11*0.25*0.73*0.73 ! --> UPDATE USING CORRECT VALUE !! RHOC
    UM=9.1717E18      !en Msol     

!*     background
    ROTE=RODO*(1.0+ZETA)**3
    RETE=RE0/(1.0+ZETA)
    FACT=1.0/(ROTE*RETE**3)
 
    WRITE(*,*) FACT, ROTE, RETE, RODO, RE0

    ALLOCATE(SDIST(NPARTT))
    ALLOCATE(INDP(NPARTT))
    SDIST(:)=0.

    ALLOCATE(NPV(0:NVOID))
    ALLOCATE(VIP(NPARTT))

    DO IVOID=0, NVOID
       NPV(IVOID)=0
    ENDDO
    VIP(:)=0
    CONTA=0

    LOW1=SUM(NPART(0:NL)) !EQUAL TO NPARTT
    LOW1=NPART(0)

    ALLOCATE(NP0(NVOID))
    ALLOCATE(MP0(NVOID))
    NP0=0
    MP0=0.
    DO IP=1, LOW1
           BAS=RXPA(IP)                                                    
           I=INT(((BAS-RADX(1))/DX)+0.49999) + 1                               
           BAS=RYPA(IP)                                                    
           J=INT(((BAS-RADY(1))/DY)+0.49999) + 1                               
           BAS=RZPA(IP)                                                    
           K=INT(((BAS-RADZ(1))/DZ)+0.49999) + 1 
           IND=MARCA(I,J,K)
           IF(MASAP(IP) .LE. 0. ) WRITE(*,*) 'WARNING!!! MASAP=0!', IP, MASAP(IP)
           IF(IND .GT. 0 ) THEN
              NP0(IND)=NP0(IND)+1
              MP0(IND)=MP0(IND)+MASAP(IP)
           ENDIF
        ENDDO

     !DO IVOID=1, NVOIDL
     !   IND=INDICEL(IVOID)
     !   WRITE(*,*) IVOID, IND, VOLNEW(IND), NP0(IND), UM*MP0(IND), UM*MP0(IND)/VOLNEW(IND)/DENSC, &
     !        MP0(IND)/(VOLNEW(IND)*((RETE)**3.))/ROTE
     !ENDDO

    LOW1=NPART(0)
    WRITE(*,*) LOW1, NPARTT,  &
         SUM(MASAP(1:LOW1))/((142.*RETE)**3.)/ROTE, SUM(MASAP(1:npartt))/((142.*RETE)**3.)/ROTE, &
         UM*SUM(MASAP(1:npartt))/(142.**3.)/DENSC


    DO IVOID=1, NVOIDL
       IND=INDICEL(IVOID)

       REQ=(VOLNEW(IND)*3./(4.*PI))**(1./3.)
       RSHELL=1.3*REQ

       write(*,*) ivoid, ind, req, uvoid(ind)
       !FIND CENTER
       IF(FLAG_C .EQ. 1) THEN !center of the void
          XV=XC(IND) !try with different centers
          YV=YC(IND)
          ZV=ZC(IND) 
          !XV=0.5*(INICIOX(IND)+FINALX(IND))
          !YV=0.5*(INICIOY(IND)+FINALY(IND))
          !ZV=0.5*(INICIOZ(IND)+FINALZ(IND))
       ELSE 
          IF(FLAG_C .EQ. 2) THEN !max diverV
          
             BAS=-1. !MAX DIV
             BASX=INICIOX(IND) !COORD OF MAX DIV
             BASY=INICIOY(IND) !COORD OF MAX DIV
             BASZ=INICIOZ(IND) !COORD OF MAX DIV
             DO IX=INICIOX(IND), FINALX(IND)
                DO JY=INICIOY(IND), FINALY(IND)
                   DO KZ=INICIOZ(IND), FINALZ(IND)
                      IF(DIVERCO(IX,JY,KZ) .GT. BAS .AND. &
                           MARCA(IX, JY, KZ) .EQ. IND) THEN !MAX DIV
                         BAS=DIVERCO(IX,JY,KZ)
                         BASX=IX
                         BASY=JY
                         BASZ=KZ
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ELSE IF(FLAG_C .EQ. 3) THEN !min density
             BAS=1000. !MIN DENS
             BASX=INICIOX(IND) !COORD OF MIN DENS
             BASY=INICIOY(IND) !COORD OF MIN DENS
             BASZ=INICIOZ(IND) !COORD OF MIN DENS
             DO IX=INICIOX(IND), FINALX(IND)
                DO JY=INICIOY(IND), FINALY(IND)
                   DO KZ=INICIOZ(IND), FINALZ(IND)
                      IF(U1CO(IX,JY,KZ) .LT. BAS .AND. &
                           MARCA(IX, JY, KZ) .EQ. IND) THEN !MAX DIV
                         BAS=U1CO(IX,JY,KZ)
                         BASX=IX
                         BASY=JY
                         BASZ=KZ
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO

             XV=RADX(BASX)
             YV=RADY(BASY)
             ZV=RADZ(BASZ)
          ENDIF
       ENDIF
     

       write(*,*) 'low1:', low1
       !locate particles in voids
! $OMP PARALLEL DO SHARED(LOW1,RXPA, RYPA, RZPA, XV, YV, ZV,  &
! $OMP     RSHELL &
! $OMP     INIX, FINX, INIY, FINY, INIZ, FINZ, IVOID, VIP), &
! $OMP PRIVATE(IP, D ), &
! $OMP REDUCTION(+: NPV, CONTA)
        DO IP=1, LOW1
           D=SQRT((RXPA(IP)-XV)**2.+(RYPA(IP)-YV)**2.+(RZPA(IP)-ZV)**2. )
           IF(D .LE. RSHELL) THEN
              CONTA=CONTA+1
              NPV(IVOID)=NPV(IVOID)+1 
              VIP(CONTA)=IP
              SDIST(NPV(IVOID))=D
           ENDIF  
        ENDDO   
        NPARTV=NPV(IVOID) 
        write(*,*) npartv, conta
        CALL INDEXX(NPARTV, SDIST, INDP)

        write(*,*) 'sorted'
        NPARTBIN=NPARTV*FRACR
        NMAXR=NPARTV/NPARTBIN+1 !max numb of radial bins
        
       WRITE(15,*) IND, REQ, XV, YV, ZV, NPARTV, INT(NPARTV/NPARTBIN), REQ
       WRITE(50,*) IND, REQ, XV, YV, ZV, NPARTV, INT(NPARTV/NPARTBIN), REQ
        !WRITE(*,*) 'Void :', IVOID, NPARTV, NPARTBIN, NMAXR
 
       ALLOCATE(URAD(NMAXR)) 
        ALLOCATE(URADC(0:NMAXR)) 
        ALLOCATE(RAD(0:NMAXR))
        
        RAD(0)=0.
 
        IRAD=0
        IP1=1
        IP2=NPARTBIN
        URADC(0)=0.
        URADC(:)=0. !cumulative
        
        DO WHILE (IP2 .LE. NPARTV)

           IPAR1=INDP(IP1) !void list
           IPAR2=INDP(IP2)
           !WRITE(*,*) SDIST(IPAR1), SDIST(IPAR2), IP1, IPAR1, IP2, IPAR2

           IRAD=IRAD+1

           IF(IRAD .GT. NMAXR) WRITE(*,*) 'Error in nmaxr!', nmaxr, irad, ind

           URAD(IRAD)=0.
          
           DO IP=IP1, IP2
              IPV=INDP(IP) !index in the void list
              IPAR=IPV+NPV(IVOID-1)
              IPP=VIP(IPAR) !index in the global list
              URAD(IRAD)=URAD(IRAD)+MASAP(IPP)             
              !if(IP .LT. 6 ) WRITE(*,*) MASAP(IPP), MASAP(IPP)*UM
           ENDDO

           URADC(IRAD)=URADC(IRAD-1)+URAD(IRAD) !SUM(URAD(1:IRAD))
           !RAD(IRAD)=0.5*(SQRT(SDIST(INDP(IPAR2)))+RAD(IRAD-1))  !o 0.5*(SQRT(SDIST(IP2) 
           RAD(IRAD)=SDIST(IPAR2)

           VSHELL=4.*PI*(SDIST(IPAR2)**3.-SDIST(IPAR1)**3.)/3.
           VSHELL2=(4.*PI*(SDIST(IPAR2)**3.-SDIST(IPAR1)**3.)/3.)*((RETE)**3.)
           VCIRC=4.*PI*(SDIST(IPAR2)**3.)/3.
           VCIRC2=4.*PI*(SDIST(IPAR2)**3.)*((RETE)**3.)/3.
           WRITE(15,*) IND, RAD(IRAD), RAD(IRAD)/REQ, VCIRC, VCIRC2,  &
                URAD(IRAD)*UM, URADC(IRAD)*UM, &
                URADC(IRAD)*UM/(VCIRC*(RETE**3.))/DENSC-1., &
                URADC(IRAD)*UM/(VCIRC)/DENSC-1., &
                FACT*URADC(IRAD)/VCIRC-1., &
                URADC(IRAD)/VCIRC2/ROTE-1.

           WRITE(50,*) IND, RAD(IRAD), RAD(IRAD)/REQ, URADC(IRAD)/VCIRC2/ROTE-1., &
               URAD(IRAD)/VSHELL2/ROTE-1., URAD(IRAD)*UM/VSHELL/DENSC-1., &
               URAD(IRAD)*UM/VSHELL2/DENSC-1., URAD(IRAD)*UM, URADC(IRAD)*UM, VCIRC, VCIRC/ROTE, & 
               URADC(IRAD)
           !WRITE(15,*) RAD(IRAD), RAD(IRAD)/REQ, VSHELL, sqrt(sdist(ipar1)), sqrt(sdist(ipar2)), &
           !     sqrt(sdist(indp(lastp))), URAD(IRAD)*UM, URADC(IRAD)*UM, &
           !     URAD(IRAD)*UM/VSHELL/DENSC-1., FACT, &      
           !     FACT*URAD(IRAD)/VSHELL-1., FACT*UM*URAD(IRAD)/VSHELL-1.
               
           IP1=IP1+NPARTBIN
           IP2=IP2+NPARTBIN

        ENDDO

       DEALLOCATE(URAD, URADC, RAD)
 
    ENDDO

    CLOSE(50)

  END SUBROUTINE PROFILES_SPH

!-------------------------------------------------------------------------------------------------


   SUBROUTINE PROFILES_SQ(NPART, NL, RODO, ZETA, NVOID, NVOIDL, XC, YC, ZC, &
       INDICE, INDICEL,UVOID, VOLNEW, REQ, EPS, UMEAN, LEV, LEVP, &
       UU, UUG, UUS, NX,NY,NZ, DDX, DDY, DDZ, RX1, RY1, RZ1, &
       NPATCH, PARE, PATCHNX, PATCHNY, & 
       PATCHNZ,PATCHX, PATCHY, PATCHZ, &
       PATCHRX, PATCHRY, PATCHRZ, FILEO) 
    USE COMMONDATA
    IMPLICIT NONE
    !input variables
    INTEGER:: NVOID, NVOIDL, NL, LEV, LEVP
    INTEGER:: NX, NY, NZ
    INTEGER NPART(0:NL) 
    INTEGER, DIMENSION(NVOID):: UVOID !PUT IN COMMONDATA ?
    INTEGER, DIMENSION(NVOID):: INDICE, INDICEL
    !INTEGER:: MARCA(NCOX, NCOY, NCOZ)
    REAL*4, DIMENSION(NX,NY,NZ):: UU, UUG, UUS
    REAL*4:: RX1, RY1, RZ1, DDX, DDY, DDZ
    REAL*4, DIMENSION(NVOID):: XC, YC, ZC, VOLNEW, EPS, UMEAN, REQ
    REAL*4:: RODO, ZETA, FRACR
    INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
    INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
    INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHZ(NPALEV)
    REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
    CHARACTER(LEN=*) :: FILEO
    !local variables
    INTEGER:: IX1, IX2, JY1, JY2, KZ1, KZ2, NRAD, IR
    INTEGER:: IVOID, NPARTBIN, LASTP, IP1, IP2, IP, IPP, IPAR, IPAR1, IPAR2
    INTEGER:: NMAXR, IND,  NPARTV, I, J, K, IRAD, IPV, CONTA, N
    INTEGER:: NCELL, NCELLV, NCELLR, NCELLE, IRAD_RE
    INTEGER:: NCOXP, NCOYP, NCOZP, LOW1
    INTEGER:: BASX, BASY, BASZ
    INTEGER:: IX, JY, KZ, FLAG_LASTP
    INTEGER:: INIX, FINX, INIY, FINY, INIZ, FINZ
    INTEGER:: IXV, JYV, KZV, IXV0, JYV0, KZV0
    INTEGER:: IRR, LEV0, IPA, IPA0, PX, PY, PZ
    REAL*4:: DXR2, DYR2, DZR2
    REAL*4:: TIME_1, TIME_2
    REAL*4:: RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ, LENG
    REAL*4:: DENSC, BAS, XV, YV, ZV, UM, RETE, ROTE, FACT
    REAL*4:: MAXX, MAXY, MAXZ, DISTX, DISTY, DISTZ, REQ0
    REAL*4:: VSHELL, VSHELL2, RSHELL, VCIRC, VCIRC2, D
    REAL*4:: SCALEF, SIZEB
    real*4:: VOL20, VOL10
    INTEGER, ALLOCATABLE:: INDP(:), NP0(:)
    REAL*4, ALLOCATABLE:: SDIST(:), MP0(:)
    !REAL*4, ALLOCATABLE:: UU(:,:,:)
    CHARACTER(LEN=2):: LIR
    CHARACTER(LEN=1):: S
    CHARACTER(LEN=50) :: FILEO1, FILEO2, FILEO3
    !output variables
    REAL*4, ALLOCATABLE:: RAD(:), URAD(:), URADC(:), URADG(:), URADS(:), URAD_DIFF(:)
    INTEGER, PARAMETER:: FLAG_C=1
    REAL*4, PARAMETER:: FCELL=0. !fraction of cells within Re belonging to the void: marca=ind

    !sort voids again

    IR=LEVP !LEVP level for profile computation
    IR=LEV !LEV level used for void search/hierarchy

    !ALLOCATE(UU(NX,NY,NZ))
    !UU(1:NX,1:NY,1:NZ)=U1CO(1:NX,1:NY,1:NZ)

!    DO IR=LEV, LEVP
!
       S='p'
       IF(IR .LT. 0) S='m'
       WRITE(LIR,'(I1)') ABS(IR)
       FILEO1=TRIM(ADJUSTL(FILEO))//'_'//TRIM(ADJUSTL(S))//TRIM(ADJUSTL(LIR))
       WRITE(*,*) 'Profiles @lev=',LIR,' written to: ',FILEO1
       OPEN(UNIT=50, FILE=FILEO1)
       
       FILEO2=TRIM(ADJUSTL(FILEO))//'_bar_'//TRIM(ADJUSTL(S))//TRIM(ADJUSTL(LIR))
       WRITE(*,*) 'Profiles (baryonic)  @lev=',LIR,' written to: ',FILEO2
       OPEN(UNIT=51, FILE=FILEO2)

       FILEO3=TRIM(ADJUSTL(FILEO))//'_st_'//TRIM(ADJUSTL(S))//TRIM(ADJUSTL(LIR))
       WRITE(*,*) 'Profiles (stellar)  @lev=',LIR,' written to: ',FILEO3
       OPEN(UNIT=52, FILE=FILEO3)

       SCALEF=(3./(4.*PI))**(1./3.)
!
       DO IVOID=1, NVOIDL

          IND=INDICEL(IVOID)
          NCELLV=COUNT(MARCA .EQ. IND) !total number of coarse cells in the void

          REQ0=REQ(IND)
          RSHELL=1.3*REQ(IND)
          
          XV=XC(IND) !try with different centers
          YV=YC(IND)
          ZV=ZC(IND) 

          !LEV 1
          IXV=INT(((XV-RX1)/DDX)+0.49999) + 1                               
          JYV=INT(((YV-RY1)/DDY)+0.49999) + 1                               
          KZV=INT(((ZV-RZ1)/DDZ)+0.49999) + 1 
          !LEV 0
          IXV0=INT(((XV-RADX0(1))/DX0)+0.49999) + 1                               
          JYV0=INT(((YV-RADY0(1))/DY0)+0.49999) + 1                               
          KZV0=INT(((ZV-RADZ0(1))/DZ0)+0.49999) + 1 


!       !find patch of the void center
          LEV0=0
          IPA=0
          DO IRR=1, NLEVELS
             DXR2=DX0/(2.**IRR)
             DYR2=DY0/(2.**IRR)
             DZR2=DZ0/(2.**IRR)
             CALL PAFIND(XV,YV,ZV, IRR-1, DXR2, DYR2, DZR2, NLEVELS, NPATCH, PARE, PATCHNX, PATCHNY, & 
                  PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
                  PATCHRX, PATCHRY, PATCHRZ, &
                  IPA0, PX, PY, PZ)
             IF(IPA0 .GT. 0) THEN
             IPA=IPA0
             LEV0=IRR
          ENDIF
          !
       ENDDO


       NRAD=INT(RSHELL/DDX+0.4999)+1

       ! test void: compute profile only if: Ncell(MARCA=IND) ge fcell*Ncell(Re)
       ALLOCATE(URAD(0:NRAD))
       ALLOCATE(URADG(0:NRAD))
       ALLOCATE(URADS(0:NRAD))
       ALLOCATE(URADC(0:NRAD))
       ALLOCATE(URAD_DIFF(0:NRAD))
       ALLOCATE(RAD(0:NRAD))
       
       URAD(:)=0.
       URADG(:)=0.
       URADS(:)=0.
       URAD(:)=0.
       URAD_DIFF(:)=0.
       RAD(:)=0.
       URADC(:)=0.
       N=1
       NCELL=0

       IRAD=1
       URAD(1)=UU(IXV,JYV,KZV)
       URADG(1)=UUG(IXV,JYV,KZV)
       URADS(1)=UUS(IXV,JYV,KZV)
       URAD_DIFF(1)=URAD(1) !DENSITY IN SHELLS
       RAD(1)=SCALEF*DDX
       VOL10=DDX**3. !OR RAD(1)**3

       DO IRAD=2, NRAD

        !IF(IXV-N .LT. 1 .OR. IXV+N .GT. NX .OR. &
        !    JYV-N .LT. 1 .OR. JYV+N .GT. NY .OR. &
        !    KZV-N .LT. 1 .OR. KZV+N .GT. NZ ) EXIT

       SIZEB=(2*N+1)*DDX
       RAD(IRAD)=SCALEF*SIZEB
       VOL20=SIZEB**3.
         
       NCELL=0
       DO IX=IXV-N, IXV+N
          DO JY=JYV-N, JYV+N
             DO KZ=KZV-N, KZV+N
                IX2=IX
                JY2=JY
                KZ2=KZ
                IF(IX .LT. 1) IX2=IX+NX
                IF(IX .GT. NX) IX2=IX-NX
                IF(JY .LT. 1) JY2=JY+NY
                IF(JY .GT. NY) JY2=JY-NY
                IF(KZ .LT. 1) KZ2=KZ+NZ
                IF(KZ .GT. NZ) KZ2=KZ-NZ

                NCELL=NCELL+1
                !URAD(IRAD)=URAD(IRAD)+UU(IX,JY,KZ)     
                !URADG(IRAD)=URADG(IRAD)+UUG(IX,JY,KZ)
                !URADS(IRAD)=URADS(IRAD)+UUS(IX,JY,KZ)
                URAD(IRAD)=URAD(IRAD)+UU(IX2,JY2,KZ2)     
                URADG(IRAD)=URADG(IRAD)+UUG(IX2,JY2,KZ2)
                URADS(IRAD)=URADS(IRAD)+UUS(IX2,JY2,KZ2)

                URAD_DIFF(IRAD)=(URAD(IRAD)*VOL20-URAD(IRAD-1)*VOL10)/(VOL20-VOL10) !/ncell ????

             ENDDO
          ENDDO
       ENDDO

       URAD(IRAD)=URAD(IRAD)/NCELL !DENS(<R)
       URADG(IRAD)=URADG(IRAD)/NCELL
       URADS(IRAD)=URADS(IRAD)/NCELL

       !SIZEB=(IX2-IX1+1)*DDX 

       VOL10=VOL20
       N=N+1
       
    ENDDO !loop on radial bins

    !dark matter + gas
      NRAD=IRAD-1
      WRITE(50,*) IND, REQ0, XV, YV, ZV, NRAD, EPS(IND), UMEAN(IND)
      DO IRAD=1, NRAD
         WRITE(50,*) IRAD, RAD(IRAD), RAD(IRAD)/REQ0, URAD(IRAD)-1.
      ENDDO

      !BARIONYC MATTER (gas)
      NRAD=IRAD-1
      WRITE(51,*) IND, REQ0, XV, YV, ZV, NRAD, EPS(IND), UMEAN(IND)
      DO IRAD=1, NRAD
         WRITE(51,*) IRAD, RAD(IRAD), RAD(IRAD)/REQ0, URADG(IRAD)-1.
      ENDDO

      !stars
      NRAD=IRAD-1
      WRITE(52,*) IND, REQ0, XV, YV, ZV, NRAD, EPS(IND), UMEAN(IND)
      DO IRAD=1, NRAD
         WRITE(52,*) IRAD, RAD(IRAD), RAD(IRAD)/REQ0, URADS(IRAD)-1.
      ENDDO
!
      DEALLOCATE(RAD, URAD, URADG, URADS, URADC, URAD_DIFF)
!      
!
    ENDDO !loop on voids
!    
    !DEALLOCATE(UU)
    CLOSE(50)
    CLOSE(51)
    CLOSE(52)

! ENDDO ! loop on levels

    RETURN

  END SUBROUTINE PROFILES_SQ

!****************************************************************************************************

  SUBROUTINE PROFILES_SQ1(IND, NPART, NL, RODO, ZETA, NVOID, NVOIDL, XC, YC, ZC, &
       INDICE, INDICEL, UVOID,  VOLNEW, EPS, UMEAN, LEV, FILEO) 
    USE COMMONDATA
    IMPLICIT NONE
    !input variables
    INTEGER:: NVOID, NVOIDL, NL, LEV
    INTEGER NPART(0:NL) 
    INTEGER, DIMENSION(NVOID):: UVOID
    INTEGER, DIMENSION(NVOID_MAX)::INDICE, INDICEL
    !INTEGER:: MARCA(NCOX, NCOY, NCOZ)
    REAL*4, DIMENSION(NVOID):: XC, YC, ZC, VOLNEW, EPS, UMEAN
    REAL*4:: RODO, ZETA, FRACR
    CHARACTER(LEN=*) :: FILEO
    !local variables
    INTEGER:: IX1, IX2, JY1, JY2, KZ1, KZ2, NRAD
    INTEGER:: IVOID, NPARTBIN, LASTP, IP1, IP2, IP, IPP, IPAR, IPAR1, IPAR2
    INTEGER:: NMAXR, IND,  NPARTV, I, J, K, IRAD, IPV, CONTA, N
    INTEGER:: NCELL, NCELLV, NCELLR, NCELLE, IRAD_RE
    INTEGER:: NCOXP, NCOYP, NCOZP, LOW1, NX, NY, NZ
    INTEGER:: BASX, BASY, BASZ
    INTEGER:: IX, JY, KZ, FLAG_LASTP
    INTEGER:: INIX, FINX, INIY, FINY, INIZ, FINZ
    INTEGER:: IXV, JYV, KZV
    REAL*4:: RX1, RY1,RZ1, DDX,DDY,DDZ
    REAL*4:: TIME_1, TIME_2
    REAL*4:: RADX1, RADXNX, RADY1, RADYNY, RADZ1, RADZNZ, LENG
    REAL*4:: DENSC, BAS, XV, YV, ZV, UM, RETE, ROTE, FACT
    REAL*4:: MAXX, MAXY, MAXZ, DISTX, DISTY, DISTZ, REQ
    REAL*4:: VSHELL, VSHELL2, RSHELL, VCIRC, VCIRC2, D
    REAL*4:: SCALEF, SIZEB
    INTEGER, ALLOCATABLE:: INDP(:), NP0(:)
    REAL*4, ALLOCATABLE:: SDIST(:), MP0(:)
    REAL*4, ALLOCATABLE:: UU(:,:,:)
    !output variables
    REAL*4, ALLOCATABLE:: RAD(:), URAD(:), URADC(:)
    INTEGER, PARAMETER:: FLAG_C=1
    REAL*4, PARAMETER:: FCELL=0. !fraction of cells within Re belonging to the void: marca=ind

    !sort voids again
    OPEN(UNIT=50, FILE=FILEO)

    IF(LEV .EQ. -1) THEN

       NX=NCOX
       NY=NCOY
       NZ=NCOZ
       RX1=RADX(1)
       RY1=RADY(1)
       RZ1=RADZ(1)
       DDX=DX
       DDY=DY
       DDZ=DZ
       ALLOCATE(UU(NX,NY,NZ))
       UU(1:NX,1:NY,1:NZ)=U1CO(1:NX,1:NY,1:NZ)

     ELSE IF(LEV .EQ. 0) THEN
        
        NX=NHYX
        NY=NHYY
        NZ=NHYZ
        RX1=RADX0(1)
        RY1=RADY0(1)
        RZ1=RADZ0(1)
        DDX=DX0
        DDY=DY0
        DDZ=DZ0
        ALLOCATE(UU(NX,NY,NZ))
        UU(1:NX,1:NY,1:NZ)=U1(1:NX,1:NY,1:NZ)

     ELSE IF(LEV .GT. 0) THEN

        NX=(2**LEV)*NHYX
        NY=(2**LEV)*NHYY
        NZ=(2**LEV)*NHYZ
        DDX=DX0/(2**LEV)
        DDY=DY0/(2**LEV)
        DDZ=DZ0/(2**LEV)
        RX1=RADX0(1)-0.5*DDX
        RY1=RADY0(1)-0.5*DDY
        RZ1=RADZ0(1)-0.5*DDZ
        ALLOCATE(UU(NX,NY,NZ))
        UU(1:NX,1:NY,1:NZ)=UR(1:NX,1:NY,1:NZ)

     ENDIF

    SCALEF=(3./(4.*PI))**(1./3.)
   


    DO IVOID=1, NVOIDL

       IND=INDICEL(IVOID)
       NCELLV=COUNT(MARCA .EQ. IND) !total number of coarse cells in the void

       REQ=(VOLNEW(IND)*3./(4.*PI))**(1./3.)
       RSHELL=1.3*REQ

       XV=XC(IND) !try with different centers
       YV=YC(IND)
       ZV=ZC(IND) 

       !LEV CO --> MARCA defined only with lev CO   
       IXV=INT(((XV-RADX(1))/DX)+0.49999) + 1                               
       JYV=INT(((YV-RADY(1))/DY)+0.49999) + 1                               
       KZV=INT(((ZV-RADZ(1))/DZ)+0.49999) + 1 

       IRAD_RE=INT(REQ/DX+0.4999) !=0 if Re is within IXV, JYV, KZV

       !IF(IXV-IRAD_RE .LT. 1 .OR. IXV+IRAD_RE .GT. NCOX .OR. &
       !     JYV-IRAD_RE .LT. 1 .OR. JYV+IRAD_RE .GT. NCOY .OR. &
       !     KZV-IRAD_RE .LT. 1 .OR. KZV+IRAD_RE .GT. NCOZ ) CYCLE

       !NCELLR=0
       !NCELLE=0      
       !DO IX= IXV-IRAD_RE, IXV+IRAD_RE
       !   DO JY= JYV-IRAD_RE, JYV+IRAD_RE          
       !      DO KZ= KZV-IRAD_RE, KZV+IRAD_RE
       !        NCELLE=NCELLE+1 
       !        IF(MARCA(IX,JY,KZ) .EQ. IND) NCELLR=NCELLR+1                  
       !     ENDDO
       !  ENDDO
      !ENDDO
 
      !NCELLE=(INT(2*IRAD_RE+1))**3 !numb of cells within Re

       ncellr=0
       ncelle=0

      IF(NCELLR .GE. FCELL*NCELLE) THEN ! void suitable for profile

       !LEV    
       IXV=INT(((XV-RX1)/DDX)+0.49999) + 1                               
       JYV=INT(((YV-RY1)/DDY)+0.49999) + 1                               
       KZV=INT(((ZV-RZ1)/DDZ)+0.49999) + 1 


       NRAD=INT(RSHELL/DDX+0.4999)+1

       ! test void: compute profile only if: Ncell(MARCA=IND) ge fcell*Ncell(Re)
       ALLOCATE(URAD(0:NRAD))
       ALLOCATE(URADC(0:NRAD))
       ALLOCATE(RAD(0:NRAD))

       URAD(:)=0.
       RAD(:)=0.
       URADC(:)=0.
       N=1
       NCELL=0

       IRAD=1
       URAD(1)=UU(IXV,JYV,KZV)
       RAD(1)=SCALEF*DDX

       DO IRAD=2, NRAD

        IF(IXV-N .LT. 1 .OR. IXV+N .GT. NX .OR. &
            JYV-N .LT. 1 .OR. JYV+N .GT. NY .OR. &
            KZV-N .LT. 1 .OR. KZV+N .GT. NZ ) EXIT
         
       NCELL=0
       DO IX=IXV-N, IXV+N
          DO JY=JYV-N, JYV+N
             DO KZ=KZV-N, KZV+N
                NCELL=NCELL+1
                URAD(IRAD)=URAD(IRAD)+UU(IX,JY,KZ)                
             ENDDO
          ENDDO
       ENDDO

       URAD(IRAD)=URAD(IRAD)/NCELL !DENS(<R)

       SIZEB=(IX2-IX1+1)*DDX 
       SIZEB=(2*N+1)*DDX
       RAD(IRAD)=SCALEF*SIZEB
       !IF(RAD(IRAD) .GT. RSHELL) WRITE(*,*) 'WARNING!! RAD(IRAD)>RSHELL', RAD(IRAD), IRAD, RSHELL
         
       N=N+1
       
    ENDDO !loop on radial bins

      NRAD=IRAD-1
      WRITE(50,*) IND, REQ, XV, YV, ZV, NRAD, EPS(IND), UMEAN(IND)
      DO IRAD=1, NRAD
         WRITE(50,*) IRAD, RAD(IRAD), RAD(IRAD)/REQ, URAD(IRAD)-1.
      ENDDO

      DEALLOCATE(RAD, URAD, URADC)
      
      ENDIF !NCELLR > FCELL*NCELLV
    ENDDO !loop on voids
    
    DEALLOCATE(UU)
    CLOSE(50)

    RETURN

  END SUBROUTINE PROFILES_SQ1
!-------------------------------------------------------------------------------------------------

  SUBROUTINE VMESH1(IND, IR,  RINIX, RFINX, RINIY, RFINY, RINIZ, RFINZ, XV, YV, ZV, &
       NX0, NY0, NZ0, &
       NPATCH, PARE, PATCHNX, PATCHNY, & 
       PATCHNZ,PATCHX, PATCHY, PATCHZ, &
       PATCHRX, PATCHRY, PATCHRZ, NX, NY, NZ)
    USE COMMONDATA  
    IMPLICIT NONE        
    !input variables:
    INTEGER:: IND, IR, NX0, NY0, NZ0
    REAL*4::  RINIX, RFINX, RINIY, RFINY, RINIZ, RFINZ,XV, YV, ZV
!    INTEGER NLEVELS
    INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
    INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
    INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHZ(NPALEV)
    REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
    !local variables   
    INTEGER:: INIX0, FINX0, INIY0, FINY0, INIZ0, FINZ0, NX, NY, NZ, NPX, NPY, NPZ
    INTEGER:: KX, KY, KZ, IX0, JY0, KZ0, IXR, JYR, KZR, IXRR, JYRR, KZRR, IXR1, JYR1, KZR1
    INTEGER:: IXCO, JYCO, KZCO
    INTEGER:: IPX, IPY, IPZ, LX, LY, LZ, NREF, IX00, JY00, KZ00, IPA0, PX0, PY0, PZ0, NN0, IX000, JY000, KZ000
    REAL*4:: RX0, RY0, RZ0, DXR, DYR, DZR, UMEANR
    REAL*4:: RX1, RY1, RZ1 
    !INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: FLAGR
    !output varibales
    INTEGER, PARAMETER:: NSEC=0
    

!IR=1 !generalizzarla a tutti i livelli

DXR=DX0/(2.**IR)
DYR=DY0/(2.**IR)
DZR=DZ0/(2.**IR)

RX1=RINIX-0.5*DXR !attenzione a non cambiare RX1, RY1, RZ1
RY1=RINIY-0.5*DYR
RZ1=RINIZ-0.5*DZR

!generalize to level IR
!RX1=RADX0(1)-0.5*(DX0/(2.**IR))*(2.**IR-1)
!RY1=RADY0(1)-0.5*(DX0/(2.**IR))*(2.**IR-1)
!RZ1=RADZ0(1)-0.5*(DX0/(2.**IR))*(2.**IR-1)


WRITE(*,*) 'RINI, RFIN LEV 0:', RINIX, RFINX, RINIY, RFINY, RINIZ, RFINZ

INIX0=MAX(1,INT((RINIX-RADX0(1))/DX0+0.49999)+1-NSEC) !level0
FINX0=MIN(NX0,INT((RFINX-RADX0(1))/DX0+0.49999)+1+NSEC)
INIY0=MAX(1,INT((RINIY-RADY0(1))/DY0+0.49999)+1-NSEC)
FINY0=MIN(NY0,INT((RFINY-RADY0(1))/DY0+0.49999)+1+NSEC)
INIZ0=MAX(1,INT((RINIZ-RADZ0(1))/DZ0+0.49999)+1-NSEC)
FINZ0=MIN(NZ0,INT((RFINZ-RADZ0(1))/DZ0+0.49999)+1+NSEC)
 

NX=0
NY=0
NZ=0

!dimensions of the fix mesh within the void at level IR
NX=(2.**(IR))*(FINX0-INIX0+1)
NY=(2.**(IR))*(FINY0-INIY0+1)
NZ=(2.**(IR))*(FINZ0-INIZ0+1)

write(*,*) 'NX0, NY0, NZ0:', NX0, NY0, NZ0
write(*,*) 'NX, NY, NZ:',nx, ny, nz, inix0, finx0, iniy0, finy0, iniz0, finz0

ALLOCATE(UR(NX, NY, NZ))
ALLOCATE(DIVR(NX, NY, NZ)) 
ALLOCATE(FLAGR(NX0,NY0,NZ0)) 

UR(:,:,:)=-2.
DIVR(:,:,:)=0.
FLAGR(:,:,:)=0.


!loop over lev0 cells
NREF=0 !number of cells refined: if NREF/NCELL0 << do not use this level, use the coarser

! $OMP PARALLEL DO SHARED(INIX0, FINX0, INIY0, FINY0, INIZ0, FINZ0, &
! $OMP          FLAGR, RADX0, RADY0, RADZ0, DX0, DY0, DZ0, DXR, DYR, DZR, IR, &
! $OMP          NLEVELS, NPATCH, PARE, PATCHNX, PATCHNY, PATCHNZ,PATCHX, PATCHY, PATCHZ, &
! $OMP          PATCHRX, PATCHRY, PATCHRZ
! $OMP PRIVATE(IX0, JY0, KZ0, RX0, RY0, RZ0, IPA0, PX0, PY0, PZ0)
! $OMP REDUCTION(+: NREF
! $OMP IF(FINX0-INIX0+1 .GE. 10)
DO IX0=INIX0, FINX0
   DO JY0=INIY0, FINY0
      DO KZ0=INIZ0, FINZ0

         IF(FLAGR(IX0,JY0,KZ0) .EQ. 0 ) THEN !cell not yet tested

         !center of the cell (IX0, JY0, KZ0) at level 0 - DXR
         RX0=RADX0(1)+(IX0-1)*DX0-0.5*DXR
         RY0=RADY0(1)+(JY0-1)*DY0-0.5*DYR
         RZ0=RADZ0(1)+(KZ0-1)*DZ0-0.5*DZR

         !RX0=RADX0(1)+(IX0-1)*DX0-0.5*(DX0/(2.**IR))*(2.**IR-1)
         !RY0=RADY0(1)+(JY0-1)*DY0-0.5*(DY0/(2.**IR))*(2.**IR-1)
         !RZ0=RADZ0(1)+(KZ0-1)*DZ0-0.5*(DZ0/(2.**IR))*(2.**IR-1)

         IPA0=0

         write(28,*) IX0, RX0, DX0, DXR, RADX0(1)+(IX0-1)*DX0, RADX0(1), JY0, RY0, DY0, DYR, RADY0(1)+(JY0-1)*DY0,RADY0(1), &
              KZ0, RZ0, DZ0, DZR, RADZ0(1)+(KZ0-1)*DZ0, RADZ0(1)
         CALL PAFIND(RX0, RY0, RZ0, IR-1, DXR, DYR, DZR, NLEVELS, NPATCH, PARE, PATCHNX, PATCHNY, & 
              PATCHNZ,PATCHX, PATCHY, PATCHZ, &
              PATCHRX, PATCHRY, PATCHRZ, &
              IPA0, PX0, PY0, PZ0)

          IF(IPA0 .GT. 0) THEN ! cell is refined


           FLAGR(IX0, JY0, KZ0)=1
 
           IF(2*INT(PX0/2) .eq. PX0) WRITE(*,*) 'WARNING: PX0 is even! (IPA0, PX0, PY0, PZ0):', IPA0, PX0, PY0, PZ0, RX0, &
                RY0, RZ0, IX0, JY0, KZ0, RADX0(1), DX0, DXR, RADX0(1)+(IX0-1)*DX0-0.5*DXR
           IF(2*INT(PY0/2) .eq. PY0) WRITE(*,*) 'WARNING: PY0 is even! (IPA0, PX0, PY0, PZ0):', IPA0, PX0, PY0, PZ0, RX0, &
                RY0, RZ0, IX0, JY0, KZ0
           IF(2*INT(PZ0/2) .eq. PZ0) WRITE(*,*) 'WARNING: PZ0 is even! (IPA0, PX0, PY0, PZ0):', IPA0, PX0, PY0, PZ0, RX0, &
                RY0, RZ0, IX0, JY0, KZ0

           NREF=NREF+1 
           IXR1=2*(IX0-INIX0+1)-1 !coord of the cell in the fix IR level (in the void mesh)
           JYR1=2*(JY0-INIY0+1)-1
           KZR1=2*(KZ0-INIZ0+1)-1
           !IXR1=(2**IR)*(IX0-INIX0+1)-1
           !JYR1=(2**IR)*(JY0-INIY0+1)-1
           !KZR1=(2**IR)*(KZ0-INIZ0+1)-1


           IXCO=INT(IX0/2.-0.4999)+1 !INT((IX0-1)/2)+1
           JYCO=INT(JY0/2.-0.4999)+1 !INT((JY0-1)/2)+1
           KZCO=INT(KZ0/2.-0.4999)+1 !INT((KZ0-1)/2)+1
           UMEANR=0

           DO IXR=IXR1, IXR1+1
              DO JYR=JYR1, JYR1+1
                 DO KZR=KZR1, KZR1+1
                    UMEANR=UMEANR+u11(PX0+IXR-IXR1, PY0+JYR-JYR1, PZ0+KZR-KZR1,ipa0)

                 ENDDO
              ENDDO
           ENDDO

           !IXR1=(2**IR)*(IX0-INIX0+1)-(2**IR-1) !all levels
           !JYR1=(2**IR)*(JY0-INIY0+1)-(2**IR-1)
           !KZR1=(2**IR)*(KZ0-INIZ0+1)-(2**IR-1)

           LX=PATCHX(IPA0)
           LY=PATCHY(IPA0)
           LZ=PATCHZ(IPA0)

           !extend to neighbour cells in the same patch
           NPX=MIN(PATCHNX(IPA0), PX0+NX-IXR1)
           NPY=MIN(PATCHNY(IPA0), PY0+NY-JYR1)
           NPZ=MIN(PATCHNZ(IPA0), PZ0+NZ-KZR1)

           IF(2*INT(NPX/2) .ne. NPX) WRITE(*,*) 'WARNING: NPX is odd!:', NPX, PX0, NX, IXR1, PATCHNX(IPA0)
           IF(2*INT(NPY/2) .ne. NPY) WRITE(*,*) 'WARNING: NPY is odd!:', NPY, PY0, NY, JYR1, PATCHNY(IPA0)
           IF(2*INT(NPZ/2) .ne. NPZ) WRITE(*,*) 'WARNING: NPZ is odd!:', NPZ, PZ0, NZ, KZR1, PATCHNZ(IPA0)

           IXR=IXR1-1
           JYR1=JYR1-1           
           KZR1=KZR1-1
            DO IPX=PX0, NPX
               IXR=IXR+1 
               IX00=INIX0+INT((IXR-1)/(2.**IR))

               JYR=JYR1
               DO IPY=PY0, NPY
                  JYR=JYR+1
                  JY00=INIY0+INT((JYR-1)/(2.**IR))
                  KZR=KZR1
                  DO IPZ=PZ0, NPZ
                     KZR=KZR+1
                     UR(IXR, JYR, KZR)=U11(IPX, IPY, IPZ, IPA0)
                     DIVR(IXR, JYR, KZR)=DIVER(IPX, IPY, IPZ, IPA0)
                     !IX00=INT((IPX+1)/2)+LX-1
                     !JY00=INT((IPY+1)/2)+LY-1
                     !KZ00=INT((IPZ+1)/2)+LZ-1

                     !IX00=INIX0+INT(IXR/2.-0.4999) !INT((IXR-1)/2.)
                     !JY00=INIY0+INT(JYR/2.-0.4999)
                     !KZ00=INIZ0+INT(KZR/2.-0.4999)

                     !IX00=INIX0+INT((IXR-1)/(2.**IR))
                     !JY00=INIY0+INT((JYR-1)/(2.**IR))
                     KZ00=INIZ0+INT((KZR-1)/(2.**IR))

                     FLAGR(IX00, JY00, KZ00)=1
                  ENDDO
               ENDDO
            ENDDO

         ELSE !if 0-level cell is not refined copy values of the cell at the 0-level in the IR-level cell

            FLAGR(IX0, JY0, KZ0)=1
            IXR=2*(IX0-INIX0+1)-1
            JYR=2*(JY0-INIY0+1)-1
            KZR=2*(KZ0-INIZ0+1)-1
            !IXR=(2**IR)*(IX0-INIX0+1)-1
            !JYR=(2**IR)*(JY0-INIY0+1)-1
            !KZR=(2**IR)*(KZ0-INIZ0+1)-1

            DO IXRR=IXR, IXR+1
               DO JYRR=JYR,JYR+1
                  DO KZRR=KZR,KZR+1
                     UR(IXRR,JYRR,KZRR)=U1(IX0, JY0, KZ0)
                     DIVR(IXRR,JYRR,KZRR)=DIVER0(IX0, JY0, KZ0)

                  ENDDO
               ENDDO
            ENDDO
                           
         ENDIF !if ipa0>0

      ELSE
         NREF=NREF+1
      ENDIF !if FLAGR=0
   ENDDO
ENDDO
ENDDO


!chech if all the lev0 celss have been tested:

NN0=0
DO IX0=INIX0, FINX0
   DO JY0=INIY0, FINY0
      DO KZ0=INIZ0, FINZ0
         IF(FLAGR(IX0,JY0, KZ0) .EQ. 0) THEN
            NN0=NN0+1
         ENDIF
         
      ENDDO
   ENDDO
ENDDO

IF(COUNT(FLAGR(INIX0:FINX0,INIY0:FINY0,INIZ0:FINZ0) .EQ. 0) .GT. 0) &
        WRITE(*,*) 'WARNING!!! Routine MESHV: cells not tested exist !!'

WRITE(*,*) 'Mean density at level 1:', SUM(UR-1.)/(8.*NHYX*NHYY*NHYZ)

DEALLOCATE(FLAGR)
RETURN

  END SUBROUTINE VMESH1

!-----------------------------------------------------------------------------------------------

  SUBROUTINE VMESH2(IND, IR,  RINIX, RFINX, RINIY, RFINY, RINIZ, RFINZ, XV, YV, ZV, &
       NX0, NY0, NZ0, &
       NPATCH, PARE, PATCHNX, PATCHNY, & 
       PATCHNZ,PATCHX, PATCHY, PATCHZ, &
       PATCHRX, PATCHRY, PATCHRZ, NX, NY, NZ)
    USE COMMONDATA  
    IMPLICIT NONE        
    !input variables:
    INTEGER:: IND, IR, NX0, NY0, NZ0
    !INTEGER:: MARCA(NCOX, NCOY, NCOZ) 
    REAL*4::  RINIX, RFINX, RINIY, RFINY, RINIZ, RFINZ,XV, YV, ZV
    INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
    INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
    INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHZ(NPALEV)
    REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
    !local variables   
    INTEGER:: INIX0, FINX0, INIY0, FINY0, INIZ0, FINZ0, NX, NY, NZ, NPX, NPY, NPZ
    INTEGER:: KX, KY, KZ, IX0, JY0, KZ0, IXR, JYR, KZR, IXRR, JYRR, KZRR, IXR1, JYR1, KZR1
    INTEGER:: IXCO, JYCO, KZCO, N, IX, JY
    INTEGER:: IPX, IPY, IPZ, LX, LY, LZ, NREF, IX00, JY00, KZ00, IPA0, PX0, PY0, PZ0, NN0
    INTEGER:: IX000, JY000, KZ000
    INTEGER:: N1, N2, N3, IPA, LOW1, LOW2 
    INTEGER:: IXRA, IXRB, JYRA, JYRB, KZRA, KZRB
    REAL*4:: RX0, RY0, RZ0, DXR, DYR, DZR, UMEANR
    REAL*4:: RX1, RY1, RZ1
    REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: UPA

    !output variables
    INTEGER, PARAMETER:: NSEC=0
    

!IR=1 !generalizzarla a tutti i livelli

DXR=DX0/(2.**IR)
DYR=DY0/(2.**IR)
DZR=DZ0/(2.**IR)

!IR=1
!RX1=RINIX-0.5*DXR !attenzione a non cambiare RX1, RY1, RZ1
!RY1=RINIY-0.5*DYR
!RZ1=RINIZ-0.5*DZR

!generic level IR
RX1=RADX0(1)-0.5*DX0*(1.-1./(2**IR))
RY1=RADY0(1)-0.5*DY0*(1.-1./(2**IR))
RZ1=RADZ0(1)-0.5*DZ0*(1.-1./(2**IR))


NX=0
NY=0
NZ=0

!dimensions of the fix mesh within the void at level IR
!NX=(2.**(IR))*(FINX0-INIX0+1)
!NY=(2.**(IR))*(FINY0-INIY0+1)
!NZ=(2.**(IR))*(FINZ0-INIZ0+1)

NX=(2.**(IR))*NX0
NY=(2.**(IR))*NY0
NZ=(2.**(IR))*NZ0

write(*,*) 'NX0, NY0, NZ0:', NX0, NY0, NZ0, IR
write(*,*) 'NX, NY, NZ:',NX, NY, NZ
write(*,*) 'RX1, RY1, RZ1:', RX1, RY1, RZ1

ALLOCATE(UR(NX, NY, NZ))
ALLOCATE(DIVR(NX, NY, NZ)) 
ALLOCATE(MARCAR(NX, NY, NZ)) 
ALLOCATE(FLAGAMR(NX, NY, NZ)) 
ALLOCATE(FLAGR(NX0,NY0,NZ0)) 

UR(:,:,:)=-2.
DIVR(:,:,:)=0.
FLAGR(:,:,:)=0.
MARCAR(:,:,:)=0
FLAGAMR(:,:,:)=0

LOW1=SUM(NPATCH(0:IR-1))+1
LOW2=SUM(NPATCH(0:IR))
write(*,*) 'LOW1, LOW1:', LOW1, LOW2
!write(*,*) 'UR(1,1,1):', UR(1,1,30), UR(1,1,1)

DO IPA=LOW1, LOW2

N1=PATCHNX(IPA)
N2=PATCHNY(IPA)
N3=PATCHNZ(IPA)


IXR1=INT((PATCHRX(IPA)-0.5*DXR-RX1)/DXR+0.5)
JYR1=INT((PATCHRY(IPA)-0.5*DYR-RY1)/DYR+0.5)
KZR1=INT((PATCHRZ(IPA)-0.5*DZR-RZ1)/DZR+0.5)

IF(IPA .EQ. 501) THEN
   IXRA=IXR1+1
   IXRB=IXR1+N1
   JYRA=JYR1+1
   JYRB=JYR1+N1
   KZRA=KZR1+1
   KZRB=KZR1+N1
ENDIF

!IR=1
!IXR1=2*INT((PATCHRX(IPA)-RADX0(1))/DX0)
!JYR1=2*INT((PATCHRY(IPA)-RADY0(1))/DX0)
!KZR1=2*INT((PATCHRZ(IPA)-RADZ0(1))/DX0)

IF(IPA .EQ. 500) THEN
   WRITE(73,*) IPA, N1, N2, N3, PATCHRX(IPA),PATCHRY(IPA), &
     PATCHRZ(IPA), DXR, DYR, DZR, RX1, RY1, RZ1, IXR1, JYR1, KZR1
   OPEN(UNIT=16,FILE='map2', FORM='UNFORMATTED')
   ALLOCATE(UPA(N1,N2,N3))
ELSE IF (IPA .EQ. 501) THEN
   WRITE(73,*) IPA, N1, N2, N3, PATCHRX(IPA),PATCHRY(IPA), &
     PATCHRZ(IPA), DXR, DYR, DZR, RX1, RY1, RZ1, IXR1, JYR1, KZR1
      OPEN(UNIT=16,FILE='map3', FORM='UNFORMATTED')
      ALLOCATE(UPA(N1,N2,N3))
ENDIF

DO IPX=1, N1
DO IPY=1, N2
DO IPZ=1, N3

IXR=IXR1+IPX
JYR=JYR1+IPY
KZR=KZR1+IPZ


!IF(FLAGAMR(IXR, JYR, KZR) .EQ. 0) THEN

UR(IXR, JYR, KZR)=U11(IPX, IPY, IPZ, IPA)
DIVR(IXR, JYR, KZR)=DIVER(IPX, IPY, IPZ, IPA)
FLAGAMR(IXR, JYR, KZR)=1

IX0=INT((IXR-1)/(2.**IR))+1
JY0=INT((JYR-1)/(2.**IR))+1
KZ0=INT((KZR-1)/(2.**IR))+1
IF( IPA .EQ. 500 .or. IPA .EQ. 501) WRITE(73,*) IPA, IPX, IPY, IPZ, U11(IPX, IPY, IPZ, IPA), &
     IXR, JYR, KZR, UR(IXR, JYR, KZR), U1(IX0, JY0, KZ0), IX0, JY0, KZ0


!IR=1
!IX0=INT((IXR-1)/2.)+1
!JY0=INT((JYR-1)/2.)+1
!KZ0=INT((KZR-1)/2.)+1

FLAGR(IX0, JY0, KZ0)=1

IXCO=INT((IX0-1)/2.)+1
JYCO=INT((JY0-1)/2.)+1
KZCO=INT((KZ0-1)/2.)+1
MARCAR(IXR, JYR, KZR)=MARCA(IXCO, JYCO, KZCO) 

!ENDIF

ENDDO
ENDDO
ENDDO


ENDDO

!FLAGR(:,:,:)=0

write(*,*) 'Number of refined cells:',count(flagr .eq. 1)

DO IX0=1, NHYX
   DO JY0=1, NHYY
      DO KZ0=1, NHYZ
         IF(FLAGR(IX0, JY0, KZ0) .EQ. 0 ) THEN
            IXCO=INT((IX0-1)/2.)+1
            JYCO=INT((JY0-1)/2.)+1
            KZCO=INT((KZ0-1)/2.)+1
            !IR=1
            !IXR1=2*(IX0-1)+1
            !JYR1=2*(JY0-1)+1
            !KZR1=2*(KZ0-1)+1
            IXR1=(2**IR)*(IX0-1)+1
            JYR1=(2**IR)*(JY0-1)+1
            KZR1=(2**IR)*(KZ0-1)+1
            N=(2**IR)-1
            DO IXR=IXR1, IXR1+N
               DO JYR=JYR1, JYR1+N
                  DO KZR=KZR1, KZR1+N
                     UR(IXR, JYR, KZR)=U1(IX0, JY0, KZ0)
                     DIVR(IXR, JYR, KZR)=DIVER0(IX0, JY0, KZ0)
                     MARCAR(IXR, JYR, KZR)=MARCA(IXCO, JYCO, KZCO) 
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO
   ENDDO
ENDDO

OPEN(UNIT=9,FILE='test',FORM='UNFORMATTED')
WRITE(9)  (((UR(IXR,JYR,KZR), KZR=1,NZ), JYR=1,NY), IXR=1, NX)  
CLOSE(9) 

WRITE(*,*) 'MARCAR>0:', COUNT(MARCAR .GT. 0), COUNT(MARCAR .GT. 0)/(512.**3.), COUNT(MARCAR .GT. 100000)

DEALLOCATE(FLAGR)
RETURN

  END SUBROUTINE VMESH2

!----------------------------------------------------------------------------------------------------------
 

SUBROUTINE VMESH(IR, NX, NY, NZ, DXR, DYR, DZR, RX1, RY1, RZ1, &
     NX0, NY0, NZ0, &
     NPATCH, PARE, PATCHNX, PATCHNY, & 
     PATCHNZ,PATCHX, PATCHY, PATCHZ, &
     PATCHRX, PATCHRY, PATCHRZ)
    USE COMMONDATA  
    IMPLICIT NONE        
    !input variables:
    INTEGER:: IR, NX0, NY0, NZ0
    INTEGER:: NX, NY, NZ, NXX, NYY, NZZ !NCOX, NCOY, NCOZ --> cells of the mesh used by the voidfinder
    REAL*4:: DXR, DYR, DZR, RX1, RY1, RZ1
    !INTEGER:: MARCA(NXX, NYY, NZZ) 
    !REAL*4::  RINIX, RFINX, RINIY, RFINY, RINIZ, RFINZ,XV, YV, ZV
    INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
    INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
    INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHZ(NPALEV)
    REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
    !local variables   
    INTEGER:: INIX0, FINX0, INIY0, FINY0, INIZ0, FINZ0,  NPX, NPY, NPZ
    INTEGER:: KX, KY, KZ, IX0, JY0, KZ0, IXR, JYR, KZR, IXRR, JYRR, KZRR, IXR1, JYR1, KZR1
    INTEGER:: IXCO, JYCO, KZCO, N, IX, JY
    INTEGER:: IPX, IPY, IPZ, LX, LY, LZ, NREF, IX00, JY00, KZ00, IPA0, PX0, PY0, PZ0, NN0
    INTEGER:: IX000, JY000, KZ000
    INTEGER:: N1, N2, N3, IPA, LOW1, LOW2 
    INTEGER:: IXRA, IXRB, JYRA, JYRB, KZRA, KZRB
    REAL*4:: RX0, RY0, RZ0, UMEANR
    REAL*4:: RX, RY, RZ
    REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: UPA
    !INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: FLAGAMR
    !output varibales
    INTEGER, PARAMETER:: NSEC=0


    WRITE(*,*) 'VMESH: IR=',IR

    !DXR=DX0/(2.**IR)
    !DYR=DY0/(2.**IR)
    !DZR=DZ0/(2.**IR)
    
    !NX=(2.**IR)*NX0
    !NY=(2.**IR)*NY0
    !NZ=(2.**IR)*NZ0

    !RX1=RADX0(1)-0.5*DXR ! valid only for IR=1
    !RY1=RADY0(1)-0.5*DYR
    !RZ1=RADZ0(1)-0.5*DZR

    WRITE(*,*) 'NX0, NY0, NZ0:', NX0, NY0, NZ0
    WRITE(*,*) 'NX, NY, NZ:',NX, NY, NZ

    ALLOCATE(UR(NX, NY, NZ))
    ALLOCATE(UGR(NX, NY, NZ))
    ALLOCATE(USR(NX, NY, NZ))
    ALLOCATE(DIVR(NX, NY, NZ))
    ALLOCATE(FLAGR(NX, NY, NZ))
    !ALLOCATE(FLAGAMR(NX, NY, NZ))
    !ALLOCATE(MARCAR(NX, NY, NZ))

    UR(:,:,:)=0.
    UGR(:,:,:)=0.
    USR(:,:,:)=0.
    DIVR(:,:,:)=0.
    FLAGR(:,:,:)=0
    FLAGAMR(:,:,:)=0
    !MARCAR(:,:,:)=0



    IF(IR .EQ. 0) THEN

       DO KZR=1, NZ
          DO JYR=1, NY
             DO IXR=1, NX

                FLAGAMR(IXR, JYR, KZR)=1
                FLAGR(IXR, JYR, KZR)=1
                UR(IXR, JYR, KZR)= U1(IXR, JYR, KZR)
                UGR(IXR, JYR, KZR)= U1G(IXR, JYR, KZR)
                USR(IXR, JYR, KZR)= U1S(IXR, JYR, KZR)
                DIVR(IXR, JYR, KZR)=DIVER0(IXR, JYR, KZR)

             ENDDO
          ENDDO
       ENDDO

    ELSE

    ALLOCATE(UPA(NX, NY, NZ))
    UPA(:,:,:)=-2.

    LOW1=SUM(NPATCH(0:IR-1))+1
    LOW2=SUM(NPATCH(0:IR))
    loop_pa: DO IPA=LOW1,LOW2

          N1= PATCHNX(IPA)
          N2= PATCHNY(IPA)
          N3= PATCHNZ(IPA)

          IXR1=INT(((PATCHRX(IPA)-0.5*DXR)-RX1)/DXR+0.5)
          JYR1=INT(((PATCHRY(IPA)-0.5*DYR)-RY1)/DYR+0.5)
          KZR1=INT(((PATCHRZ(IPA)-0.5*DZR)-RZ1)/DZR+0.5)


          DO IPZ=1, N3
             DO IPY=1, N2
                DO IPX=1, N1

                   IXR=IXR1+IPX
                   JYR=JYR1+IPY
                   KZR=KZR1+IPZ

                   IF(IXR .LT.1 .OR. IXR .GT. NX) WRITE(*,*) 'WARNING ON IXR!', IXR, IXR1, IPX
                   IF(JYR .LT.1 .OR. JYR .GT. NY) WRITE(*,*) 'WARNING ON JYR!', JYR, JYR1, IPY
                   IF(KZR .LT.1 .OR. KZR .GT. NZ) WRITE(*,*) 'WARNING ON KZR!', KZR, KZR1, IPZ


                   !MEAN VALUE
                   FLAGAMR(IXR, JYR, KZR)=FLAGAMR(IXR, JYR, KZR)+1
                   UR(IXR, JYR, KZR)=UR(IXR, JYR, KZR)+U11(IPX, IPY, IPZ, IPA)
                   UGR(IXR, JYR, KZR)=UGR(IXR, JYR, KZR)+U11G(IPX, IPY, IPZ, IPA)
                   USR(IXR, JYR, KZR)=USR(IXR, JYR, KZR)+U11S(IPX, IPY, IPZ, IPA)
                   DIVR(IXR, JYR, KZR)=DIVR(IXR, JYR, KZR)+DIVER(IPX, IPY, IPZ, IPA)
                   !IF(DIVR(IXR, JYR, KZR) .EQ. 0 ) WRITE(*,*) 'DIVR=0 (1)',DIVER(IPX, IPY, IPZ, IPA), ipa, ipx, ipy, ipz
                   IX0=INT((IXR-1)/2.)+1
                   JY0=INT((JYR-1)/2.)+1
                   KZ0=INT((KZR-1)/2.)+1

                   FLAGR(IX0, JY0, KZ0)=1 ! --> put in commondata

                      !IXCO=INT((IX0-1)/2.)+1
                      !JYCO=INT((JY0-1)/2.)+1
                      !KZCO=INT((KZ0-1)/2.)+1

                      !MARCAR(IXR, JYR, KZR)=MARCA(IXCO, JYCO, KZCO)                   

                     
                ENDDO
             ENDDO
          ENDDO

       ENDDO loop_pa !loop on patches

 
       DO KZR=1, NZ
          DO JYR=1, NY
             DO IXR=1, NX
                IF(FLAGAMR(IXR,JYR,KZR) .GE. 1) THEN
                   UR(IXR, JYR, KZR)=UR(IXR, JYR, KZR)/REAL(FLAGAMR(IXR,JYR,KZR))
                   UGR(IXR, JYR, KZR)=UGR(IXR, JYR, KZR)/REAL(FLAGAMR(IXR,JYR,KZR))
                   USR(IXR, JYR, KZR)=USR(IXR, JYR, KZR)/REAL(FLAGAMR(IXR,JYR,KZR))
                   DIVR(IXR, JYR, KZR)=DIVR(IXR, JYR, KZR)/REAL(FLAGAMR(IXR,JYR,KZR))
                   IF(DIVR(IXR, JYR, KZR) .EQ. 0 ) WRITE(*,*) 'DIVR=0 (2)', IXR, JYR, KZR, FLAGAMR(IXR,JYR,KZR)
                ENDIF
             ENDDO
          ENDDO
       ENDDO


    
       DO KZ0=1, NZ0
          DO JY0=1, NY0
             DO IX0=1, NX0

                IXCO=INT((IX0-1)/2.)+1
                JYCO=INT((JY0-1)/2.)+1
                KZCO=INT((KZ0-1)/2.)+1

                IF(FLAGR(IX0, JY0, KZ0) .EQ. 0 ) THEN
                   
                   IXR1=2*IX0-1
                   JYR1=2*JY0-1
                   KZR1=2*KZ0-1
                   DO KZR=KZR1, KZR1+1
                      DO JYR=JYR1, JYR1+1
                         DO IXR=IXR1, IXR1+1
                            UR(IXR, JYR, KZR)= U1(IX0, JY0, KZ0)
                            UGR(IXR, JYR, KZR)= U1G(IX0, JY0, KZ0)
                            USR(IXR, JYR, KZR)= U1S(IX0, JY0, KZ0)
                            DIVR(IXR, JYR, KZR)=DIVER0(IX0, JY0, KZ0)
                            IF(DIVER0(IX0, JY0, KZ0) .EQ. 0) WRITE(67,*) 'WARNING!!! DIVER0=0', IX0, JY0, KZ0
                            !MARCAR(IXR, JYR, KZR)=MARCA(IXCO, JYCO, KZCO)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO


    ENDIF


  END SUBROUTINE VMESH


!******************************************************************************** 
     SUBROUTINE VOIDFIND_PAR(NUM, LEV, NX, NY, NZ, DDX, DDY, DDZ, RX1, RY1, RZ1,  &
          !MARCA, 
          NVOID1, REQP, DENS_THRE, DENS_THRE2, GRAD_THRE, FLAG_DIV, VMIN0, RMIN0, NVOID, &
          NPATCH, PARE, PATCHNX, PATCHNY, &
          PATCHNZ,PATCHX, PATCHY, PATCHZ, &  
          PATCHRX, PATCHRY, PATCHRZ)          
!********************************************************************************
     USE COMMONDATA!, ONLY:  DIVERCO, FLAGV, U1CO
     IMPLICIT NONE
!input variables
     INTEGER:: NUM, LEV, NX,NY,NZ, FLAG_DIV, NVOID1 
     REAL*4:: DDX, DDY, DDZ, RX1, RY1, RZ1, VMIN0, RMIN0
     REAL*4:: DENS_THRE, DENS_THRE2, GRAD_THRE
     !INTEGER, DIMENSION(NX,NY,NZ):: MARCA
     INTEGER NPATCH(0:NLEVELS), PARE(NPALEV)        
     INTEGER PATCHNX(NPALEV), PATCHNY(NPALEV), PATCHNZ(NPALEV)
     INTEGER PATCHX(NPALEV),  PATCHY(NPALEV),  PATCHZ(NPALEV)
     REAL*4  PATCHRX(NPALEV), PATCHRY(NPALEV), PATCHRZ(NPALEV)
     REAL*4, DIMENSION(NVOID1):: REQP
!local variables
     INTEGER:: NXYZ, II, I1, I, J, K, IVOID, II2, II0
     INTEGER:: L1, IX, JY, KZ, IXX, JYY, KZZ
     INTEGER:: FLAG1, FLAG, FLAGX1, FLAGX2,FLAGY1, FLAGY2, FLAGZ1, FLAGZ2
     INTEGER:: INICIO, FINAL, CENTRO
     REAL*4:: DIJK, VV, REQ0
     REAL*4 DDENSX(NX), DDENSY(NY), DDENSZ(NZ)
     REAL*4 DDENSX0(NHYX), DDENSY0(NHYY), DDENSZ0(NHYZ)
     REAL*4 DDENSXR(NMAX), DDENSYR(NMAY), DDENSZR(NMAZ)
     REAL*4 DDENSX2, DDENSY2, DDENSZ2, DDENS2
     REAL*4, ALLOCATABLE:: DDD(:), DDD2(:)
     REAL*4 UDENSX(NX), DIVX(NX)
     REAL*4 UDENSY(NY), DIVY(NY)
     REAL*4 UDENSZ(NZ), DIVZ(NZ)
     REAL*4 UDENSX0(NHYX), DIVX0(NHYX)
     REAL*4 UDENSY0(NHYY), DIVY0(NHYY)
     REAL*4 UDENSZ0(NHYZ), DIVZ0(NHYZ)
     REAL*4 A,B,C
     REAL*4, DIMENSION(NX,NY,NZ):: DIVERCO2, UCO_2
     INTEGER, ALLOCATABLE:: DDDX(:), DDDY(:), DDDZ(:), &
           DDDX2(:), DDDY2(:), DDDZ2(:), INDICE(:), INDICE2(:)
     INTEGER:: CENTRO0(NX,NY,NZ) !common?
     INTEGER:: FLAG2
     INTEGER:: IX0, JY0, KZ0, IR, N1
     INTEGER:: INICIO0, FINAL0, INICIO1, FINAL1
     INTEGER:: IPA0, IR0, PX, PY, PZ
     REAL*4 DIV0, RX, RY, RZ, DXR, DYR, DZR, DDENS    
     INTEGER:: FLAG_DIV_P, FLAG_GAL_P, NSEC_P, NVOID_MAX_P, NUM_MAX
     INTEGER:: INDV
!parameters
     INTEGER, PARAMETER:: NSEC=1
!output variables
     INTEGER:: NVOID
     INTEGER:: IPRO, OMP_GET_THREAD_NUM, NN
     INTEGER, DIMENSION(NUM):: FLAGX1_P, FLAGX2_P, FLAGY1_P, FLAGY2_P, FLAGZ1_P, FLAGZ2_P

     

     !DDX=DX0*(2.**LEV)
     !DDY=DY0*(2.**LEV)
     !DDZ=DZ0*(2.**LEV)
     !RX1=RADX0(1)-0.5*DX0*(2.**LEV-1)/(2**LEV) 
     !RY1=RADY0(1)-0.5*DY0*(2.**LEV-1)/(2**LEV)
     !RZ1=RADZ0(1)-0.5*DZ0*(2.**LEV-1)/(2**LEV)

     IF(FLAG_GAS .EQ. 1) THEN !USE DIVERG AS DIVERCO
        DO IX=1, NX
           DO JY=1, NY
              DO KZ=1, NZ
                 DIVERCO(IX, JY, KZ)=DIVERGCO(IX, JY, KZ)
              ENDDO
           ENDDO
        ENDDO
     ENDIF

     !max divergence
     NXYZ=0
     DIVERCO2(1:NX,1:NY,1:NZ)= DIVERCO(1:NX,1:NY,1:NZ)*REAL(FLAGV(1:NX,1:NY,1:NZ)) 
 
    IF(FLAG_VEL == 0) THEN !if velocity field is not given I rank cells according to the inverse density (low density cells come sfirst)
        DIVERCO(1:NX,1:NY,1:NZ)=1./U1CO(1:NX,1:NY,1:NZ)
        DIVERCO2(1:NX,1:NY,1:NZ)= DIVERCO(1:NX,1:NY,1:NZ)*REAL(FLAGV(1:NX,1:NY,1:NZ)) 
     ENDIF


     !new diverco is =diverco if diverco>0 and 0 elsewhere
     NXYZ=COUNT(DIVERCO2(1:NX,1:NY,1:NZ).GT.0)


     ALLOCATE(DDD(NXYZ))
     ALLOCATE(DDD2(NXYZ))
     ALLOCATE(DDDX2(NXYZ))
     ALLOCATE(DDDY2(NXYZ))
     ALLOCATE(DDDZ2(NXYZ))
     ALLOCATE(DDDX(NXYZ))
     ALLOCATE(DDDY(NXYZ))
     ALLOCATE(DDDZ(NXYZ))
     ALLOCATE(INDICE2(NXYZ))
     ALLOCATE(INDICE(NXYZ)) 
      
!$OMP PARALLEL DO SHARED(NXYZ,DDD, DDD2,DDDX2,DDDZ2,DDDX,DDDY,DDDZ,INDICE, INDICE2), &
!$OMP     PRIVATE(I)
     DO I=1, NXYZ
        DDD(I)=0.0 !diverco collapsed
        DDD2(I)=0.0 !DDD ordered from largest to smallest, temp
        DDDX2(I)=0 !x index temp
        DDDY2(I)=0 !y index temp
        DDDZ2(I)=0 !z index temp
        DDDX(I)=0 !x index
        DDDY(I)=0 !y index
        DDDZ(I)=0 !z index
        INDICE2(I)=0 !indices of sorted diver (1=smallest, NXYZ=largest)
        INDICE(I)=0 !indices of sorted diver (1=largest, NXYZ=smallest)
     ENDDO


!*parallel 
!variables for parallelization:
     FLAG_GAL_P=FLAG_GAL
     FLAG_DIV_P=FLAG_DIV
     IF(FLAG_VEL == 0) FLAG_DIV_P=0

     NSEC_P=NSEC
     NVOID_MAX_P=NVOID_MAX

!$OMP PARALLEL DO SHARED(NVOID_MAX_P,INICIOX, FINALX, INICIOY, FINALY, &
!$OMP     INICIOZ, FINALZ, ICX, ICY, ICZ), &
!$OMP     PRIVATE(I)
     DO I=1, NVOID_MAX_P
        INICIOX(I)=0
        FINALX(I)=0
        INICIOY(I)=0
        FINALY(I)=0
        INICIOZ(I)=0
        FINALZ(I)=0
        ICX(I)=0
        ICY(I)=0
        ICZ(I)=0
     ENDDO

     DDENSX(1:NX)=0.0
     DDENSY(1:NY)=0.0
     DDENSZ(1:NZ)=0.0
     DIVX(1:NX)=0.0
     DIVY(1:NY)=0.0
     DIVZ(1:NZ)=0.0
     DDENSX0(1:NHYX)=0.0
     DDENSY0(1:NHYY)=0.0
     DDENSZ0(1:NHYZ)=0.0
     DIVX0(1:NHYX)=0.0
     DIVY0(1:NHYY)=0.0
     DIVZ0(1:NHYZ)=0.0

     !collapse DIVERCO in DD
     II=0
     DO K=1,NZ                                                        
     DO J=1,NY                                                        
     DO I=1,NX
       DIJK=0.0
       DIJK=DIVERCO2(I,J,K) !max diver
       !DIJK=UCO_2(I,J,K) !min dens
       IF(DIJK .GT. 0.0) THEN ! div
        !IF(DIJK .LT. DENS_THRE+1.) THEN ! density   
        II=II+1
        DDD(II)=DIJK
        DDDX2(II)=I
        DDDY2(II)=J
        DDDZ2(II)=K
       ENDIF
     ENDDO
     ENDDO
     ENDDO
     IF (II.NE.NXYZ) THEN
        WRITE(*,*) 'WARNING,NXYZ', II, NXYZ
        return !STOP
     END IF

     !!!  ORDER DIVERGENCE : from largest to smallest: I want centers in cells of max diver
     CALL INDEXX(NXYZ,DDD(1:NXYZ),INDICE2(1:NXYZ)) !sort diverco from smallest to largest


!max diver: reverse order
       I1=NXYZ

 !*$OMP PARALLEL DO SHARED(NXYZ,DDD, DDD2,DDDX2,INDICE, INDICE2), &
 !*$OMP     PRIVATE(I), &
 !*$OMP REDUCTION(-:I1)
       DO I=1,NXYZ
        DDD2(I1)=DDD(INDICE2(I)) !diverco ordered from largest to smallest
        INDICE(I1)=INDICE2(I)
        I1=I1-1
       END DO 

!min dens
!density ordered from smallest to largest
!!        I1=1
!! !*$OMP PARALLEL DO SHARED(NXYZ,DDD, DDD2,DDDX2,INDICE, INDICE2), &
!! !*$OMP     PRIVATE(I), &
!! !*$OMP REDUCTION(-:I1)
!!      DO I=1,NXYZ
!!        DDD2(I1)=DDD(INDICE2(I)) !diverco ordered from largest to smallest
!!        INDICE(I1)=INDICE2(I)
!!        I1=I1+1
!!       END DO 


!$OMP PARALLEL DO SHARED(NXYZ,DDD, DDDX,DDDY,DDDZ), &
!$OMP     PRIVATE(I)
       DO I=1,NXYZ
        DDD(I)=0.0
        DDDX(I)=0
        DDDY(I)=0
        DDDZ(I)=0
       END DO
 
 
!$OMP PARALLEL DO SHARED(NXYZ,DDD, DDD2,DDDX2,DDDY2,DDDZ2,&
!$OMP     DDDX,DDDY,DDDZ,INDICE), &
!$OMP     PRIVATE(I)
       DO I=1,NXYZ
        DDD(I)=DDD2(I)
        DDDX(I)=DDDX2(INDICE(I))
        DDDY(I)=DDDY2(INDICE(I))
        DDDZ(I)=DDDZ2(INDICE(I))
       END DO

       DEALLOCATE(DDD2, DDDX2, DDDY2, DDDZ2, INDICE2)


!%%%%%%%%%%%   loop on potential centers, starting from cells with max diver    %%%%%%%%%%%%

       NVOID=0

       DO L1=1, NXYZ
       
       IX=DDDX(L1)    !!!celda centro
       JY=DDDY(L1)
       KZ=DDDZ(L1)


       FLAG2=0

       !If cell (IX, JY, KZ) falls within another void I exclude it 
       FLAG1=0 !

       DO IVOID=1, NVOID 
          IF(IX .GE. INICIOX(IVOID) .AND. IX .LE. FINALX(IVOID) .AND. &
               JY .GE. INICIOY(IVOID) .AND. JY .LE. FINALY(IVOID) .AND. &
               KZ .GE. INICIOZ(IVOID) .AND. KZ .LE. FINALZ(IVOID)) FLAG1=1         

          IF(FLAG1==1) EXIT
       ENDDO
       


       !IF(IX==1 .OR. JY==1 .OR. KZ==1 .OR. &
       !     IX==NX .OR. JY==NY .OR. KZ==NZ) FLAG1=1 !use security border: NSEC
       !new 24/09
       IF(IX .LE. NSEC+1 .OR. JY .LE. NSEC+1 .OR. KZ .LE. NSEC+1 .OR. &
            IX .GE. NX-NSEC .OR. JY .GE. NY-NSEC .OR. KZ .GE. NZ-NSEC) FLAG1=1 !use security border: NSEC

       !qua dovrei usare i bordi del vuoto...
       !ind=marcav(ix, jy, kz)
       !if(ix .eq. inicioxv(ind) ...)


       IF(FLAG1==0) THEN !cell of the center does not fall in other voids and not at the edge


       !CENTRO0(IX,JY,KZ)=1
       NVOID=NVOID+1
       IF(NVOID>NVOID_MAX) THEN
          WRITE(*,*) 'NVOID > NVOID_MAX!! Increase NVOID_MAX '
          RETURN
       ENDIF


!*      VECTORIZANDO...       
       UDENSX(:)=U1CO(:,JY,KZ)  
       UDENSY(:)=U1CO(IX,:,KZ) 
       UDENSZ(:)=U1CO(IX,JY,:) 
       DIVX(:)=DIVERCO(:,JY,KZ) !non e' DIVERCO2 ??
       DIVY(:)=DIVERCO(IX,:,KZ)
       DIVZ(:)=DIVERCO(IX,JY,:)

       IF(IX==1 .OR. JY==1 .OR. KZ==1) WRITE(*,*) 'IX, JY, KZ:', IX, JY, KZ, NVOID


       FLAGX1=0
       FLAGX2=0
       FLAGY1=0
       FLAGY2=0
       FLAGZ1=0
       FLAGZ2=0
       FLAG=0
       INICIOX(NVOID)=IX
       FINALX(NVOID)=IX
       INICIOY(NVOID)=JY
       FINALY(NVOID)=JY
       INICIOZ(NVOID)=KZ
       FINALZ(NVOID)=KZ
       
       FLAGX1_P(:)=0
       FLAGX2_P(:)=0
       FLAGY1_P(:)=0
       FLAGY2_P(:)=0
       FLAGZ1_P(:)=0
       FLAGZ2_P(:)=0


       INDV=MARCAP(IX, JY, KZ) !ID of the parent void when looking for subvoids
                              ! =0 when looking for voids

       REQ0=0.
       IF(INDV .GT. 0) THEN ! SUBVOID SEARCH ONLY: INDV ALWAYS =0 FOR VOIDS
          REQ0=0.73*REQP(INDV)
          IF(REQ0 .LT. RMIN0) THEN
             FLAG=1 !SELECT ONLY LARGE VOIDS
             !CYCLE
          ENDIF
       ENDIF
       !WRITE(68,*) IX,JY,KZ,INDV, REQ0, FLAG
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%   start expanding the protovoids in all the directions    %%%%%%%%%%%%%%%%%%%

       DO WHILE(FLAG==0) !if FLAG=1 means that void can not be expanded, 
                         ! because does not satisfy the void conditions in at least one face
                          
!!!!!    expand along  X  

!*parallel 


          IPRO=1
          FLAGX1_P(:)=FLAGX1
          FLAGX2_P(:)=FLAGX2

          NN=FINALZ(NVOID)-INICIOZ(NVOID)+1

          NUM_MAX=1 !! serial !!
          IF(NN .GE. NUM) NUM_MAX=NUM !! parallel !!

!$OMP PARALLEL DO SHARED(NUM, NVOID, INICIOZ, FINALZ, INICIOY, FINALY, INICIOX, FINALX, &
!$OMP            U1CO, DIVERCO, FLAG_SUB, FLAGX1_P, FLAGX2_P, INDV, MARCAP, &
!$OMP            NX, NSEC_P, DDX, GRAD_THRE, FLAG_GAL_P, FLAG_DIV_P, DENS_THRE2, INICIO, NN), &
!$OMP            PRIVATE(KZZ, JYY, UDENSX, DIVX, II, DDENSX, II2, DDENS2, &
!$OMP            FLAGX2, FLAGX1, IPRO), & 
!$OMP          IF (NN .GE. NUM)
          OUTX:DO KZZ=INICIOZ(NVOID), FINALZ(NVOID)

!*parallel 
             IF(NUM .EQ. 1) THEN
                IPRO=1
             ELSE
               IPRO=1 !! serial !!
               !IPRO=OMP_GET_THREAD_NUM()+1 !! parallel !!
             ENDIF

             FLAGX1=FLAGX1_P(IPRO)
             FLAGX2=FLAGX2_P(IPRO)


          DO JYY=INICIOY(NVOID), FINALY(NVOID)


!*parallel 
             IF(FLAGX1==0 .OR. FLAGX2==0) THEN

         !*      VECTORIZANDO...       
                   UDENSX(:)=U1CO(:,JYY,KZZ)  
                   DIVX(:)=DIVERCO(:,JYY,KZZ)

                !+X
                IF(FINALX(NVOID) .GE. NX-NSEC_P)  THEN             
                   FLAGX2=1 !USE NSEC
                   FLAG2=1 ! void will be excluded
                ELSE
                   II=FINALX(NVOID)+1
                   IF(FLAG_SUB(II, JYY, KZZ) .EQ. 0.OR. MARCAP(II,JYY,KZZ) .NE. INDV) FLAGX2=1 !useful only for subvoid search    
                   ! if looking for voids FLAG_SUB should be always =1 
                ENDIF

                IF(FLAGX2 .EQ. 0 ) THEN
                   
                   !II=FINALX(NVOID)+1

                   DDENSX(II)=(UDENSX(II+1)-UDENSX(II-1))/(2.D0*DDX) !use scalar

                   IF(DDENSX(II)  .GE. GRAD_THRE) FLAGX2=1
                   IF(FLAG_GAL_P==1 .AND. FLAGX2==1 .AND. II .LT. NX-1) THEN !check gradient in the next cell: use NSEC?
                      II2=II+1
                      DDENS2=(UDENSX(II2+1)-UDENSX(II2-1))/(2.D0*DDX)
                      IF(DDENS2 .LT. GRAD_THRE) FLAGX2=0 
                   ENDIF
                   IF(FLAG_DIV_P==1 .AND. DIVX(II) .LT. 0.) FLAGX2=1
                   IF(UDENSX(II) .GE. DENS_THRE2+1.) FLAGX2=1

                ENDIF

          
                !-X        
                IF(INICIOX(NVOID) .LE. 1+NSEC_P)  THEN
                   FLAGX1=1
                   FLAG2=1
                ELSE
                   II=INICIOX(NVOID)-1
                   IF(FLAG_SUB(II, JYY, KZZ) .EQ. 0 .OR. MARCAP(II,JYY,KZZ) .NE. INDV) FLAGX1=1 
                ENDIF

                IF(FLAGX1 .EQ. 0 ) THEN

                   !II=INICIOX(NVOID)-1

                   DDENSX(II)=(UDENSX(II-1)-UDENSX(II+1))/(2.D0*DDX)

                   IF(DDENSX(II) .GE. GRAD_THRE) FLAGX1=1

                   IF(FLAG_GAL_P==1 .AND. FLAGX1==1 .AND. II .GT. 2) THEN !check gradient in the next cell
                      II2=II-1
                      DDENS2=(UDENSX(II2-1)-UDENSX(II2+1))/(2.D0*DDX)
                      IF(DDENS2 .LT. GRAD_THRE) FLAGX1=0 
                   ENDIF
                   IF(FLAG_DIV_P==1 .AND. DIVX(II) .LT. 0.) FLAGX1=1
                   IF(UDENSX(II) .GE. DENS_THRE2+1.) FLAGX1=1
                   !IF(FLAGX1==0) INICIO=II
                ENDIF

                !*IF(FLAGX1==1 .AND. FLAGX2==1) EXIT OUTX !*parallel: remove this line for parallel run

!*parallel 
             ENDIF

    ENDDO
!*parallel 
    FLAGX1_P(IPRO)=FLAGX1
!*parallel 
    FLAGX2_P(IPRO)=FLAGX2
 ENDDO OUTX

!*parallel 
 FLAGX1=MAXVAL(FLAGX1_P(1:NUM_MAX))
!*parallel 
 FLAGX2=MAXVAL(FLAGX2_P(1:NUM_MAX))

 
 
!!!!!    expand along  Y
          IPRO=1
          FLAGY1_P(:)=FLAGY1
          FLAGY2_P(:)=FLAGY2

          NN=FINALZ(NVOID)-INICIOZ(NVOID)+1
          NUM_MAX=1 !! serial !!
          IF(NN .GE. NUM) NUM_MAX=NUM !! parallel !!

!$OMP PARALLEL DO SHARED(NUM, NVOID, INICIOZ, FINALZ, INICIOY, FINALY, INICIOX, FINALX, &
!$OMP            U1CO, DIVERCO, FLAG_SUB, FLAGY1_P, FLAGY2_P, INDV, MARCAP, &
!$OMP            NX, NSEC_P, DDY, GRAD_THRE, FLAG_GAL_P, FLAG_DIV_P, DENS_THRE2, INICIO, NN), &
!$OMP            PRIVATE(KZZ, IXX, UDENSY, DIVY, II, DDENSY, II2, DDENS2, &
!$OMP            FLAGY2, FLAGY1, IPRO) , & 
!$OMP          IF (NN .GE. NUM)
       OUTY: DO KZZ=INICIOZ(NVOID), FINALZ(NVOID)

!*parallel 
             IF(NUM .EQ. 1) THEN
                IPRO=1
             ELSE
                IPRO=1 !! serial !!
                !IPRO=OMP_GET_THREAD_NUM()+1 !! parallel !!
             ENDIF

             FLAGY1=FLAGY1_P(IPRO)
             FLAGY2=FLAGY2_P(IPRO)

          DO IXX=INICIOX(NVOID), FINALX(NVOID)
             
!*parallel 
             IF(FLAGY1==0 .OR. FLAGY2==0) THEN

!*      VECTORIZANDO...       
          UDENSY(:)=U1CO(IXX,:,KZZ) 
          DIVY(:)=DIVERCO(IXX,:,KZZ)
                    

          !+Y
          IF(FINALY(NVOID) .GE. NY-NSEC_P)  THEN
             FLAGY2=1
             FLAG2=1
          ELSE
             II=FINALY(NVOID)+1
             IF(FLAG_SUB(IXX, II, KZZ) .EQ. 0.OR. MARCAP(IXX, II,KZZ) .NE. INDV) FLAGY2=1 
          ENDIF

          IF(FLAGY2 .EQ. 0 ) THEN

          !II=FINALY(NVOID)+1

          DDENSY(II)=(UDENSY(II+1)-UDENSY(II-1))/(2.D0*DDY)

          IF(DDENSY(II) .GE. GRAD_THRE) FLAGY2=1

          IF(FLAG_GAL_P==1 .AND. FLAGY2==1 .AND. II .LT. NY-1) THEN !check gradient in the next cell
                II2=II+1
                DDENS2=(UDENSY(II2+1)-UDENSY(II2-1))/(2.D0*DDY)
                IF(DDENS2 .LT. GRAD_THRE) FLAGY2=0 
          ENDIF
          IF(FLAG_DIV_P==1 .AND. DIVY(II) .LT. 0.) FLAGY2=1
          IF(UDENSY(II) .GE. DENS_THRE2+1.) FLAGY2=1
          !IF(FLAGY2==0) FINAL=II   
          ENDIF

          !-Y
          IF(INICIOY(NVOID) .LE. 1+NSEC_P)  THEN
             FLAGY1=1
             FLAG2=1
          ELSE
             II=INICIOY(NVOID)-1
             IF(FLAG_SUB(IXX, II, KZZ) .EQ. 0.OR. MARCAP(IXX,II,KZZ) .NE. INDV) FLAGY1=1 
          ENDIF
          
          IF(FLAGY1 .EQ. 0 ) THEN

          !II=INICIOY(NVOID)-1
 
          DDENSY(II)=(UDENSY(II-1)-UDENSY(II+1))/(2.D0*DDY)

          IF(DDENSY(II) .GE. GRAD_THRE) FLAGY1=1

          IF(FLAG_GAL_P==1 .AND. FLAGY1==1 .AND. II .GT. 2) THEN !check gradient in the next cell
                II2=II-1
                DDENS2=(UDENSY(II2-1)-UDENSY(II2+1))/(2.D0*DDY)
                IF(DDENS2 .LT. GRAD_THRE) FLAGY1=0 
          ENDIF
          IF(FLAG_DIV_P==1 .AND. DIVY(II) .LT. 0.) FLAGY1=1
          IF(UDENSY(II) .GE. DENS_THRE2+1.) FLAGY1=1
          !IF(FLAGY1==0) INICIO=II  

          ENDIF
          
          !IF(FLAGY1==1 .AND. FLAGY2==1) EXIT OUTY
       ENDIF !IF (FLAGY1=0 .OR. FLAGY2-0)

       ENDDO

!*parallel 
    FLAGY1_P(IPRO)=FLAGY1
!*parallel 
    FLAGY2_P(IPRO)=FLAGY2

    ENDDO OUTY

!*parallel 
 FLAGY1=MAXVAL(FLAGY1_P(1:NUM_MAX))
!*parallel 
 FLAGY2=MAXVAL(FLAGY2_P(1:NUM_MAX))

!!!!!    expand along  Z
          IPRO=1
          FLAGZ1_P(:)=FLAGZ1
          FLAGZ2_P(:)=FLAGZ2

          NN=FINALY(NVOID)-INICIOY(NVOID)+1
          NUM_MAX=1 !! serial !!
          IF(NN .GE. NUM) NUM_MAX=NUM !! parallel !!

          IF(FLAGZ1.EQ.0 .OR. FLAGZ2 .EQ.0) THEN

!$OMP PARALLEL DO SHARED(NUM, NVOID, INICIOZ, FINALZ, INICIOY, FINALY, INICIOX, FINALX, &
!$OMP            U1CO, DIVERCO, FLAG_SUB, FLAGZ1_P, FLAGZ2_P,INDV, MARCAP,&
!$OMP            NZ, NSEC_P, DDZ, GRAD_THRE, FLAG_GAL_P, FLAG_DIV_P, DENS_THRE2, INICIO, NN), &
!$OMP            PRIVATE(JYY, IXX, UDENSZ, DIVZ, II, DDENSZ, II2, DDENS2, &
!$OMP            FLAGZ2, FLAGZ1, IPRO), & 
!$OMP          IF (NN .GE. NUM)
       OUTZ: DO JYY=INICIOY(NVOID), FINALY(NVOID)
!*parallel 
             IF(NUM .EQ. 1) THEN
                IPRO=1
             ELSE
                IPRO=1 !! serial !!
                !IPRO=OMP_GET_THREAD_NUM()+1 !! parallel !!
             ENDIF

             FLAGZ1=FLAGZ1_P(IPRO)
             FLAGZ2=FLAGZ2_P(IPRO)
             

             DO IXX=INICIOX(NVOID), FINALX(NVOID)
!*parallel 
                IF(FLAGZ1==0 .OR. FLAGZ2==0) THEN

!*      VECTORIZANDO...       
                   UDENSZ(:)=U1CO(IXX,JYY,:) 
                   DIVZ(:)=DIVERCO(IXX,JYY,:)
                   
 
          !+Z : FLAGZ2
                   IF(FINALZ(NVOID) .GE. NZ-NSEC_P)  THEN
                      FLAGZ2=1
                      FLAG2=1
                   ELSE
                      II=FINALZ(NVOID)+1
                      IF(FLAG_SUB(IXX, JYY, II) .EQ. 0.OR. MARCAP(IXX, JYY,II) .NE. INDV) FLAGZ2=1
                   ENDIF


                   IF(FLAGZ2 .EQ. 0 ) THEN

                      !II=FINALZ(NVOID)+1

                      DDENSZ(II)=(UDENSZ(II+1)-UDENSZ(II-1))/(2.D0*DDZ)

                      IF(DDENSZ(II) .GE. GRAD_THRE) FLAGZ2=1

                      IF(FLAG_GAL_P==1 .AND. FLAGZ2==1 .AND. II .LT. NZ-1) THEN !check gradient in the next cell
                         II2=II+1
                         DDENS2=(UDENSZ(II2+1)-UDENSZ(II2-1))/(2.D0*DDZ)
                         IF(DDENS2 .LT. GRAD_THRE) FLAGZ2=0 
                      ENDIF

                      IF(FLAG_DIV_P==1 .AND. DIVZ(II) .LT. 0.) FLAGZ2=1

                      IF(UDENSZ(II) .GE. DENS_THRE2+1.) FLAGZ2=1

                      !IF(FLAGZ2==0) FINAL=II  

                   ENDIF


            !-Z : FLAGZ1
                   IF(INICIOZ(NVOID) .LE. 1+NSEC_P)  THEN
                      FLAGZ1=1
                      FLAG2=1
                   ELSE
                      II=INICIOZ(NVOID)-1
                      IF(FLAG_SUB(IXX, JYY, II) .EQ. 0.OR. MARCAP(IXX,JYY,II) .NE. INDV) FLAGZ1=1 
                   ENDIF

                   IF(FLAGZ1 .EQ. 0 ) THEN

                      !II=INICIOZ(NVOID)-1
                      
                      DDENSZ(II)=(UDENSZ(II-1)-UDENSZ(II+1))/(2.D0*DDZ)

                      IF(DDENSZ(II) .GE. GRAD_THRE) FLAGZ1=1

                      IF(FLAG_GAL_P==1 .AND. FLAGZ1==1 .AND. II .GT. 2) THEN !check gradient in the next cell
                         II2=II-1
                         DDENS2=(UDENSZ(II2-1)-UDENSZ(II2+1))/(2.D0*DDZ)
                         IF(DDENS2 .LT. GRAD_THRE) FLAGZ1=0 
                      ENDIF

                      IF(FLAG_DIV_P==1 .AND. DIVZ(II) .LT. 0.) FLAGZ1=1

                      IF(UDENSZ(II) .GE. DENS_THRE2+1.) FLAGZ1=1

                      !IF(FLAGZ2==0) INICIO=II  
                   ENDIF

                   !IF(FLAGZ1==1 .AND. FLAGZ2==1) EXIT OUTZ

                ENDIF !IF (FLAGZ1=0 .OR. FLAGZ2-0)

             ENDDO

!*parallel 
             FLAGZ1_P(IPRO)=FLAGZ1
!*parallel 
             FLAGZ2_P(IPRO)=FLAGZ2

          ENDDO OUTZ

!*parallel 
          FLAGZ1=MAXVAL(FLAGZ1_P(1:NUM_MAX))
!*parallel 
          FLAGZ2=MAXVAL(FLAGZ2_P(1:NUM_MAX))

       ENDIF


!%%%  FLAGXYZ=0: face of the parallelepiped can be  expanded --> move by one row

    IF(FLAGX1 == 0) INICIOX(NVOID)=INICIOX(NVOID)-1 
    IF(FLAGX2 == 0) FINALX(NVOID)=FINALX(NVOID)+1
    IF(FLAGY1 == 0) INICIOY(NVOID)=INICIOY(NVOID)-1
    IF(FLAGY2 == 0) FINALY(NVOID)=FINALY(NVOID)+1
    IF(FLAGZ1 == 0) INICIOZ(NVOID)=INICIOZ(NVOID)-1
    IF(FLAGZ2 == 0) FINALZ(NVOID)=FINALZ(NVOID)+1

!%%% if edges close to the box boundaries stop expansion
! remove void ?
    IF(INICIOX(NVOID) .LE. 1+NSEC) FLAGX1=1   !use NSEC
    IF(FINALX(NVOID) .GE. NX-NSEC) FLAGX2=1 
    IF(INICIOY(NVOID) .LE. 1+NSEC) FLAGY1=1
    IF(FINALY(NVOID) .GE. NY-NSEC) FLAGY2=1
    IF(INICIOZ(NVOID) .LE. 1+NSEC) FLAGZ1=1
    IF(FINALZ(NVOID) .GE. NZ-NSEC) FLAGZ2=1

!%%% flag=1 in all the directions: no more expansion is possible
 IF(FLAGX1==1 .AND. FLAGX2==1 .AND. FLAGY1==1 .AND. FLAGY2==1 .AND. &
      FLAGZ1==1 .AND. FLAGZ2==1) FLAG=1 !no faces to be expanded --> exit loop
   
ENDDO !void expansion: do while (FLAG==0)



   !    IF(INICIOX(NVOID) .GE. FINALX(NVOID) .OR. INICIOY(NVOID) .GE. FINALY(NVOID) .OR. &
   !      INICIOZ(NVOID) .GE. FINALZ(NVOID)) FLAG2=1 !exclude  thin voids
       IF(INICIOX(NVOID) .GE. FINALX(NVOID) .AND. INICIOY(NVOID) .GE. FINALY(NVOID) .AND. &
         INICIOZ(NVOID) .GE. FINALZ(NVOID)) FLAG2=1 !exclude one-cell side voids
       IF(INICIOX(NVOID) .GT. FINALX(NVOID) .OR. INICIOY(NVOID) .GT. FINALY(NVOID) .OR. &
         INICIOZ(NVOID) .GT. FINALZ(NVOID)) FLAG2=1 !include one-cell size voids

!%%% save void quantities


       ICX(NVOID)=IX
       ICY(NVOID)=JY
       ICZ(NVOID)=KZ
       RINIXCO(NVOID)=RX1+(INICIOX(NVOID)-1)*DDX-0.5*DDX
       RFINXCO(NVOID)=RX1+(FINALX(NVOID)-1)*DDX+0.5*DDX
       RINIYCO(NVOID)=RY1+(INICIOY(NVOID)-1)*DDY-0.5*DDX
       RFINYCO(NVOID)=RY1+(FINALY(NVOID)-1)*DDY+0.5*DDX
       RINIZCO(NVOID)=RZ1+(INICIOZ(NVOID)-1)*DDZ-0.5*DDX
       RFINZCO(NVOID)=RZ1+(FINALZ(NVOID)-1)*DDZ+0.5*DDX

       VV=(RFINXCO(NVOID)-RINIXCO(NVOID))*(RFINYCO(NVOID)-RINIYCO(NVOID))*&
            (RFINZCO(NVOID)-RINIZCO(NVOID))

       IF(VV .LT. VMIN0) FLAG2=1 !use minimum void when memory is limited

       !remove voids at the borders and too thin voids (one cell wide in one direction)
       IF( INICIOX(NVOID)==1 .OR. FINALX(NVOID)==NX .OR. &
            INICIOY(NVOID)==1 .OR. FINALY(NVOID)==NY  .OR. &
            INICIOZ(NVOID)==1 .OR. FINALZ(NVOID)==NZ .OR. FLAG2==1) THEN !use 2 and NX-1 that have been used to flag the borders

!%%% remove void from the list   
          INICIOX(NVOID)=0       
          FINALX(NVOID)=0
          INICIOY(NVOID)=0
          FINALY(NVOID)=0
          INICIOZ(NVOID)=0
          FINALZ(NVOID)=0
          RINIXCO(NVOID)=0.
          RFINXCO(NVOID)=0.
          RINIYCO(NVOID)=0.
          RFINYCO(NVOID)=0.
          RINIZCO(NVOID)=0.
          RFINZCO(NVOID)=0.
          ICX(NVOID)=0
          ICY(NVOID)=0
          ICZ(NVOID)=0
          NVOID=NVOID-1  
 
       ENDIF


          IF(NVOID .GE. 1) THEN
             A=(FINALX(NVOID)-INICIOX(NVOID)+1)*DDX
             B=(FINALY(NVOID)-INICIOY(NVOID)+1)*DDY
             C=(FINALZ(NVOID)-INICIOZ(NVOID)+1)*DDZ 
             VOL(NVOID)=A*B*C

       !look for refinements at the edges
          END IF
       ENDIF !if flag1=0

       !write(67,*) NVOID, INICIOX(NVOID), FINALX(NVOID), INICIOY(NVOID), FINALY(NVOID), &
       !     INICIOZ(NVOID), FINALZ(NVOID)

    ENDDO !loop on center cells: L!=1, NXYZ

  END SUBROUTINE VOIDFIND_PAR

!-------------------------------------------------------------------------
SUBROUTINE HALOES(ITER, FILEO, NX, NY, NZ, NVOID, XC, YC, ZC, REQ)
USE COMMONDATA!, ONLY: RADX, RADY, RADZ, DX, DY, DZ, NGAL_MAX, FR, NCELL_GAL !put in parameter file                                                
IMPLICIT NONE
!input variables                                                                                                                                  
INTEGER:: NX, NY, NZ, ITER, NVOID, LOW1, LOW2
!INTEGER, DIMENSION(NX, NY, NZ):: MARCA
!INTEGER, DIMENSION(NVOID):: INDICE                                                                                                               
REAL*4, DIMENSION(NVOID):: XC, YC, ZC, REQ
CHARACTER(LEN=*) :: FILEO
!local variables                                                                                                                                  
INTEGER:: ITT, NHAL, NPARTTOT, II, IH, CONTA, I, J, K, IND, NVOIDH, FLAG, INDW, &
     JJ, KK, IV, INDW1, IN, ICONTA, IC, KMIN, KMAX, JMIN, JMAX, IMIN, IMAX
REAL*4:: XHAL, YHAL, ZHAL, X, DR, D, DMIN
REAL*4, ALLOCATABLE, DIMENSION(:):: MH, XCH, YCH, ZCH, DIST, AGEM, METM, RVIRH 
INTEGER, ALLOCATABLE, DIMENSION(:):: FLAGW, IDHALV, IDH, NPARTH, LEVH, REALCLUS
!REAL*4, DIMENSION(NGAL_MAX) :: MH, XCH, YCH, ZCH, DIST, AGEM, METM 
!INTEGER, DIMENSION(NGAL_MAX):: FLAGW, IDHALV
!INTEGER, DIMENSION(NGAL_MAX):: IDH, NPARTH, LEVH
INTEGER, ALLOCATABLE, DIMENSION(:):: IND2, NHALV
!output variables               
CHARACTER(LEN=50):: FILEH


CALL NOMFILE2('.', ITER, FILEH)
WRITE(*,*) 'File with DM haloes:', FILEH

OPEN(UNIT=21, FILE=FILEH, ACTION='READ')
READ(21,*)
READ(21,*) ITT, X, NHAL
READ(21,*)

ALLOCATE(IDH(NHAL))
ALLOCATE(MH(NHAL))
ALLOCATE(XCH(NHAL))
ALLOCATE(YCH(NHAL))
ALLOCATE(ZCH(NHAL))
ALLOCATE(DIST(NHAL))
ALLOCATE(RVIRH(NHAL))
ALLOCATE(FLAGW(NHAL))
ALLOCATE(IDHALV(NHAL))
ALLOCATE(NPARTH(NHAL))
ALLOCATE(REALCLUS(NHAL))

DO IH=1, NHAL
   READ(21,*) IDH(IH), XCH(IH), YCH(IH), ZCH(IH), MH(IH), RVIRH(IH), NPARTH(IH),X, REALCLUS(IH) !coord are in Mpc
ENDDO

CLOSE(21)


       OPEN(UNIT=30, FILE=FILEO)
       ALLOCATE(IND2(NVOID))
       ALLOCATE(NHALV(NVOID))

       !look for galaxies in voids and void walls                                                                                       
          CONTA=0 !tot numb of gal in voids                                                                                            
          NVOIDH=0 !num of voids with gal                                                                                                         
          NHALV(:)=0
          IDHALV(:)=0
          IND2(:)=0
          FLAGW(:)=0
          DIST(:)=0.
          DO IH=1, NHAL
             IF(REALCLUS(IH) .GE. 0) CYCLE
             IND=0
             INDW=0
             INDW1=0
             FLAG=0 !1=within void                                                                                                                
             DMIN=1.E10

             XHAL=XCH(IH)
             YHAL=YCH(IH)
             ZHAL=ZCH(IH)

             !locate galaxy in the coarse grid                                                                     
             I=INT(((XHAL-RADX(1))/DX+0.49999)+1)
             J=INT(((YHAL-RADY(1))/DY+0.49999)+1)
             K=INT(((ZHAL-RADZ(1))/DZ+0.49999)+1)
             IND=MARCA(I,J,K)
             IF(IND .LT. 0) WRITE(*,*) 'WARNING!! IN GALAXIES: IND < 0!, CASE NOT HANDLED'

             IF(IND .GT. 0) THEN
                IF(IND .GT. NVOID) WRITE(*,*) 'WARNING!! IND > NVOID', IND, NVOID
                FLAG=1
             ELSE
                !LOOK IN WALLS                                                                                                            
                KMIN=MAX(K-NCELL_GAL, 1)
                KMAX=MIN(K+NCELL_GAL, NZ)
                JMIN=MAX(J-NCELL_GAL, 1)
                JMAX=MIN(J+NCELL_GAL, NY)
                IMIN=MAX(I-NCELL_GAL, 1)
                IMAX=MIN(I+NCELL_GAL, NX)

                DO KK=KMIN, KMAX
                   DO JJ=JMIN, JMAX
                      DO II=IMIN, IMAX
                         INDW=MARCA(II,JJ,KK)
                         IF(INDW .GT. 0) THEN
                            IF(INDW .GT. NVOID) WRITE(*,*) 'WARNING!! INDW > NVOID', INDW, NVOID
                            !check if radial interval satisfies: FR*REW                    
                            DR=REQ(INDW)*FR
                            D=SQRT((RADX(I)-RADX(II))**2.+(RADY(J)-RADY(JJ))**2.+&
                                 (RADZ(K)-RADZ(KK))**2.) ! galaxy distance to the closest void edge                               
                            IF(D .LE. DR) THEN ! galaxies lie within 1.25Re                                             
                               FLAG=1 !CHECK ID INDW /= IND                                                                      
                               INDW1=INDW
                               DMIN=MIN(D,DMIN) ! if galaxy falls in several cells of the wall i consider that with the min dist to the gal       
                               !break loop?                                                                                                       
                            ENDIF
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF(FLAG .EQ. 1) THEN
                CONTA=CONTA+1 !numb of gal in walls/voids (counting double occurrances)                                         
                IF(IND .GT. 0) THEN
                   IND2(CONTA)=IND !REMOVE?                                                                                              
                   NHALV(IND)=NHALV(IND)+1
                   IDHALV(CONTA)=IH
                ELSE
                   IND2(CONTA)=INDW1
                   NHALV(INDW1)=NHALV(INDW1)+1
                   IDHALV(CONTA)=IH
                   FLAGW(CONTA)=1
                   DIST(CONTA)=DMIN
                ENDIF
             ENDIF
          ENDDO

          NVOIDH=SUM(NHALV(1:NVOID))

          WRITE(*,*) '************************************************'
          WRITE(*,*) 'Total number of galaxies in voids:', CONTA
          WRITE(*,*) 'Number of voids with galaxies:', NVOIDH


          WRITE(30,*) NVOIDH, CONTA

           DO IND=1, NVOID
             IF(NHALV(IND) .GT. 0 )THEN
                WRITE(30, *) IND, NHALV(IND), REQ(IND)
                DO IC=1, CONTA
                   IF(IND2(IC) .EQ. IND) THEN
                      IH=IDHALV(IC)
                      WRITE(30,*) IDH(IH), NPARTH(IH), MH(IH), XCH(IH), YCH(IH), ZCH(IH), &
                           FLAGW(IC), DIST(IC) !FLAGW=0: within void; =1: within wall           
                   ENDIF
                ENDDO
             ENDIF
          ENDDO

          CLOSE(30)

          WRITE(*,*) 'void galaxies written to:', FILEO
          WRITE(*,*) '************************************************'

          DEALLOCATE(IND2, NHALV)      
          DEALLOCATE(MH, XCH, YCH, ZCH, DIST, RVIRH)
          DEALLOCATE(FLAGW, IDHALV, IDH, NPARTH, REALCLUS)

END SUBROUTINE HALOES
!-------------------------------------------------------------------------

SUBROUTINE GALAXIES(ITER, FILEO, NX, NY, NZ, NVOID, XC, YC, ZC, REQ)
USE COMMONDATA!, ONLY: RADX, RADY, RADZ, DX, DY, DZ, NGAL_MAX, FR, NCELL_GAL !put in parameter file                                                
IMPLICIT NONE
!input variables                                                                                                                                  
INTEGER:: NX, NY, NZ, ITER, NVOID, LOW1, LOW2
!INTEGER, DIMENSION(NX, NY, NZ):: MARCA
!INTEGER, DIMENSION(NVOID):: INDICE                                                                                                               
REAL*4, DIMENSION(NVOID):: XC, YC, ZC, REQ
CHARACTER(LEN=*) :: FILEO
!local variables                                                                                                                                  
INTEGER:: ITT, NHAL, NPARTTOT, II, IH, CONTA, I, J, K, IND, NVOIDH, FLAG, INDW, &
     JJ, KK, IV, INDW1, IN, ICONTA, IC, KMIN, KMAX, JMIN, JMAX, IMIN, IMAX
REAL*4:: XHAL, YHAL, ZHAL, X, DR, D, DMIN
REAL*4, ALLOCATABLE, DIMENSION(:):: MH, XCH, YCH, ZCH, DIST, AGEM, METM, MGAS 
INTEGER, ALLOCATABLE, DIMENSION(:):: FLAGW, IDHALV, IDH, NPARTH, LEVH, MAINPRO
!REAL*4, DIMENSION(NGAL_MAX) :: MH, XCH, YCH, ZCH, DIST, AGEM, METM 
!INTEGER, DIMENSION(NGAL_MAX):: FLAGW, IDHALV
!INTEGER, DIMENSION(NGAL_MAX):: IDH, NPARTH, LEVH
INTEGER, ALLOCATABLE, DIMENSION(:):: IND2, NHALV
!output variables               


write(*,*) 'in galaxies...'
 !read Halma of the corresponding iteration

       ITT=0
       OPEN(UNIT=21, FILE='../halo_results/halma_halo_stars_rp.res500_1100', ACTION='READ')
          READ(21,*)
       DO WHILE(ITT .NE. ITER)
          READ(21,*)
          READ(21,*) NHAL, NPARTTOT, II, ITT
          READ(21,*)
          READ(21,*)
          READ(21,*)
          WRITE(*,*) 'Galaxies at iteration:', ITER, ':', NHAL
          IF(ITT .NE. ITER) THEN
             DO IH=1, NHAL
                READ(21,*)
             ENDDO
          ENDIF
       ENDDO

       IF(ITT .NE. ITER) THEN
          WRITE(*,*) 'WARNING!! HALMA FILE DOES NOT CONTAIN THE CORRECT ITERATION'
       ELSE
          ALLOCATE(MH(NHAL))
          ALLOCATE(XCH(NHAL))
          ALLOCATE(YCH(NHAL))
          ALLOCATE(ZCH(NHAL))
          ALLOCATE(DIST(NHAL))
          ALLOCATE(AGEM(NHAL))
          ALLOCATE(METM(NHAL))
          ALLOCATE(FLAGW(NHAL))
          ALLOCATE(IDHALV(NHAL))
          ALLOCATE(IDH(NHAL))
          ALLOCATE(NPARTH(NHAL))
          ALLOCATE(LEVH(NHAL))
          ALLOCATE(MAINPRO(NHAL))
          ALLOCATE(MGAS(NHAL))
          DO IH=1, NHAL
             !READ(21,*) IDH(IH), NPARTH(IH), MH(IH), X,X,X,X,X,XCH(IH), YCH(IH), ZCH(IH), &
             !     X,X,AGEM(IH),X,METM(IH),X,X,X,X,X,X,X,X, LEVH(IH)
             READ(21,*) IDH(IH), NPARTH(IH), MH(IH), MGAS(IH),X,X,X,X,X,X,X,X,X,X,XCH(IH), YCH(IH), ZCH(IH), &
                  X,X,X,MAINPRO(IH),X,AGEM(IH),X,X,METM(IH),X!, LEVH(IH) ! new halma 
             XCH(IH)=XCH(IH)*1.E-3 !Kpc2Mpc
             YCH(IH)=YCH(IH)*1.E-3
             ZCH(IH)=ZCH(IH)*1.E-3
          ENDDO
       ENDIF
       CLOSE(21)

!
       OPEN(UNIT=30, FILE=FILEO)
       ALLOCATE(IND2(NVOID))
       ALLOCATE(NHALV(NVOID))

       !look for galaxies in voids and void walls                                                                                       
          CONTA=0 !tot numb of gal in voids                                                                                            
          NVOIDH=0 !num of voids with gal                                                                                                         
          NHALV(:)=0
          IDHALV(:)=0
          IND2(:)=0
          FLAGW(:)=0
          DIST(:)=0.
          DO IH=1, NHAL
             IND=0
             INDW=0
             INDW1=0
             FLAG=0 !1=within void                                                                                                                
             DMIN=1.E10

             XHAL=XCH(IH)
             YHAL=YCH(IH)
             ZHAL=ZCH(IH)

             !locate galaxy in the coarse grid                                                                     
             I=INT(((XHAL-RADX(1))/DX+0.49999)+1)
             J=INT(((YHAL-RADY(1))/DY+0.49999)+1)
             K=INT(((ZHAL-RADZ(1))/DZ+0.49999)+1)
             IND=MARCA(I,J,K)
             IF(IND .LT. 0) WRITE(*,*) 'WARNING!! IN GALAXIES: IND < 0!, CASE NOT HANDLED'

             IF(IND .GT. 0) THEN
                IF(IND .GT. NVOID) WRITE(*,*) 'WARNING!! IND > NVOID', IND, NVOID
                FLAG=1
             ELSE
                !LOOK IN WALLS                                                                                                            
                KMIN=MAX(K-NCELL_GAL, 1)
                KMAX=MIN(K+NCELL_GAL, NZ)
                JMIN=MAX(J-NCELL_GAL, 1)
                JMAX=MIN(J+NCELL_GAL, NY)
                IMIN=MAX(I-NCELL_GAL, 1)
                IMAX=MIN(I+NCELL_GAL, NX)

                DO KK=KMIN, KMAX
                   DO JJ=JMIN, JMAX
                      DO II=IMIN, IMAX
                         INDW=MARCA(II,JJ,KK)
                         IF(INDW .GT. 0) THEN
                            IF(INDW .GT. NVOID) WRITE(*,*) 'WARNING!! INDW > NVOID', INDW, NVOID
                            !check if radial interval satisfies: FR*REW                    
                            DR=REQ(INDW)*FR
                            D=SQRT((RADX(I)-RADX(II))**2.+(RADY(J)-RADY(JJ))**2.+&
                                 (RADZ(K)-RADZ(KK))**2.) ! galaxy distance to the closest void edge                               
                            IF(D .LE. DR) THEN ! galaxies lie within 1.25Re                                             
                               FLAG=1 !CHECK ID INDW /= IND                                                                      
                               INDW1=INDW
                               DMIN=MIN(D,DMIN) ! if galaxy falls in several cells of the wall i consider that with the min dist to the gal       
                               !break loop?                                                                                                       
                            ENDIF
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF


             IF(FLAG .EQ. 1) THEN
                CONTA=CONTA+1 !numb of gal in walls/voids (counting double occurrances)                                         
                IF(IND .GT. 0) THEN
                   IND2(CONTA)=IND !REMOVE?                                                                                              
                   NHALV(IND)=NHALV(IND)+1
                   IDHALV(CONTA)=IH
                ELSE
                   IND2(CONTA)=INDW1
                   NHALV(INDW1)=NHALV(INDW1)+1
                   IDHALV(CONTA)=IH
                   FLAGW(CONTA)=1
                   DIST(CONTA)=DMIN
                ENDIF
             ENDIF
          ENDDO

          NVOIDH=SUM(NHALV(1:NVOID))

          WRITE(*,*) '************************************************'
          WRITE(*,*) 'Total number of galaxies in voids:', CONTA
          WRITE(*,*) 'Number of voids with galaxies:', NVOIDH


          WRITE(30,*) NVOIDH, CONTA

          DO IND=1, NVOID
             IF(NHALV(IND) .GT. 0 )THEN
                WRITE(30, *) IND, NHALV(IND), REQ(IND), XC(IND), YC(IND), ZC(IND)
                DO IC=1, CONTA
                   IF(IND2(IC) .EQ. IND) THEN
                      IH=IDHALV(IC)
                      IF(METM(IH) .EQ. 0) METM(IH)=0.00009 !rounded at .res precision                                                                 
                      WRITE(30,*) IDH(IH), NPARTH(IH), MH(IH), XCH(IH), YCH(IH), ZCH(IH), &
                           AGEM(IH), log10(METM(IH)/0.019), LEVH(IH), FLAGW(IC), DIST(IC),MGAS(IC),  MAINPRO(IH) !FLAGW=0: within void; =1: within wall           
                   ENDIF
                ENDDO
             ENDIF
          ENDDO

          CLOSE(30)

          WRITE(*,*) 'void galaxies written to:', FILEO
          WRITE(*,*) '************************************************'

          DEALLOCATE(IND2, NHALV)      
          DEALLOCATE(MH, XCH, YCH, ZCH, DIST, AGEM, METM, MGAS)
          DEALLOCATE(FLAGW, IDHALV, IDH, NPARTH, LEVH, MAINPRO)

END SUBROUTINE GALAXIES

!---------------------------------------------------------------------------------------------------------------------------

SUBROUTINE SHAPE(NX, NY, NZ, NVOID, INDICE, NCELLV, UVOID, XC, YC, ZC, EPS, IP)
USE COMMONDATA
IMPLICIT NONE
!input variables
INTEGER:: NX, NY, NZ, NVOID
!INTEGER, DIMENSION(NX,NY,NZ):: MARCA
INTEGER, DIMENSION(NVOID):: INDICE, UVOID, NCELLV
REAL*4,DIMENSION(NVOID):: XC, YC, ZC
!local variables
INTEGER:: IX,JY,KZ, IND0, I1,I2, NROT, IV, IND
REAL*4:: RR, DELTA, AA, BB, CC, VOLM, VELL, DXX, DYY, DZZ
REAL*4,DIMENSION(3):: RK, AXIS, BASEIGENVAL
REAL, DIMENSION(NVOID, 3,3):: INERTIA
REAL*4, ALLOCATABLE, DIMENSION(:):: RADXX, RADYY, RADZZ
!output variables
REAL*4, DIMENSION(NVOID) :: EPS, IP
INTEGER, PARAMETER:: NCELLV_MIN=0

EPS(:)=0.
IP(:)=0.
INERTIA(:,:,:)=0.


ALLOCATE(RADXX(0:NX+1))
ALLOCATE(RADYY(0:NY+1))
ALLOCATE(RADZZ(0:NZ+1))
   RADXX(:)=RADX(:)
   RADYY(:)=RADY(:)
   RADZZ(:)=RADZ(:)
   DXX=DX
   DYY=DY
   DZZ=DZ

!IF(NX .EQ. NCOX) THEN
!   RADXX(:)=RADX(:)
!   RADYY(:)=RADY(:)
!   RADZZ(:)=RADZ(:)
!   DXX=DX
!   DYY=DY
!   DZZ=DZ
!ELSE IF(NX .EQ. NHYX) THEN
!   RADXX(:)=RADX0(:)
!   RADYY(:)=RADY0(:)
!   RADZZ(:)=RADZ0(:)
!   DXX=DX0
!   DYY=DY0
!   DZZ=DZ0
!ELSE IF(NX .EQ. 2*NHYX) THEN
!   RADXX(:)=RADXR(:)
!   RADYY(:)=RADYR(:)
!   RADZZ(:)=RADZR(:)
!   DXX=0.5*DX0
!   DYY=0.5*DY0
!   DZZ=0.5*DZ0
!ELSE
!   WRITE(*,*) 'IR>1, CASE NOT HANDLED YET'
!ENDIF

!inertia tensor for each void
DO KZ=1, NZ
   DO JY=1, NY
      DO IX=1, NX
         IND0=MARCA(IX,JY,KZ)
         IF(IND0 .GT. 0) THEN !filter void cell
            RK(1)=RADXX(IX)-XC(IND0)
            RK(2)=RADYY(JY)-YC(IND0)
            RK(3)=RADZZ(KZ)-ZC(IND0)

            RR=SQRT(RK(1)**2.+RK(2)**2.+RK(3)**2.)
            DO I1=1, 3
               DO I2=1,3
                  IF(I1 .EQ. I2) DELTA=1 !KHRONAKER DELTA
                  IF(I1 .NE. I2) DELTA=0
                  INERTIA(IND0, I1, I2)=INERTIA(IND0,I1,I2)+DELTA*RR*RR-RK(I1)*RK(I2)
               ENDDO
            ENDDO

         ENDIF
      ENDDO
   ENDDO
ENDDO

DO IV=1,NVOID
   IND=INDICE(IV)
   IF(UVOID(IND) .NE. -1 .OR. NCELLV(IND) .LT. NCELLV_MIN) CYCLE !definire NCELLV_MIN per poter calcl inertia

!normalize inertia
   DO I1=1,3
      DO I2=1,3

         INERTIA(IND,I1,I2)=INERTIA(IND,I1,I2)/REAL(NCELLV(IND))
         !IF(INERTIA(IND, I1,I2) .EQ. 0) WRITE(98,*) IND, NCELLV(IND), INERTIA(IND,I1,I2), I1,I2
         !IF(INERTIA(IND,I1,I2) .LE. 0.) THEN
            !IF(NCELLV(IND) .GT. 8) WRITE(*,*) 'WARNING! INERTIA=0', INERTIA(IND,I1,I2), IND, I1, I2, NCELLV(IND)
            !INERTIA(IND,I1,I2)=1. !??? check
         !ENDIF
      ENDDO
   ENDDO


!diagonalize tensor
   BASEIGENVAL(:)=0.
   CALL JACOBI(INERTIA(IND,:,:),3,BASEIGENVAL,NROT)

!principal axes   
   AXIS(1)=SQRT(2.5*(ABS(BASEIGENVAL(2)+BASEIGENVAL(3)-BASEIGENVAL(1))))
   AXIS(2)=SQRT(2.5*(ABS(BASEIGENVAL(3)+BASEIGENVAL(1)-BASEIGENVAL(2))))
   AXIS(3)=SQRT(2.5*(ABS(BASEIGENVAL(1)+BASEIGENVAL(2)-BASEIGENVAL(3))))

   CALL SORT(AXIS,3,3) !from largest to smallest
   AA=AXIS(1)
   BB=AXIS(2)
   CC=AXIS(3)

   EPS(IND)=1.-CC/AA

   VOLM=NCELLV(IND)*DXX*DYY*DZZ ! ACTUAL VOLUME OF THE VOID 

   VELL=(4./3)*PI*AA*BB*CC ! VOLUME OF THE ELLIPSOIDAL FIT

   IP(IND)=VOLM/VELL
   IF(IP(IND) .GT. 1.) THEN
      IP(IND)=-99
      EPS(IND)=-99.
   ENDIF
   !WRITE(99,*) IND, NCELLV(IND), EPS(IND), AA,BB,CC, IP(IND), VOLM, VELL

ENDDO

DEALLOCATE(RADXX, RADYY, RADZZ)
END SUBROUTINE SHAPE



