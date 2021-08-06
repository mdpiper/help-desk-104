! An example of the heat equation.
module heatf

  implicit none

         INTEGER, PARAMETER :: NN=5000,NB=30000,NL=100,NURM=1408,NSPECTRUM=5000

  ! Define the attributes of the model.
  type :: heat_model

           double precision :: TIME, TP,WKPO,ANGLE,WT(NN),FREQMIN, &
      FREQMAX, FREQNUM, TIMEBC(NB),TPBC(NB),HRMSBC(NB), &
      WSETBC(NB), SWLBC(NB),WANGBC(NB),FREQMINBC(NB),FREQMAXBC(NB), &
      FREQNUMBC(NB), &
      HRMS(NN),SIGMA(NN),H(NN),WSETUP(NN),SIGSTA(NN),  &
      XBINP(NN,NL),ZBINP(NN,NL),FBINP(NN,NL),XS(NL),  &
      YLINE(NL),DYLINE(NL),AGLINE(NL), &
      DXD2,DXDX,DX2,XB(NN),ZB(NN,NL),FB2(NN,NL),SWLDEP(NN,NL), &
      BSLOPE(NN,NL), &
      GRAV,SQR2,SQR8,PI,TWOPI,SQRG1,SQRG2, &
      WKP,CP(NN),WN(NN),WKPSIN,STHETA(NN),CTHETA(NN), &
      FSX,FSY,FE,QWX,QWY, GBX(NN),GBY(NN),GF(NN), &
      GAMMA,QBREAK(NN),DBSTA(NN),SISMAX,ABREAK(NN), &
      DVEGSTA(NN), &
      SXXSTA(NN),TBXSTA(NN), &
      SXYSTA(NN),TBYSTA(NN), &
      EFSTA(NN),DFSTA(NN), &
      XR,ZR,SSP, &
      UMEAN(NN),USTD(NN),USTA(NN),VMEAN(NN),VSTD(NN),VSTA(NN), &
      WF,SG,SPORO1,WFSGM1,GSGM1,TANPHI,BSLOP1,BSLOP2, &
      EFFB,EFFF,D50,SHIELD,GSD50S,BLP,SLP,BLD,BEDLM,CSTABN,CSEDIA, &
      PS(NN),VS(NN),QSX(NN),QSY(NN), &
      PB(NN),GSLOPE(NN),QBX(NN),QBY(NN),Q(NN), &
      VBX(NN,NL),VSX(NN,NL),VBY(NN,NL),VSY(NN,NL), &
      VY(NN,NL),DZX(NN,NL), &
      DELT,DELZB(NN,NL), &
      RBZERO,RBETA(NN),RQ(NN),RX(NN),RY(NN),RE(NN), &
      XPINP(NN,NL),ZPINP(NN,NL),ZP(NN,NL),HP(NN,NL), &
      WNU,SNP,SDP,BETA1,BETA2,ALSTA,BESTA1,BESTA2,UPMEAN(NN), &
      UPSTD(NN),DPSTA(NN),QP(NN),UPMWD(NN), &
      RWH,RCREST(NL),QO(NL),QOTF,SPRATE,SLPOT, &
      W10(NB),WANGLE(NB),WINDCD(NB),TWXSTA(NB),TWYSTA(NB), &
      AWD,WDN,EWD,CWD,AQWD,BWD,AGWD,AUWD,WPM,ALSTA2,BE2,BE4, &
      PWET(NN),USWD(NN),HWD(NN),SIGWD(NN),UMEAWD(NN),USTDWD(NN), &
      VMEAWD(NN),VSTDWD(NN),HEWD(NN),UEWD(NN),QEWD(NN),H1, &
      SWLAND(NB),HWDMIN,ZW,QD,QM,DETADY(NB),DSWLDT(NB), &
      TSQO(NL),TSQBX(NL),TSQSX(NL), &
      VEGCD,VEGCDM,VEGN(NN,NL),VEGB(NN,NL),VEGD(NN,NL), &
      VEGINP(NN,NL),VEGH(NN,NL),VEGFB(NN,NL),VEGRD(NN,NL), &
      VEGRH(NN,NL),VEGZD(NN,NL),VEGZR(NN,NL),UPROOT(NN,NL), &
      EDIKE(NN,NL),ZB0(NN,NL),DSTA(NN),DSUM(NN), &
      GDINP(NN,NL),GRINP(NN,NL),GRDINP(NN,NL),GRSD(NN,NL), &
      GRSR(NN,NL),GRSRD(NN,NL), DEEB, DEEF, &
      WMINP(NN,NL),WMNODE(NN,NL),ZMESH(NN,NL), &
      ZBSTON(NN,NL),ZPSTON(NN,NL),HPSTON(NN,NL), &
      VDSAND(NN),CPSTON,EPCLAY(NN,NL),ZP0(NN,NL),RCINP(NN,NL), &
      FCINP(NN,NL),RCLAY(NN,NL),FCLAY(NN,NL), &
      DIKETOE, TZ, RUNUPKAPPA, RUNUPPHI, &
      VMEASOMEG(NSPECTRUM),VMEASSE(NSPECTRUM),VMEASWNUM(NSPECTRUM)

          integer :: IPROFL,IANGLE,IROLL,IWIND,IPERM,IOVER,IWCINT, &
      ISEDAV,IWTRAN,IVWALL(NL),ILAB,INFILT,IPOND,ITIDE,ILINE,IQYDY, &
      IVEG,ICLAY,ISMOOTH,IDISS,IFV,IWEIBULL, &
      NWAVE,NSURG,NWIND,NTIME,NBINP(NL),NPINP(NL), &
      NPT,NPE, NMEASSPEC, JMAX(NL),JSWL(NL), JR,JCREST(NL), &
      JWD,JDRY, &
      ISWLSL,JSL,JSL1,IOFLOW,JXW,JX2,NOPOND, ISTSAN, &
      mainloop_itime
     integer :: id

     real :: dt
     real :: t
     real :: t_end

     real :: alpha

     integer :: n_x
     integer :: n_y

     real :: dx
     real :: dy

     real, pointer :: temperature(:,:)
     real, pointer :: temperature_tmp(:,:)

     character(len=90) :: basename
     logical :: enable_cshore_stdout, enable_cshore_outputfiles, &
       enable_varscheck_t0tend, enable_varscheck_alltime


  end type heat_model

  private :: initialize, set_boundary_conditions

contains

!lzhu added
      subroutine compute (cg)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      type (heat_model), intent (inout) :: cg
      character (len=64) :: varscheck_label
      integer :: dump_vars_status
      REAL YVAL
      CHARACTER FINMIN*100, VER*70, BASENAME*90
      DIMENSION DUMVEC(NN),QTIDE(NB),SMDEDY(NB)
      DOUBLE PRECISION  KC, WKZ, WKMEAN, TMEAN
      DOUBLE PRECISION URSELL,HS2H,HV2H,HV2HTOM
!      DATA EPS1, EPS2, MAXITE /1.D-3, 1.D-6, 20/
      END SUBROUTINE compute

  ! Initializes the model with values read from a file.
  subroutine initialize_from_file(model, config_file)
    character (len=*), intent (in) :: config_file
    type (heat_model), intent (out) :: model

    open(15, file=config_file)
    read(15, *) model%alpha, model%t_end, model%n_x, model%n_y
    close(15)
    call initialize(model)
  end subroutine initialize_from_file

  ! Initializes the model with default hardcoded values.
  subroutine initialize_from_defaults(model)
    type (heat_model), intent (out) :: model

    model%alpha = 0.75
    model%t_end = 20.
    model%n_x = 10
    model%n_y = 20
    call initialize(model)
  end subroutine initialize_from_defaults

  ! Allocates memory and sets values for either initialization technique.
  subroutine initialize(model)
    type (heat_model), intent (inout) :: model

    model%id = 0
    model%t = 0.
    model%dt = 1.
    model%dx = 1.
    model%dy = 1.

    allocate(model%temperature(model%n_y, model%n_x))
    allocate(model%temperature_tmp(model%n_y, model%n_x))

    model%temperature = 0.
    model%temperature_tmp = 0.

    call set_boundary_conditions(model%temperature)
    call set_boundary_conditions(model%temperature_tmp)
  end subroutine initialize

  ! Sets boundary conditions on values array.
  subroutine set_boundary_conditions(z)
    implicit none
    real, dimension (:,:), intent (out) :: z
    integer :: i, top_x

    top_x = size(z, 2)-1

    do i = 0, top_x
       z(1,i+1) = 0.25*top_x**2 - (i - 0.5*top_x)**2
    end do
  end subroutine set_boundary_conditions

  ! Frees memory when program completes.
  subroutine cleanup(model)
    type (heat_model), intent (inout) :: model

    deallocate (model%temperature)
    deallocate (model%temperature_tmp)
  end subroutine cleanup

  ! Steps the heat model forward in time.
  subroutine advance_in_time(model)
    type (heat_model), intent (inout) :: model

    call solve_2d(model)
    model%temperature = model%temperature_tmp
    model%t = model%t + model%dt
  end subroutine advance_in_time

  ! The solver for the two-dimensional heat equation.
  subroutine solve_2d(model)
    type (heat_model), intent (inout) :: model

    real :: dx2
    real :: dy2
    real :: coef
    integer :: i, j

    dx2 = model%dx**2
    dy2 = model%dy**2
    coef = model%alpha * model%dt / (2. * (dx2 + dy2))

    do j = 2, model%n_x-1
       do i = 2, model%n_y-1
          model%temperature_tmp(i,j) = &
               model%temperature(i,j) + coef * ( &
               dx2*(model%temperature(i-1,j) + model%temperature(i+1,j)) + &
               dy2*(model%temperature(i,j-1) + model%temperature(i,j+1)) - &
               2.*(dx2 + dy2)*model%temperature(i,j) )
       end do
    end do
  end subroutine solve_2d

  ! A helper routine for displaying model parameters.
  subroutine print_info(model)
    type (heat_model), intent (in) :: model

    write(*,"(a10, i8)") "n_x:", model%n_x
    write(*,"(a10, i8)") "n_y:", model%n_y
    write(*,"(a10, f8.2)") "dx:", model%dx
    write(*,"(a10, f8.2)") "dy:", model%dy
    write(*,"(a10, f8.2)") "alpha:", model%alpha
    write(*,"(a10, f8.2)") "dt:", model%dt
    write(*,"(a10, f8.2)") "t:", model%t
    write(*,"(a10, f8.2)") "t_end:", model%t_end
  end subroutine print_info

  ! A helper routine that prints the current state of the model.
  subroutine print_values(model)
    type (heat_model), intent (in) :: model
    integer :: i, j
    character(len=30) :: rowfmt

    write(rowfmt,'(a,i4,a)') '(', model%n_x, '(1x,f6.1))'
    do i = 1, model%n_y
       write(*,fmt=rowfmt) model%temperature(i,:)
    end do
  end subroutine print_values

end module heatf
