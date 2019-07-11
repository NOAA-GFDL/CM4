module vegn_photosynthesis_mod

#include "../shared/debug.inc"

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif
use fms_mod, only: error_mesg, FATAL, file_exist, close_file, check_nml_error, stdlog, &
      mpp_pe, mpp_root_pe, lowercase
use constants_mod,      only : TFREEZE
use sphum_mod,          only : qscomp

use land_constants_mod, only : BAND_VIS, Rugas,seconds_per_year, mol_h2o, mol_air
use land_debug_mod,     only : is_watch_point
use vegn_data_mod,      only : MSPECIES, PT_C4, spdata
use vegn_tile_mod,      only : vegn_tile_type
use vegn_cohort_mod,    only : vegn_cohort_type, get_vegn_wet_frac
use land_data_mod,      only : log_version

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_photosynthesis_init
public :: vegn_photosynthesis
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'vegn_photosynthesis_mod'
#include "../shared/version_variable.inc"

! values for selector of CO2 option used for photosynthesis
integer, public, parameter :: &
    VEGN_PHOT_CO2_PRESCRIBED  = 1, &
    VEGN_PHOT_CO2_INTERACTIVE = 2

! values for internal vegetation photosynthesis option selector
integer, parameter :: &
    VEGN_PHOT_SIMPLE  = 1, & ! zero photosynthesis
    VEGN_PHOT_LEUNING = 2    ! photosynthesis according to simplified Leuning model

! values for internal vegetation respiration option selector
integer, parameter :: &
    VEGN_RESP_LM3      = 1, &
    VEGN_RESP_GAUTHIER = 2

! ==== module variables ======================================================
integer :: vegn_phot_option = -1 ! selector of the photosynthesis option
integer :: vegn_resp_option = -1 ! selector of the photosynthesis option

character(32) :: photosynthesis_to_use = 'simple' ! or 'leuning'
character(32) :: respiration_to_use = 'lm3' ! or 'gauthier'

logical       :: Kok_effect  = .FALSE. ! if TRUE, Kok effect is taken in photosynthesis
real          :: light_kok   = 0.00004 !mol_of_quanta/(m^2s) PAR
real          :: Inib_factor = 0.5

real :: TmaxP=45.0, ToptP=35.0,  tshrP=0.6, tshlP=1.4 ! Parameters of T-response for photosynthesis
real :: TmaxR=65.0, ToptR=47.0,  tshrR=1.4, tshlR=1.0 ! Parameters of T-response for respiration

character(32) :: co2_to_use_for_photosynthesis = 'prescribed' ! or 'interactive'
   ! specifies what co2 concentration to use for photosynthesis calculations:
   ! 'prescribed'  : a prescribed value is used, equal to co2_for_photosynthesis
   !      specified below.
   ! 'interactive' : concentration of co2 in canopy air is used
real, public, protected :: co2_for_photosynthesis = 350.0e-6 ! concentration of co2 for
   ! photosynthesis calculations, mol/mol. Ignored if co2_to_use_for_photosynthesis is
   ! not 'prescribed'

real :: lai_eps = 0.0 ! threshold for switching to linear approximation for Ag_l

namelist /photosynthesis_nml/ &
    photosynthesis_to_use, respiration_to_use, &
    Kok_effect, light_kok, Inib_factor, &
    TmaxP, ToptP,  tshrP, tshlP, &
    TmaxR, ToptR,  tshrR, tshlR, &
    co2_to_use_for_photosynthesis, co2_for_photosynthesis, &
    lai_eps

integer, public, protected :: vegn_phot_co2_option = -1 ! selector of co2 option used for photosynthesis

contains

! ============================================================================
subroutine vegn_photosynthesis_init()
  ! ---- local vars
  integer :: unit         ! unit for namelist i/o
  integer :: io           ! i/o status for the namelist
  integer :: ierr         ! error code, returned by i/o routines

  call log_version(version, module_name, &
  __FILE__)
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=photosynthesis_nml, iostat=io)
    ierr = check_nml_error(io, 'photosynthesis_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=photosynthesis_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'photosynthesis_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif

  unit=stdlog()
  if (mpp_pe() == mpp_root_pe()) then
     write(unit, nml=photosynthesis_nml)
  endif

  ! convert symbolic names of photosynthesis options into numeric IDs to
  ! speed up selection during run-time
  if (trim(lowercase(photosynthesis_to_use))=='simple') then
     vegn_phot_option = VEGN_PHOT_SIMPLE
  else if (trim(lowercase(photosynthesis_to_use))=='leuning') then
     vegn_phot_option = VEGN_PHOT_LEUNING
  else
     call error_mesg('vegn_photosynthesis_init',&
          'vegetation photosynthesis option photosynthesis_to_use="'//&
          trim(photosynthesis_to_use)//'" is invalid, use "simple" or "leuning"',&
          FATAL)
  endif

  ! convert symbolic names of photosynthesis CO2 options into numeric IDs to
  ! speed up selection during run-time
  if (trim(lowercase(co2_to_use_for_photosynthesis))=='prescribed') then
     vegn_phot_co2_option = VEGN_PHOT_CO2_PRESCRIBED
  else if (trim(lowercase(co2_to_use_for_photosynthesis))=='interactive') then
     vegn_phot_co2_option = VEGN_PHOT_CO2_INTERACTIVE
  else
     call error_mesg('vegn_photosynthesis_init',&
          'vegetation photosynthesis option co2_to_use_for_photosynthesis="'//&
          trim(co2_to_use_for_photosynthesis)//'" is invalid, use "prescribed" or "interactive"',&
          FATAL)
  endif

  if (trim(lowercase(respiration_to_use))=='lm3') then
     vegn_resp_option = VEGN_RESP_LM3
  else if (trim(lowercase(respiration_to_use))=='gauthier') then
     vegn_resp_option = VEGN_RESP_GAUTHIER
  else
     call error_mesg('vegn_photosynthesis_init',&
          'vegetation photosynthesis option respiration_to_use="'//&
          trim(co2_to_use_for_photosynthesis)//'" is invalid, use "lm3" or "gauthier"',&
          FATAL)
  endif
end subroutine vegn_photosynthesis_init


! ============================================================================
! compute stomatal conductance, photosynthesis and respiration
subroutine vegn_photosynthesis ( vegn, &
     PAR_dn, PAR_net, cana_q, cana_co2, p_surf, drag_q, &
     soil_beta, soil_water_supply, &
     evap_demand, stomatal_cond, psyn, resp, &
     lai_kok, Anlayer,lai_light)
  type(vegn_tile_type), intent(in) :: vegn
  real, intent(in)  :: PAR_dn   ! downward PAR at the top of the canopy, W/m2
  real, intent(in)  :: PAR_net  ! net PAR absorbed by the canopy, W/m2
  real, intent(in)  :: cana_q   ! specific humidity in canopy air space, kg/kg
  real, intent(in)  :: cana_co2 ! co2 concentration in canopy air space, mol CO2/mol dry air
  real, intent(in)  :: p_surf   ! surface pressure
  real, intent(in)  :: drag_q   ! drag coefficient for specific humidity
  real, intent(in)  :: soil_beta
  real, intent(in)  :: soil_water_supply ! max supply of water to roots per unit
                                ! active root biomass per second, kg/(m2 s)
  real, intent(out) :: evap_demand ! evaporative water demand, kg/(m2 s)
  real, intent(out) :: stomatal_cond ! stomatal conductance, m/s(?)
  real, intent(out) :: psyn     ! net photosynthesis, mol C/(m2 s)
  real, intent(out) :: resp     ! leaf respiration, mol C/(m2 s)
  real, intent(out) :: lai_kok  ! LAI value for light inhibition m2/m2
  real, intent(out) :: Anlayer
  real, intent(out) :: lai_light ! LAI at which Ag=Resp


  ! ---- local constants
  real, parameter :: res_scaler = 20.0    ! scaling factor for water supply

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cohort
  integer :: sp ! shorthand for vegetation species
  real    :: water_supply ! water supply, mol H2O per m2 of leaves per second
  real    :: Ed ! evaporative demand, mol H2O per m2 of leaves per second
  real    :: fw, fs ! wet and snow-covered fraction of leaves


  ! set the default values for outgoing parameters, overriden by the calculations
  ! in gs_leuning
  lai_kok   = 0.0
  Anlayer   = 0.0
  lai_light = 0.0

  ! get the pointer to the first (and, currently, the only) cohort
  cohort => vegn%cohorts(1)

  select case (vegn_phot_option)

  case(VEGN_PHOT_SIMPLE)
     ! beta non-unity only for "beta" models
     stomatal_cond = soil_beta / (cohort%rs_min  + (1-soil_beta)/drag_q)
     cohort%An_op  = 0
     cohort%An_cl  = 0
     psyn = 0
     resp = 0
     evap_demand   = 0

  case(VEGN_PHOT_LEUNING)
     if(cohort%lai > 0) then
        ! assign species type to local var, purely for convenience
        sp = cohort%species
        ! recalculate the water supply to mol H20 per m2 of leaf per second
        water_supply = soil_water_supply/(mol_h2o*cohort%lai)

        call get_vegn_wet_frac (cohort, fw=fw, fs=fs)
        call gs_Leuning(PAR_dn, PAR_net, cohort%Tv, cana_q, cohort%lai, &
             cohort%leaf_age, p_surf, water_supply, sp, cohort%pt, cana_co2, &
             cohort%extinct, fs+fw, stomatal_cond, psyn, resp, Ed, &
             lai_kok, Anlayer,lai_light)
        ! store the calculated photosynthesis and fotorespiration for future use
        ! in carbon_int
        cohort%An_op  = psyn * seconds_per_year
        cohort%An_cl  = resp * seconds_per_year
        ! convert stomatal conductance, photosynthesis and leaf respiration from units
        ! per unit area of leaf to the units per unit area of land
        stomatal_cond = stomatal_cond*cohort%lai
        psyn          = psyn         *cohort%lai
        resp          = resp         *cohort%lai
        ! convert evaporative demand from mol H2O/(m2_of_leaf s) to
        ! kg/(m2_of_land s)
        evap_demand   = Ed*mol_h2o   *cohort%lai
     else
        ! no leaves means no photosynthesis and no stomatal conductance either
        cohort%An_op  = 0
        cohort%An_cl  = 0
        stomatal_cond = 0
        psyn          = 0
        resp          = 0
        evap_demand   = 0
     endif

  case default
     call error_mesg('vegn_stomatal_cond', &
          'invalid vegetation photosynthesis option', FATAL)
  end select

end subroutine vegn_photosynthesis


! ============================================================================
subroutine gs_Leuning(rad_top, rad_net, tl, ea, lai, leaf_age, &
                   p_surf, ws, pft, pt, ca, &
                   kappa, leaf_wet,  &
                   gs, apot, acl, Ed, &
                   lai_kok, Anlayer, lai_light)
  real,    intent(in)    :: rad_top ! PAR dn on top of the canopy, W/m2
  real,    intent(in)    :: rad_net ! Net canopy PAR to the canopy, W/m2
  real,    intent(in)    :: tl   ! leaf temperature, degK
  real,    intent(in)    :: ea   ! specific humidity in the canopy air (?), kg/kg
  real,    intent(in)    :: lai  ! leaf area index
  real,    intent(in)    :: leaf_age ! age of leaf since budburst (deciduous), days
  real,    intent(in)    :: p_surf ! surface pressure, Pa
  real,    intent(in)    :: ws   ! water supply, mol H2O/(m2 of leaf s)
  integer, intent(in)    :: pft  ! species
  integer, intent(in)    :: pt   ! physiology type (C3 or C4)
  real,    intent(in)    :: ca   ! concentration of CO2 in the canopy air space, mol CO2/mol dry air
  real,    intent(in)    :: kappa! canopy extinction coefficient (move inside f(pft))
  real,    intent(in)    :: leaf_wet ! fraction of leaf that's wet or snow-covered
  ! note that the output is per area of leaf; to get the quantities per area of
  ! land, multiply them by LAI
  real,    intent(out)   :: gs   ! stomatal conductance, m/s
  real,    intent(out)   :: apot ! net photosynthesis, mol C/(m2 s)
  real,    intent(out)   :: acl  ! leaf respiration, mol C/(m2 s)
  real,    intent(out)   :: Ed   ! evaporative demand, mol H2O/(m2 s)
  real,    intent(out)   :: lai_kok ! Lai at which kok effect is considered
  !#### Modified by PPG 2016-12-01
  real,    intent(out)   :: Anlayer
  real,    intent(out)   :: lai_light ! Lai at which Ag=Resp

  ! ---- local vars
  ! photosynthesis
  real :: vm;
  real :: kc,ko; ! Michaelis-Menten constants for CO2 and O2, respectively
  real :: ci;
  real :: capgam; ! CO2 compensation point
  real :: f2,f3;
  real :: coef0,coef1;

  real :: Resp;

  ! conductance related
  real :: b;
  real :: ds;  ! humidity deficit, kg/kg
  real :: hl;  ! saturated specific humidity at the leaf temperature, kg/kg
  real :: do1;

  ! miscellaneous
  real :: dum2;
  real, parameter :: light_crit = 0;
  real, parameter :: gs_lim = 0.25;
  real, parameter :: dLAI = 0.01 ! small LAI increment for Anlayer calculations

  !#### MODIFIED BY PPG 2016-12-01
  real :: Ag_layer
  real :: layer_light

  ! new average computations
  real :: lai_eq;
  real, parameter :: rad_phot = 0.0000046 ! PAR conversion factor of J -> mol of quanta
  real :: light_top;
  real :: par_net;
  real :: Ag;
  real :: An;
  real :: Ag_l;
  real :: Ag_rb;
  real :: anbar;
  real :: gsbar;
  real :: layer
  real :: w_scale;
  real, parameter :: p_sea = 1.0e5 ! sea level pressure, Pa
  ! soil water stress
  real :: an_w,gs_w;

  !########MODIFIED BY PPG 2016-12-05
  real :: TempFactP
  real :: TempFactR
  real :: TempFuncP
  real :: TempFuncR

  if (is_watch_point()) then
     write(*,*) '####### gs_leuning input #######'
     __DEBUG2__(rad_top, rad_net)
     __DEBUG1__(tl)
     __DEBUG1__(ea)
     __DEBUG1__(lai)
     __DEBUG1__(leaf_age)
     __DEBUG1__(p_surf)
     __DEBUG1__(ws)
     __DEBUG1__(pft)
     __DEBUG1__(ca)
     __DEBUG1__(kappa)
     __DEBUG1__(leaf_wet)
     __DEBUG1__(pt)
     write(*,*) '####### end of ### gs_leuning input #######'
  endif
  layer = 0.05
  b=0.01;
  do1=0.09 ; ! kg/kg
  if (pft < 2) do1=0.15;

  ! Convert Solar influx from W/(m^2s) to mol_of_quanta/(m^2s) PAR,
  ! empirical relationship from McCree is light=rn*0.0000046
  light_top = rad_top*rad_phot;
  par_net   = rad_net*rad_phot;

  ! calculate humidity deficit, kg/kg
  call qscomp(tl, p_surf, hl)
  ds = max(hl-ea,0.0)

  ! capgam=0.209/(9000.0*exp(-5000.0*(1.0/288.2-1.0/tl))); - Foley formulation, 1986

  ko=0.25   *exp(1400.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
  kc=0.00015*exp(6000.0*(1.0/288.2-1.0/tl))*p_sea/p_surf;
  vm=spdata(pft)%Vmax*exp(3000.0*(1.0/288.2-1.0/tl));

  !decrease Vmax due to aging of temperate deciduous leaves
  !(based on Wilson, Baldocchi and Hanson (2001)."Plant,Cell, and Environment", vol 24, 571-583)
  if (spdata(pft)%leaf_age_tau>0 .and. leaf_age>spdata(pft)%leaf_age_onset) then
     vm=vm*exp(-(leaf_age-spdata(pft)%leaf_age_onset)/spdata(pft)%leaf_age_tau)
  endif

  capgam=0.5*kc/ko*0.21*0.209; ! Farquhar & Caemmerer 1982

  ! Find respiration for the whole canopy layer
  if (light_top>light_kok) then
     lai_kok=min(log(light_top/light_kok)/kappa,lai)
  else
     lai_kok = 0.0
  endif
  if (Kok_effect) then
     ! modify vm for Vmax later and add a temperature function to it.
     Resp=(1-Inib_factor)*spdata(pft)%gamma_resp*vm*lai_kok+spdata(pft)%gamma_resp*vm*(lai-lai_kok)
  else
     Resp=spdata(pft)%gamma_resp*vm*lai
  endif
  select case (vegn_resp_option)
  case(VEGN_RESP_LM3)
     TempFuncR=(1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE)))
     TempFuncP=(1.0+exp(0.4*(5.0-tl+TFREEZE)))*(1.0+exp(0.4*(tl-45.0-TFREEZE)))
  case(VEGN_RESP_GAUTHIER)
     TempFactP=(TmaxP-(tl-TFREEZE))/(TmaxP-ToptP)
     if (TempFactP < 0.) TempFactP=0.

     TempFactR=(TmaxR-(tl-TFREEZE))/(TmaxR-ToptR)
     if (TempFactR < 0.) TempFactR=0.
     TempFuncR=1/((TempFactR**tshrR)*exp((tshrR/tshlR)*(1.-(TempFactR**tshlR))))
     TempFuncP=1/((TempFactP**tshrP)*exp((tshrP/tshlP)*(1.-(TempFactP**tshlP))))
  end select
  Resp=Resp/TempFuncR
  Anlayer=(-spdata(pft)%gamma_resp*vm*(layer)/TempFuncR)/dLAI

  ! ignore the difference in concentrations of CO2 near
  !  the leaf and in the canopy air, rb=0.

  Ag_l=0.;
  Ag_rb=0.;
  Ag=0.;
  anbar=-Resp/lai;
  gsbar=b;
  lai_light = 0.0

  ! find the LAI level at which gross photosynthesis rates are equal
  ! only if PAR is positive
  if ( light_top > light_crit .and. par_net > 0) then
     if (pt==PT_C4) then ! C4 species
        coef0=(1+ds/do1)/spdata(pft)%m_cond;
        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
        if (ci>capgam) then
           f2=vm;
           f3=18000.0*vm*ci;

           dum2=min(f2,f3)

           ! find LAI level at which rubisco limited rate is equal to light limited rate
           lai_eq = -log(dum2/(kappa*spdata(pft)%alpha_phot*light_top))/kappa;
           lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

           ! gross photosynthesis for light-limited part of the canopy
           if (lai>lai_eps) then
              Ag_l = spdata(pft)%alpha_phot * par_net &
                   * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))
           else
              ! approximation for very small LAI: needed because the general formula above
              ! produces division by zero if LAI is very close to zero
              Ag_l = spdata(pft)%alpha_phot * par_net &
                   * (lai-lai_eq)/lai
           endif

           !#### MODIFIED BY PPG 2017-06-10
           ! find the LAI at which Ag_l equal Resp - This part is optional
           lai_light= log((light_top*kappa*spdata(pft)%alpha_phot * (ci-capgam))/ &
                            (spdata(pft)%gamma_resp*vm*(ci+2*capgam)))/kappa
           lai_light = min(max(0.0,lai_light),lai)

           ! gross photosynthesis for rubisco-limited part of the canopy
           Ag_rb  = dum2*lai_eq

           !#### MODIFIED BY PPG 2017-06-13
           !Correct gross photosynthesis for Temperature Response
           Ag=(Ag_l+Ag_rb)/TempFuncP

           An=Ag-Resp;
           anbar=An/lai;

           if(anbar>0.0) then
               gsbar=anbar/(ci-capgam)/coef0;
           endif

           !#### MODIFIED BY PPG 2016-12-01
           !Calculate the amount of light reaching a fine extra layer of leaves
           layer_light=par_net*(exp(-kappa*lai)-exp(-kappa*(lai+layer)))/(1-exp(-(lai+layer)*kappa))
           !Calculate photosynthesis for the specific layer
           Ag_layer= spdata(pft)%alpha_phot * (ci-capgam)/(ci+2.*capgam) * layer_light

           !Calculate Net Photosynthesis in the layer
           if (Kok_effect .and. layer_light > light_kok) then
              Anlayer=(Ag_layer/TempFuncP-(1-Inib_factor)*spdata(pft)%gamma_resp*vm*layer/TempFuncR)/dLAI
           else
              Anlayer=(Ag_layer/TempFuncP-spdata(pft)%gamma_resp*vm*layer/TempFuncR)/dLAI
           endif
        endif ! ci>capgam
     else ! C3 species
        coef0=(1+ds/do1)/spdata(pft)%m_cond;
        coef1=kc*(1.0+0.209/ko);
        ci=(ca+1.6*coef0*capgam)/(1+1.6*coef0);
        f2=vm*(ci-capgam)/(ci+coef1);
        f3=vm/2.;
        dum2=min(f2,f3);
        if (ci>capgam) then
                      ! find LAI level at which rubisco limited rate is equal to light limited rate
           lai_eq=-log(dum2*(ci+2.*capgam)/(ci-capgam)/ &
                       (spdata(pft)%alpha_phot*light_top*kappa))/kappa;
           lai_eq = min(max(0.0,lai_eq),lai) ! limit lai_eq to physically possible range

           !#### MODIFIED BY PPG 2017-06-10
           ! find the LAI at which Ag_l equal Resp - This part is optional
           lai_light= log((light_top*kappa*spdata(pft)%alpha_phot * (ci-capgam))/ &
                            (spdata(pft)%gamma_resp*vm*(ci+2*capgam)))/kappa
           lai_light = min(max(0.0,lai_light),lai)

           ! gross photosynthesis for light-limited part of the canopy
           if (lai>lai_eps) then
              Ag_l = spdata(pft)%alpha_phot * (ci-capgam)/(ci+2.*capgam) * par_net &
                   * (exp(-lai_eq*kappa)-exp(-lai*kappa))/(1-exp(-lai*kappa))
           else
              ! approximation for very small LAI: needed because the general formula above
              ! produces division by zero if LAI is very close to zero
              Ag_l = spdata(pft)%alpha_phot * (ci-capgam)/(ci+2.*capgam) * par_net &
                   * (lai-lai_eq)/lai
           endif

           ! gross photosynthesis for rubisco-limited part of the canopy
           Ag_rb  = dum2*lai_eq

           !#### MODIFIED BY PPG 2017-06-13
           !Correct gross photosynthesis for Temperature Response
           Ag=(Ag_l+Ag_rb)/TempFuncP

           An=Ag-Resp;
           anbar=An/lai;

           if(anbar>0.0) then
               gsbar=anbar/(ci-capgam)/coef0;
           endif

           !#### MODIFIED BY PPG 2016-12-01
           ! Calculate the amount of light reaching a fine extra layer of leaves
           layer_light=par_net*(exp(-kappa*lai)-exp(-kappa*(lai+layer)))/(1-exp(-(lai+layer)*kappa))
           ! Calculate photosynthesis for the specific layer
           Ag_layer= spdata(pft)%alpha_phot * (ci-capgam)/(ci+2.*capgam) * layer_light

           ! Calculate Net Photosynthesis in the layer
           if (Kok_effect .and. layer_light > light_kok) then
              Anlayer=(Ag_layer/TempFuncP-(1-Inib_factor)*spdata(pft)%gamma_resp*vm*layer/TempFuncR)/dLAI
           else
              Anlayer=(Ag_layer/TempFuncP-spdata(pft)%gamma_resp*vm*layer/TempFuncR)/dLAI
           endif
        endif ! ci>capgam
     endif
  endif ! light is available for photosynthesis


  an_w=anbar;
  if (an_w > 0.) then
     an_w=an_w*(1-spdata(pft)%wet_leaf_dreg*leaf_wet);
  endif

  gs_w=gsbar*(1-spdata(pft)%wet_leaf_dreg*leaf_wet);

  if (gs_w > gs_lim) then
      if(an_w > 0.) an_w = an_w*gs_lim/gs_w;
      gs_w = gs_lim;
  endif

#if 1
  ! find water availability
  ! diagnostic demand

  Ed=gs_w*ds*mol_air/mol_h2o;
  ! the factor mol_air/mol_h2o makes units of gs_w and humidity deficit ds compatible:
  ! ds*mol_air/mol_h2o is the humidity deficit in [mol_h2o/mol_air]

  if (Ed>ws) then
     w_scale=ws/Ed;
     gs_w=w_scale*gs_w;
     if(an_w > 0.0) an_w = an_w*w_scale;
     if(an_w < 0.0.and.gs_w >b) gs_w=b;
     if (is_watch_point()) then
        write(*,*)'#### gs is water-limited'
        __DEBUG1__(w_scale)
        __DEBUG3__(gs_w, an_w, b)
     endif
  endif
  gs=gs_w;
  apot=an_w;
  acl=-Resp/lai;


#else
! no water limitation on stomata
   gs=gsbar;
   apot=anbar;
   acl=-Resp/lai;
#endif

   ! finally, convert units of stomatal conductance to m/s from mol/(m2 s) by
   ! multiplying it by a volume of a mole of gas
   gs = gs * Rugas * Tl / p_surf

   if (is_watch_point()) then
      __DEBUG3__(gs, apot, acl)
   endif
end subroutine gs_Leuning

end module vegn_photosynthesis_mod
