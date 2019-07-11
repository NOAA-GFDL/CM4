module land_io_mod

use mpp_domains_mod, only : mpp_pass_sg_to_ug

use mpp_io_mod, only : fieldtype, mpp_get_info, mpp_get_fields, mpp_get_axis_data, &
     mpp_read, validtype, mpp_get_atts, MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE, &
     axistype, mpp_open, mpp_close, mpp_is_valid, mpp_get_file_name, mpp_get_field_index

use axis_utils_mod, only : get_axis_bounds

use constants_mod,     only : PI
use fms_mod, only : file_exist, error_mesg, FATAL, stdlog, mpp_pe, &
     mpp_root_pe, string, check_nml_error, close_file

use mpp_mod, only: mpp_sync
#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif


use horiz_interp_mod,  only : horiz_interp_type, &
     horiz_interp_new, horiz_interp_del, horiz_interp

use land_numerics_mod, only : nearest, bisect
use nf_utils_mod,      only : nfu_validtype, nfu_get_dim, nfu_get_dim_bounds, &
     nfu_get_valid_range, nfu_is_valid, nfu_inq_var, nfu_get_var
use land_data_mod, only : log_version, lnd, horiz_interp_ug

implicit none
private

! ==== public interface ======================================================
public :: init_cover_field
public :: read_field
public :: read_land_io_namelist

public :: print_netcdf_error

public :: input_buf_size
public :: new_land_io
! ==== end of public interface ===============================================

interface read_field
   module procedure read_field_N_2D, read_field_N_3D
   module procedure read_field_I_2D, read_field_I_3D
   module procedure read_field_N_2D_int, read_field_N_3D_int
   module procedure read_field_I_2D_int, read_field_I_3D_int
end interface

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),__FILE__,__LINE__)

! ==== module constants ======================================================
character(len=*), parameter :: module_name = 'land_io_mod'
#include "../shared/version_variable.inc"

real, parameter :: DEFAULT_FILL_INT  = -HUGE(1)
real, parameter :: DEFAULT_FILL_REAL = -HUGE(1.0)

logical :: module_is_initialized = .false.
character(len=64)  :: interp_method = "conservative"
integer :: input_buf_size = 65536 ! input buffer size for tile and cohort reading
logical :: new_land_io = .true.
namelist /land_io_nml/ interp_method, input_buf_size, new_land_io

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine read_land_io_namelist()
  integer :: io, ierr, unit


  module_is_initialized = .TRUE.

  call log_version (version, module_name, &
  __FILE__)

#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=land_io_nml, iostat=io)
     ierr = check_nml_error(io, 'land_io_nml')
#else
  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;
     do while (ierr /= 0)
        read (unit, nml=land_io_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_io_nml')
     enddo
10   continue
     call close_file (unit)
  endif
#endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=land_io_nml)
     call close_file (unit)
  endif

  if(trim(interp_method) .NE. "conservative" .AND. trim(interp_method) .NE. "conserve_great_circle") then
     call error_mesg ( module_name,'interp_method should be "conservative" or "conserve_great_circle"', FATAL)
  endif

  if (input_buf_size <= 0) then
     call error_mesg ( module_name,'input_buf_size must be larger than zero', FATAL)
  endif

end subroutine read_land_io_namelist


! ============================================================================
! This procedure creates and initializes a field of fractional coverage.
subroutine init_cover_field( &
     cover_to_use, filename, cover_field_name, frac_field_name, &
     lonb, latb, uniform_cover, input_cover_types, frac)
  character(len=*), intent(in) :: cover_to_use
  character(len=*), intent(in) :: filename
  character(len=*), intent(in) :: cover_field_name, frac_field_name
  real            , intent(in) :: lonb(:,:), latb(:,:) ! boundaries of the grid cells
  integer         , intent(in) :: uniform_cover
  integer         , intent(in) :: input_cover_types(:)
  real            , intent(out):: frac(:,:) ! output-global map of soil fractional coverage

  ! ---- local vars ---------------------------------------------------------
  integer :: l, k     ! iterators
  integer :: cover_id
  real    :: maxfrac, total

  if( .not. module_is_initialized ) &
       call error_mesg(module_name,'land_io_init is not called', FATAL)

  frac = 0

  if (cover_to_use == 'multi-tile') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
  else if (cover_to_use=='single-tile') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
     do l = 1,size(frac,1)
        total = sum(frac(l,:))
        if (total <= 0) cycle ! operate on valid input data points only
        maxfrac=0 ; cover_id=1
        do k = 1,size(frac,2)
           if(frac(l,k).gt.maxfrac) then
              maxfrac=frac(l,k)
              cover_id=k
           endif
        enddo
        ! set all fractions except dominant fraction to zero
        frac(l,:) = 0.0
        frac(l,cover_id) = total
     enddo
  else if (cover_to_use == 'uniform') then
     call read_cover_field(filename,cover_field_name,frac_field_name,lonb,latb,input_cover_types,frac)
     do l = 1,size(frac,1)
        total = sum(frac(l,:))
        if (total <= 0) cycle ! operate on valid input data points only
        ! set all fractions except dominant fraction to zero
        frac(l,:) = 0.0
        frac(l,uniform_cover) = total
     enddo
  else
     call error_mesg ( module_name,'illegal value of cover_to_use '//cover_to_use, FATAL )
  endif

end subroutine init_cover_field


! ============================================================================
subroutine read_cover_field(file, cover_field_name, frac_field_name,&
     lonb, latb, input_cover_types, frac)
  character(len=*)  , intent(in)  :: file            ! file to read from
  character(len=*)  , intent(in)  :: cover_field_name, frac_field_name
  real              , intent(in)  :: lonb(:,:),latb(:,:) ! boundaries of the model grid
  real              , intent(out) :: frac(:,:)     ! resulting fractions
  integer, optional , intent(in)  :: input_cover_types(:)

  ! --- local vars
! integer :: ncid, varid
  integer :: input_unit , ndim , nvar , natt , nrec , iret
  type(fieldtype), allocatable, dimension(:) :: fields
  type(fieldtype) :: field

  if (.not.file_exist(file)) call error_mesg(module_name,'input file "'//trim(file)//'" does not exist',FATAL)

! If field named 'cover' does not exist in file then read field named 'frac'
! 'cover' does not exist in either ground_type.nc or cover_type.nc
! The extent of the third dimension is 10 in ground_type.nc and 11 in cover_type.nc
  call mpp_open(input_unit, trim(file), action=MPP_RDONLY, form=MPP_NETCDF, &
       threading=MPP_MULTI, fileset=MPP_SINGLE, iostat=iret)
  call mpp_get_info(input_unit,ndim,nvar,natt,nrec)
  allocate(fields(nvar))
  call mpp_get_fields(input_unit,fields)

! if(nf_inq_varid(ncid,cover_field_name,varid)==NF_NOERR) then
!    call do_read_cover_field(ncid,varid,lonb,latb,input_cover_types,frac)
! else if ( nf_inq_varid(ncid,frac_field_name,varid)==NF_NOERR) then
!     call do_read_fraction_field(ncid,varid,lonb,latb,input_cover_types,frac)

  if(get_field(fields,cover_field_name,field)==0) then
     call do_read_cover_field(input_unit,field,lonb,latb,input_cover_types,frac)
  else if ( get_field(fields,frac_field_name,field)==0) then
     call do_read_fraction_field(input_unit,field,lonb,latb,input_cover_types,frac)
  else
     call error_mesg(module_name,&
          'neither "'//trim(cover_field_name)//'" nor "'//&
          frac_field_name//'" is present in input file "'//trim(file)//'"' ,&
          FATAL)
  endif
  call mpp_close(input_unit)

end subroutine read_cover_field

! ============================================================================
function get_field(fields,field_name,field)
  type(fieldtype), intent(in) :: fields(:)
  character(len=*), intent(in) :: field_name
  type(fieldtype), intent(out) :: field
  integer :: get_field, n
  character(len=256) :: name_out

  n =  mpp_get_field_index(fields,trim(field_name))
  if ( n > 0 ) then
    get_field = 0
    field = fields(n)
  else
    get_field = 1
  endif

end function get_field


! ============================================================================
subroutine do_read_cover_field(input_unit, field, lonb, latb, input_cover_types, frac)
  integer, intent(in)  :: input_unit
  type(fieldtype), intent(in) :: field
  real   , intent(in)  :: lonb(:,:),latb(:,:)
  integer, intent(in)  :: input_cover_types(:)
  real   , intent(out) :: frac(:,:)

  ! ---- local vars
  integer :: nlon, nlat, k
  integer, allocatable :: in_cover(:,:)
  real, allocatable    :: in_lonb(:), in_latb(:), x(:,:), r_in_cover(:,:)
  type(horiz_interp_type) :: interp
  type(validtype) :: v
  integer :: in_j_start, in_j_end, in_j_count ! limits of the latitude belt we read
  integer :: ndim, dimlens(2)
  type(axistype) :: axes_centers(2), axis_bounds
  integer :: start(4), count(4)
  character(len=256) :: name
  real :: min_in_latb, max_in_latb, y

  ! check the field dimensions
  call mpp_get_atts(field, ndim=ndim, name=name)
  if (ndim.ne.2) call error_mesg('do_read_cover_field',&
          'cover field "'//trim(name)//'" in file "'//trim(mpp_get_file_name(input_unit))// &
          '" must be two-dimensional (lon,lat)', FATAL)

  ! get size of the longitude and latitude axes
  call mpp_get_atts(field, name=name, siz=dimlens, axes=axes_centers)
  nlon = dimlens(1); nlat = dimlens(2)
  allocate ( in_lonb(nlon+1), in_latb(nlat+1) )

  call get_axis_bounds(axes_centers(1), axis_bounds, axes_centers)
  call mpp_get_axis_data(axis_bounds, in_lonb)
  call get_axis_bounds(axes_centers(2), axis_bounds, axes_centers)
  call mpp_get_axis_data(axis_bounds, in_latb)
  in_lonb = in_lonb*PI/180; in_latb = in_latb*PI/180

  ! to minimize the i/o and work done by horiz_interp, find the boundaries
  ! of latitude belt in input data that covers the entire latb array
  min_in_latb = minval(in_latb); max_in_latb = maxval(in_latb)
  y = minval(latb)
  if (y<min_in_latb) then
     in_j_start = 1
  else if (y>max_in_latb) then
     in_j_start = nlat
  else
     in_j_start=bisect(in_latb, y)
  endif

  y = maxval(latb)
  if (y<min_in_latb) then
     in_j_end = 1
  else if (y>max_in_latb) then
     in_j_end = nlat
  else
     in_j_end = bisect(in_latb, y)
  endif
  in_j_count = in_j_end - in_j_start + 1

  ! check for unreasonable values
  if (in_j_start<1) &
     call error_mesg('do_read_cover_field','reading field "'//trim(name)//'" from file "'&
                     //trim(mpp_get_file_name(input_unit))//'" input latitude start index ('&
                     //trim(string(in_j_start))//') is out of bounds', FATAL)
  if (in_j_count<1) &
     call error_mesg('do_read_cover_field','reading field "'//trim(name)//'" from file "'&
                     //trim(mpp_get_file_name(input_unit))//'" computed input latitude count for domain'&
                     //' is not positive, perhaps input data do not cover entire globe', FATAL)
  if (in_j_start+in_j_count-1>nlat) &
     call error_mesg('do_read_cover_field','reading field "'//trim(name)//'" from file "'&
                     //trim(mpp_get_file_name(input_unit))//'input latitude count ('&
                     //trim(string(in_j_count))//') is too large (start index='&
                     //trim(string(in_j_start))//')', FATAL)

  ! allocate input data buffers
  allocate ( x(nlon,in_j_count), in_cover(nlon,in_j_count), r_in_cover(nlon,in_j_count) )

  ! read input data
! iret = nf_get_vara_int(ncid,varid, (/1,in_j_start/), (/nlon,in_j_count/), in_cover)
  start = (/1,in_j_start,1,1/)
  count(1:2) = shape(in_cover)
  count(3:4) = 1
  call mpp_read( input_unit, field, r_in_cover, start, count )
  in_cover = r_in_cover
  call mpp_get_atts(field,valid=v)

  call horiz_interp_new(interp, in_lonb,in_latb(in_j_start:in_j_start+in_j_count), &
       lonb,latb, interp_method=trim(interp_method))
  frac=0
  do k = 1,size(input_cover_types(:))
     x=0
     where(mpp_is_valid(r_in_cover,v).and.in_cover==input_cover_types(k)) x = 1
     call horiz_interp_ug(interp,x,frac(:,k))
  enddo

  call horiz_interp_del(interp)

  ! clean up memory
  deallocate(in_lonb, in_latb, in_cover, x)

end subroutine do_read_cover_field


! ============================================================================
 subroutine do_read_fraction_field(input_unit,field,lonb,latb,input_cover_types,frac)
  integer, intent(in)  :: input_unit
  type(fieldtype), intent(in) :: field
  real   , intent(in)  :: lonb(:,:),latb(:,:)
  integer, intent(in)  :: input_cover_types(:)
  real   , intent(out) :: frac(:,:)

  ! ---- local vars
  integer :: nlon, nlat, ntypes, k, cover
  real, allocatable :: in_frac(:,:,:)
  real, allocatable :: in_lonb(:), in_latb(:)
  real, allocatable :: in_mask(:,:)
  type(horiz_interp_type) :: interp
  type(validtype) :: v
  integer :: in_j_end, in_j_start, in_j_count ! limits of the latitude belt we read
  integer :: ndim, dimlens(3)
  type(axistype) :: axes_centers(3), axis_bounds
  integer :: start(4), count(4)
  character(len=256) :: name
  real :: min_in_latb, max_in_latb, y

  ! check the field dimensions
  call mpp_get_atts(field, ndim=ndim, name=name)
  if (ndim.ne.3) call error_mesg('do_read_fraction_field', &
         'fraction field "'//trim(name)//'" in file "'//trim(mpp_get_file_name(input_unit))// &
         '" must be three-dimensional (lon,lat,_)', FATAL)

  ! get size of the longitude and latitude axes
  call mpp_get_atts(field, name=name, siz=dimlens, axes=axes_centers)
  nlon = dimlens(1); nlat = dimlens(2) ; ntypes = dimlens(3)
  allocate ( in_lonb(nlon+1), in_latb(nlat+1) )

  call get_axis_bounds(axes_centers(1), axis_bounds, axes_centers)
  call mpp_get_axis_data(axis_bounds, in_lonb)
  call get_axis_bounds(axes_centers(2), axis_bounds, axes_centers)
  call mpp_get_axis_data(axis_bounds, in_latb)
  in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0
  ! find the boundaries of latitude belt in input data that covers the
  ! entire latb array
  min_in_latb = minval(in_latb); max_in_latb = maxval(in_latb)
  y = minval(latb)
  if (y<min_in_latb) then
     in_j_start = 1
  else if (y>max_in_latb) then
     in_j_start = nlat
  else
     in_j_start=bisect(in_latb, y)
  endif

  y = maxval(latb)
  if (y<min_in_latb) then
     in_j_end = 1
  else if (y>max_in_latb) then
     in_j_end = size(in_latb)-1
  else
     in_j_end = bisect(in_latb, y)
  endif
  in_j_count = in_j_end - in_j_start + 1

  ! check for unreasonable values
  if (in_j_start<1) &
     call error_mesg('do_read_fraction_field','reading field "'//trim(name)//'" from file "'&
                     //trim(mpp_get_file_name(input_unit))//'input latitude start index ('&
                     //trim(string(in_j_start))//') is out of bounds', FATAL)
  if (in_j_count<1) &
     call error_mesg('do_read_fraction_field','reading field "'//trim(name)//'" from file "'&
                     //trim(mpp_get_file_name(input_unit))//'computed input latitude count for domain'&
                     //' is not positive, perhaps input data do not cover entire globe', FATAL)
  if (in_j_start+in_j_count-1>nlat) &
     call error_mesg('do_read_fraction_field','reading field "'//trim(name)//'" from file "'&
                     //trim(mpp_get_file_name(input_unit))//'input latitude count ('&
                     //trim(string(in_j_count))//') is too large (start index='&
                     //trim(string(in_j_start))//')', FATAL)

  allocate( in_mask(nlon,in_j_count), in_frac(nlon,in_j_count,ntypes) )

  ! read input data
! iret = nf_get_vara_double(ncid, varid, (/1,in_j_start,1/), (/nlon,in_j_count,ntypes/), in_frac)
  start = (/1,in_j_start,1,1/)
  count(1:3) = shape(in_frac)
  count(4) = 1
  call mpp_read( input_unit, field, in_frac, start, count ) ! interface called here is mpp_read_region_r3D
  call mpp_get_atts(field,valid=v)

  ! Initialize horizontal interpolator; we assume that the valid data mask is
  ! the same for all levels in input frac array. This is probably a good assumption
  ! in all cases.
  where(mpp_is_valid(in_frac(:,:,1),v))
     in_mask = 1.0
  elsewhere
     in_mask = 0.0
  end where
  call horiz_interp_new(interp, &
       in_lonb,in_latb(in_j_start:in_j_start+in_j_count), lonb,latb,&
       interp_method=trim(interp_method), mask_in=in_mask)

  frac = 0
  do k = 1,size(input_cover_types)
     cover = input_cover_types(k)
     if (cover<1.or.cover>ntypes) then
        cycle ! skip all invalid indices in the array of input cover types
     endif

     call horiz_interp_ug(interp,in_frac(:,:,cover),frac(:,k))
  enddo

  ! clean up memory
  call horiz_interp_del(interp)
  deallocate(in_lonb, in_latb, in_frac, in_mask)

end subroutine do_read_fraction_field

! ============================================================================
subroutine read_field_N_2D_int(filename, varname, data_ug, interp, fill)
  character(*), intent(in)  :: filename
  character(*), intent(in)  :: varname
  integer,      intent(out) :: data_ug(:)
  character(*), intent(in), optional :: interp  ! kind of interpolation
  integer,      intent(in), optional :: fill    ! fill value for missing pints
  ! ---- local vars
  real :: data3(size(data_ug,1),1)
  real :: fill_

  fill_ = DEFAULT_FILL_INT
  if (present(fill)) fill_ = fill

  call read_field_N_3D(filename, varname, data3, interp, fill_)
  data_ug = nint(data3(:,1))
end subroutine read_field_N_2D_int

! ============================================================================
subroutine read_field_N_3D_int(filename, varname, data_ug, interp, fill)
  character(*), intent(in)  :: filename
  character(*), intent(in)  :: varname
  integer,      intent(out) :: data_ug(:,:)
  character(*), intent(in), optional :: interp
  integer,      intent(in), optional :: fill
  ! ---- local vars
  real :: data3(size(data_ug,1),size(data_ug,2))
  real :: fill_

  fill_ = DEFAULT_FILL_INT
  if (present(fill)) fill_ = fill

  call read_field_N_3D(filename, varname, data3, interp, fill_)
  data_ug = nint(data3(:,:))
end subroutine read_field_N_3D_int

! ============================================================================
subroutine read_field_I_2D_int(ncid, varname, data_ug, interp, fill)
  integer,      intent(in)  :: ncid
  character(*), intent(in)  :: varname
  integer,      intent(out) :: data_ug(:)
  character(*), intent(in), optional :: interp
  integer,      intent(in), optional :: fill
  ! ---- local vars
  real :: data3(size(data_ug,1),1)
  real :: fill_

  fill_ = DEFAULT_FILL_INT
  if (present(fill)) fill_ = fill

  call read_field_I_3D(ncid, varname, data3, interp, fill_)
  data_ug = nint(data3(:,1))
end subroutine read_field_I_2D_int

! ============================================================================
subroutine read_field_I_3D_int(ncid, varname, data_ug, interp, fill)
  integer,      intent(in)  :: ncid
  character(*), intent(in)  :: varname
  integer,      intent(out) :: data_ug(:,:)
  character(*), intent(in), optional :: interp
  integer,      intent(in), optional :: fill
  ! ---- local vars
  real    :: data3(size(data_ug,1),size(data_ug,2))
  real :: fill_

  fill_ = DEFAULT_FILL_INT
  if (present(fill)) fill_ = fill

  call read_field_I_3D(ncid, varname, data3, interp, fill_)
  data_ug = nint(data3(:,:))
end subroutine read_field_I_3D_int

! ============================================================================
subroutine read_field_N_2D(filename, varname, data_ug, interp, fill)
  character(*), intent(in)  :: filename
  character(*), intent(in)  :: varname
  real,         intent(out) :: data_ug(:)
  character(*), intent(in),  optional :: interp
  real,         intent(in), optional :: fill
  ! ---- local vars
  real    :: data3(size(data_ug,1),1)

  call read_field_N_3D(filename, varname, data3, interp, fill)
  data_ug = data3(:,1)
end subroutine read_field_N_2D

! ============================================================================
subroutine read_field_N_3D(filename, varname, data_ug, interp, fill)
  character(*), intent(in)  :: filename
  character(*), intent(in)  :: varname
  real,         intent(out) :: data_ug(:,:)
  character(*), intent(in), optional :: interp
  real,         intent(in), optional :: fill
  ! ---- local vars
  integer :: ierr, input_unit

  ! Files read: biodata.nc, geohydrology.nc, soil_brdf.nc
  input_unit = -9999
  call mpp_open(input_unit, trim(filename), action=MPP_RDONLY, form=MPP_NETCDF, &
           threading=MPP_MULTI, fileset=MPP_SINGLE, iostat=ierr)
  call read_field_I_3D(input_unit, varname, data_ug, interp, fill)
  call mpp_sync()
  call mpp_close(input_unit)
end subroutine read_field_N_3D

! ============================================================================
subroutine read_field_I_2D(ncid, varname, data_ug, interp, fill)
  integer,      intent(in)  :: ncid
  character(*), intent(in)  :: varname
  real,         intent(out) :: data_ug(:)
  character(*), intent(in), optional :: interp
  real,         intent(in), optional :: fill
  ! ---- local vars
  real    :: data3(size(data_ug,1),1)
  logical :: mask3(size(data_ug,1),1)

  call read_field_I_3D(ncid, varname, data3, interp, fill)
  data_ug = data3(:,1)
end subroutine read_field_I_2D

! ============================================================================
subroutine read_field_I_3D(input_unit, varname, data_ug, interp, fill)
  integer,      intent(in)  :: input_unit
  character(*), intent(in)  :: varname
  real,         intent(out) :: data_ug(lnd%ls:,:)
  character(*), intent(in), optional :: interp
  real,         intent(in), optional :: fill

  ! TODO: possibly check the size of the output array

  ! ---- local vars
  integer :: nlon, nlat, nlev ! size of input grid
  integer :: varndims ! number of variable dimension
  integer :: dimlens(1024) ! sizes of respective dimensions
  real,    allocatable :: in_lonb(:), in_latb(:), in_lon(:), in_lat(:)
  real,    allocatable :: in_data(:,:,:) ! input buffer
  logical, allocatable :: lmask(:,:,:) ! mask of valid input values
  real,    allocatable :: rmask(:,:,:) ! real mask for interpolator
  real,    allocatable :: data_sg(:,:,:) ! data on structured grid
  real,    allocatable :: omask(:,:) ! mask of valid output data
  real,    allocatable :: data2(:,:)
  integer :: k,imap,jmap,l
  type(validtype) :: v
  type(horiz_interp_type) :: hinterp
  integer :: ndim,nvar,natt,nrec
  type(axistype), allocatable ::  varaxes(:)
  type(axistype):: axis_bnd
  type(fieldtype), allocatable :: fields(:)
  type(fieldtype) :: fld
  character(len=256) :: file_name
  integer :: jstart, jend
  character(len=20) :: interp_
  real    :: fill_

  interp_ = 'bilinear'
  if(present(interp)) interp_ = interp

  fill_ = DEFAULT_FILL_REAL
  if (present(fill)) fill_=fill

  ! find the field in the file
  call mpp_get_info(input_unit,ndim,nvar,natt,nrec)
  file_name = mpp_get_file_name(input_unit)
  allocate(fields(nvar))
  call mpp_get_fields(input_unit, fields)
  k = mpp_get_field_index(fields,trim(varname))
  if(k > 0) then
    fld = fields(k)
  else
    call error_mesg('read_field','variable "'//trim(varname)//'" not found in file "'//trim(file_name)//'"',FATAL)
  endif

  ! get the dimensions of our variable
  call mpp_get_atts(fld, ndim=varndims, siz=dimlens, valid=v)
  if(varndims<2.or.varndims>3) then
     call error_mesg('read_field','variable "'//trim(varname)//'" in file "'//trim(file_name)//&
          '" is '//string(varndims)//'D, but only reading 2D or 3D variables is supported', FATAL)
  endif
  allocate(varaxes(varndims))
  call mpp_get_atts(fld, axes=varaxes)
  nlon = dimlens(1) ; nlat = dimlens(2)
  nlev = 1;
  if (varndims==3) nlev=dimlens(3)
  if(nlev/=size(data_ug,2)) then
     call error_mesg('read_field','3rd dimension length of the variable "'&
          //trim(varname)//'" ('//trim(string(nlev))//') in file "'//trim(file_name)//&
          '" is different from the expected size of data ('// trim(string(size(data_ug,2)))//')', &
          FATAL)
  endif

  allocate (in_lon  (nlon),   in_lat  (nlat),  &
            in_lonb (nlon+1), in_latb (nlat+1) )

  ! read boundaries of the grid cells in longitudinal direction
  call mpp_get_axis_data(varaxes(1), in_lon)
  call mpp_get_axis_data(varaxes(2), in_lat)
  in_lon = in_lon*PI/180.0; in_lat = in_lat*PI/180.0
  call get_axis_bounds(varaxes(1), axis_bnd, varaxes)
  call mpp_get_axis_data(axis_bnd, in_lonb)
  call get_axis_bounds(varaxes(2), axis_bnd, varaxes)
  call mpp_get_axis_data(axis_bnd, in_latb)
  in_lonb = in_lonb*PI/180.0; in_latb = in_latb*PI/180.0

  select case(trim(interp_))
  case('nearest')
     allocate (in_data(nlon, nlat, nlev), lmask(nlon, nlat, nlev))
     ! read input data. In case of nearest interpolation we need global fields.
     call mpp_read(input_unit, fld, in_data)
     lmask = mpp_is_valid(in_data,v)
     do k = 1,size(data_ug,2)
        do l = lnd%ls,lnd%le
           call nearest (lmask(:,:,k), in_lon, in_lat, lnd%ug_lon(l), lnd%ug_lat(l), imap, jmap)
           if(imap<0 .or. jmap<0) then
              call error_mesg('read_field', 'imap or jamp is negative' ,FATAL)
           endif
           data_ug(l,k) = in_data(imap,jmap,k)
        enddo
     enddo
     deallocate(in_data,lmask)
  case('bilinear')
     ! we do bilinear interpolation directly on unstructured grid
     call jlimits(minval(lnd%ug_lat),maxval(lnd%ug_lat),in_lat,jstart,jend)
     call read_input()
     call horiz_interp_new(hinterp, in_lonb, in_latb(jstart:jend+1), &
            reshape(lnd%ug_lon,[lnd%le-lnd%ls+1,1]), & ! reshape converts 1D array (N) to array of shape (N,1)
            reshape(lnd%ug_lat,[lnd%le-lnd%ls+1,1]), &
            interp_method='bilinear')
     allocate(omask(size(data_ug,1),1), data2(size(data_ug,1),1))
     do k = 1,size(data_ug,2)
        call horiz_interp(hinterp, in_data(:,:,k), data2(:,:), mask_in=rmask(:,:,k), mask_out=omask(:,:))
        data_ug(:,k) = data2(:,1)
        where (omask(:,1)==0.0) data_ug(:,k) = fill_
     enddo
     call horiz_interp_del(hinterp)
     deallocate(in_data, rmask, omask, data2)
  case('conservative')
     ! conservative interpolation is done on structured grid, and then interpolated values
     ! are passed to unstructured grid
     call jlimits(minval(lnd%sg_latb),maxval(lnd%sg_latb),in_lat,jstart,jend)
     call read_input() ! allocates and fills in_data, rmask
     ! we create horiz interpolator inside the loop, because data masks may be different for
     ! different levels
     allocate (data_sg(lnd%is:lnd%ie,lnd%js:lnd%je,size(data_ug,2))  ) ! data on structured grid
     do k = 1,size(data_ug,2)
        call horiz_interp_new(hinterp, in_lonb, in_latb(jstart:jend+1), lnd%sg_lonb, lnd%sg_latb, &
             mask_in=rmask(:,:,k), interp_method='conservative')
        data_sg(:,:,k) = fill_
        call horiz_interp(hinterp,in_data(:,:,k),data_sg(:,:,k))
        call horiz_interp_del(hinterp)
     enddo
     call mpp_pass_sg_to_ug(lnd%ug_domain, data_sg, data_ug)
     deallocate(in_data, rmask, data_sg)
  case default
     call error_mesg('read_field','Unknown interpolation method "'//trim(interp_)//'". use "nearest", "bilinear", or "conservative"', FATAL)
  end select

  deallocate(in_lonb, in_latb, in_lon, in_lat)
  deallocate(varaxes, fields)

  contains ! internal subroutines

  subroutine read_input
    ! Note that it changes data in host subroutine, so it must be internal
    ! ---- local vars
    integer :: start(4), count(4)

    allocate(in_data(nlon,jstart:jend,nlev), rmask(nlon,jstart:jend,nlev))
    start(1)  = 1;       count(1)  = nlon
    start(2)  = jstart;  count(2)  = jend-jstart+1
    start(3)  = 1;       count(3)  = nlev
    start(4:) = 1;       count(4:) = 1
    ! read input data
    call mpp_read(input_unit, fld, in_data, start, count)
    where (mpp_is_valid(in_data,v))
       rmask = 1.0
    elsewhere
       rmask = 0.0
    end where
  end subroutine read_input

  subroutine jlimits(minlat, maxlat, in_lat, jstart, jend)
    ! to minimize memory footprint, find the latitudinal boundaries in input
    ! data grid that cover our domain.

    ! This subroutine doesn't have to be internal, but it is not used (and perhaps not
    ! useful) anywhere else
    real,    intent(in)  :: minlat, maxlat
    real,    intent(in)  :: in_lat(:)
    integer, intent(out) :: jstart,jend

    integer :: nlat, j
    nlat = size(in_lat)

    jstart = 1; jend = nlat
    do j = 1, nlat
       if(minlat < in_lat(j)) then
          jstart = j-1
          exit
       endif
    enddo
    jstart = max(jstart-1,1)

    do j = 1, nlat
       if(maxlat < in_lat(j)) then
          jend = j
          exit
       endif
    enddo
    jend = min(jend+1,nlat)
  end subroutine jlimits

end subroutine read_field_I_3D

! ============================================================================
subroutine print_netcdf_error(ierr, file, line)
  ! prints out NetCDF library error message, including file name and line number
  integer,          intent(in) :: ierr ! error code
  character(len=*), intent(in) :: file ! name of the file
  integer,          intent(in) :: line ! number of line in the file

  ! ---- local vars
  character(len=1024) :: mesg

  if (ierr.ne.NF_NOERR) then
     write(mesg, "('File ',a,' Line ',i4.4,' :: ',a)") &
          trim(file),line,trim(NF_STRERROR(ierr))
     call error_mesg('NetCDF', mesg, FATAL)
  endif
end subroutine print_netcdf_error

end module
