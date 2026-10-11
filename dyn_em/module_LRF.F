! AUTHOR: Yi-Ling Hwong (CCRC), with modifications by Tim Raupach, 2020.
!
! Functions for forcing vertical temperature or moisture tendencies to implement
! Herman & Kuang (2013) Linear Response Function (LRF).

module module_LRF
  
CONTAINS

  subroutine zero_tendencies(        &
       threlaxten, qvrelaxten,       &
       urelaxten, vrelaxten,         &
       ids, ide, jds, jde, kds, kde, &
       ims, ime, jms, jme, kms, kme, &
       ips, ipe, jps, jpe, kps, kpe, &
       its, ite, jts, jte, kts, kte)
    ! Zero tendencies so other tendency functions can be additive.
    !
    ! threlaxten: Potential temperature (theta) relaxation tendency field [K s-1].
    ! qvrelaxten: Water vapour relaxation tendency field [kg kg-1 s-1].
    ! urelaxten: U wind relaxation tendency field [m s-2].
    ! vrelaxten: V wind relaxation tendency field [m s-2].
    ! ids, ide, jds, jde: kds, kde: Horizontal domain dimensions.
    ! ims, ime, jms, jme, kms, kme: Memory dimensions.
    ! ips, ipe, jps, jpe, kps, kpe: Patch dimensions.
    ! its, ite, jts, jte, kts, kte: Tile dimensions.

    implicit none
    
    integer, intent(in) :: ids,ide,jds,jde,kds,kde, &
                           ims,ime,jms,jme,kms,kme, &  
                           ips,ipe,jps,jpe,kps,kpe, &
                           its,ite,jts,jte,kts,kte  

    ! Input/output variables.
    real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: &
         threlaxten, qvrelaxten, urelaxten, vrelaxten

    ! Local variables.
    integer :: i,j,k

    ! Zero all values at mass points.
    do i = its, min(ite, ide-1)
       do k = kts, min(kte, kde-1)
          do j = jts, min(jte, jde-1)    
             threlaxten(i,k,j) = 0
             qvrelaxten(i,k,j) = 0
             urelaxten(i,k,j) = 0
             vrelaxten(i,k,j) = 0
          end do
       end do
    end do
    
  end subroutine zero_tendencies

subroutine relax_var_to_profile(    &
       ten, var, p, tau,              &
       target_profile, is_temp, pi,      &
       ids, ide, jds, jde, kds, kde,   &
       ims, ime, jms, jme, kms, kme,   &
       ips, ipe, jps, jpe, kps, kpe,   &
       its, ite, jts, jte, kts, kte)
    ! Relax given variable toward a profile read from file. 
    !
    ! ---
    ! ten: The variable tendency field. 
    ! var: The variable itself. 
    ! p: Hydrostatic pressure [Pa].
    ! tau: Relaxation time, over which a full relaxation is made [s].
    ! target_profile: The target profile [m s-1].
    ! ids, ide, jds, jde: kds, kde: Horizontal domain dimensions.
    ! ims, ime, jms, jme, kms, kme: Memory dimensions.
    ! ips, ipe, jps, jpe, kps, kpe: Patch dimensions.
    ! its, ite, jts, jte, kts, kte: Tile dimensions.
    ! 
    ! Optional Variables for temperature conversion
    ! is_temp: Check whether temp is being relaxed. 
    ! pi: The Exner function (pi_phy in WRF).

    implicit none

    real, dimension(ims:ime, kms:kme, jms:jme), intent(in)    :: var, p
    real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: ten
    real, intent(in)                 :: tau
    real, dimension(kte), intent(in) :: target_profile
    integer, intent(in) :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         ips, ipe, jps, jpe, kps, kpe, &
         its, ite, jts, jte, kts, kte
    logical, optional, intent(in)    :: is_temp
    real, optional, dimension(ims:ime,kms:kme,jms:jme), intent(in)    :: pi
    integer :: i,j,k
    real :: inv_tau
    real, dimension(kms:kme) :: means
    real :: relax_value

    call horizontal_means(v=var, out=means, &
       ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
       ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
       ips=ips, ipe=ipe, jps=jps, jpe=jpe, kps=kps, kpe=kpe, &
       its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

    ! Relax variable. 
    do j=jts,min(jte, jde-1)
       do k=kts,min(kte, kde-1)
          do i=its,min(ite, ide-1)
             inv_tau = 1/tau

             ! Apply relaxation to shift horizontal means to reference point.
             relax_value = (inv_tau * -(means(k) - target_profile(k)))

             if (present(is_temp)) then
                if (is_temp) then
                   relax_value = relax_value / pi(i,k,j)
                endif
             endif

             ten(i,k,j) = ten(i,k,j) + relax_value
          enddo
       enddo
    enddo

  end subroutine relax_var_to_profile
  
  subroutine horizontal_means(v, out, &
       ids, ide, jds, jde, kds, kde, &
       ims, ime, jms, jme, kms, kme, &
       ips, ipe, jps, jpe, kps, kpe, &
       its, ite, jts, jte, kts, kte)
    ! Calculate horizontal mean values per level.
    !
    ! v: The variable for which to find horizontal means.
    ! out: Output array, one value per level in v.
    ! ids, ide, jds, jde: kds, kde: Horizontal domain dimensions.
    ! ims, ime, jms, jme, kms, kme: Memory dimensions.
    ! ips, ipe, jps, jpe, kps, kpe: Patch dimensions.
    ! its, ite, jts, jte, kts, kte: Tile dimensions.

    USE module_dm, ONLY: wrf_dm_sum_real
    
    implicit none
 
    real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: v
    real, dimension(kms:kme), intent(inout) :: out
    integer, intent(in) :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         ips, ipe, jps, jpe, kps, kpe, &
         its, ite, jts, jte, kts, kte
    
    integer :: i, j, k
    integer :: stag_x, stag_y
    integer :: i_from, i_to, j_from, j_to, k_from, k_to
    real :: sum, count

    do k = kts, min(kte, kde-1)
       sum = 0
       count = 0
       
       do j = jts, min(jte, jde-1)
          do i = its, min(ite, ide-1)
             sum = sum + v(i,k,j)
             count = count + 1
          enddo
       enddo
       
       ! Compute count and sum over whole layer in parallel.
       sum = wrf_dm_sum_real(sum)
       count = wrf_dm_sum_real(count)
       
       ! Check that all domain points have been counted.
       if(.not. count .eq. (ide-ids) * (jde-jds)) then
          call wrf_error_fatal("Wind nudging grid mean domain size mismatch.")
       endif

       out(k) = sum / count
    enddo
  end subroutine horizontal_means
  
  subroutine read_profile_file(num_levels, profile_file, vals)
    ! Read a profile from file. Expect one value per line.
    ! ---
    ! num_levels: The number of levels (values) to read.
    ! profile_file: The file to read from.
    ! vals: The read values (output).

    implicit none
    integer, intent(in) :: num_levels
    character(len=*), intent(in) :: profile_file
    real, dimension(num_levels), intent(inout) :: vals

    integer :: k
    
    open(unit=10, file=profile_file, form='formatted', status='old')
    rewind(10)

    k = 1
    do while(k .le. num_levels)
       read(10,*) vals(k)
       k = k + 1
    enddo

    close(unit=10, status='keep')
  end subroutine read_profile_file
  
  subroutine read_target_profiles(num_levels, profile_file, val)
    ! Read profile from a file.
    ! ---
    ! num_levels: The number of levels (values) to read, one per line.
    ! profile_file: The file to read from.
    ! val: The read values (output).
    
    implicit none
    integer, intent(in) :: num_levels
    character(len=*), intent(in) :: profile_file
    real, dimension(num_levels), intent(inout) :: val
    write(6,*) 'Loading target profiles for nudging.'

    call read_profile_file(num_levels, profile_file, val)
  end subroutine read_target_profiles

end module module_LRF