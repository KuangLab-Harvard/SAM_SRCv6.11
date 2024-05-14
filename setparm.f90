	
subroutine setparm
	
!       initialize parameters:

use vars
!use micro_params
use params
use microphysics, only: micro_setparm
use sgs, only: sgs_setparm
use movies, only : irecc
use instrument_diagnostics, only: zero_instr_diag
implicit none
	
integer icondavg, ierr, ios, ios_missing_namelist, place_holder

NAMELIST /PARAMETERS/ dodamping, doupperbound, docloud, doprecip, &
                dolongwave, doshortwave, dosgs, dz, doconstdz, &
                docoriolis, docoriolisz, dosurface, dolargescale, doradforcing, &
		            fluxt0,fluxq0,tau0,tabs_s,z0,nelapse, nelapsemin, dt, dx, dy,  &
                fcor, ug, vg, nstop, caseid, case_restart,caseid_restart, &
		            nstat, nstatfrq, nprint, nrestart, doradsimple, &
		            nsave3D, nsave3Dstart, nsave3Dend, dosfcforcing, &
		            donudging_uv, donudging_tq, &
                donudging_t, donudging_q, tauls,tautqls,&
                nudging_uv_z1, nudging_uv_z2, nudging_t_z1, nudging_t_z2, &
                nudging_q_z1, nudging_q_z2, dofplane, &
		            timelargescale, longitude0, latitude0, day0, nrad, &
		            OCEAN,LAND,SFC_FLX_FXD,SFC_TAU_FXD, soil_wetness, &
                doensemble, nensemble, dowallx, dowally, &
                nsave2D, nsave2Dstart, nsave2Dend, qnsave3D, & 
                docolumn, save2Dbin, save2Davg, save3Dbin, &
                save2Dsep, save3Dsep, dogzip2D, dogzip3D, restart_sep, &
	              doseasons, doperpetual, doradhomo, dosfchomo, doisccp, &
                domodis, domisr, dodynamicocean, ocean_type, delta_sst, &
                depth_slab_ocean, Szero, deltaS, timesimpleocean, &
                dosolarconstant, solar_constant, zenith_angle, rundatadir, &
                dotracers, output_sep, perturb_type, &
                doSAMconditionals, dosatupdnconditionals, &
                doscamiopdata, iopfile, dozero_out_day0, &
                nstatmom, nstatmomstart, nstatmomend, savemomsep, savemombin, &
                nmovie, nmoviestart, nmovieend, nrestart_skip, &
                bubble_x0,bubble_y0,bubble_z0,bubble_radius_hor, &
                bubble_radius_ver,bubble_dtemp,bubble_dq, dosmoke, dossthomo, &
                rad3Dout, nxco2, dosimfilesout, notracegases, &
                doradlat, doradlon, ncycle_max, doseawater, SLM, LES_S

! Parameters added by Kuang Lab at Harvard
NAMELIST /KUANG_PARAMS/ dompiensemble, &
                dolayerperturb, tperturbi, qperturbi, tperturbA, qperturbA, & ! linear response perturbation: layer by layer (Song Qiyu, 2022)
                nstartperturb, nperturbstep, dorandmultisine, &
                dorandmultisineoddonly, douniformintime, dowhitenoiseforcing, &
                delt_perturbt, delt_perturbq, &
                increase_delt, nstep_increase_start, nstep_increase_end, &
                delt_perturbt_end, delt_perturbq_end, &
                iensemble, nT, nQ, maxperturbperiod, & ! perturbation signal: randmultisine (Qiyu, 2022)
                doidealizedrad, dobulksfc, icopy, &
                nstep_separate_statfile, &
                wavenumber_factor,nstartlinearwave,nsteplinearwavebg,nsteplinearwave,&
                wavedampingtime,wavetqdampingtime,&
                doadvectbg,doparameterizedwave,dointernalnoise
                

!bloss: Create dummy namelist, so that we can figure out error code
!       for a mising namelist.  This lets us differentiate between
!       missing namelists and those with an error within the namelist.
NAMELIST /BNCUIODSBJCB/ place_holder

!----------------------------------
!  Read namelist variables from the standard input:
!------------

open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
read (55,PARAMETERS,IOSTAT=ierr)
if (ierr.ne.0) then
     !namelist error checking
        write(*,*) '****** ERROR: bad specification in PARAMETERS namelist'
        call task_abort()
end if
close(55)

!----------------------------------
!  Read namelist for Kuang_Lab options from same prm file:
!------------
open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

!bloss: get error code for missing namelist (by giving the name for
!       a namelist that doesn't exist in the prm file).
read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
rewind(55) !note that one must rewind before searching for new namelists

!bloss: read in UWOPTIONS namelist
read (UNIT=55,NML=KUANG_PARAMS,IOSTAT=ios)
if (ios.ne.0) then
  if(masterproc) write(*,*) 'ios_missing_namelist = ', ios_missing_namelist
  if(masterproc) write(*,*) 'ios for KUANG_PARAMS = ', ios
   !namelist error checking
   if(ios.ne.ios_missing_namelist) then
     rewind(55) !note that one must rewind before searching for new namelists
     read (UNIT=55,NML=KUANG_PARAMS)
     if(masterproc) then
       write(*,*) '****** ERROR: bad specification in KUANG_PARAMS namelist'
     end if
      call task_abort()
   elseif(masterproc) then
      write(*,*) '****************************************************'
      write(*,*) '******* No KUANG_PARAMS namelist in prm file *******'
      write(*,*) '****************************************************'
   end if
end if
close(55)

! write namelist values out to file for documentation
if(masterproc) then
      open(unit=55,file='./'//trim(case)//'/'//trim(case)//'_'//trim(caseid)//'.nml',&
            form='formatted')
      write (55,nml=PARAMETERS)
      write (55,nml=KUANG_PARAMS)
      write(55,*) 
      close(55)
end if

!------------------------------------
!  Set parameters 


        ! Allow only special cases for separate output:

        output_sep = output_sep.and.RUN3D
        if(output_sep)  save2Dsep = .true.

	if(RUN2D) dy=dx

	if(RUN2D.and.YES3D.eq.1) then
	  print*,'Error: 2D run and YES3D is set to 1. Exitting...'
	  call task_abort()
	endif
	if(RUN3D.and.YES3D.eq.0) then
	  print*,'Error: 3D run and YES3D is set to 0. Exitting...'
	  call task_abort()
	endif

        if(docoriolis.and..not.dofplane.or.doradlat) dowally=.true.

	if(ny.eq.1) dy=dx
	dtn = dt

	notopened2D = .true.
	notopened3D = .true.

        call zero_instr_diag() ! initialize instruments output 
        call sgs_setparm() ! read in SGS options from prm file.
        call micro_setparm() ! read in microphysical options from prm file.

        if(dosmoke) then
           epsv=0.
        else    
           epsv=0.61
        endif   

        if(navgmom_x.lt.0.or.navgmom_y.lt.0) then  
            nstatmom        = 1
            nstatmomstart    = 999999999
            nstatmomend      = 999999999
        end if

        if(doseawater) then
          salt_factor = 0.981
        else
          salt_factor = 1.
        end if

        if(tautqls.eq.99999999.) tautqls = tauls
        
        dtfactor = 1.

        !===============================================================
        ! KUANG_LAB ADDITION

        if(dompiensemble.AND.dompi) then
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Using the Kuang_Lab Ensemble Run Method'
            write(*,*) '  This will turn off MPI in the model run, such that'
            write(*,*) '  each subdomain is run independently of each other.'
            write(*,*) '  However, MPI is turned on for saving of output and'
            write(*,*) '  restart files.'
            write(*,*) '*********************************************************'
          end if
        else if(dompiensemble.AND.(.NOT.dompi)) then
          dompiensemble = .false.
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Do not use the Kuang_Lab Ensemble Run Method'
            write(*,*) '  MPI is not called because number of processors = 1'
            write(*,*) '  Setting dompiensemble to FALSE'
            write(*,*) '*********************************************************'
          end if
        end if
          
        if(dolayerperturb) then
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Using the Kuang_Lab layer-by-layer perturbation'
            write(*,*) '  This is meant to calculate linear response function.'
            if(tperturbi.gt.0.and.tperturbA.ne.0.) then
              write(*,*) '  Add temperature perturbation to layer ', tperturbi
            end if
            if(qperturbi.gt.0.and.qperturbA.ne.0.) then
              write(*,*) '  Add water vapor perturbation to layer ', qperturbi
            end if
            write(*,*) '*********************************************************'
          end if
        end if
        
        if(dorandmultisine) then
          if((nT.le.0).or.(nQ.le.0).or.(nT.gt.nzm).or.(nQ.gt.nzm)) then
            print*,'Check nT and nQ in randmultisine!'
            call task_abort()
          end if
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Using Random Phase Multi-Sine perturbations as input'
            write(*,*) '  signal of the system. These perturbations are added to'
            write(*,*) '  field every nstat steps, while statistic field is'
            write(*,*) '  written in output at the same steps.'
            if(dorandmultisineoddonly) then
              write(*,*) '  - Only odd frequency perturbations are included.'
            else
              write(*,*) '  - All frequency perturbations are included.'
            end if
            if(increase_delt) then
              write(*,*) '  - delt_t/q are: ', delt_perturbt, delt_perturbq
              write(*,*) '  - before step number: ', nstep_increase_start
              write(*,*) '  - then linearly increase to:', delt_perturbt_end, delt_perturbq_end
              write(*,*) '  - by step number: ', nstep_increase_end
              write(*,*) '  - then stay constant.'
            else
              write(*,*) '  - delt_t/q are: ', delt_perturbt, delt_perturbq
            end if
            write(*,*) '*********************************************************'
          end if
        end if

        if(doidealizedrad) then
          if(masterproc) then
            write(*,*) '  Doing idealized radiation.'
            if(dolongwave.or.doshortwave.or.doradforcing) then
              dolongwave = .false.
              doshortwave = .false.
              doradforcing = .false.
              write(*,*) '  Setting dolongwave, doshortwave, doradforcing as false.'
            end if
          end if
        end if

        !===============================================================
        ! UW ADDITION

        !bloss: set up conditional averages
        ncondavg = 1 ! always output CLD conditional average
        if(doSAMconditionals) ncondavg = ncondavg + 2
        if(dosatupdnconditionals) ncondavg = ncondavg + 3
        if(allocated(condavg_factor)) then ! avoid double allocation when nrestart=2
          DEALLOCATE(condavg_factor,condavg_mask,condavgname,condavglongname)
        end if
        ALLOCATE(condavg_factor(nzm,ncondavg), & ! replaces old cloud_factor, core_factor
             condavg_mask(nx,ny,nzm,ncondavg), & ! nx x ny x nzm indicator arrays
             condavgname(ncondavg), & ! short names (e.g. CLD, COR, SATUP)
             condavglongname(ncondavg), & ! long names (e.g. cloud, core, saturated updraft)
             STAT=ierr)
        if(ierr.ne.0) then
             write(*,*) '**************************************************************************'
             write(*,*) 'ERROR: Could not allocate arrays for conditional statistics in setparm.f90'
             call task_abort()
        end if
        
        ! indicators that can be used to tell whether a particular average
        !   is present.  If >0, these give the index into the condavg arrays
        !   where this particular conditional average appears.
        icondavg_cld = -1
        icondavg_cor = -1
        icondavg_cordn = -1
        icondavg_satup = -1
        icondavg_satdn= -1
        icondavg_env = -1

        icondavg = 0
        icondavg = icondavg + 1
        condavgname(icondavg) = 'CLD'
        condavglongname(icondavg) = 'cloud'
        icondavg_cld = icondavg

        if(doSAMconditionals) then
           icondavg = icondavg + 1
           condavgname(icondavg) = 'COR'
           condavglongname(icondavg) = 'core'
           icondavg_cor = icondavg

           icondavg = icondavg + 1
           condavgname(icondavg) = 'CDN'
           condavglongname(icondavg) = 'downdraft core'
           icondavg_cordn = icondavg
        end if
           
        if(dosatupdnconditionals) then
           icondavg = icondavg + 1
           condavgname(icondavg) = 'SUP'
           condavglongname(icondavg) = 'saturated updrafts'
           icondavg_satup = icondavg

           icondavg = icondavg + 1
           condavgname(icondavg) = 'SDN'
           condavglongname(icondavg) = 'saturated downdrafts'
           icondavg_satdn = icondavg

           icondavg = icondavg + 1
           condavgname(icondavg) = 'ENV'
           condavglongname(icondavg) = 'unsaturated environment'
           icondavg_env = icondavg
        end if
           
        ! END UW ADDITIONS
        !===============================================================

        irecc = 1


end
