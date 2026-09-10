! ==============================================================================
! apply_qc_check_template.inc
! ------------------------------------------------------------------------------
! C-Preprocessor (CPP) template for precision and rank-agnostic QC checks.
! ==============================================================================

subroutine ROUTINE_NAME(var_name, var_array, action_config, rc, mask)
  character(len=*),      intent(in)    :: var_name
  real(VAR_KIND),        intent(inout) :: var_array(VAR_RANK)
  integer,               intent(in)    :: action_config
  integer,               intent(out)   :: rc
  logical,     optional, intent(in)    :: mask(VAR_RANK)
  
  real(ESMF_KIND_R8) :: limit_min_r8, limit_max_r8
  real(VAR_KIND)     :: limit_min, limit_max
  logical            :: is_known
  integer            :: num_bad_cells
  character(len=256) :: log_msg

  rc = ESMF_SUCCESS
  if (action_config == QC_ACTION_IGNORE) return
  
  ! Find limits via O(log N) Binary Search
  call get_limits(var_name, limit_min_r8, limit_max_r8, is_known)
  if (.not. is_known) return

  ! Safely cast 64-bit registry limits to the array's native precision
  limit_min = real(limit_min_r8, VAR_KIND)
  limit_max = real(limit_max_r8, VAR_KIND)

  ! Evaluate violations (Safely bypassing fill-values if masked)
  if (present(mask)) then
    num_bad_cells = count((var_array < limit_min .or. var_array > limit_max) .and. mask)
  else
    num_bad_cells = count(var_array < limit_min .or. var_array > limit_max)
  end if

  ! Take Configured Action
  if (num_bad_cells > 0) then
    select case (action_config)
      case (QC_ACTION_WARN)
        write(log_msg, '(A,A,A,I8,A)') "QC WARNING: ", trim(var_name), " has ", num_bad_cells, " out-of-bounds cells."
        call ESMF_LogWrite(trim(log_msg), ESMF_LOGMSG_WARNING, rc=rc)
        
      case (QC_ACTION_CLAMP)
        write(log_msg, '(A,A,A,I8,A)') "QC CLAMPING: ", trim(var_name), " limits exceeded. Clamping ", num_bad_cells, " cells."
        call ESMF_LogWrite(trim(log_msg), ESMF_LOGMSG_INFO, rc=rc)
        
        if (present(mask)) then
          where (mask .and. var_array < limit_min) var_array = limit_min
          where (mask .and. var_array > limit_max) var_array = limit_max
        else
          var_array = max(limit_min, min(var_array, limit_max))
        end if
        
      case (QC_ACTION_ABORT)
        write(log_msg, '(A,A,A)') "QC FATAL: ", trim(var_name), " violated physical limits! Aborting integration."
        call ESMF_LogWrite(trim(log_msg), ESMF_LOGMSG_ERROR, rc=rc)
        rc = ESMF_RC_VAL_OUTOFRANGE
    end select
  end if
end subroutine ROUTINE_NAME
