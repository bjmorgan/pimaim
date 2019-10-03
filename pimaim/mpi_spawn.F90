MODULE MPI_SPAWN

! Contains routines for use if pimaim is
! spawned by the MPI Spawn command, used by
! ppfit

IMPLICIT NONE

CONTAINS


    SUBROUTINE read_command_line_args
        !
        ! Subroutine for handling the set-up of a ppfit run of pimaim.
        ! Reads in command line arguments when pimaim is spawned
        ! Sets the values of ppfit_{dipoles,forces,stresses} so we
        ! know what to return to ppfit for fitting against
        !
        ! Args:
        !   None
        ! Returns:
        !   None


        use commondata, only: ppfit_dipoles, ppfit_forces, ppfit_stresses, exit_code
        use mpipara!, only: ierr, iam
        IMPLICIT NONE

        INTEGER:: nargs, i, mpi_comm_parent
        CHARACTER(len=255):: spawndir
        CHARACTER(len=32)::forces, dipoles, stresses, arg, prev_arg
                 
        nargs = command_argument_count()
        call mpi_COMM_get_parent(mpi_comm_parent,ierr)
        DO i=1,nargs
            call getarg(i,arg)
            SELECT CASE(arg)
            CASE('dipoles')
                ppfit_dipoles = .TRUE.
            CASE('forces')
                ppfit_forces = .TRUE.
            CASE('stresses')
                ppfit_stresses = .TRUE.
            CASE('-spawndir')
                CALL getarg(i+1,spawndir)
                CALL chdir(TRIM(spawndir))
!            case('-parentdir')
!                CALL getarg(i+1,parentdir)
!                CALL chdir(TRIM(parentdir))
                
            CASE default
                call getarg(i-1,prev_arg)
                IF (prev_arg .ne. '-spawndir') THEN
                    IF (iam ==0) THEN
                        write(6,*) ' *********************** ERROR **************************'
                        write(6,*) ' Unknown command line argument ', arg, 'options are:'
                        write(6,*) '   dipoles'
                        write(6,*) '   forces'
                        write(6,*) '   stresses'
                        write(6,*) '   -spawndir DIRECTORY_TO_SPAWN_IN'
                        write(6,*) ' *********************** EXITING ************************'
                    END IF
                    call close_down_mpi()
                    stop
                END IF
            END SELECT
        END DO


    END SUBROUTINE read_command_line_args

    SUBROUTINE send_data_to_parent()

        ! Sends data from spawned pimaim process back to ppfit
        ! for fitting.
        
        use mpipara
        use commondata, only: xmu, ymu, zmu, frrx, frry, frrz,stress_tensor,&
                              ppfit_dipoles, ppfit_forces, ppfit_stresses, exit_code
        IMPLICIT NONE

        DOUBLE PRECISION, ALLOCATABLE, dimension(:):: trans_dat
        INTEGER:: no_ele, mpi_comm_parent, length, dum_err
        CHARACTER(MPI_MAX_ERROR_STRING):: message

        no_ele = size(frrx)
        call mpi_COMM_get_parent(mpi_comm_parent,ierr)
        print*, "Sending data in send" 
        if (exit_code == 0) then
            allocate(trans_dat(3*no_ele))
   
            if (ppfit_dipoles) then
                trans_dat(1:no_ele) = xmu(:)
                trans_dat(no_ele+1:no_ele*2) = ymu(:)
                trans_dat((2*no_ele)+1:no_ele*3) = zmu(:)
                call MPI_SEND(trans_dat(:),3*no_ele,mpi_double,0,77,mpi_comm_parent,ierr)
            end if

 
            if (ppfit_forces) then   
            
                trans_dat(1:no_ele) = frrx(:)
                trans_dat(no_ele+1:no_ele*2) = frry(:)
                trans_dat((2*no_ele)+1:no_ele*3) = frrz(:)
                call MPI_SEND(trans_dat(:),3*no_ele,mpi_double,0,66,mpi_comm_parent,ierr)
            end if

            if (ppfit_stresses) then
                call MPI_SEND(dble(stress_tensor(:)),6,mpi_double,0,88,mpi_comm_parent,ierr)
            end if

        end if

        if (ierr .ne. MPI_success) THEN
                write(6,*) 'Error detected in send_data_to__parent(). Exiting'
            call mpi_error_string(ierr, message, length, dum_err)
            write(6,*) message
            call mpi_abort(MPI_COMM_WORLD, 1, ierr)
        ELSE
           deallocate(trans_dat)
           print*,"sent and deallocated"
        END IF
      

    END SUBROUTINE send_data_to_parent

    SUBROUTINE send_exit_code_to_parent
        !send exit code to parent - used by ppfit to check for failed convergence
        ! and to set chi values to inf if it doesnt converge
        use mpipara
        use commondata, only: exit_code

        INTEGER :: mpi_comm_parent

        call mpi_COMM_get_parent(mpi_comm_parent,ierr)
    
        IF (iam .eq. 0) THEN
              call MPI_SEND(exit_code,1,MPI_INTEGER,0,55,mpi_comm_parent,ierr)
!        ELSE
!             CALL MPI_BCast(exit_code, 1, MPI_INTEGER, MPI_PROC_NULL, mpi_comm_parent, ierr)
        END IF

    END SUBROUTINE send_exit_code_to_parent

    SUBROUTINE check_optimisation_status(mpi_comm_parent)

        ! Check whether optimisation has finished and whether to 
        ! run pimaim. Recieves ppfit_stop from parent, which will break loop
        
        use mpipara
        use commondata, only: ppfit_stop
        IMPLICIT NONE

        INTEGER,INTENT(IN):: mpi_comm_parent
        DOUBLE PRECISION, ALLOCATABLE, dimension(:):: trans_dat
        INTEGER::dum_err,length
        CHARACTER(MPI_MAX_ERROR_STRING):: message


        call MPI_BCast(ppfit_stop, 1, MPI_INTEGER, 0, mpi_comm_parent, ierr)
        !call MPI_BARRIER(mpi_comm_parent,ierr)

        if (ierr .ne. MPI_success) THEN
            write(6,*) 'Error detected in check_optimisation_status(). Exiting'
            call mpi_error_string(ierr, message, length, dum_err)
            write(6,*) message
            call mpi_abort(MPI_COMM_WORLD, 1, ierr)
        END IF


    END SUBROUTINE check_optimisation_status


    SUBROUTINE disconnect_from_parent
        
        ! Handles the disconnect from parent when pimaim is spawned
        ! from ppfit: Broadcasts the exit_code back to ppfit to let
        ! it know if the run was succesful, holds the ppfit until
        ! pimaim completes and then disconnects the intracommunicator
        !
        ! Args:
        !   None
        !
        ! Returns:
        !   None

        use mpipara!, only: ierr, iam
        use commondata, only: exit_code
        IMPLICIT NONE


        INTEGER:: mpi_comm_parent, length, dum_err
        CHARACTER(MPI_MAX_ERROR_STRING):: message


        call mpi_COMM_get_parent(mpi_comm_parent,ierr)
            
        IF (mpi_comm_parent .ne. MPI_COMM_NULL) THEN
            !IF (iam .eq. 0) THEN
            !    CALL MPI_BCast(exit_code, 1, MPI_INTEGER, MPI_ROOT, mpi_comm_parent, ierr)
            !ELSE
            !    CALL MPI_BCast(exit_code, 1, MPI_INTEGER, MPI_PROC_NULL, mpi_comm_parent, ierr)
            !END IF
            write(6,*) "About to call barrier"
            CALL MPI_BARRIER(mpi_comm_parent,ierr)
            write(6,*) "about to call disconnect"
            CALL MPI_COMM_DISCONNECT(mpi_comm_parent,ierr)
        end if

        if (ierr .ne. MPI_success) THEN
            write(6,*) 'Error detected in disconnect_from_parent(). Exiting'
            call mpi_error_string(ierr, message, length, dum_err)
            write(6,*) message
            call mpi_abort(MPI_COMM_WORLD, 1, ierr)
        END IF

    END SUBROUTINE disconnect_from_parent

    SUBROUTINE close_down_mpi()
        ! Closes down mpi and stops pimaim
        ! Args:
        !   None
        ! Returns:
        !   None

        use mpipara

        IMPLICIT NONE

        INTEGER:: length, dum_err
        CHARACTER(MPI_MAX_ERROR_STRING):: message


        call disconnect_from_parent()
        call mpi_finalize(ierr)

        if (ierr .ne. MPI_success) THEN
            write(6,*) 'Error detected in close_down_mpi(). Exiting'
            call mpi_error_string(ierr, message, length, dum_err)
            write(6,*) message
            call mpi_abort(MPI_COMM_WORLD, 1, ierr)
        END IF

    END SUBROUTINE close_down_mpi


END MODULE MPI_SPAWN
