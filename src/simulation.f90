program Planet_Simulation
    use constants
    use vectors
    implicit none
    integer(ik) :: N_obj, units, center, duration, N_step, print_step, save_step, curr_step=0, save_count=0, i, j, ios
    real(rk) :: h, t=0, G_param, sts, cput1, cput2
    real(rk),allocatable :: mass(:)
    type(vector),allocatable :: x0(:), v0(:), x(:), v(:), a(:), F(:), a_prev(:)

    call cpu_time(cput1)

    !!! READ USER PROVIDED CONFIG FILE !!!
    open(unit=9, file="config.txt", action='read', status='old', iostat=ios)
        if (ios /= 0) then
            print *, "ERROR: Cannot open/find config file."      ! Error message if file cannot be opened/doesn't exist
            stop
        end if

        read(9,*,iostat=ios) N_obj, units, center, duration, N_step, print_step, save_step
        if (ios < 0) then
            print *, "ERROR: Amount of input settings is incorrect."    ! Error message if end of file is reached before all the parameters are read
            stop
        end if
        call verify_config()        !Subroutine that checks if all the parameters are in correct format

        h = real(duration,rk)/real(N_step,rk)   ! time step length


        !========================== CONFIG PARAMETERS ==============================!
        !                                                                           !
        ! N_obj: number of objects used in simulation                               !
        ! units: (1: m|s|kg)  (2: au|day|kg)  (3: pc|year|M_sol)                    !
        ! center: Which object the simulation revolves around (0: Center of mass)   !
        !         (1 - N_obj: Object of given index) (else: No center object)       !
        ! duration: simulation time in selected units                               !
        ! N_step: number of steps used for the simulation, more                     !
        ! print_step: amount of steps between each console print out                !
        ! save_step: amount of steps between each position file save                !
        ! h: length of time step in selected units                                  !
        !                                                                           !
        !===========================================================================!

    close(9)

    !!! READ USER PROVIDED DATA INPUT FILE !!!
    open(unit=10,file="input.dat",action='read',status='old', iostat=ios)
        if (ios /= 0) then
            print *, "ERROR: Cannot open input data file."
            stop
        end if

        allocate(mass(N_obj), x0(N_obj), v0(N_obj), x(N_obj), v(N_obj), a(N_obj), F(N_obj), a_prev(N_obj), stat=ios)  ! allocate arrays to given object amount
        if (ios /= 0) then                          ! Allocation status check.
            print *, "ERROR: Allocation failed."
            stop
        end if

        !========================== ALLOCATED ARRAYS ===============================!
        !                                                                           !
        ! mass: real array containing masses of the objects                         !
        ! x0: vector array containing initial positions of the objects              !
        ! v0: vector array containing initial velocities of the objects             !
        ! x: vector array containing current positions of the objects               !
        ! v: vector array containing current velocities of the objects              !
        ! a: vector array containing current accelerations of the objects           !
        ! F: vector array containing current gravitational forces of the objects    !  
        ! a_prev: acceleration vector array of the previous step,                   !
        !         used in velocity calculation                                      !
        !                                                                           !
        !===========================================================================!

        read(10,*,iostat=ios) mass, x0%x, x0%y, x0%z, v0%x, v0%y, v0%z  ! read masses, positions and velocities from file

        if (ios < 0) then                           ! Check whether there are insufficient data points.
            print *, "ERROR: End of file. Insufficient data points or number of objects is too large."
            stop
        end if

        read(10,*,iostat=ios) sts       ! Check if all given data has been read.
        if (ios == 0) then
            print *, "ERROR: End of file not reached. Too many data points or number of objects is too small."
            stop
        end if
    close(10)

    !!! SIMULATION MOTION PARAMETER INITIALIZATION !!!
    x = x0              ! Set initial positions
    v = v0              ! Set initial velocities
    do i=1, N_obj       ! Calculate starting values of gravitational force and acceleration
        F(i)=vector(0,0,0)
        do j=1, N_obj       ! calculates gravitational pull of each object
            if (j/=i) then
                F(i) = F(i) + F_g(mass(i),mass(j),x(i),x(j))
            end if
        end do
        a(i)=F(i) / mass(i)     ! Newtons second law.
    end do

    call sim_stat()         ! print initials into the console
    !call status_test()      ! debugging function, prints out all the motion parameters


    !!! BEGINS SIMULATION !!!
    open(unit=11,file="output.dat",action='write',status='replace',iostat=ios)
        if (ios /= 0) then
            print *, "ERROR: Cannot create output.dat file."
            stop
        end if

        call centrify()         ! puts the given center object at origin
        write(11,*) t, x        ! write the initial time and position into a file
        save_count = save_count + 1

        do curr_step=1, N_step      ! postition, velocity and acceleration calculated with velocity velvet algorithm
            a_prev = a

            t = t + h                       ! time
            x = x + v*h + 0.5_rk * a * h**2 ! position

            do i=1, N_obj                   ! gravitational force
                F(i)=vector(0,0,0)
                do j=1, N_obj
                    if (j/=i) then
                        F(i) = F(i) + F_g(mass(i),mass(j),x(i),x(j))
                    end if
                end do
            end do

            a = F / mass                    ! acceleration
            v = v + 0.5_rk*(a_prev + a)*h   ! velocity

            if (mod(curr_step,save_step)==0) then   ! save time and position into file
                call centrify()
                write(11,*) t, x   
                save_count = save_count + 1
            end if

            if (mod(curr_step,print_step)==0) call sim_stat()   ! prints status into console
        
            !call status_test()     ! debugging function, prints out all the motion parameters of current step
        end do
    close(11)

    call sim_stat()         ! end of simulation status


    !!! FUNCTIONS AND SUBROUTINES !!!
    contains
    subroutine verify_config() ! verify the given config parameters to be correct
        implicit none
        character(len=1) :: proc

        if (N_obj < 1) then
            print *,"ERROR: Number of objects is less than 1."
            stop
        end if

        if (duration <= 0) then
            print *,"ERROR: Simulation duration is not a positive value."
            stop
        end if

        if (N_step < 1) then
            print *,"ERROR: Number of steps is less than 1."
            stop
        end if

        if (print_step < 1) then
            print *,"ERROR: print_step is less than 1."
            stop
        end if

        if (save_step < 1) then
            print *,"ERROR: save_step is less than 1."
            stop
        end if

        !!! DEFINE UNITS !!!
        if (units == 1) then
            G_param = G_SI
        else if (units == 2) then
            G_param = G_aud
        else if (units == 3) then
            G_param = G_sol
        else
            print *, "ERROR: Wrong unit type given."
            stop
        end if

        if (N_obj * N_step / save_step > 1.3e7) then                                ! file size safe check
            print *,"WARNING: Simulation output file size will be in gigabytes. Proceed? (y/n)"
            read (5,*) proc
            if (proc /= "y") stop
        end if

    end subroutine

    subroutine sim_stat() ! prints current simulation status into console
        print '(A,G3.0,A)', "Simulation progress: ", 100*curr_step/N_step,"%"       !simulation progress in %
        call cpu_time(cput2)
        print '(A,G0)', "Time elapsed (s): ", cput2-cput1                           !time elapsed irl in seconds
        print '(A,G0)', "Number of objects: ", N_obj

        if (units==1) print '(A,G0.5)', "Time in simulation (s): ",t
        if (units==2) print '(A,G0.5)', "Time in simulation (day): ",t
        if (units==3) print '(A,G0.5)', "Time in simulation (year): ",t

        print '(A,G0)', "Current step: ",curr_step
        print '(A,G0)', "Number of steps written to file: ",save_count
        print '(A)', "Current positions: "

        do i=1, N_obj
            print '(A,G0,A,3(G0.5,X))', "Obj.No. ",i, ": ", x(i)
        end do
        !print '(A,G0.5)', "Center of mass position: ", vecsum(mass*x)/sum(mass)
        print *
    end subroutine

    subroutine status_test() ! prints current t, x, v, a, F into console, used for debugging
        print '(A,G0)',"Time: ",t
        do i=1, N_obj
            print '(A,G0)', "Obj.No. ",i
            print '(A,3(G0.5,X))', "x: ", x(i)
            print '(A,3(G0.5,X))', "v: ", v(i)
            print '(A,3(G0.5,X))', "a: ", a(i)
            print '(A,3(G0.5,X))', "F: ", F(i)
        end do
        print *
    end subroutine

    subroutine centrify()   ! centers the coordinate system so that the given object is at origin
        if (center == 0) then
            x = x - vecsum(mass*x)/sum(mass)    ! subtract the position of center of mass
        else if (center > 0 .and. center <= N_obj) then
            x = x - x(center)                   ! subtract the position of center object
        else
            return
        end if
    end subroutine

    elemental function F_g(m1, m2, x1, x2) result(res)    !gravitational force, direction away from m1, towards m2
        implicit none
        real(rk),intent(in) :: m1,m2
        type(vector),intent(in) :: x1,x2
        type(vector) :: res
        
        res = G_param * m1 * m2 / vecdist(x1, x2)**2 * unitvec(x2-x1)
    end function

end program Planet_Simulation