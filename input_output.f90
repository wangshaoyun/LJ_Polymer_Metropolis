module input_output
  implicit none
  
  save

  contains

subroutine initialize_parameters
  !------------------------------------!
  !Input and initialize system, timing and histogram parameters.
  !Allocate arrays.
  !------------------------------------!
  use global_variables
  implicit none
  logical alive 

  !
  !Judge whether restart or continue
  Inquire(file='start_time.txt',exist=alive)
  if (alive) then
    open(11,file='./start_time.txt')
      read(11,*) restart_or_continue
    close(11)
  else
    restart_or_continue=0
  end if

  !
  !Input parameters
  call read_data

  Lx = (NN/rho)**(1./3)
  Ly = Lx
  Lz = Lx

  !Write data
  call write_data_to_screen

  !
  !Allocate arrays and initialize them
  call allocatte_arrays_and_initialize

end subroutine initialize_parameters


subroutine read_data
  use global_variables
  implicit none

  open(unit=100, file='system_data.txt')
    read(100,*) rho
    read(100,*) Beta
    read(100,*) NN
    read(100,*) StepNum0
    read(100,*) StepNum
    read(100,*) DeltaStep1
    read(100,*) DeltaStep2
    read(100,*) dr
    read(100,*) best_accpt_ratio
    read(100,*) delta_dr
  close(100)

end subroutine read_data


subroutine write_data_to_screen
  use global_variables
  implicit none

  write(*,*) '******************system_data***********************'
  write(*,*) 'Total particles,          NN :',    NN
  write(*,*) 'Length of the box,        Lx :',    Lx
  write(*,*) 'Width of the box,         Ly :',    Ly
  write(*,*) 'Height of the box,        Lz :',    Lz
  write(*,*) '****************************************************'

  write(*,*) '******************running_steps*********************'
  write(*,*) 'Preheating steps             :', StepNum0
  write(*,*) 'Running steps                :', StepNum
  write(*,*) 'Total steps                  :', (StepNum0+StepNum)
  write(*,*) 'DeltaStep1                   :', DeltaStep1
  write(*,*) 'DeltaStep2                   :', DeltaStep2
  write(*,*) 'Distance of each move        :', dr
  write(*,*) '****************************************************'
end subroutine write_data_to_screen


subroutine allocatte_arrays_and_initialize
  use global_variables
  implicit none

  allocate( pos(NN, 3) )

end subroutine allocatte_arrays_and_initialize


subroutine continue_read_data(l)
  !------------------------------------!
  !
  !------------------------------------!
  use global_variables
  implicit none
  integer, intent(out) :: l
  integer :: i, j 

  open(20,file='./data/pos1.txt')
    read(20,*) ((pos(i,j),j=1,3),i=1,NN)
  close(20)
  open(19,file='./start_time.txt')
    read(19,*)
    read(19,*) l
    read(19,*) dr
    read(19,*) total_time
  close(19)
  
end subroutine continue_read_data


subroutine compute_physical_quantities
  !----------------------------------------!
  ! compute pressure
  !input:
  !  pos
  !output:
  !  pressure
  !External Variables:
  !  Ngl, Nml, Npe, NN,
  !Reference:
  !Frenkel, Smit, 'Understanding molecular simulation: from
  !algorithm to applications', Elsevier, 2002, pp.52, Eq. (3.4.1).
  !----------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  integer i,j,k
  real*8 :: rr, pressure, Rg
  real*8, dimension(3) :: rij
  
  Rg     = 0
  do i = 1, NN-1
    do j = i, NN
        call rij_and_rr(rij, rr, j, k)
        Rg  = Rg + rr
    end do
  end do
  Rg     = Rg / NN / (NN-1)
  !
  !Calculate Pressure
  call compute_pressure(pressure)
  !
  !Output Pressure
  open(37, position='append', file='./data/pressure.txt')
    write(37,370) 1.*j, Rg, pressure
    370 format(3F15.6)
  close(37)
  
end subroutine compute_physical_quantities


subroutine write_pos
  !----------------------------------------!
  !write position to pos.txt
  !input:
  !  pos
  !External Variants:
  !  NN
  !----------------------------------------!
  use global_variables
  implicit none
  integer :: i

  open(30,file='./data/pos.txt')
    do i=1, NN
      write(30,300) pos(i,1), pos(i,2), pos(i,3)
      300 format(3F15.6)
    end do
  close(30)

end subroutine write_pos


subroutine write_pos1(j)
  !----------------------------------------!
  !write position to pos1.txt
  !input:
  !  pos
  !External Variants:
  !  NN
  !----------------------------------------!
  use global_variables
  implicit none
  integer, intent(in) :: j
  integer :: i

  open(30,file='./data/pos1.txt')
    do i=1, NN
      write(30,300) pos(i,1), pos(i,2), pos(i,3)
      300 format(3F15.6)
    end do
  close(30)

  open(32,file='./start_time.txt')
    write(32,*) 1
    write(32,*) j
    write(32,*) dr
    call cpu_time(finished)
    total_time=total_time+finished-started
    call cpu_time(started)
    write(32,*) total_time
  close(32)

end subroutine write_pos1


subroutine write_physical_quantities(j, EE, EE1, DeltaE)
  !----------------------------------------!
  !
  !----------------------------------------!
  use global_variables
  implicit none
  integer, intent(in) :: j
  real*8,  intent(in) :: EE
  real*8,  intent(in) :: EE1
  real*8,  intent(in) :: DeltaE

  open(37,position='append', file='./data/energy_and_time.txt')
    write(37,370) 1.*j, EE, EE1, DeltaE, accpt_ratio
    370 format(5F15.6)
  close(37)

end subroutine write_physical_quantities


subroutine write_time(time)
  !------------------------------------!
  !
  !------------------------------------!
  use global_variables
  implicit none

  real*8, intent(in) :: time
  open(10,file='./data/time.txt')
    write(10,*) 'time:(seconds)', real(total_time)
    write(10,*) 'time:(hours)  ', real(total_time/3600)
    write(10,*) 'time:(days)   ', real(total_time/86400)
    write(10,*) 'Lx:           ', real(Lx)
    write(10,*) 'Ly:           ', real(Ly)
    write(10,*) 'Lz:           ', real(Lz)
    write(10,*) 'NN:           ', NN
  close(10)

end subroutine write_time


end module input_output
