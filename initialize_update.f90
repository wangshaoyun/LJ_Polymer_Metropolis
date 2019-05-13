module initialize_update
  !------------------------------------!
  !Use system parameters to initialize
  !the positions and update the positions.
  !------------------------------------!
  implicit none

  contains


subroutine Initialize_position
  !--------------------------------------!
  !Initialize linear chains.
  !
  !Input
  !  pos
  !Output
  !  pos
  !External Variables
  !  
  !Routine Referenced:
  !   rij_and_rr, period_condition_rij
  !Reference:
  !   Spherical coordinate on Wiki to generate uniform distribution
  !   on the sphere surface.
  !--------------------------------------!
  use global_variables
  implicit none
  integer :: i, j, k, m, n, x, y, p
  integer, intent(inout) :: l
  real*8 :: theta, rnd1, rnd2, rnd3, rsqr
  real*8, dimension(3) :: rij


  do i=1, Ngl
    l = (i-1)*Ngl + 1
    x=(i-1)/nint(sqrt(Ngl*1.))+1
    y=mod(i-1,nint(sqrt(Ngl*1.)))+1
    pos(l,1)=Lx/nint(sqrt(Ngl*1.))*(x-0.5)-Lx/2
    pos(l,2)=Ly/nint(sqrt(Ngl*1.))*(y-0.5)-Ly/2
    pos(l,3)=Lz/nint(sqrt(Ngl*1.))-Lz/2
    do k=2, Nml
      l=l+1
      m=1
      p=0
      do while (m==1)
        m=0
        call random_number(rnd1)
        if (p<10) then
          rnd1=rnd1/5
        else
          rnd1=rnd1**2
        end if
        call random_number(rnd2)
        pos(l,1)=pos(l-1,1)+R_bond*sin(pi*rnd1)*cos(2*pi*rnd2)
        pos(l,2)=pos(l-1,2)+R_bond*sin(pi*rnd1)*sin(2*pi*rnd2)
        pos(l,3)=pos(l-1,3)+R_bond*cos(pi*rnd1)
        !periodic condition
        call period_condition_rij(pos(l,1:3))
        !
        !Judge whether the particle is close the former paritcle
        !too much.
        do n=1,l-1
          call rij_and_rr(rij,rsqr,n,l)
          if (rsqr<0.7 .or. pos(l,3)<1) then
            m=1
            p=p+1
            cycle
          end if
        end do
      end do
    end do
  end do

end subroutine Initialize_position


! subroutine Initialize_position
!   !------------------------------------!
!   !Initialize position
!   !   This program is used to initialize the position of
!   !   Polyelectrolytes and ions, and parameters, energy of
!   !   the potential.
!   !Input
!   !   pos, random_or_uniform
!   !Output
!   !   pos
!   !External Variables
!   !   pos, random_or_uniform
!   !Routine Referenced:
!   !1.subroutine random_grafted
!   !   initialize chains by randomly grafting on the plate
!   !2.subroutine uniform_grafted
!   !   initialize chains by uniformly grafting on the plate
!   !3.subroutine initialize_ions
!   !   initialize ions in the system
!   !------------------------------------!
!   use global_variables
!   implicit none
!   integer :: i,j,k,l,m
!   real*8 :: rnd1, rnd2, rnd3, rr
!   real*8, dimension(3) :: rij

!   pos=0

!   do i = 1, Ngl

!     !
!     ! the start monomer of each polymer
!     m = 1
!     k = (i-1)*Nml+1
!     do while (m==1)
!       m = 0
!       call random_number(rnd1)
!       call random_number(rnd2)
!       call random_number(rnd3)
!       pos(k,1) = rand1*Lx-Lx/2
!       pos(k,2) = rand2*Ly-Ly/2
!       pos(k,3) = rand3*Lz-Lz/2
!       !
!       !Jugde whether the particle is close the former paritcle
!       !too much.
!       do l = 1, k-1
!         call rij_and_rr(rij,rsqr,l,k)
!         if (rsqr<0.8) then
!           m=1
!           cycle
!         end if
!       end do
!     end do

!     !
!     !the remaining monomers of that polymer
!     do j = 2, Nml
!       k = k + 1
!       m = 1
!       do while (m == 1)
!         m = 0
!         call random_number(rnd1)
!         call random_number(rnd2)
!         pos(k,1) = pos(k-1,1) + R_bond*cos(2*pi*rnd2)*sin(pi*rnd1)
!         pos(k,2) = pos(k-1,2) + R_bond*sin(2*pi*rnd2)*sin(pi*rnd1)
!         pos(k,3) = pos(k-1,3) + R_bond*cos(pi*rnd1)
!         !periodic condition
!         call period_condition_rij(pos(k,1:3))
!         !
!         !Jugde whether the particle is close the former paritcle
!         !too much.
!         do l = 1, k-1
!           call rij_and_rr(rij,rsqr,l,k)
!           if (rsqr<0.8) then
!             m=1
!             cycle
!           end if
!         end do
!       end do
!     end do 

!   end do

! end subroutine Initialize_position


subroutine Monte_Carlo_Move( EE, DeltaE )
  !------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8, intent(inout) :: EE
  real*8, intent(out)   :: DeltaE
  integer :: j
  real*8 :: EE1, EE2

  do j = 1, NN-Ngl
!     call total_energy(EE1)

    call Choose_Particle
    call New_Position
    call Delta_Energy(DeltaE)
    call Move_or_not(EE, DeltaE)

    !
    !test EE2-EE1 = DeltaE
    ! call total_energy(EE2)
    ! write(*,*) EE2 - EE1, DeltaE, EE2, EE1  
  end do

end subroutine Monte_Carlo_Move


subroutine choose_particle
  !------------------------------------!
  !This subroutine is used to choose a particle ip to move.
  !   
  !Input
  !   
  !Output
  !   ip
  !External Variables
  !   NN, Nm, Npe, ip
  !Routine Referenced:
  !1.
  !------------------------------------!
  use global_variables
  implicit none
  real*8 :: rnd                 

  call random_number(rnd)
  ip = int(rnd*NN) + 1
  !
  !The monomer anchored on the plate can't move, so we need to choose again.
  do while( mod(ip,Nml) == 1 .and. ip <= Npe )
    call random_number(rnd)
    ip = int(rnd*NN) + 1
  end do

end subroutine choose_particle


subroutine New_Position
  !--------------------------------------!
  !This program is used to generate new position.
  !   
  !Input
  !   ip
  !Output
  !   pos1
  !External Variables
  !   pos, pos1, dr
  !Routine Referenced:
  !1. Periodic_condition( rr(2) )
  !--------------------------------------!
  use global_variables
  implicit none
  real*8 :: rnd(3)

  pos_ip0 = pos(ip,:)
  call random_number(rnd)
  pos_ip1(1:3) = pos_ip0(1:3) + (rnd - 0.5D0) * dr
  pos_ip1(4)   = pos_ip0(4)
  call periodic_condition( pos_ip1(1:2) )

end subroutine New_Position


subroutine Move_or_not(EE, DeltaE)
  !--------------------------------------!
  !
  !   
  !Input
  !   
  !Output
  !   
  !External Variables
  !   pos, pos_ip0, pos_ip1, ip, Beta
  !Routine Referenced:
  !1.
  !Reference:
  !
  !--------------------------------------!
  use global_variables
  use compute_energy
  implicit none
  real*8,  intent(in)   :: DeltaE
  real*8,  intent(inout) :: EE
  real*8  :: rnd
  !
  !Judge whether move or not
  if ( DeltaE < 0 ) then
    pos(ip,1:3) = pos(ip,1:3) + pos_ip1(1:3) - pos_ip0(1:3)
    EE = EE + DeltaE
  else 
    call random_number(rnd)
    if ( rnd < Exp(-Beta*DeltaE) ) then
      pos(ip,1:3) = pos(ip,1:3) + pos_ip1(1:3) - pos_ip0(1:3)
      EE = EE + DeltaE
    endif
  endif
end subroutine Move_or_not

end module initialize_update



