program dhdl_final
implicit none
real*8 :: grid_min,grid_max,grid_width,x,dum,norm,value,beta_cv
real*8 :: value_dhdl_z,dhdl_current,dhdl_lambda
real*8 :: dhdl_z,num_phi_z,dhdl_prob,dhdl,int_num_phi_z,sum_int
integer :: nbin,bin,i,steps,xi,yj,zk,x_len,y_len,z_len
REAL*8, PARAMETER :: kj_to_kcal=1.d0/4.214d0
!open(2,file='dhdl.xvg',status='old')
!open(3,file='interpol.dat',status='old')
open(2,file='Pu.dat',status='old')
open(1,FILE='2_input.dat',STATUS='OLD')


read(1,*) x_len,y_len,z_len !No. of data points 


call step_count(2,steps)
!print *, 'steps=', steps
rewind(2)

norm=0.d0
dhdl_z=0.d0
num_phi_z=0.d0
dhdl_prob=0.d0
dhdl=0.d0
int_num_phi_z=0.d0
sum_int=0.d0
do xi=1,x_len
   do yj=1,y_len
      do zk=1,z_len
         read(2,*) dum,dum,dum,num_phi_z,dhdl_z
      
         dhdl_prob = dhdl_z * num_phi_z
         dhdl = dhdl + dhdl_prob
         norm = norm + num_phi_z      
      end do
      read(2,*)
   end do
   read(2,*)
end do
!print*, "norm =", norm
rewind(2)
do xi=1,x_len
   do yj=1,y_len
      do zk=1,z_len
      read(2,*) dum,dum,dum,num_phi_z,dhdl_z
      int_num_phi_z = num_phi_z/norm
      sum_int=sum_int+int_num_phi_z
      end do
      read(2,*)
   end do
   read(2,*)
end do
!print*, "Integration of A(z) =", sum_int

open(10,file='dhdl.dat',status='replace')
write(10,'("XXXX"X(1f18.9))')dhdl/norm
print*, "XXXX", dhdl/(norm)

end program

subroutine step_count(file_number,steps)
integer :: file_number, steps, ios,i
steps=0
do
 read(file_number,*,iostat=ios)
 if(ios.ne.0) exit
  steps=steps+1
end do
end subroutine

SUBROUTINE grid_min_max(num,grid_min,grid_max)
implicit none
integer :: num, ios
real*8 :: grid_min, grid_max
real*8 :: cv, dumm

rewind(num)
read(num,*,iostat=ios) dumm,cv
if (ios.ne.0) stop 'error reading colvar file'
grid_min=cv
grid_max=cv

rloop : DO
     read(num,*,iostat=ios)dumm,cv 
        if(ios.ne.0) exit rloop
     grid_min=MIN(cv,grid_min)
     grid_max=MAX(cv,grid_max)
end do rloop
end subroutine

