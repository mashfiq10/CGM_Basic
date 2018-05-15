!Sk. Mashfiqur Rahman
!Oklahoma State University
!conjugate gradient method basic

program CGM_example
implicit none
integer::NA,k,i,nx,ny
double precision,allocatable::b(:),x(:),r(:),p(:),d(:)
double precision::alpha,beta,eps
double precision::m,z,y

open(7,file='input.txt')
read(7,*)nx 	!number of rows
read(7,*)NA 	!total no. of iterations
close(7)

ny=nx

allocate(b(0:nx),x(0:nx),r(0:nx),p(0:nx),d(0:nx))

do i=0,nx
  write(*,*)'Enter the constants:'
  write(*,*)
  read(*,*) b(i)
end do

do i=0,nx
p(i)=b(i)
r(i)=b(i)
end do

do i=0,nx
  x(i)=0.0d0
end do

eps=1.0d-10

do k=1,NA

call matrix(nx,ny,p,d) 

z=0.0d0
y=0.0d0 
do i=0,nx
z= z+ (r(i) * r(i))
y= y+ (d(i) * p(i))
end do

alpha = z/y

do i=0,nx
x(i) = x(i) + alpha * p(i)
end do

m=0.0d0
do i=0,nx
r(i) = r(i)-(alpha * d(i))
m= m + r(i) * r(i)
end do
beta = m/z

do i=0,nx
p(i) = r(i)+ (beta * p(i))
end do

!stopping criteria- the iteration will close whenever r(1:2) will be 0
if (r(0).lt.eps .and. r(1).lt.eps .and. r(2).lt.eps .and. r(3).lt.eps) then
   exit
end if
end do
 
write(*,*) "solution is: ", x

end


subroutine matrix(nx,ny,wf,w)
implicit none
double precision,dimension(0:nx,0:ny)::A
double precision,allocatable:: u(:)
double precision,dimension(0:nx)::w,wf
double precision::p
integer::i,j,nx,ny

allocate(u(0:nx))
do i=0,nx
 u(i) = wf(i)
end do

do i=0,nx
  do j=0,ny
    if (i.eq.j) then
      A(i,j)= 4.0d0
    else if (j.eq.i+1) then
      A(i,j)= 1.0d0
    else if (i.eq.j+1) then
      A(i,j)= 2.0d0
    else 
      A(i,j)=0.0d0
    end if
  end do
end do

p=0.0d0
do i=0,nx
  do j=0,ny
    p=p+A(i,j)*u(j)
  end do
  w(i)=p
  p=0.0d0
end do

deallocate(u)
return
end