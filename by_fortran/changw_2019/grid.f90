program main
parameter (nx=71,ny=41)
real lat(ny),lon(nx)
real s(nx,ny)
open(1,file='D:\YEWEI\practice6\grid.grd',form='binary')
lat(1)=15.0
lon(1)=70.0
do j=1,ny-1
lat(j+1)=lat(j)+1.0
end do
do i=1,nx-1
lon(i+1)=lon(i)+1.0
end do
do i=1,nx
do j=1,ny
s(i,j)=1
end do
end do
write(1)s
close(1)
end
