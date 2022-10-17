      program yewei
      parameter(nx=37,ny=17,mo=12,yr=4)
      real var(nx,ny,mo,yr),cli(nx,ny,mo),ano(nx,ny,mo,yr),mean(nx,ny,mo)
      integer i,j,m,y,irec
      open(1,file='D:\YEWEI\practice1\h500.dat')
       do  y=1,4
       do  m=1,12
       read(1,1000)   
       read(1,3000) ((var(i,j,m,y),i=1,nx),j=1,ny)
       enddo
       enddo
1000   format(2i7)
3000   format(37f8.1)
      close(1)
       do i=1,nx
         do j=1,ny
          do m=1,mo
             cli(i,j,m)=0
           do y=1,yr
             cli(i,j,m)=cli(i,j,m)+ var(i,j,m,y)
           end do
             cli(i,j,m)=cli(i,j,m)/yr
          end do
       end do 
       end do
        open(2,file='D:\YEWEI\practice1\varcli.grd',form='unformatted',access='direct',recl=nx*ny*4 )
              irec=0
         do m=1,mo
              irec=irec+1
         write(2,rec=irec)((cli(i,j,m),i=1,nx),j=1,ny)
         enddo
         close(2)
        do i=1,nx
       do j=1,ny
       do m=1,mo
        mean(i,j,m)=0
       do y=1,yr
       ano(i,j,m,y)=var(i,j,m,y)-cli(i,j,m)
       mean(i,j,m)=mean(i,j,m)+ano(i,j,m,y)**2
       end do
        mean(i,j,m)=sqrt(1/yr*mean(i,j,m))
       end do
       end do
       end do
      open(3,file='D:\YEWEI\practice1\ano.grd',form='unformatted',access='direct',recl=nx*ny*4 )
              irec=0
         do  y=1,yr
         do m=1,mo
              irec=irec+1
         write(3,rec=irec)((ano(i,j,m,y),i=1,nx),j=1,ny)
         end do
         end do
         close(3)
         open(4,file='D:\YEWEI\practice1\mean.grd',form='unformatted',access='direct',recl=nx*ny*4)
              irec=0
         do m=1,mo
              irec=irec+1
         write(4,rec=irec)((mean(i,j,m),i=1,nx),j=1,ny)
         enddo
         close(4)
       end
