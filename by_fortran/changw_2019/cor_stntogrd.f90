program corstntogrd
real cor(160)
open(1,file="D:\YEWEI\practice6\cor.dat",form='formatted')
do i=1,160
read(1,*)cor(i)
end do
close(1)
call stntogrd(cor)
end
subroutine stntogrd(x)
real lat(160),lon(160),x(160)
character *8,stid(160)
open(2,file='D:\YEWEI\practice6\china.dat')
do 20 k=1,160
20 read(2,'(f5.2,2x,f6.2)')lat(k),lon(k)
close(2)
do 2 i=1,160
2 stid(i)=char(i)
open(3,file='D:\YEWEI\practice6\cor_gr.grd',form='binary')
tim=0.0
nlev=1
nflag=1
do 40 i=1,160
write(3) stid(i),lat(i),lon(i),tim,nlev,nflag,x(i)
40 continue
nlev=0
write(3) stid(i-1),lat(i-1),lon(i-1),tim,nlev,nflag
close(3)
return
end
