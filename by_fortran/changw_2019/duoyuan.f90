  program main 
  integer N,k,MM
      parameter (N=5) 
      parameter (k=3) 
      parameter (MM=4) 
c******** Note that  MM=K+1  ***************** 
      dimension y(N),x(k,N),a(MM),b(mm,mm),v(k) 
c     double precision	x,y,a,b,v,q,s,r,u 
 
      data x/1.1,2.0,3.2,1.0,2.0,3.2,1.2,1.8,3.0, 
     *	       1.1,1.9,2.9,0.9,2.1,2.9/ !3个预报因子，抽样5次
      data y/10.1,10.2,10.0,10.1,10.0/ !抽样5次 
	  call isqt2(x,y,k,mm,n,a,q,s,r,v,u,b,dyy) 
	  write(*,88) a(1) !把第一个系数b0输出
88	format(/1x,'b 0=',f9.5) 
	   do 89 j=2,mm !89代表把4个系数按顺序输出
89	    write(*,100) j-1,a(j) 
100	   format(1x,'b',i2,'=',f9.5) !100代表定义一个输出格式b=
cccccccccccccccccccccccccccccccccccccccccccccccccccccc 
	   write(*,20)q,s,r 
 20	   format(1x,'q=',f13.6,3x,'s=',f13.6,3x,'r=',f13.6) !20,22,30分别在定义变量输出格式
	   write(*,22)q,u,dyy 
 22	   format(1x,'q=',f13.6,3x,'u=',f13.6,3x,'dyy=',f13.6) 
	   write(*,30)(i,v(i),i=1,k) 
30	  format(1x,'v(',i2,')=',f13.6) 
	   write(*,40)u 
 40	   format(1x,'u=',f13.6) 
	    open(6,file='table') !打开了一个文件
	  write(6,180) 
 180	  format(/2x,'regression coefficients:') !斜杠代表下一行输出
	  write(6,88) a(1) !按顺序输入系数
	   do 189 j=2,mm 
189	     write(6,100) j-1,a(j) 
	  write(6,200) 
200	format(/1x,'Generic Analysis of Variance Table for the Multiple 
     * Linear Regression') 
	  write(6,202) 
202	 format(/1x,'----------------------------------------------------- 
     *---------------') 
	  write(6,204) 
204	format(/3x,'Source       df       SS           MS') 
	  write(6,202) 
	  write(6,206)	n-1,dyy 
206	format(/1x,'Total       n-1=',i2,'  SST=',f8.4) 
	  u2=u/real(K) 
	  write(6,208)	k,u,u2 
208   format(/1x,'Regression    K=',i2,'  SSR=',f8.4,'   MSR=SSR/K=' 
     *,f8.4) 
	   q2=q/real(n-k-1) 
	  write(6,209)	n-k-1,q,q2 
209   format(/1x,'Residual  n-k-1=',i2,'  SSE=',f8.4,'   MSE=SSE/(n-k-1) 
     *=',f8.4) 
	    f=(u/real(k))/(q/real(n-k-1)) 
	    write(6,220) f 
220	format(/1x,'                   F=MSR/MSE=',f9.4) 
	  write(6,202) 
	  close(6) 
	   stop 
	  end 
 
	  subroutine isqt2(x,y,m,mm,n,a,q,s,r,v,u,b,dyy) 
	   dimension x(m,n),y(n),a(mm),b(mm,mm),v(m) 
c	  real	x,y,a,b,v,q,s,r,u,yy,dyy,p,pp 
c	  double precision  x,y,a,b,v,q,s,r,u,yy,dyy,p,pp 
	  b(1,1)=n 
	  do 20 j=2,mm 
	    b(1,j)=0.0 
	    do 10 i=1,n 
10	    b(1,j)=b(1,j)+x(j-1,i) 
	    b(j,1)=b(1,j) 
20	    continue 
	  do 50 i=2,mm 
	   do 40 j=i,mm 
	     b(i,j)=0.0 
	     do 30 k=1,n 
30	     b(i,j)=b(i,j)+x(i-1,k)*x(j-1,k) 
	     b(j,i)=b(i,j) 
40	   continue 
50	     continue 
	    a(1)=0.0 
	   do 60 i=1,n 
60	    a(1)=a(1)+y(i) 
	   do 80 i=2,mm 
	    a(i)=0.0 
	    do 70 j=1,n 
70	    a(I)=a(i)+x(i-1,j)*y(j) 
80	  continue 
	  call achol(b,mm,1,a,l) 
	  yy=0.0 
	  do 90 i=1,n 
90	  yy=yy+y(i)/n 
	   q=0.0 
	   dyy=0.0 
	   u=0.0 
cccccccccccccccccccccccccccccccccc 
	   do 110 i=1,n 
	   p=a(1) 
	   do 100 j=1,m 
100	   p=p+a(j+1)*x(j,i) 
	   q=q+(y(i)-p)*(y(i)-p) 
	   dyy=dyy+(y(i)-yy)*(y(i)-yy) 
	   u=u+(yy-p)*(yy-p) 
 
110	   continue 
 
cccccccccccccccccccccccccccccccccc 
	   s=sqrt(q/n) 
	   r=sqrt(1.0-q/dyy) 
	   do 150 j=1,m 
	   p=0.0 
	   do 140 i=1,n 
	   pp=a(1) 
	   do 130 k=1,m 
	   if(k.ne.j)pp=pp+a(k+1)*x(k,i) 
130	   continue 
	   p=p+(y(i)-pp)*(y(i)-pp) 
140	   continue 
	   v(j)=sqrt(1.0-q/p) 
150	   continue 
	   return 
	   end 
 
	 subroutine achol(a,n,m,d,l) 
	 dimension a(n,n),d(n,m) 
c	   real a,d 
c	  double precision a,d 
	 l=1 
	 if(a(1,1)+1.0.eq.1.0)then 
	    l=0 
	    write(*,30) 
	    return 
	  endif 
	    a(1,1)=sqrt(a(1,1)) 
	  do 10 j=2,n 
10	      a(1,j)=a(1,j)/a(1,1) 
	  do 100 i=2,n 
	    do 20 j=2,i 
20	     a(i,i)=a(i,i)-a(j-1,i)*a(j-1,i) 
	     if(a(i,i)+1.0.eq.1.0)then 
	     l=0 
	     write(*,30) 
	     return 
	     endif 
30	   format(1x,'fail') 
	   a(i,i)=sqrt(a(i,i)) 
	   if(i.ne.n)then 
	   do 50 j=I+1,n 
	   do 40 k=2,i 
40	   a(i,j)=a(i,j)-a(k-1,i)*a(k-1,j) 
50	   a(i,j)=a(i,j)/a(I,i) 
	   endif 
100	   continue 
	   do 130 j=1,m 
	   d(1,j)=d(1,j)/a(1,1) 
	   do 120 i=2,n 
	   do 110 k=2,i 
110	      d(i,j)=d(i,j)-a(k-1,i)*d(k-1,j) 
	    d(i,j)=d(i,j)/a(i,i) 
120	       continue 
130	       continue 
	    do 160 j=1,m 
	    d(n,j)=d(n,j)/a(n,n) 
	    do 150 k=n,2,-1 
	    do 140 i=k,n 
140	       d(k-1,j)=d(k-1,j)-a(k-1,i)*d(i,j) 
	     d(k-1,j)=d(k-1,j)/a(k-1,k-1) 
150		continue 
160		continue 
	    return 
	     end 