; 4.10 读取 ERA5数据
;    目的：绘制0℃距离地面的高度随时间变化图
; 4.15 将时次文件处理之后的时间剖面图
;    通过 restick属性设置XB的time label
;      参考函数：https://www.ncl.ucar.edu/Document/Functions/User_contributed/time_axis_labels.shtml
;    通过 grib_stime2time 将字符串时间转化为YYYYCCDDHHMMSS形式 用整除分段 cd_inv_cdalendar转化为Julian历
;      参考函数: https://www.ncl.ucar.edu/Document/Functions/Built-in/cd_inv_calendar.shtml
;    通过 gsn_csm_pres_hgt 绘制不同高度的剖面图
;      参考绘图集: https://www.ncl.ucar.edu/Applications/height_time.shtml 
;9.17 pressure插值函数，int2p http://www.ncl.ucar.edu/Document/Functions/Built-in/int2p.shtml
;9.18 将气压换算成高度 http://www.ncl.ucar.edu/Document/Functions/Built-in/stdatmus_p2tdz.shtml
;  ncl /mnt/h/d/ncl_related/ERA5_code/time_h_0_label.ncl


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

;1. 数据读入
  pic_dir  ="/mnt/h/d/ncl_related/pictures/tmp/ERA5/"
  data_dir = "/mnt/h/d/data/ERA5/"    
  data1    = addfile(data_dir+"ERA5.nc", "r") 
  tmp      = short2flt(data1 ->t(0:48:6,{10:1000},{10:35},{85:105}))
  pres1    = data1 ->level({10:1000})
  time     = data1 ->time(0:48:6)
  tmp      =tmp-273.15

;3. 选定站点
  tmp_pointcha = tmp(time|:,level|:,{latitude|28.75},{longitude|97.5})
  tmp_cha = tmp_pointcha(level|:,time|:)

  tmp_pointcuo = tmp(time|:,level|:,{latitude|28.0},{longitude|92.0})
  tmp_cuo = tmp_pointcuo(level|:,time|:)

  tmp_pointmi = tmp(time|:,level|:,{latitude|29.75},{longitude|95.75})
  tmp_mi = tmp_pointmi(level|:,time|:)

;4. 等压面插值，间隔0.25hPa
  pres2=fspan(10,1000,3961)
  pres2@units="hPa"
  linlog=2
  t_out_cha = int2p_n_Wrap(pres1,tmp_cha,pres2,linlog,0)
  t_out_cuo = int2p_n_Wrap(pres1,tmp_cuo,pres2,linlog,0)
  t_out_mi = int2p_n_Wrap(pres1,tmp_mi,pres2,linlog,0)

  ; printVarSummary(t_out_cha) ;[LV | 1981] x [time | 9]
  ; exit

;5. 找到0度所在海拔高度(altitude)
 pres_cha=new(9,float)
 pres_cuo=new(9,float)
 pres_mi =new(9,float)

  do i=0,8
   x1D  = t_out_cha(:,i); convert to 1D array
   x2D  = t_out_cuo(:,i); convert to 1D array
   x3D  = t_out_mi(:,i); convert to 1D array

   ; printVarSummary(t_out_cha(:,0)) ;输出1维 [LV | 1981]
   ; exit
   ; printVarSummary(x1D) ;输出1维 [1981]
   ; exit
   copy_VarCoords_not_n(t_out_cha,pres2,1);除了time以外的坐标信息
   copy_VarCoords_not_n(t_out_cha,x1D,1);除了time以外的坐标信息
   copy_VarCoords_not_n(t_out_cha,x2D,1);除了time以外的坐标信息
   copy_VarCoords_not_n(t_out_cha,x3D,1);除了time以外的坐标信息
   ;printVarSummary(pres2) ;输出1维 [LV | 1981]
   ;exit
   ;printVarSummary(x1D) ;输出1维[LV | 1981]
   ;exit
   ;print(ind((x1D.lt.0.01).and.(x1D.gt.(-0.01)))) ;输出（0）1188
   ;exit
   ;print(x3D)
   ;exit
   icha=ind((x1D.lt.0.02).and.(x1D.gt.(-0.02))) ;0线，0.02--0.02
   icuo=ind((x2D.lt.0.02).and.(x2D.gt.(-0.02)))  
   imi=ind((x3D.lt.0.02).and.(x3D.gt.(-0.02)))

   pres_cha(i)=pres2(icha(0))
   pres_cuo(i)=pres2(icuo(0));(1214)  -0.02356793  (1215)  0.01488853
   pres_mi(i)=pres2(imi(0));同理，没有0.01范围内的
   delete(icha)
   delete(icuo)
   delete(imi)
   delete(x1D)
   delete(x2D)
   delete(x3D)
   ;print(pres_cha(0)) ;60400pa
   ;exit
  end do

;6. 换算成那个高度
 pres_cha@units="hPa"
 pres_cuo@units="hPa"
 pres_mi@units="hPa"
 tdz_cha = stdatmus_p2tdz(pres_cha)    ; ===> tdz(3,9) ;f要求hPa
 tdz_cuo = stdatmus_p2tdz(pres_cuo)    ; ===> tdz(3,9) 
 tdz_mi = stdatmus_p2tdz(pres_mi)    ; ===> tdz(3,9) 

 ;print(tdz_cha)
 ;exit

;7. 合并数据，并设置时间属性
  data = new((/3,9/),float)
  cuo  = new(9,float)
  cha  = new(9,float)
  mi   = new(9,float)
  cuo(4:6)  =tdz_cuo(2,4:6)-4281
  mi   =tdz_mi(2,:)-2737
  cha  =tdz_cha(2,:)-2331

  data(0,:) = cuo
  data(1,:) = cha
  data(2,:) = mi
  data!1="time"
  data&time=int2dble(time)

  ;将数据分段,画点,注意，0度线，错那只有一段
  data_split  = new((/6,9/),float)

  data_split(0,4:6)  = cuo(4:6) ;花
  data_split(1,0)    = cha(0) ;圈
  data_split(2,1:8)  = cha(1:8) ;点
  data_split(3,0:1)  = mi(0:1) ;圈
  data_split(4,2)    = mi(2) ;点
  data_split(4,5:8)  = mi(5:8) 
  data_split(5,3:4)  = mi(3:4) ;花
 
wks  = gsn_open_wks ("ps", pic_dir+"time_h_0_label")
gsn_define_colormap(wks,"wgne15")
colors=(/14,8,4/) ;红 绿 蓝

;属性设置
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize=True
  res@tmXTOn = False ;关刻度
  res@tmYROn = False
  res@tiYAxisString="H (m)"
  res@tmBorderThicknessF  =3   ; border 6x as thick,图框粗细
  res@tmXBMajorThicknessF =3
  res@tmXBMinorThicknessF =3
  res@tmYLMajorThicknessF =3
  res@tmYLMinorThicknessF =3  
  res@tmXBLabelAngleF = 48.   ;角度
  res@tmEqualizeXYSizes =True
  res@tmXBLabelFontHeightF = 0.015 ; x轴标签字号

  res@gsnLeftStringFontHeightF=0.015
  res@gsnLeftString="(a)"
  res@gsnLeftStringOrthogonalPosF=-0.1
  res@gsnLeftStringParallelPosF =0.01
  res@txPerimOn=True
  res@txPerimThicknessF=2.5

  res@trYMaxF=2000;扩大y轴便签显示的最大值
  res@trYMinF=0

  ;设置画图区域大小
  res@vpWidthF=0.5
  res@vpHeightF=0.4

  Tres=res
  Tres@xyMarkLineMode = "Markers"               ; Markers *and* lines
  Tres@xyMarkers =(/3,4,16,4,16,3/)
  Tres@xyMarkerThicknessF=5
  Tres@xyMarkerSizeF=10

  ;时间轴设置
  trestick = True
  trestick@ttmFormat = "%D/%H" ; 时间轴的格式
  trestick@ttmAxis   = "XB"
  trestick@ttmNumTicks = 9 ; 显示所有9个点


;绘点
time_axis_labels(data&time,Tres,trestick) ; call the formatting procedure 
plot  = gsn_csm_xy (wks,data&time,data_split,Tres) ; create plot

;绘线
;把察隅和波密合在一起画，错那单独画
  nline=2
  nrow=9
  dumLine = new((/nline,nrow/), graphic) 
  x= new(nrow, typeof(data&time))
  y= new(nrow, typeof(data&time))

  ;属性设置
  resLine = True
  resLine@gsLineThicknessF = 8
  resLine@tfPolyDrawOrder="PreDraw"

  ;开始画线
  x=data&time
  do j=0,nline-1
    y=data(j+1,:)
    if (j.eq.0) then
      resLine@gsLineColor = colors(2)
    end if
    if (j.eq.1) then
      resLine@gsLineColor = colors(1)
    end if 
  do i = 0, nrow-2
    xx = (/ x(i), x(i+1)/)
    yy = (/ y(i), y(i+1)/)
    dumLine(j,i) = gsn_add_polyline(wks, plot, xx, yy, resLine)
  end do
  delete(y)
  end do

;单独画错那
  delete(x)
  nrow_c=2
  dumLine_c = new(nrow_c, graphic) 
  x= new(nrow_c, typeof(data&time))
  y= new(nrow_c, typeof(data&time))
  x=data&time(5:6)
  y=data(0,5:6)

  resLine_c = True
  resLine_c@gsLineThicknessF = 8
  resLine_c@gsLineColor=colors(0)
  resLine_c@tfPolyDrawOrder="PreDraw"
  do i = 0, nrow_c-2
    xx = (/ x(i), x(i+1)/)
    yy = (/ y(i), y(i+1)/)
    dumLine_c(i) = gsn_add_polyline(wks, plot, xx, yy, resLine_c)
  end do

 
; 绘制图例
  lgcolours = (/0,0,0,14,8,4/)

  resLg = True
  resLg@lgMonoItemType=False
  resLg@lgItemTypes = (/"Markers","Markers","Markers","Lines","Lines","Lines"/)
  resLg@lgMarkerColors=(/"black","black","black","black","black","black"/)
  resLg@lgMarkerIndexes=(/4,3,16,0,0,0/)
  resLg@lgMarkerSizeF = 0.02  
  resLg@lgMarkerThicknessF=5
  resLg@lgMonoDashIndex=True ;一定要把这个打开！不然线型不能一致
  ;resLg@lgLineDashIndexes=(/0,0,0,0,0,0/) ;这个好像打开也没用
  resLg@lgLineColors = lgcolours
  resLg@lgLineThicknessF = 8 
  resLg@vpWidthF = 0.25
  resLg@vpHeightF = 0.2
  resLg@lgLabelFontHeightF = 0.15
  resLg@lgPerimOn=False
 
  lbid = gsn_create_legend(wks, 6, (/"No precipitation",\
  "Snowfall","Rainfall","Cona","Bomi","Zayu"/), resLg)
 
; 将图例放置在图中
  amres = True
  amres@amParallelPosF = -0.32 ;左右，正右左负
  amres@amOrthogonalPosF = 0.32 ;上下，正下负上
  dumLg = gsn_add_annotation(plot, lbid, amres)


draw(plot)
frame(wks)
end
