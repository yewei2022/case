; 4.10 读取 所有 时次文件的TMP_P0_L100_GLL0变量
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
;9.24 http://www.ncl.ucar.edu/Document/Functions/Built-in/ind.shtml

;load "~/offical_sample_script/download/pvs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

;1. 数据读入
  path_rawdata = "/cygdrive/D/ncl_related/data/tmp/"
  path_picture="/cygdrive/D/ncl_related/pictures/tmp/" 
  nc_files = systemfunc("ls "+path_rawdata+"fnl_*_*_*.nc"); 
  f = addfiles(nc_files,"r") ; f[0-8],9个文件
  ListSetType(f,"join"); 打开多个文件并且返回文件句柄(是列表),指定如何返回文件变量
  f_time =  f[:]->LAND_P0_L1_GLL0@initial_time; 注意原始数据中有些变量如 LAND_P0_L1_GLL0 有时间属性
  ; print(f_time) ; 是每个文件的具体时间,具有单位, initial_time : 10/27/2008 (00:00)

  tmpdata=f[:]->$"TMP_P0_L100_GLL0"$ ;第一个维度是 join 格式（列表），下面将其替换为f_time
  tmpdata!0 = "time" ;修改维度名
  tmpdata!1 = "lv"
  tmp    =tmpdata-273.15
  copy_VarCoords(tmpdata, tmp)

  pres11  =f[:]->$"lv_ISBL0"$  ;Pa,26层

  ;printVarSummary(pres1) ;输出9*26
  ;exit
  pres1 =pres11(0,:)
  ;printVarSummary(pres1) ;输出是一维[lv_ISBL0 | 26]
  ;exit
  ;pres   =pres1/100
  ;pres@units="hPa"


;2. 时间格式的修改
  ; tmp&time = f_time
  time_set = grib_stime2itime(f_time) ; 将grib格式时间由字符串转化为整型
  time_set_yy = time_set/1000000
  time_set_mm = time_set/10000%100
  time_set_dd = time_set/100%100
  time_set_hh = time_set%100
  time_set_min = time_set_hh
  time_set_sec = time_set_hh
  time_set_min = 0
  time_set_sec = 0
  ; print(time_set_hh)
  ; print(time_set)
  ; tmp的时间格式重新写为 ju格式
  tmp&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)
  ; print(tmp&time) ; [time | 9] x [lv_ISBL0 | 26] x [lat_0 | 181] x [lon_0 | 360]

;3. 选定站点
  tmp_pointcha = tmp(time|:,lv|:,{lat_0|28},{lon_0|97})
  tmp_cha = tmp_pointcha(lv|:,time|:)

  tmp_pointcuo = tmp(time|:,lv|:,{lat_0|28},{lon_0|92})
  tmp_cuo = tmp_pointcuo(lv|:,time|:)

  ;tmp_pointmi = tmp(time|:,lv|:,{lat_0|29},{lon_0|94}) ;米林
  ;tmp_pointmi_plot = tmp_pointmi(lv|:,time|:)

  ;波密的经纬度
  tmp_pointmi = tmp(time|:,lv|:,{lat_0|30},{lon_0|96})
  tmp_mi = tmp_pointmi(lv|:,time|:)




;4. 等压面插值
  pres2=fspan(1000,100000,3961)
  pres2@units="Pa"
  linlog=2
  t_out_cha = int2p_n_Wrap(pres1,tmp_cha,pres2,linlog,0)
  t_out_cuo = int2p_n_Wrap(pres1,tmp_cuo,pres2,linlog,0)
  t_out_mi = int2p_n_Wrap(pres1,tmp_mi,pres2,linlog,0)

  ;printVarSummary(t_out_cha) ;[LV | 1981] x [time | 9]
  ;exit

;5. 找到0度所在海拔高度(altitude)
 pres_cha=new(9,float)
 pres_cuo=new(9,float)
 pres_mi =new(9,float)

  do i=0,8
   x1D  = t_out_cha(:,i); convert to 1D array
   x2D  = t_out_cuo(:,i); convert to 1D array
   x3D  = t_out_mi(:,i); convert to 1D array

   ;printVarSummary(t_out_cha(:,0)) ;输出1维 [LV | 1981]
   ;exit
   ;printVarSummary(x1D) ;输出1维 [1981]
   ;exit
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
 pres_cha=pres_cha*0.01
 pres_cha@units="hPa"
 pres_cuo=pres_cuo*0.01
 pres_cuo@units="hPa"
 pres_mi=pres_mi*0.01
 pres_mi@units="hPa"

 tdz_cha = stdatmus_p2tdz(pres_cha)    ; ===> tdz(3,9) 
 tdz_cuo = stdatmus_p2tdz(pres_cuo)    ; ===> tdz(3,9) 
 tdz_mi = stdatmus_p2tdz(pres_mi)    ; ===> tdz(3,9) 

 ;print(tdz_cha)
 ;exit




;7. 合并数据，并设置时间属性
  data      = new((/3,9/),float)
  data(0,4:7) = tdz_cuo(2,4:7)-4281
  data(1,:) = tdz_cha(2,:)-2331
  data(2,:) = tdz_mi(2,:)-2737
  data!1="time"
  data&time=tmp&time


;8. 
  wks  = gsn_open_wks ("png", path_picture+"time_h_0")          ; send graphics to PNG file


 ;------------------------------绘图属性设置------------------------------------------
res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize=True
  res@tmXTOn = False ;关刻度
  res@tmYROn = False
  res@tiYAxisString="H/m"
  res@tiYAxisFont=22 ;y轴标题字体
  res@tmXBLabelFont=22
  res@tmXBLabelAngleF = 48.   ;角度
  res@tmYLLabelFont=22
  res@tmXBLabelFontHeightF = 0.015 ; x轴标签字号
  res@tmYLLabelFontHeightF =0.015

  res@trYMaxF=max(data)+1000;扩大y轴便签显示的最大值
  res@trYMinF=0

  ;设置画图区域大小
  res@vpWidthF=0.5
  res@vpHeightF=0.4

  ;时间轴设置
  trestick = True
  trestick@ttmFormat = "%D/%H" ; 时间轴的格式
  trestick@ttmAxis   = "XB"
  trestick@ttmNumTicks = 9 ; 显示所有9个点


  Tres=res
  ;Tres@xyDashPatterns  = (/1,2,3/)                        ; Make curves all solid
  Tres@xyDashPattern =0
  Tres@xyLineThicknessF= 8
  Tres@xyMarkLineMode = "MarkLines"               ; Markers *and* lines
  Tres@xyMarker =16
  Tres@xyMarkerThicknessF=8
  Tres@xyLineColors =(/"blue","red","purple"/) ;cuo,cha,bo

  ;图例
  Tres@pmLegendDisplayMode    = "Always"                              ;add lengend
  Tres@xyExplicitLegendLabels = (/"Cuona","Chayu","Bomi"/)          ;lengend title by english
  Tres@pmLegendWidthF         = 0.11                                   ;lengend width
  Tres@pmLegendHeightF        = 0.1                                   ;lengend height

  Tres@lgLabelFontHeightF = 0.015

  ;放在图中
  Tres@pmLegendOrthogonalPosF  = -1.21                              ;
  Tres@pmLegendParallelPosF    = 0.865


;绘图
time_axis_labels(data&time,Tres,trestick) ; call the formatting procedure 
plot  = gsn_csm_xy (wks,data&time,data,Tres) ; create plot


draw(plot)
frame(wks)
end
