; 4.10 读取 所有 时次文件,计算Fr和理查森数Ri
;    目的：绘制time变化图
; 4.15 将时次文件处理平均之后的时间剖面图
;    通过 restick属性设置XB的time label
;      参考函数：https://www.ncl.ucar.edu/Document/Functions/User_contributed/time_axis_labels.shtml
;    通过 grib_stime2time 将字符串时间转化为YYYYCCDDHHMMSS形式 用整除分段 cd_inv_cdalendar转化为Julian历
;      参考函数: https://www.ncl.ucar.edu/Document/Functions/Built-in/cd_inv_calendar.shtml
;    通过 gsn_csm_pres_hgt 绘制不同高度的剖面图
;      参考绘图集: https://www.ncl.ucar.edu/Applications/height_time.shtml 
; 6.15 插值 linint2要求经纬度单增
;参考https://www.ncl.ucar.edu/Document/Functions/Built-in/linint2.shtml
;代码路径  ncl /mnt/d/ncl_related/code/froude_number_xy.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


begin

;-------------1.变量数据读入--------------------

  path_rawdata = "/mnt/d/ncl_related/data/fnl/"
  path_picture ="/mnt/d/ncl_related/pictures/terrain/"  
  nc_files     = systemfunc("ls "+path_rawdata+"fnl_*_*_*.nc")
  f            = addfiles(nc_files,"r") ; f[0-12]
  ListSetType(f,"join")    ;打开多个文件并且返回文件句柄(是列表),指定如何返回文件变量
  f_time       =  f[:]->LAND_P0_L1_GLL0@initial_time; 注意原始数据中有些变量如 LAND_P0_L1_GLL0 有时间属性

  u_in         =f[:]->$"UGRD_P0_L100_GLL0"$  ;单位：m/s lv_ISBL0=26层
  v_in         =f[:]->$"VGRD_P0_L100_GLL0"$ ;26
  tmp_in       =f[:]->$"TMP_P0_L100_GLL0"$
  rh_in        =f[:]->$"RH_P0_L100_GLL0"$ ;lv_ISBL4=21层 %
  pres1        =f[:]->$"lv_ISBL4"$
  pres         =pres1(0,:) ;Pa
  pres_hpa     =pres*0.01
  pres_hpa@units="hPa"
  lon_in       =u_in&lon_0
  lat_in       =u_in&lat_0

  ; 修改维度名
  tmp_in!0 = "time"
  rh_in!0="time"

;---------------------2.时间格式的修改-----------------------------
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
  ; 时间格式重新写为 ju格式

  tmp_in&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)
  rh_in&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)



;----------------------------------------------网格点插值----------------------------------------
lon_out=fspan(80,110,301)
lat_out=fspan(0, 40, 401)
lon_out@units="degrees_east"
lat_out@units="degrees_north"
   
  u     =linint2_Wrap(lon_in, lat_in(::-1), u_in(:,5:25,::-1,:), True, lon_out, lat_out, 0)
  delete(u_in)
  v     =linint2_Wrap(lon_in, lat_in(::-1), v_in(:,5:25,::-1,:), True, lon_out, lat_out, 0)
  delete(v_in) 
  tmp   =linint2_Wrap(lon_in, lat_in(::-1), tmp_in(:,5:25,::-1,:), True, lon_out, lat_out, 0)
  delete(tmp_in)
  rh    =linint2_Wrap(lon_in, lat_in(::-1), rh_in(:,:,::-1,:), True, lon_out, lat_out, 0)
  delete(rh_in)

;----------------------------数据处理：湿Froude数----------------------------------------
  ; latstr=26
  ; latend=30
  ; lonstr=90
  ; lonend=98

  w = mixhum_ptrh((conform(tmp,pres_hpa,1)),tmp, rh, 1) ; 1代表混合比，mixing ratio (kg/kg),要求气压mb
  tvw= tmp*(1.0+0.61*w)  ;虚温
  copy_VarCoords(tmp, tvw)
  theta_v=pot_temp(pres, tvw, 1, False)  ; 虚位温，Pa,K

 ;接下来湿静力稳定度
  tdz_pres1= stdatmus_p2tdz(pres_hpa)    ; 将气压(要求hPa)转换成高度 tdz(3,9) ,tdz(2,:)-height
  tdz_pres =tdz_pres1(2,:)
  copy_VarCoords(pres, tdz_pres)


;法一直接差分各层，再取700-900平均--------------------------------------------------------
  dtdz     = center_finite_diff_n (theta_v,tdz_pres,False,0,1)
  g=9.8
  Nw2=g*dtdz/theta_v
  Nw=sqrt(Nw2) 
  copy_VarCoords(tmp, Nw)
  Nw_pres= Nw(time|:,{lv_ISBL0|70000:90000},LAT_0|:,LON_0|:) ;确定范围
  Nw_pres_ave= dim_avg_n_Wrap(Nw_pres, 1) ; 对lv_ISBL0做平均
; ;法二700-900的位温和高度再相除结果差不多----------------------------------------------------------
;   delta_theta_v= theta_v(:,{70000},:,:)-theta_v(:,{90000},:,:)
;   delta_h      = tdz_pres({70000})-tdz_pres({90000})
;   dtdz         = delta_theta_v/delta_h
;   copy_VarCoords_not_n(theta_v,dtdz,1)
;   g=9.8
;   theta_v_pres  = theta_v(time|:,{lv_ISBL0|70000:90000},LAT_0|:,LON_0|:) ;确定范围
;   theta_v_pres_ave= dim_avg_n_Wrap(theta_v_pres, 1) ; 对lv_ISBL0做平均
;   Nw2=g*dtdz/theta_v_pres_ave
;   Nw=sqrt(Nw2) 
;   copy_VarCoords(dtdz, Nw)
;   Nw_pres_ave= Nw(time|:,LAT_0|:,LON_0|:) ;确定范围

;700-900平均风速
  copy_VarCoords(tmp, u)
  copy_VarCoords(tmp, v)
  ; u_pres= u(time|:,{lv_ISBL0|70000:90000},LAT_0|:,LON_0|:) ;确定范围
  ; u_pres_ave= dim_avg_n_Wrap(u_pres, 1) ; 对lv_ISBL0做平均
  v_pres= v(time|:,{lv_ISBL0|70000:90000},LAT_0|:,LON_0|:) ;确定范围
  v_pres_ave= dim_avg_n_Wrap(v_pres, 1) ; 对lv_ISBL0做平均
  ; U=sqrt(u_pres_ave^2+v_pres_ave^2)
  ; copy_VarCoords_not_n(tmp, U,1)
  h=3000.
  Fr=v_pres_ave/Nw_pres_ave/h
  copy_VarCoords_not_n(tmp, Fr,1)

  ; Fr_area= Fr(time|:,{LAT_0|20:25},{LON_0|90:98}) ;确定范围
  ; Fr_ave= dim_avg_n_Wrap(Fr_area, (/1,2/)) ; 对lv_ISBL0做平均
  Fr_ave=Fr(:,{26.0},{96.0}) ;2700-0.3356853 2706-0.3860504
  ; print(Fr_ave(4)) 
  ; exit

; ;-----------------------计算理查森数---------------------
;   dudz= center_finite_diff_n (u,tdz_pres,False,0,1)
;   dvdz= center_finite_diff_n (v,tdz_pres,False,0,1)
;   copy_VarCoords(tmp, dudz)
;   copy_VarCoords(tmp, dvdz)
;   dUdz2=dudz^2+dvdz^2
;   dUdz2_inv = 1. / where(dUdz2.ne.0,dUdz2,dUdz2@_FillValue)
;   Ri=Nw2*dUdz2_inv
;   copy_VarCoords(tmp, Ri)
;    Ri_area= Ri(time|:,lv_ISBL0|:,{LAT_0|25:35},{LON_0|90:98}) ;确定范围
;   ; Ri_ave= dim_avg_n_Wrap(Ri_area, (/2,3/)) ; 对lv_ISBL0做平均 
;   print(min(Ri_area(4,:,:,:)))
;   exit






; ;将变量整合到一个数组
;  data  = new((/2,12/),float)

;  data(0,:) = vvel
;  data(1,:) = w_area_ave
;  data!1="time"
;  data&time=u_sig&time
;  ;print(data&time)
;  ;exit 



wks  = gsn_open_wks ("png",path_picture+"froude_number_time_v")          ; send graphics to PNG file



;---------------------------------------绘图属性设置------------------------------------------
res = True
   res@gsnDraw  = False
   res@gsnFrame = False

 
  ;res@gsnMaximize=True
  res@tmXTOn = False ;关刻度
  res@tmYROn = False
  ;res@gsnLeftString = "hPa" ;左上角字符
  res@tiYAxisString="Moist Froude number"
  res@tiYAxisFont=22
  ;设置坐标轴
  res@tmXBLabelFont   =22
  res@tmXBLabelAngleF = 48.   ;角度
  res@tmYLLabelFont   =22
  res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
  res@tmYLLabelFontHeightF = 0.02
  ; ;设置画图区域大小
  res@vpWidthF =0.7
  res@vpHeightF=0.4

  ;>>>时间轴设置
  trestick = True
  trestick@ttmFormat   = "%D/%H" ; 时间轴的格式
  trestick@ttmAxis    = "XB"
  trestick@ttmNumTicks = 12    ; 显示所有12个点


frouderes=res
  frouderes@xyDashPattern  = 0                         ; Make curves all solid
  frouderes@xyLineThicknessF= 8
  frouderes@xyMarkLineMode = "MarkLines"               ; Markers *and* lines
  frouderes@xyMarker       = 16               ; 2 different markers
  frouderes@xyMarkerColors ="black"; 2 different colors
  frouderes@xyLineColors   ="black"
  ; frouderes@trYMinF=min(Fr_ave)-0.5



;------------------------------------------绘图-----------------------------------------------------


time_axis_labels(Fr_ave&time,frouderes,trestick) ; call the formatting procedure 
plot  = gsn_csm_xy (wks,Fr_ave&time,Fr_ave,frouderes) ; create plot

; ;-------------------------------------------------------------------------------------
; ; 绘制图例
;   resLg = True
;   colours=(/"blue","red"/)
 
;   resLg@lgItemType = "MarkLines"

;   resLg@lgMonoMarkerIndex = True
;   resLg@lgMarkerColors = colours
;   resLg@lgMarkerIndex = 1
;   resLg@lgMarkerSizeF = 0.04
 
 
;   resLg@lgMonoDashIndex = True
;   resLg@lgDashIndex = 0
;   resLg@lgLineColors = colours
;   resLg@lgLineThicknessF = 3
 
;   resLg@vpWidthF = 0.14
;   resLg@vpHeightF = 0.1
 
;   resLg@lgPerimFill = 0
;   resLg@lgPerimFillColor = "Background"
 
;   resLg@lgLabelFontHeightF = 0.08
 
;   ;resLg@lgTitleFontHeightF = 0.015
;   ;resLg@lgTitleString = ""
 
;   lbid = gsn_create_legend(wks, 2, (/"W","W~B~f~N~"/), resLg)
 
; ; 将图例放置在图中
;   amres = True
;   amres@amParallelPosF = 0.398
;   amres@amOrthogonalPosF = -0.37
;   dumLg = gsn_add_annotation(plot, lbid, amres)

draw(plot)
frame(wks)
end