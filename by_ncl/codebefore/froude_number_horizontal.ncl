; 湿Froude数，先画平面图，看哪块最大
; 4.10 读取 所有 时次文件
;    目的：绘制time变化图， lat 26~32 lon 92~100 区域平均
; 4.15 将时次文件处理平均之后的时间剖面图
;    通过 restick属性设置XB的time label
;      参考函数：https://www.ncl.ucar.edu/Document/Functions/User_contributed/time_axis_labels.shtml
;    通过 grib_stime2time 将字符串时间转化为YYYYCCDDHHMMSS形式 用整除分段 cd_inv_cdalendar转化为Julian历
;      参考函数: https://www.ncl.ucar.edu/Document/Functions/Built-in/cd_inv_calendar.shtml
;    通过 gsn_csm_pres_hgt 绘制不同高度的剖面图
;      参考绘图集: https://www.ncl.ucar.edu/Applications/height_time.shtml 
; 6.15 插值 linint2要求经纬度单增
;参考https://www.ncl.ucar.edu/Document/Functions/Built-in/linint2.shtml
;load "~/offical_sample_script/download/pvs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


begin

;-------------1.变量数据读入--------------------

  path_rawdata = "/cygdrive/D/ncl_related/data/fnl/"
  path_picture ="/cygdrive/D/ncl_related/pictures/terrain/"  
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

  w = mixhum_ptrh((conform(tmp,pres_hpa,1)),tmp, rh, 1) ; 1代表混合比，mixing ratio (kg/kg),要求气压mb
  tvw= tmp*(1.0+0.61*w)  ;虚温
  copy_VarCoords(tmp, tvw)
  theta_v=pot_temp(pres, tvw, 1, False)  ; 虚位温，Pa,K
 ;接下来湿静力稳定度
  tdz_pres = stdatmus_p2tdz(pres_hpa)    ; 将气压(要求hPa)转换成高度 tdz(3,9) ,tdz(2,:)-height
  dtdz     = center_finite_diff_n (theta_v,tdz_pres(2,:),False,0,1)
  g=9.80665
  Nw=g/theta_v*dtdz
  Nw=sqrt(Nw) ; [time | 12] x [lv_ISBL0 | 21] x [LAT_0 | 401] x [LON_0 | 301]
  copy_VarCoords(tmp, Nw)
  Nw_pres= Nw(time|:,{lv_ISBL0|70000:90000},LAT_0|:,LON_0|:) ;确定范围
  Nw_pres_ave= dim_avg_n_Wrap(Nw_pres, 1) ; 对lv_ISBL0做平均
;700-900平均风速
  copy_VarCoords(tmp, u)
  copy_VarCoords(tmp, v)
  u_pres= u(time|:,{lv_ISBL0|70000:90000},LAT_0|:,LON_0|:) ;确定范围
  u_pres_ave= dim_avg_n_Wrap(u_pres, 1) ; 对lv_ISBL0做平均
  v_pres= v(time|:,{lv_ISBL0|70000:90000},LAT_0|:,LON_0|:) ;确定范围
  v_pres_ave= dim_avg_n_Wrap(v_pres, 1) ; 对lv_ISBL0做平均
  U=sqrt(u_pres_ave^2+v_pres_ave^2)
  copy_VarCoords_not_n(tmp, U,1)
  h=3000.
  Fr=U/Nw_pres_ave/h
  copy_VarCoords(U,Fr)
  Fr_area= Fr(time|:,{LAT_0|20:30},{LON_0|88:98}) ;确定范围
  ; printVarSummary(Fr_area)
  ; exit






; ;将变量整合到一个数组
;  data  = new((/2,12/),float)

;  data(0,:) = vvel
;  data(1,:) = w_area_ave
;  data!1="time"
;  data&time=u_sig&time
;  ;print(data&time)
;  ;exit 



wks  = gsn_open_wks ("png",path_picture+"froude_number_panel")          ; send graphics to PNG file



;---------------------------------------绘图属性设置------------------------------------------
res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnAddCyclic=False

 
  ;res@gsnMaximize=True
  res@tmXTOn = False ;关刻度
  res@tmYROn = False
  ;设置坐标轴
  res@tmXBLabelFont=22
  res@tmYLLabelFont=22
  res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
  res@tmYLLabelFontHeightF =0.02



frouderes=res
  ;设置地图
  frouderes = res
  frouderes@mpFillOn                = False
  frouderes@mpMinLatF               = 18.
  frouderes@mpMaxLatF               = 32.
  frouderes@mpMinLonF               = 85.
  frouderes@mpMaxLonF               = 100.
  frouderes@gsnMajorLonSpacing      =2. 
  frouderes@gsnMinorLonSpacing      =1.    
  frouderes@gsnMajorLatSpacing      =2. 
  frouderes@gsnMinorLatSpacing      =1
  frouderes@mpDataBaseVersion       = "MediumRes"
  frouderes@mpDataSetName           = "Earth..4"
  frouderes@mpOutlineSpecifiers     = (/"China:states"/)

  frouderes@cnFillOn = True ;cn means ContourPlot
  frouderes@cnLinesOn = False   
  frouderes@cnFillPalette = "WhiteGreen"
  frouderes@cnLineLabelsOn = False  
  frouderes@cnInfoLabelOn = False ; set false will not set an information label 
  ; frouderes@cnLevelSelectionMode = "ManualLevels" 
  ; frouderes@cnMinLevelValF       = 0.5
  ; frouderes@cnLevelSpacingF      = 0.1              ; set contour spacing
  frouderes@lbLabelBarOn         = False
  ; frouderes@lbBoxLinesOn         = True ; Label bar
  ; frouderes@lbLabelAutoStride    = True            ; let NCL determine label spacing
  ; frouderes@lbLabelFontHeightF   = 0.02
  ; frouderes@lbOrientation        = "Vertical"         ; vertical label bar
  ; frouderes@pmLabelBarWidthF     = 0.03     
  ; frouderes@lbBoxEndCapStyle     = "TriangleBothEnds"

  panel_res                       = True
  panel_res@gsnPanelMainString    = "moist froud number" ; use "txString" for older versions
  panel_res@gsnPanelLabelBar      =True
  ;panel_res@gsnPanelRowSpec = True 
  ; If the resource gsnPanelRowSpec is set to True,
  ; then dims will represent the number of plots per row.



;------------------------------------------绘图-----------------------------------------------------
 plot = new(8,graphic)
 do i=0,7
    plot(i)=gsn_csm_contour_map(wks, Fr_area(i,:,:), frouderes)
 end do

gsn_panel(wks,plot,(/2,4/),panel_res)

;draw(plot)
;frame(wks)
end