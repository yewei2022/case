; 4.10 读取 所有 时次文件
;    目的：绘制time剖面图， lat 26~30 lon 90~98 区域平均
; 4.15 将时次文件处理平均之后的时间剖面图
;    通过 restick属性设置XB的time label
;      参考函数：https://www.ncl.ucar.edu/Document/Functions/User_contributed/time_axis_labels.shtml
;    通过 grib_stime2time 将字符串时间转化为YYYYCCDDHHMMSS形式 用整除分段 cd_inv_cdalendar转化为Julian历
;      参考函数: https://www.ncl.ucar.edu/Document/Functions/Built-in/cd_inv_calendar.shtml
;    通过 gsn_csm_pres_hgt 绘制不同高度的剖面图
;      参考绘图集: https://www.ncl.ucar.edu/Applications/height_time.shtml 
;load "~/offical_sample_script/download/pvs.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
;;;>>>数据读入
  ;场数据
  ; 一个示范路径：indows_path = "/mnt/D/DATAFILE/NCLDATA/bc_data/RegCM_PRE_daily_025x025.nc"
  path_rawdata = "/mnt/D/DATAFILE/NCLDATA/lx_ex_data/ye/" ; 修改路径
  nc_files = systemfunc("ls "+path_rawdata+"fnl_*_*_*.nc"); warning?通配符？
  f = addfiles(nc_files,"r") ; f[0-12]
  ListSetType(f,"join"); 打开多个文件并且返回文件句柄(是列表)

  rh= f[:]->$"RH_P0_L100_GLL0"$ 
  vvel= f[:]->$"VVEL_P0_L100_GLL0"$ 
  f_time =  f[:]->LAND_P0_L1_GLL0@initial_time; 注意原始数据中有些变量如 LAND_P0_L1_GLL0 有时间属性
  ;其中第一个维度是 join 格式（列表），下面将其替换为f_time
  ; 修改维度名
  rh!0 = "time"
  rh!1 = "lv"
  vvel!0 = "time"
  vvel!1 = "lv"  
;;;>>>时间格式的修改
  ; rh&time = f_time
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
  ; rh的时间格式重新写为 ju格式
  rh&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)

  vvel&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)
;;;
  ; print(rh&time) ; [time | 13] x [lv_ISBL0 | 26] x [lat_0 | 181] x [lon_0 | 360]
;;;>>>数据处理 ： 做区域平均： lat 26-30 lon 90-98

  rh_areas = rh(time|:,lv|:,{lat_0|26:32},{lon_0|92:97})
  vvel_areas = vvel(time|:,lv|:,{lat_0|26:32},{lon_0|92:97})
 
  ; pvs(rh_areas)

  rh_area_ave = dim_avg_n(rh_areas, (/2,3/)) ; 对lat 和 lon 做区域平均
  vvel_area_ave = dim_avg_n(vvel_areas, (/2,3/)) ; 对lat 和 lon 做区域平均

  ; printMinMax(rh_area_ave, 1)
  ; pvs(rh_area_ave) ; 区域平均后变量属性消失
  copy_VarCoords_not_n(rh_areas,rh_area_ave,(/2,3/)) ; 拷贝除了lat和lon以外的坐标信息
  rh_area_ave_plot = rh_area_ave(lv|:,time|:)

  copy_VarCoords_not_n(vvel_areas,vvel_area_ave,(/2,3/)) ; 拷贝除了lat和lon以外的坐标信息
  vvel_area_ave_plot = vvel_area_ave(lv|:,time|:)
  printVarSummary(rh_area_ave_plot&time)
  printVarSummary(vvel_area_ave_plot&time)
  print(rh_area_ave_plot&time)


;;;>>>进行作图 ： 时间剖面图(未做调整，检查用)
  wks  = gsn_open_wks ("x11", "vvel_rh")          ; send graphics to PNG file



;---------------------------------------绘图属性设置------------------------------------------
   res = True
   res@gsnDraw = False
   res@gsnFrame = False
   res@tmXTOn = False ;关刻度
   res@tmYROn = False
   res@gsnLeftString = "hPa" ;左上角字符
   res@tiYAxisString="P/hPa"
  ;---------相对湿度-------------
   rhres=res
   rhres@cnFillOn             = True                  ; turn on color
   rhres@cnFillPalette        = "MPL_GnBu"              ; set color map
   rhres@cnLinesOn            = False                 ; no contour lines
   rhres@lbBoxLinesOn = True ; Label bar
   rhres@lbOrientation        = "vertical"            ; 垂直坐标
   rhres@pmLabelBarOrthogonalPosF  = 0.08
   rhres@lbBoxEndCapStyle = "TriangleBothEnds"
   rhres@pmLabelBarWidthF = 0.05     ;色标变窄
   rhres@tmYRMode             = "Automatic"          ; turn off special labels on right axis   


   ;---------垂直速度-----------------
   vvelres=res
   vvelres@cnFillOn = False  ;cn means ContourPlot
   vvelres@cnLinesOn = True   
   vvelres@cnLineColor="black"
   vvelres@cnLineThicknessF =4  
   vvelres@cnLineLabelsOn = True 
   vvelres@cnLineLabelFormat="0@;*.1f" 
   vvelres@cnLineLabelBackgroundColor = -1  
   vvelres@cnLabelMasking = True
   vvelres@cnLevelSpacingF=0.1
   vvelres@cnInfoLabelOn = False ; set false will not set an information label
   vvelres@gsnContourZeroLineThicknessF = 5            ; doubles thickness of zero contour
   vvelres@gsnContourNegLineDashPattern = 1         ; sets negative contours to dash pattern 1
   vvelres@tmYRMode             = "Automatic"          ; turn off special labels on right axis 



  ;>>>时间轴设置
  ; res=True
  rhres@tmXBLabelFontHeightF = 0.004 ; 旋转角度
  rhres@tmXBLabelAngleF      = 45.   ;字号
  vvelres@tmXBLabelFontHeightF = 0.004 ; 旋转角度
  vvelres@tmXBLabelAngleF      = 45.   ;字号
  trestick = True
  trestick@ttmFormat = "%D/%H" ; 时间轴的格式
  trestick@ttmAxis   = "XB"
  trestick@ttmNumTicks = 13    ; 显示所有12个点

;------------------------------------------绘图-----------------------------------------------------
  time_axis_labels(rh_area_ave_plot&time,rhres,trestick) ; call the formatting procedure
  time_axis_labels(vvel_area_ave_plot&time,vvelres,trestick) ; call the formatting procedure
  plot_rh = gsn_csm_pres_hgt(wks,  rh_area_ave_plot, rhres)         ; create plot    
  plot_vvel = gsn_csm_pres_hgt(wks,  vvel_area_ave_plot, vvelres)         ; create plot   
  overlay(plot_rh, plot_vvel)
  draw(plot_rh)
  frame(wks)
end