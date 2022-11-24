; 4.10 读取 所有 时次文件
;    目的：相对湿度和垂直速度随时间的变化，插值取波密
; 4.15 将时次文件处理平均之后的时间剖面图
;    通过 restick属性设置XB的time label
;      参考函数：https://www.ncl.ucar.edu/Document/Functions/User_contributed/time_axis_labels.shtml
;    通过 grib_stime2time 将字符串时间转化为YYYYCCDDHHMMSS形式 用整除分段 cd_inv_cdalendar转化为Julian历
;      参考函数: https://www.ncl.ucar.edu/Document/Functions/Built-in/cd_inv_calendar.shtml
;    通过 gsn_csm_pres_hgt 绘制不同高度的剖面图
;      参考绘图集: https://www.ncl.ucar.edu/Applications/height_time.shtml 
; 代码路径  ncl /mnt/d/ncl_related/code/time_vvel_rh.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;读取数据
  fnl_dir  = "/mnt/d/data/fnl/"  
  pic_dir  ="/mnt/d/ncl_related/pictures/vapour/"
  nc_files = systemfunc("ls "+fnl_dir+"fnl_*_*_*.nc"); warning?通配符？
  f = addfiles(nc_files,"r") ; f[0-12]
  ListSetType(f,"join"); 打开多个文件并且返回文件句柄(是列表),指定如何返回文件变量
  ;printVarSummary(f); f的连接类型位join，要求数据维数一致
  ;listsettype进行了修改；默认为cat 直接连接不生成新的维
  f_time =  f[:]->LAND_P0_L1_GLL0@initial_time; 注意原始数据中有些变量如 LAND_P0_L1_GLL0 有时间属性
  ;print(f_time) ; 是每个文件的具体时间,具有单位, initial_time : 10/27/2008 (00:00)
  rh_in  = f[:]->$"RH_P0_L100_GLL0"$ 
  vvel_in= f[:]->$"VVEL_P0_L100_GLL0"$ 
  ;其中第一个维度是 join 格式（列表），下面将其替换为f_time
  ; 设置时间维度名
  rh_in!0 = "time" 
  vvel_in!0 = "time" 
;时间格式的修改
  ; rh_in&time = f_time
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
  ; rh_in的时间格式重新写为 ju格式
  rh_in&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)
  vvel_in&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)

  ; print(rh_in&time) ; [time | 12] x [lv_ISBL4 | 21] x [lat_0 | 181] x [lon_0 | 360] vvel也是
;数据处理 ： 做区域平均： lat 26-32 lon 92-98
  rh_areas = rh_in(time|:,lv_ISBL4|:,{lat_0|28:32},{lon_0|89:98})
  vvel_areas = vvel_in(time|:,lv_ISBL4|:,{lat_0|28:32},{lon_0|89:98})
  rh_area_ave = dim_avg_n(rh_areas, (/2,3/)) ; 对lat 和 lon 做区域平均
  vvel_area_ave = dim_avg_n(vvel_areas, (/2,3/)) ; 对lat 和 lon 做区域平均

  ; pvs(rh_area_ave) ; 区域平均后变量属性消失
  copy_VarCoords_not_n(rh_areas,rh_area_ave,(/2,3/)) ; 拷贝除了lat和lon以外的坐标信息
  rh_area_ave_plot = rh_area_ave(lv_ISBL4|:,time|:)
  copy_VarCoords_not_n(vvel_areas,vvel_area_ave,(/2,3/)) ; 拷贝除了lat和lon以外的坐标信息
  vvel_area_ave_plot = vvel_area_ave(lv_ISBL4|:,time|:)

wks  = gsn_open_wks ("png", pic_dir+"vvel_rh_time")          ; send graphics to PNG file

;---------------------------------------绘图属性设置------------------------------------------
  res = True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@gsnMaximize = True
  res@tmXTOn      = False ;关刻度
  res@tmYROn      = False
  ;res@gsnLeftString = "hPa" ;左上角字符
  res@tiYAxisString="P(hPa)"
  ; res@tmYLMode="Manual"
  ; res@tmYLTickSpacingF =5000
  ; res@trYReverse            = True   ; Reverse Y axis
  ; ; res@trYLog                = True  ;1.对数气压坐标,默认False,这里用不了
  ; ; res@trYAxisType          = "LogAxis" ;2.对数气压坐标,这里用不了
  ; res@gsnYAxisIrregular2Log = True    ;3.Log Convert Y axis to logarithmic;Linear, 

  res@tmXBLabelAngleF     = 45.   ;角度
  res@tmBorderThicknessF  =4    ; border 6x as thick,图框粗细
  res@tmXBMajorThicknessF =4
  res@tmXBMinorThicknessF =4
  res@tmYLMajorThicknessF =4 
  res@tmYLMinorThicknessF =4
  ; res@tmEqualizeXYSizes   =True 
  ; res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
  ; ;设置画图区域大小
  res@vpWidthF =0.6
  res@vpHeightF=0.4



  ;---------相对湿度-------------
   rhres                      =res
   rhres@cnFillOn             = True                  ; turn on color
   rhres@cnFillPalette        = "MPL_GnBu"              ; set color map
   rhres@cnLinesOn            = False                 ; no contour lines
   rhres@lbBoxLinesOn         = True ; Label bar
   ;rhres@lbOrientation       = "vertical"            ; 垂直坐标
   rhres@pmLabelBarOrthogonalPosF  = 0.08
   rhres@lbBoxEndCapStyle     = "TriangleBothEnds"
   ;rhres@pmLabelBarWidthF    = 0.05     ;色标变细
   rhres@pmLabelBarHeightF    =0.05     ;色标变窄
   ; turn off special labels on right axis ，用于gsn_csm_pres_hgt
   rhres@tmYRMode             = "Automatic"          
   rhres@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
 


   ;---------垂直速度-----------------
   vvelres           =res
   vvelres@cnFillOn  = False  ;cn means ContourPlot
   vvelres@cnLinesOn = True   
   vvelres@cnLineColor="black"
   vvelres@cnLineThicknessF =7  
   vvelres@cnLineLabelsOn = True 
   vvelres@cnLineLabelFormat="0@;*.1f" 
   vvelres@cnLineLabelBackgroundColor = -1  
   vvelres@cnLabelMasking = True
   vvelres@cnLevelSpacingF=0.1
   vvelres@cnInfoLabelOn = False ; set false will not set an information label
   vvelres@gsnContourZeroLineThicknessF = 8            ; doubles thickness of zero contour
   vvelres@gsnContourNegLineDashPattern = 1         ; sets negative contours to dash pattern 1
   vvelres@tmYRMode             = "Automatic"          ; turn off special labels on right axis 


;时间轴设置
  trestick = True
  trestick@ttmFormat = "%D/%H" ; 时间轴的格式
  trestick@ttmAxis   = "XB"
  trestick@ttmNumTicks = 12    ; 显示所有12个点

  ; pres11  =f[:]->$"lv_ISBL4"$  ;Pa,21层,12*21
  ; pres1 =pres11(0,:)
  ; pres_hpa=pres1*0.01
  ; ; pres_hpa@units="hPa"
  ; rh_area_ave_plot!0="lv"
  ; rh_area_ave_plot&lv=pres_hpa
  ; vvel_area_ave_plot!0="lv"
  ; vvel_area_ave_plot&lv=pres_hpa
  ; res@sfYArray         = pres_hpa                   ; data for y-axis tm's


;------------------------------------------绘图-----------------------------------------------------
  time_axis_labels(rh_area_ave_plot&time,rhres,trestick) ; call the formatting procedure
  time_axis_labels(vvel_area_ave_plot&time,vvelres,trestick) ; call the formatting procedure
  ; plot_rh   = gsn_csm_contour(wks,rh_area_ave_plot({100:600},:), rhres)
  ; plot_vvel = gsn_csm_contour(wks,vvel_area_ave_plot({100:600},:), vvelres)
  plot_rh = gsn_csm_pres_hgt(wks,  rh_area_ave_plot({10000:70000},:), rhres)         ; create plot    
  plot_vvel = gsn_csm_pres_hgt(wks, vvel_area_ave_plot({10000:70000},:), vvelres)     ; create plot   
  overlay(plot_rh, plot_vvel)
  draw(plot_rh)
  frame(wks)
end