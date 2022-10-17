; 4.10 读取 所有 时次文件
;    目的：绘制time剖面图， lat 26~30 lon 90~98 区域平均
; 4.15 将时次文件处理平均之后的时间剖面图
;    通过 restick属性设置XB的time label
;      参考函数：https://www.ncl.ucar.edu/Document/Functions/User_contributed/time_axis_labels.shtml
;    通过 grib_stime2time 将字符串时间转化为YYYYCCDDHHMMSS形式 用整除分段 cd_inv_cdalendar转化为Julian历
;      参考函数: https://www.ncl.ucar.edu/Document/Functions/Built-in/cd_inv_calendar.shtml
;    通过 gsn_csm_pres_hgt 绘制不同高度的剖面图
;      参考绘图集: https://www.ncl.ucar.edu/Applications/height_time.shtml 
load "~/offical_sample_script/download/pvs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
begin
;;;>>>数据读入
  ;场数据
  ; 一个示范路径：indows_path = "/mnt/D/DATAFILE/NCLDATA/bc_data/RegCM_PRE_daily_025x025.nc"
  path_rawdata = "/home/YEWEI/draw/" ; 修改路径
  nc_files = systemfunc("ls "+path_rawdata+"fnl_*_*_*.nc"); warning?通配符？
  f = addfiles(nc_files,"r") ; f[0-12]
  ListSetType(f,"join"); 打开多个文件并且返回文件句柄(是列表),指定如何返回文件变量
  ; printVarSummary(f); f的连接类型位join，listsettype进行了修改
  f_time =  f[:]->LAND_P0_L1_GLL0@initial_time; 注意原始数据中有些变量如 LAND_P0_L1_GLL0 有时间属性
  ; print(f_time) ; 是每个文件的具体时间,具有单位, initial_time : 10/27/2008 (00:00)
  ; test_levels = f[:]->$"HGT_P0_L100_GLL0"$ ; 下面用一个变量test_levels做示范(HGT数据)
  ; test_levels = f[:]->$"TMP_P0_L100_GLL0"$ ; (TMP数据)
  test_levels = f[:]->$"VVEL_P0_L100_GLL0"$ ; 下面用一个变量test_levels做示范(vvel数据)
  ; pvs(test_levels)
  ; pvs(test_levels) ; 其中第一个维度是 join 格式（列表），下面将其替换为f_time
  ; 修改维度名
  test_levels!0 = "time"
  test_levels!1 = "lv"
;;;>>>时间格式的修改
  ; test_levels&time = f_time
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
  ; test_levels的时间格式重新写为 ju格式
  test_levels&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)
;;;
  ; print(test_levels&time) ; [time | 13] x [lv_ISBL0 | 26] x [lat_0 | 181] x [lon_0 | 360]
;;;>>>数据处理 ： 做区域平均： lat 26-30 lon 90-98
  test_levels_areas = test_levels(time|:,lv|:,{lat_0|26:30},{lon_0|90:98})
  ; pvs(test_levels_areas)

  test_areas_ave = dim_avg_n(test_levels_areas, (/2,3/)) ; 对lat 和 lon 做区域平均
  ; printMinMax(test_areas_ave, 1)
  ; pvs(test_areas_ave) ; 区域平均后变量属性消失
  copy_VarCoords_not_n(test_levels_areas,test_areas_ave,(/2,3/)) ; 拷贝除了lat和lon以外的坐标信息
  test_areas_ave_plot = test_areas_ave(lv|:,time|:)
    ; pvs(test_areas_ave_plot) ; 
  ; print(test_areas_ave_plot&time)
   
    ; 修改高度属性单位 （不必要
      ; test_areas_ave&lv_ISBL0 = test_areas_ave&lv_ISBL0/100.0;
      ; test_areas_ave&lv_ISBL0@unit = "hPa"
      ; print(test_areas_ave&lv_ISBL0)

;;;>>>进行作图 ： 时间剖面图(未做调整，检查用)
  wks  = gsn_open_wks ("x11", "h_lat" )          ; send graphics to PNG file
  res                      = True                  ; plot mods desired
  
  res@cnFillOn             = True                  ; turn on color
  res@cnFillPalette        = "BlWhRe"              ; set color map
  res@cnLinesOn            = False                 ; no contour lines
  res@lbOrientation        = "vertical"            ; 垂直坐标
  res@pmLabelBarOrthogonalPosF  = 0.08

  ;>>>时间轴设置
  res@tmXBLabelFontHeightF = 0.004 ; 旋转角度
  res@tmXBLabelAngleF      = 45.   ;字号
  restick = True
  restick@ttmFormat = "%D/%H" ; 时间轴的格式
  restick@ttmAxis   = "XB"
  restick@ttmNumTicks = 13    ; 显示所有13个点
  time_axis_labels(test_areas_ave_plot&time,res,restick) ; call the formatting procedure
  plot = gsn_csm_pres_hgt(wks, test_areas_ave_plot, res)         ; create plot    
end