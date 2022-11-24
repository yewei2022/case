; 4.10 读取 ERA5数据
;    目的：相对湿度和垂直速度(Pa/s)随时间的变化，Wf随时间的变化,两张图自己设置位置分布
; 4.15 将时次文件处理平均之后的时间剖面图
;    通过 restick属性设置XB的time label
;    参考函数：https://www.ncl.ucar.edu/Document/Functions/User_contributed/time_axis_labels.shtml
; ncl /mnt/e/d/ncl_related/ERA5_code/time_vvel_rh_wf.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

undef("read_elev_data") ;读地形文件的函数
  function read_elev_data(topo_file)
  local nlat, nlon, topo_file, lat, lon
  begin
  ;---Read data as a straight binary file
  nlat = 2160
  nlon = 4320
  setfileoption("bin","ReadByteOrder","BigEndian")
  elev= cbinread(topo_file,(/nlat,nlon/),"short")

  ;---Create 1D coordinate arrays
  lat       = fspan(90,-90,nlat)
  lon       = fspan(0,360,nlon)
  lat!0     = "lat"
  lon!0     = "lon"
  lat@units = "degrees_north"
  lon@units = "degrees_east"
  lat&lat   = lat
  lon&lon   = lon

  ;---Attach the coordinate arrays
  elev!0    = "lat"
  elev!1    = "lon"
  elev&lat  = lat
  elev&lon  = lon
 
  return(elev)
  end


begin

;1.读取数据
  ;1.1 变量数据
  pic_dir  ="/mnt/e/d/ncl_related/pictures/vapour/ERA5/"
  data_dir = "/mnt/e/d/data/ERA5/"    
  data1    = addfile(data_dir+"ERA5.nc", "r") 
  rh_in    =short2flt(data1 ->r(0:48:6,:,{10:35},{85:105}))
  vvel_in  =short2flt(data1 ->w(0:48:6,:,{10:35},{85:105}))
  data2    = addfile(data_dir+"ERA5_wind.nc", "r") 
  u_in     =short2flt(data2 ->u100(0:48:6,{10:35},{85:105}))
  v_in     =short2flt(data2 ->v100(0:48:6,{10:35},{85:105}))

  ; ; examine time series;two methods-two var type
  ; time    =data1->time
  ; YYYYMMDDHH_f=cd_calendar(time,-3)

  ; timeset_ymdh=data1->time
  ; path_format = "%Y%N%D%H"   ; YYYYMMMDDHH
  ; path_time  = cd_string(timeset_ymdh,path_format)

  ;1.2 地形数据
  setfileoption("nc","MissingToFillValue",False)    ;设置缺省值 
  ;参考网址 http://www.ncl.ucar.edu/Document/Functions/Built-in/setfileoption.shtml
  elev_old  =read_elev_data("/mnt/e/d/data/other/ETOPO5.DAT")
  elev_float=short2flt(elev_old({10:35},{85:105}))
  ; copy_VarCoords(elev_old, elev_float)
  elev=where(elev_float.ge.0., elev_float, elev_float@_FillValue)
  copy_VarCoords(elev_float, elev)
  print(elev_float({27.98},{92.0}))
  print(elev_float({29.85},{95.77}))
  exit


  ;将水平速度插值成和地形相同分辨率
  lon_in =u_in&longitude
  lat_in =u_in&latitude
  lon_out=fspan(85,105,241)
  lat_out=fspan(10, 35, 241)
  lon_out@units="degrees_east"
  lat_out@units="degrees_north"   
  u     =linint2_Wrap(lon_in, lat_in, u_in, False, lon_out, lat_out, 0)
  delete(u_in)
  v     =linint2_Wrap(lon_in, lat_in, v_in, False, lon_out, lat_out, 0)
  delete(v_in) 

  ; 将地形高分辨率插值到水平风速低分辨率
  ; lon_out =u&longitude
  ; lat_out =u&latitude
  ; elev_lon=elev_float&lon
  ; elev_lat=elev_float&lat
  ; elev = area_hi2lores_Wrap (elev_lon,elev_lat,elev_float, False, 1, lon_out,lat_out, False)
  ; delete(elev_old)
  ; delete(elev_float)
  ; printVarSummary(u)
  ; exit


;2 数据处理 
  ;2.1  rh/vvel 做区域平均： lat 28-33 lon 89-98
  rh_areas = rh_in(time|:,level|:,{latitude|28:33},{longitude|89:98})
  vvel_areas = vvel_in(time|:,level|:,{latitude|28:33},{longitude|89:98})
  rh_area_ave = dim_avg_n(rh_areas, (/2,3/)) 
  vvel_area_ave = dim_avg_n(vvel_areas, (/2,3/)) 
  copy_VarCoords_not_n(rh_areas,rh_area_ave,(/2,3/)) 
  rh_area_ave_plot = rh_area_ave(level|:,time|:)
  copy_VarCoords_not_n(vvel_areas,vvel_area_ave,(/2,3/))
  vvel_area_ave_plot = vvel_area_ave(level|:,time|:)


  ;2.2 calculate terrain forcing velocity
  ;地形梯度
  nlat=dimsizes(elev&lat)
  mlon=dimsizes(elev&lon)
  lat    =elev&lat
  lon    =elev&lon

  dlon   = (lon(2)-lon(1))*0.0174533 ; convert to radians 0.017=3.14/180
  dHdX   = new ( dimsizes(elev),typeof(elev)); pre-allocate space
  copy_VarCoords(elev, dHdX)
  do nl=0,nlat-1                      ; loop over each latitude
  dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
  dHdX(nl:nl,:) = center_finite_diff_n (elev(nl:nl,:), dX , False,0,1) ; 1指的经度那一维
  end do

  dY   =(lat(2)-lat(1))*6378388.*0.0174533 ;标量
  ;一个纬度的长度，如果格点是1*1做差分时可以直接用这个dY，但如果不是1*1需要变一下
  dHdY = new ( dimsizes(elev),typeof(elev)); pre-allocate space
  copy_VarCoords(elev, dHdY)
  dHdY= center_finite_diff_n (elev, dY , False,0,0)

  ;算速度,会自动输出离设定位置最近的点的值
  dHdY_cuo=conform_dims(dimsizes(v(:,{27.98},{91.95})), dHdY({27.98},{91.95}), 0) ;0.03286
  dHdX_cuo= conform_dims(dimsizes(u(:,{27.98},{91.95})), dHdX({27.98},{91.95}), 0)
  wf_cuo   = v(:,{27.98},{91.95})*dHdY_cuo ;+u(:,{27.98},{91.95})*dHdX_cuo 

  dHdY_bomi= conform_dims(dimsizes(v(:,{29.87},{95.77})), dHdY({29.87},{95.77}), 0) ;0.0246
  dHdX_bomi= conform_dims(dimsizes(u(:,{29.87},{95.77})), dHdX({29.87},{95.77}), 0) 
  wf_bomi   = v(:,{29.87},{95.77})*dHdY_bomi ;+u(:,{29.87},{95.77})*dHdX_cuo 

  print(max(wf_cuo))
  print(max(wf_bomi))
  ; exit

  ;将变量整合到一个数组
  data  = new((/2,9/),float)
  data(0,:) = wf_cuo
  data(1,:) = wf_bomi
  data!1="time"
  data&time=rh_in&time


wks  = gsn_open_wks ("x11", pic_dir+"vvel_rh_wf_time_u100")          ; send graphics to PNG file

;3 绘图属性设置
  res = True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@tmXTOn      = False ;关刻度
  res@tmYROn      = False
  res@tiYAxisString="P (hPa)"
  res@tiXAxisFontHeightF =0.015
  res@tmXBMinorOn =False
  ; res@tmYLMode="Manual"
  ; res@tmYLTickSpacingF =5000
  ; res@trYReverse            = True   ; Reverse Y axis
  ; ; res@trYLog                = True  ;1.对数气压坐标,默认False,这里用不了
  ; ; res@trYAxisType          = "LogAxis" ;2.对数气压坐标,这里用不了
  ; res@gsnYAxisIrregular2Log = True    ;3.Log Convert Y axis to logarithmic;Linear, 
  res@tmXBLabelAngleF     = 45.   ;角度

  ; ;设置画图区域大小
  res@vpWidthF =0.55
  res@vpHeightF=0.3

  ;相对湿度
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
  rhres@pmLabelBarOrthogonalPosF  = 0.012 ;色标离图的距离
  rhres@cnFillDrawOrder    = "Predraw"  ;调整画图的顺序 
  rhres@lbLabelFontHeightF =0.0125
  rhres@vpXF      = 0.213092
  rhres@vpYF      = 0.888514 ;0.488514

  ;垂直速度
  vvelres           =res
  vvelres@cnFillOn  = False  ;cn means ContourPlot
  vvelres@cnLinesOn = True   
  vvelres@cnLineColor="black"
  vvelres@cnLineThicknessF =3  
  vvelres@cnLineLabelsOn = True 
  vvelres@cnLineLabelFormat="0@;*.1f" 
  vvelres@cnLineLabelBackgroundColor = -1  
  vvelres@cnLabelMasking = True
  vvelres@cnLevelSpacingF=0.1
  vvelres@cnInfoLabelOn = False ; set false will not set an information label
  vvelres@gsnContourZeroLineThicknessF = 4            ; doubles thickness of zero contour
  vvelres@gsnContourNegLineDashPattern = 1         ; sets negative contours to dash pattern 1
  vvelres@tmYRMode             = "Automatic"          ; turn off special labels on right axis 

  ;xy图属性
  wfres = True
  wfres@gsnDraw  = False
  wfres@gsnFrame = False
  wfres@tmXTOn = False ;关刻度
  wfres@tmYROn = False
  wfres@tmXBMinorOn =False
  ; wfres@tiYAxisFontHeightF=0.02
  wfres@tmXBLabelAngleF  = 45.   ;角度
  wfres@tiYAxisString="W~B~f~N~ (m/s)"
  wfres@tiXAxisFontHeightF =0.015

  wfres@xyDashPattern  = 0                         ; Make curves all solid
  wfres@xyLineThicknessF= 3
  wfres@xyMarkLineMode = "MarkLines"               ; Markers *and* lines
  wfres@xyMarker      = 16               ; 2 different markers
  wfres@xyMarkerColors = (/"blue","purple"/) ; 2 different colors
  wfres@xyLineColors =(/"blue","purple"/)
  wfres@trXMaxF=953928

  wfres@tmXMajorGrid                = True          ; implement x grid 
  wfres@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
  wfres@tmXMajorGridLineDashPattern = 2             ; select short dash lines  
  wfres@tmYMajorGrid                = True          ; implement y grid
  wfres@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
  wfres@tmYMajorGridLineDashPattern = 2             ; select short dash lines
  ; add a legend
  wfres@pmLegendDisplayMode    = "Always"              ; turn on legend
  wfres@pmLegendSide           = "Top"                 ; Change location of 
  wfres@pmLegendParallelPosF   = 0.85                   ; move units right
  wfres@pmLegendOrthogonalPosF = -0.9                  ; more neg = down
  wfres@pmLegendWidthF         = 0.07                  ; Change width and
  wfres@pmLegendHeightF        = 0.035                  ; height of legend.
  wfres@lgLabelFontHeightF     = 0.0125                   ; change font height
  wfres@lgPerimOn              = False                 ; no box around
  ; labels for the legend
  wfres@xyExplicitLegendLabels = (/"Cuona","Bomi"/)

  ; ;添加辅助线
  ; wfres@gsnYRefLine=(/0.0277,0.0753/)
  ; wfres@gsnYRefLineDashPattern =4
  ; wfres@gsnYRefLineThicknessF =4
  ; wfres@gsnYRefLineColors  =(/"purple","blue"/)

  ;设置画图区域大小和位置
  wfres@vpWidthF =0.55
  wfres@vpHeightF=0.2
  wfres@vpXF      = 0.213092  
  wfres@vpYF      = 0.480000

  ;时间轴
  trestick = True
  trestick@ttmFormat = "%D/%H" ; 时间轴的格式
  trestick@ttmAxis   = "XB"
  ; trestick@ttmNumTicks = 9    ; 显示9个点
  trestick@ttmMajorStride =1

;绘图
  time_axis_labels(rh_area_ave_plot&time,rhres,trestick) 
  time_axis_labels(vvel_area_ave_plot&time,vvelres,trestick)  
  plot_rh = gsn_csm_pres_hgt(wks,  rh_area_ave_plot({100:700},:), rhres)           
  plot_vvel = gsn_csm_pres_hgt(wks, vvel_area_ave_plot({100:700},:), vvelres) 
  overlay(plot_rh, plot_vvel)

  time_axis_labels(data&time,wfres,trestick) ; call the formatting procedure 
  plot_wf  = gsn_csm_xy (wks,data&time,data,wfres) 

  draw(plot_rh)
  draw(plot_wf)
  frame(wks)
end