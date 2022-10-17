; 是一审前的图！
; 0.1*0.1插值
;过各时刻的台风中心经向剖面，目前，环流，温度距平，高度距平
;换时次时需改动49，1处
;注意，风暴的经向移动较大，不同时刻横坐标范围需要变化，
;2600-2606,8-28N;2612-2700,10-30N;2706 15-35
;因为计算的是10个经纬度范围内的距平,不能超过把
;   ncl /mnt/h/d/ncl_related/ERA5_code/seclon_storm_topo2.ncl


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"


begin

;-1 读地形文件

  topo        = addfile("/mnt/h/d/data/other/ETOPO2v2c_f4.nc", "r")
  elev        = topo->z({15:50},{70:110}) 
  ; printVarSummary(elev)
  ; exit 
  lat_topo    = fspan(15,50,1050)
  lon_topo    = fspan(70,110,1200)
  lat_topo!0     = "lat"
  lon_topo!0     = "lon"
  lat_topo@units = "degrees_north"
  lon_topo@units = "degrees_east"
  lat_topo&lat   = lat_topo
  lon_topo&lon   = lon_topo

  ;---Attach the coordinate arrays
  elev!0    = "lat"
  elev!1    = "lon"
  elev&lat  = lat_topo
  elev&lon  = lon_topo
  ;; Convert terrain data from units "m" to "hPa", it is described as a high pressure formula
  elev        = 1013.25*(1-elev*0.0065/288.15)^5.25145
  lat1        = lat_topo  ;用于gsn_xy 地形图

;2 读变量
  pic_dir ="/mnt/h/d/ncl_related/pictures/storm/ERA5/topo2/"
  data_dir= "/mnt/h/d/data/ERA5/"   
  data    = addfile(data_dir+"ERA5.nc", "r")  
  time    =data->time
  YYYYMMDDHH_f=cd_calendar(time,-3)
count=1
  need_time=(/2008102612,2008102706/)
  plot_str=(/"(b)26/12","(d)27/06"/)

  time_locat = ind(YYYYMMDDHH_f.eq.need_time(count)) ;找该时间所处位置
  path_time  = tostring(need_time(count))
  u_short  = data ->u(time_locat,{200:1000},{10:40},{75:105}) 
  v_short  = data ->v(time_locat,{200:1000},{10:40},{75:105})
  tmp_short= data ->t(time_locat,{200:1000},{10:40},{75:105})
  u_in     =short2flt(u_short)
  v_in     =short2flt(v_short)
  tmp_in   =short2flt(tmp_short)
  vvel_in  =short2flt(data->w(time_locat,{200:1000},{10:40},{75:105}))  
  hgtData  =short2flt(data->z(time_locat,{200:1000},{10:40},{75:105}))  
  hgt_in   = hgtData/10/9.8
  copy_VarCoords(hgtData, hgt_in) ;dgpm

  lv      =data->level({200:1000})  ;hPa
  lon_in    =data->longitude({75:105})
  lat_in    =data->latitude({10:40})
  lv1     =lv*100
  lv1@units="Pa"

;3 读取风暴中心位置
  fiTY = "/mnt/h/d/data/TC/bio042008.txt"
  ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
  nrow = numAsciiRow(fiTY)
  YYYYMMDDHH = new(nrow, "string")
  lat = new(nrow, "float")
  lon = new(nrow, "float")
  cdata = asciiread(fiTY, -1, "string")
  YYYYMMDDHH = str_get_field(cdata, 3, ",")
  lat = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
  lon = stringtofloat(str_get_field(cdata, 8, ",")) *0.1
  DateChar = stringtochar(YYYYMMDDHH)
  MM = chartostring(DateChar(:,5:6))
  DD = chartostring(DateChar(:,7:8))
  HH = chartostring(DateChar(:,9:10))
  tctime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
  stdlon=0.
  stdlat=0.
  do i = 0, nrow-1
    if (i.eq.tctime_ind) then ;最方便
      stdlon=lon(i)
      stdlat=lat(i)
    end if
  end do

; lon0=floattoint(stdlon);用于绘图范围
; lat0=floattoint(stdlat)
;4 固定topo经纬度
  elev1=elev(:,{stdlon})       ;固定某一经度的地形数据    
  elev1@_FillValue= -9.96921e+36


;5 网格点插值
  lon_out=fspan(75,105,301)
  lat_out=fspan(10,40, 301)
  u      =linint2_Wrap(lon_in, lat_in, u_in, True, lon_out, lat_out, 0)
  delete(u_in)
  v      =linint2_Wrap(lon_in, lat_in, v_in, True, lon_out, lat_out, 0)
  delete(v_in)
  vvel   =linint2_Wrap(lon_in, lat_in, vvel_in, True, lon_out, lat_out, 0)
  delete(vvel_in)
  tmp    =linint2_Wrap(lon_in, lat_in, tmp_in, True, lon_out, lat_out, 0)
  delete(tmp_in)
  hgt      =linint2_Wrap(lon_in, lat_in, hgt_in, True, lon_out, lat_out, 0)
  delete(hgt_in)

  ;6 设定垂直速度矢量
  ;将Pa/s转换成m/s   
   w1     = omega_to_w(vvel, conform(vvel,lv1, 0), tmp)  ; w[*][*][*] (m/s)，要求温度K，气压Pa
   w      =w1*100.
   copy_VarCoords(w1, w)
   sec_w  = w(:,:,{stdlon})

  ;设定经向风速v
   sec_v   = v(:,:,{stdlon})

;10 计算距平
  min_lon=stdlon-10
  max_lon=stdlon+10
  min_lat=stdlat-10
  max_lat=stdlat+10
  ;温度距平
  tmp_area= tmp(:,{min_lat:max_lat},{min_lon:max_lon}) 
  tmp_avg=dim_avg_n_Wrap(tmp_area, (/1,2/)) 
  tmp_dev=tmp-conform(tmp,tmp_avg,0)
  copy_VarCoords(tmp, tmp_dev)
  sec_tmp_dev  =tmp_dev(:,:,{stdlon})

  ;高度距平
  hgt_area= hgt(:,{min_lat:max_lat},{min_lon:max_lon}) 
  hgt_avg=dim_avg_n_Wrap(hgt_area, (/1,2/)) 
  hgt_dev=hgt-conform(hgt,hgt_avg,0)
  copy_VarCoords(hgt, hgt_dev)
  sec_hgt_dev  =hgt_dev(:,:,{stdlon})


;11 打开绘图窗口
wks = gsn_open_wks("ps",pic_dir+stdlon+"E_dev_tmp_hgt_"+path_time+"_topo2") 


;12 绘图属性设置
  res = True
  res@gsnDraw       = False
  res@gsnFrame      = False
  ; res@gsnMaximize   = True
  res@tmXTOn        = False ;关刻度
  res@tmYROn        = False
  res@gsnLeftString = "" ;左上角字符
  res@gsnRightString=""
  res@tiYAxisString ="P(hPa)"
  ; res@tiYAxisFontHeightF =0.02
  res@tmYRMode      = "Automatic" ;左右Y轴均显示气压层
  res@tmBorderThicknessF  =2.5    ; border 6x as thick,图框粗细
  res@tmXBMajorThicknessF =2.5
  res@tmXBMinorThicknessF =2.5
  res@tmYLMajorThicknessF =2.5 
  res@tmYLMinorThicknessF =2.5   

  ;设置坐标轴
  res@tmXBMode       ="Explicit"

  res@tmXBValues     =(/18,22,26,30/)
  res@tmXBLabels     =(/"18~S~o~N~N","22~S~o~N~N","26~S~o~N~N","30~S~o~N~N"/)
  res@tmXBMinorValues=(/19,20,21,23,24,25,27,28,29,31,32,33/)

  ; res@tmXBValues     =(/10,15,20,25/)
  ; res@tmXBLabels     =(/"10N","15N","20N","25N"/)
  ; res@tmXBMinorValues=(/8,9,11,12,13,14,16,17,18,19,21,22,23,24,26,27,28/)

  tmpres=res
  tmpres@cnFillOn      = True  ;cn means ContourPlot
  tmpres@cnLinesOn     = False  
  tmpres@cnInfoLabelOn = False ; set false will not set an information label 
  tmpres@cnFillPalette = "BlueWhiteOrangeRed"
  tmpres@cnLevelSelectionMode ="ExplicitLevels" 
  tmpres@cnLevels      =(/-4,-3,-2,-1,0,1,2,3,4/)
  tmpres@cnFillColors  =(/25,45,65,75,110,130,178,188,200,230/)
  ; tmpres@cnLevelSelectionMode = "ManualLevels"
  ; tmpres@cnMinLevelValF  =-4
  ; tmpres@cnMaxLevelValF  =4
  ; tmpres@cnLevelSpacingF = 1
  tmpres@lbBoxLinesOn           = True ; Label bar
  tmpres@lbLabelAutoStride      = True            ; let NCL determine label spacing
  ;res@lbOrientation           = "Vertical"         ; vertical label bar
  tmpres@pmLabelBarWidthF       = 0.5     ;色标长度
  tmpres@pmLabelBarHeightF      = 0.07 
  tmpres@pmLabelBarOrthogonalPosF  = 0. ;色标离图的距离,打开有效，即便为0
  tmpres@lbBoxEndCapStyle          = "TriangleBothEnds"
  tmpres@cnFillDrawOrder           = "Predraw"  ;调整画图的顺序   


  ;垂直速度矢量叠加在contour上，属性设置合并
  tmpres@vcMapDirection        = False ;绘制垂直剖面图要设置为False
  tmpres@vcGlyphStyle          = "LineArrow"
  tmpres@vcLineArrowThicknessF = 3
  tmpres@vcLineArrowHeadMaxSizeF = 0.01
  tmpres@vcLineArrowHeadMinSizeF = 0.01

  tmpres@vcMinDistanceF        = 0.02
  tmpres@vcRefMagnitudeF       = 10. ;设置单位量级
  tmpres@vcRefLengthF          = 0.02
  tmpres@vcRefAnnoOn           = True ;anno代表矢量注释
  tmpres@vcRefAnnoFontHeightF  = 0.015
  tmpres@vcRefAnnoZone         = 1
  tmpres@vcRefAnnoSide         = "Top"
  tmpres@vcRefAnnoJust         = "TopRight"
  ;tmpres@vcRefAnnoString1 = "5 g cm~S~-1~N~hPa~S~-1~Ns~S~-1~N)"
  tmpres@vcRefAnnoString1      = "10"
  tmpres@vcRefAnnoString2On    = False
  tmpres@vcVectorDrawOrder     ="PreDraw"    

  hgtres=res
  hgtres@cnFillOn      = False  ;cn means ContourPlot
  hgtres@cnLinesOn     = True 
  ;res@cnMonoLineColor = False
  hgtres@cnLineLabelFontHeightF=0.018
  ; hgtres@cnLineLabelPlacementMode   = "Constant" ;将标签绘制为虚线的一部分；默认Randomized
  hgtres@cnLineDashSegLenF          = 0.1 ;调整每根虚线中短线密集程度
  hgtres@cnLineLabelDensityF         =0.8
  ; hgtres@cnLabelMasking             = True ;标签设置成白色这个就关掉，不然实线是断开的
  hgtres@cnLineLabelBackgroundColor = "white" ;设置标签背景透明-1
  hgtres@cnInfoLabelOn              = False ; set false will not set an information label
  hgtres@cnLineThicknessF           = 3
  hgtres@gsnContourNegLineDashPattern=15
  hgtres@cnLevelSelectionMode       = "ManualLevels" 
  hgtres@cnLevelSpacingF            = 1
  hgtres@cnLineDrawOrder            ="PreDraw"
  hgtres@cnLabelDrawOrder           ="PreDraw"

hgtres@gsnLeftStringFontHeightF=0.02
hgtres@gsnLeftString=plot_str(count)
hgtres@gsnLeftStringOrthogonalPosF=-0.1
hgtres@gsnLeftStringParallelPosF =0.01
hgtres@txPerimOn=True
hgtres@txPerimThicknessF=2.5
hgtres@txBackgroundFillColor ="white"

  etores=True
  etores@gsnMaximize      = True
  etores@gsnDraw          = False           ; don't draw the plots yet
  etores@gsnFrame         = False           ; don't advance the frame
  etores@trYMaxF          = 1000.
  etores@gsnYRefLine      = 1000.          ; create a reference line 
  etores@gsnBelowYRefLineColor = "gray60"  ; above ref line fill red

strlat=18
endlat=33
plot_base= gsn_csm_pres_hgt_vector(wks,sec_tmp_dev(:,{strlat:endlat}),\
                                    sec_v(:,{strlat:endlat}),\
                                    sec_w(:,{strlat:endlat}),tmpres)
plot_hgt   = gsn_csm_pres_hgt(wks, sec_hgt_dev(:,{strlat:endlat}), hgtres ) 
plot_topo  = gsn_csm_xy(wks,lat1,elev1,etores) 

overlay(plot_base,plot_hgt)
overlay(plot_base,plot_topo)

;某时刻风暴位置标记
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.03                     ; font smaller. default big
  txres@txFontThicknessF=3
  txres@txFontColor   ="red"
  dum = gsn_add_text(wks,plot_base,"~F37~p",stdlat,1000,txres) 
draw(plot_base)
frame(wks)
end
;2008102600,182N,881E;2008102606,187N,884E;
;2008102612,198N,888E;2008102618,210N,896E;
;2008102700,234N,904E;2008102706,250N,906E,