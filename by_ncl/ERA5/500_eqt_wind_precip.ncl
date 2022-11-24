;500等压面上相当位温，风矢量，整层可降水量-------------------
;尝试mask掉我定义的框之外的数据
;尝试mask掉高原以外的数据
;时间改1处
;  ncl /mnt/h/d/ncl_related/ERA5_code/500_eqt_wind_precip.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "/mnt/h/d/ncl_related/code/shapefile_utils.ncl"

begin

;读变量-----------------------------------------------------

  data_dir= "/mnt/h/d/data/ERA5/"   
  pic_dir ="/mnt/h/d/ncl_related/pictures/eqt/ERA5/"
  data1    = addfile(data_dir+"ERA5.nc", "r") 
  time     =data1->time
  YYYYMMDDHH_f=cd_calendar(time,-3)

count=3
  need_time=(/2008102612,2008102700,2008102706,2008102718/)
  plot_str=(/"(a)26/12","(b)27/00","(c)27/06","(d)27/18"/)

  time_locat = ind(YYYYMMDDHH_f.eq.need_time(count)) ;找该时间所处位置
  path_time =tostring(need_time(count))
  u_short  = data1 ->u(time_locat,{400:600},{20:40},{80:105}) 
  v_short  = data1 ->v(time_locat,{400:600},{20:40},{80:105})
  tmp_short= data1 ->t(time_locat,{400:600},{20:40},{80:105})
  rh_short = data1 ->r(time_locat,{400:600},{20:40},{80:105})
  u     =short2flt(u_short)
  v     =short2flt(v_short)
  tmp   =short2flt(tmp_short)
  rh    =short2flt(rh_short)

  data2    = addfile(data_dir+"ERA5_single_level.nc", "r") 
  pw_short = data2->tcwv(time_locat,{20:40},{80:105});units:kgm-2 Total column water vapour
  pw_in    =short2flt(pw_short)

  pres  = data1->level({400:600})  ;hPa
  pres1 = pres*100
  pres1@units="Pa"
  delete([/data1,data2,u_short,v_short,tmp_short,rh_short,pw_short/])

  ;将pw插值成0.1分辨率,因为mask后边界有空白
  lon_in =u&longitude
  lat_in =u&latitude
  lon_out=fspan(80,105,251)
  lat_out=fspan(20,40,201)
  lon_out@units="degrees_east"
  lat_out@units="degrees_north"   
  pw     =linint2_Wrap(lon_in, lat_in, pw_in, False, lon_out, lat_out, 0)
  delete(pw_in)


;计算相当位温-----------------------------------
 eqt= pot_temp_equiv(pres1,tmp, rh, 0, "rh")  ;可直接使用这个函数计算,用Pa来算
 copy_VarCoords(tmp, eqt)

; mask掉高原以外的pw
  shp_filename  ="/mnt/h/d/data/TPshp/DBATP/DBATP_Polygon.shp"
  pw_mask1 =new(dimsizes(pw),typeof(pw),pw@_FillValue)
  pw_mask1 =shapefile_mask_data(pw,shp_filename,True) 
;尝试mask掉我定义的框之外的数据
  min_lat  =28.            ; min of current latitude polygon
  max_lat  =33.            ; max       "
  min_lon  =89.            ; min       "    longitude polygon
  max_lon  =98.
  pw_mask = pw_mask1({min_lat:max_lat},{min_lon:max_lon}) ;[LAT_0 | 200] x [LON_0 | 219]

; 读取风暴位置
  TC = "/mnt/h/d/data/TC/bio042008.txt"
  ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
  nrow = numAsciiRow(TC) 
  YYYYMMDDHH = new(nrow, "string")
  lat_tc = new(nrow, "float")
  lon_tc = new(nrow, "float") 
  tc_data = asciiread(TC, -1, "string")
  YYYYMMDDHH = str_get_field(tc_data, 3, ",")
  lat_tc = stringtofloat(str_get_field(tc_data, 7, ",")) *0.1
  lon_tc = stringtofloat(str_get_field(tc_data, 8, ",")) *0.1
  DateChar = stringtochar(YYYYMMDDHH)
  MM = chartostring(DateChar(:,5:6))
  DD = chartostring(DateChar(:,7:8))
  HH = chartostring(DateChar(:,9:10))

;-------------------------------------------------------------------------------------------
wks = gsn_open_wks("ps", pic_dir+path_time)

;--------------------属性设置--------------------------------
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize=True
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@tmYROn=False
  res@tmXTOn=False
  ;res@pmTickMarkDisplayMode   = "Always" ;不显示次刻度线
  res@tmBorderThicknessF  =2.5   ; border 6x as thick,图框粗细
  res@tmXBMajorThicknessF =2.5
  res@tmXBMinorThicknessF =2.5
  res@tmYLMajorThicknessF =2.5
  res@tmYLMinorThicknessF =2.5  

  mpres=res
  mpres@mpFillOn                = False
  mpres@mpMinLatF               = 24.
  mpres@mpMaxLatF               = 34.
  mpres@mpMinLonF               = 85.
  mpres@mpMaxLonF               = 100.  
  mpres@gsnMajorLatSpacing    =2. 
  mpres@gsnMinorLatSpacing    =1.
  mpres@gsnMajorLonSpacing    =5. 
  mpres@gsnMinorLonSpacing    =1.  
  mpres@mpDataBaseVersion       = "Mediumres"
  mpres@mpDataSetName           = "/mnt/d/ChinaMap/Earth..4"
  mpres@mpOutlineSpecifiers     = (/"China:states"/)
  ; mpres@mpNationalLineColor            = "black"  ;我国边界
  ; mpres@mpProvincialLineColor          = "black"  ;省界
  ; mpres@mpGeophysicalLineColor         = "black" ;海陆边界
  ; mpres@mpNationalLineThicknessF       = 1  ;边界加粗，1不至于加粗
  ; mpres@mpProvincialLineThicknessF     = 1
  ; mpres@mpGeophysicalLineThicknessF    = 1

  mpres@tmXBLabels=(/"85~S~o~N~E","90~S~o~N~E","95~S~o~N~E","100~S~o~N~E"/)
  mpres@tmYLLabels=(/"24~S~o~N~N","26~S~o~N~N","28~S~o~N~N","30~S~o~N~N",\
                    "32~S~o~N~N","34~S~o~N~N"/)

  pwres=res
  pwres@cnFillOn = True  ;cn means ContourPlot
  pwres@cnLinesOn = False   
  pwres@cnFillPalette ="precip3_16lev"
  ;pwres@cnLineLabelFormat="0@;*.1f" 
  ;pwres@cnLineLabelInterval=0.05 
  ; pwres@cnLevelSelectionMode = "ExplicitLevels" 
  pwres@cnLevelSelectionMode = "ManualLevels" 
  pwres@cnMinLevelValF      = 2.5
  pwres@cnMaxLevelValF      = 20
  pwres@cnLevelSpacingF     = 2.5
  pwres@cnInfoLabelOn       = False ; set false will not set an information label
  pwres@lbLabelFontHeightF  = 0.015
  pwres@lbBoxLinesOn        = True ; Label bar
  ; pwres@lbLabelStride =2 ;每两个一个标签
  ; pwres@lbLabelAutoStride    = True               ; let NCL determine label spacing
  ; pwres@lbOrientation        = "Vertical"         ; vertical label bar
  ; pwres@pmLabelBarWidthF     = 0.05    
  pwres@pmLabelBarHeightF    = 0.1
  pwres@pmLabelBarOrthogonalPosF  = 0.1 ;色标离图的距离
  pwres@lbBoxEndCapStyle = "TriangleBothEnds"

  eqtres=res
  eqtres@cnFillOn = False  ;cn means ContourPlot
  eqtres@cnLinesOn = True 
  eqtres@cnLineColor="black"
  eqtres@cnLineLabelsOn=True
  eqtres@cnLineLabelFontHeightF=0.016
  eqtres@cnLineLabelBackgroundColor="white" ;-1是透明
  ;res@cnMonoLineColor = False
  ; eqtres@cnLabelMasking = True ;???打开就看不到label了
  ; eqtres@cnLineLabelPlacementMode = "Constant"
  eqtres@cnLineDashSegLenF    = 0.3
  eqtres@cnInfoLabelOn = False ; set false will not set an information label
  eqtres@cnLineThicknessF = 3
  eqtres@cnLevelSelectionMode = "ManualLevels" 
  eqtres@cnLevelSpacingF = 2  
  eqtres@cnLineLabelInterval=2 
  eqtres@cnLineLabelFormat ="3.3f"


  eqtres@gsnLeftStringFontHeightF=0.025
  eqtres@gsnLeftString=plot_str(count)
  eqtres@gsnLeftStringOrthogonalPosF=-0.11
  eqtres@gsnLeftStringParallelPosF =0.01
  eqtres@txPerimOn=True
  eqtres@txPerimThicknessF=2.5
  eqtres@txBackgroundFillColor ="white"

  ;矢量
  windRes = res
  ;windRes@vcGlyphStyle            = "FillArrow"   ; turn on curly vectors
  windRes@vcGlyphStyle    = "LineArrow"
  windRes@vcLineArrowThicknessF = 3
  windRes@vcMinDistanceF  = 0.04
  windRes@vcRefMagnitudeF = 10.
  windRes@vcRefLengthF    = 0.02
  windRes@vcRefAnnoOn     = True
  windRes@vcRefAnnoFontHeightF =0.015
  windRes@vcRefAnnoZone = 1
  windRes@vcRefAnnoSide = "Top"
  windRes@vcRefAnnoJust = "TopRight"
  windRes@vcRefAnnoString1   = "10 m/s"
  windRes@vcRefAnnoString2On = False
  windRes@vcLineArrowHeadMaxSizeF = 0.01
  windRes@vcLineArrowHeadMinSizeF = 0.01
  ; windRes@vcRefAnnoOrthogonalPosF=0
  ; windRes@vcRefAnnoPerimSpaceF=0.5
  windRes@vcRefAnnoArrowSpaceF=1

;----------------绘图----------------------
eqt = smth9_Wrap(eqt,0.5,0.25,False)
map       = gsn_csm_map(wks, mpres)
plot_eqt  = gsn_csm_contour(wks,eqt({500},:,:),eqtres)
plot_pw  = gsn_csm_contour(wks,pw_mask,pwres)
plot_wind = gsn_csm_vector(wks, u({500},:,:) , v({500},:,:) , windRes)

overlay(map, plot_pw)
overlay(map, plot_wind)
overlay(map, plot_eqt)
; overlay(map, plot_hgt)

; ;==================某时刻风暴位置标记============================
;   tctime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
;   stdlon=0.
;   stdlat=0.
;   do i = 0, nrow-1
;    if (i.eq.tctime_ind) then ;最方便
;      stdlon=lon_tc(i)
;      stdlat=lat_tc(i)
;    end if
;   end do
  
;   txres               = True                     ; text mods desired
;   txres@txFontHeightF = 0.035                     ; font smaller. default big
;   txres@txFontThicknessF=3
;   txres@txFontColor="red"
;   dum_tc = gsn_add_text(wks,map,"~F37~p",stdlon,stdlat,txres) 

;----------------------绘制波密站点位置-------------------------------
  lat_sta=new(3,float)
  lon_sta=new(3,float)
  lat_sta!0="lat"
  lat_sta@units="degrees_north"
  lon_sta!0="lon"
  lon_sta@units="degrees_east"
  lat_sta(0)=27.98 ;错那
  lon_sta(0)=91.95
  lat_sta(1)=29.87;波密
  lon_sta(1)=95.77
  lat_sta(2)=28.65 ;察隅
  lon_sta(2)=97.47
  ;站点标记
  mkres               = True                     ; text mods desired
  mkres@txFontColor   ="red"
  mkres@txFontHeightF = 0.045
  mkres@txFontThicknessF=3
  dum_sta = gsn_add_text(wks,map,"~F35~z",lon_sta,lat_sta,mkres) 

;-------高原轮廓线-----------------
  shp="/mnt/h/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 6.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,map,shp,lnres) 

;贴上强降水方框
  ypts = (/33.0,  33.0,  28.0, 28.0, 33.0/)
  xpts = (/89.0,  98.0,  98.0, 89.0, 89.0/)
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 3.0                       ; thickness of lines
  dum_box = new(4,graphic)  
  ; draw each line separately. Each line must contain two points.
  do i = 0 , 3
    dum_box(i)=gsn_add_polyline(wks,map,xpts(i:i+1),ypts(i:i+1),resp)      
  end do 

draw(map)
frame(wks)
end