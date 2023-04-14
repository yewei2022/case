; 500等压面上相当位温，风矢量，整层可降水量-------------------
;尝试mask掉我定义的框之外的数据
;尝试mask掉高原以外的数据
;时间改1处
;代码路径  ncl /mnt/d/ncl_related/code/500_eqt_wind_precip.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "/mnt/d/ncl_related/code/shapefile_utils.ncl"

begin

;读变量-----------------------------------------------------
  ;;;>>>手动设置一个时间 用各个数据生成一个带unit的 timeset_ymdh 并且转换为一个字符串 方便文件中读入
  year   = 2008
  month  = 10
  day    = 28
  hour   = 00
  minute = 00
  second = 0
;设置时间格式读取数据####################################
  units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
  opt = 0
  timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
  file_format  = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
  file_time    = cd_string(timeset_ymdh,file_format)
  path_format = "%Y%N%D%H"   ; 台风路径文件中的时间形式
  path_time  = cd_string(timeset_ymdh,path_format)

  data_dir="/mnt/d/data/fnl/"   
  pic_dir="/mnt/d/ncl_related/pictures/eqt/500_eqt/"
  data   = addfile(data_dir+"fnl_"+file_time+".nc", "r")                          
  u_in   = data->UGRD_P0_L100_GLL0 (5:25,:,:)  ;单位：m/s lv_ISBL0=26层
  v_in   = data->VGRD_P0_L100_GLL0 (5:25,:,:) ;26
  tmp_in = data->TMP_P0_L100_GLL0 (5:25,:,: ) ;单位：K 26层
  rh_in  = data->RH_P0_L100_GLL0    ;lv_ISBL4/lat_0/lon_0，21层
  pw_in  = data->PWAT_P0_L200_GLL0  ; (lat_0, lon_0) units:kgm-2  Precipitable water
  ; hgtData = data->HGT_P0_L100_GLL0(13,:,:)  ;13={50000}  hgt 的单位：gpm
  ; hgt_in = hgtData/10.
  ; copy_VarCoords(hgtData, hgt_in)
  pres1  = data->lv_ISBL4  ;Pa
  lon_in = data->lon_0
  lat_in = data->lat_0

;风暴数据,用于标记风暴位置或路径
 TC = "/mnt/d/data/TC/bio042008.txt"
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
 tctime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
 ; stdlon=0.
 ; stdlat=0.
 ; do i = 0, nrow-1
 ;   if (i.eq.tctime_ind) then ;最方便
 ;     stdlon=lon_tc(i)
 ;     stdlat=lat_tc(i)
 ;   end if
 ; end do

;---------------------------网格点插值--------------------
  lon_out=fspan(55,125,701)
  lat_out=fspan(0, 60, 601)
  u      =linint2_Wrap(lon_in, lat_in(::-1), u_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(u_in)
  v      =linint2_Wrap(lon_in, lat_in(::-1), v_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(v_in)
  tmp    =linint2_Wrap(lon_in, lat_in(::-1), tmp_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(tmp_in)
  rh     =linint2_Wrap(lon_in, lat_in(::-1), rh_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(rh_in)
  pw     =linint2_Wrap(lon_in, lat_in(::-1), pw_in(::-1,:), True, lon_out, lat_out, 0)
  delete(pw_in)
  ;[LAT_0 | 601] x [LON_0 | 701]    LAT_0: [ 0..60] LON_0: [55..125]
  ; hgt      =linint2_Wrap(lon_in, lat_in(::-1), hgt_in(::-1,:), True, lon_out, lat_out, 0)
  ; delete(hgt_in)


;计算相当位温-----------------------------------
 eqt= pot_temp_equiv(pres1,tmp, rh, 0, "rh")  ;可直接使用这个函数计算,用Pa来算
 copy_VarCoords(tmp, eqt)

; ;mask掉高原以外的pw
  shp_filename  ="/mnt/d/data/TPshp/DBATP/DBATP_Polygon.shp"
  pw_mask1 =new(dimsizes(pw),typeof(pw),pw@_FillValue)
  pw_mask1 =shapefile_mask_data(pw,shp_filename,True) 
;尝试mask掉我定义的框之外的数据
  min_lat  =28.            ; min of current latitude polygon
  max_lat  =32.            ; max       "
  min_lon  =89.            ; min       "    longitude polygon
  max_lon  =98.
  pw_mask = pw_mask1({min_lat:max_lat},{min_lon:max_lon}) ;[LAT_0 | 200] x [LON_0 | 219]

;-------------------------------------------------------------------------------------------
wks = gsn_open_wks("eps", pic_dir+file_time)

;--------------------属性设置--------------------------------
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  ; res@gsnMaximize=True
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@tmYROn=False
  res@tmXTOn=False
  ;res@pmTickMarkDisplayMode   = "Always" ;不显示次刻度线
  res@tmBorderThicknessF  =3   ; border 6x as thick,图框粗细
  res@tmXBMajorThicknessF =3
  res@tmXBMinorThicknessF =3
  res@tmYLMajorThicknessF =3
  res@tmYLMinorThicknessF =3  

  mpres=res
  mpres@mpFillOn                = False
  mpres@mpMinLatF               = 22.
  mpres@mpMaxLatF               = 36.
  mpres@mpMinLonF               = 85.
  mpres@mpMaxLonF               = 100.  
  mpres@gsnMajorLatSpacing    =4. 
  mpres@gsnMinorLatSpacing    =1.
  mpres@gsnMajorLonSpacing    =5. 
  mpres@gsnMinorLonSpacing    =1.  
  mpres@mpDataBaseVersion       = "Mediumres"
  mpres@mpDataSetName           = "Earth..4"
  mpres@mpOutlineSpecifiers     = (/"China:states"/)
  mpres@mpNationalLineColor            = "black"  ;我国边界
  mpres@mpProvincialLineColor          = "black"  ;省界
  mpres@mpGeophysicalLineColor         = "black" ;海陆边界
  mpres@mpNationalLineThicknessF       = 1  ;边界加粗，1不至于加粗
  mpres@mpProvincialLineThicknessF     = 1
  mpres@mpGeophysicalLineThicknessF    = 1

  pwres=res
  pwres@cnFillOn = True  ;cn means ContourPlot
  pwres@cnLinesOn = False   
  pwres@cnFillPalette ="precip3_16lev"
  ;pwres@cnLineLabelFormat="0@;*.1f" 
  ;pwres@cnLineLabelInterval=0.05 
  ; pwres@cnLevelSelectionMode = "ExplicitLevels" 
  pwres@cnLevelSelectionMode = "ManualLevels" 
  pwres@cnMinLevelValF      = 2.5
  pwres@cnMaxLevelValF      = 25
  pwres@cnLevelSpacingF     = 2.5
  pwres@cnInfoLabelOn       = False ; set false will not set an information label
  pwres@lbLabelFontHeightF  = 0.02
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
  eqtres@cnLineLabelFontHeightF=0.015
  eqtres@cnLineLabelBackgroundColor="white" ;-1是透明
  ;res@cnMonoLineColor = False
  ; eqtres@cnLabelMasking = True ;???打开就看不到label了
  ; eqtres@cnLineLabelPlacementMode = "Constant"
  eqtres@cnLineDashSegLenF    = 0.3
  eqtres@cnInfoLabelOn = False ; set false will not set an information label
  eqtres@cnLineThicknessF = 4
  eqtres@cnLevelSelectionMode = "ManualLevels" 
  eqtres@cnLevelSpacingF = 2   

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

  ; hgtRes = res
  ; hgtRes@cnFillOn  = False  ;cn means ContourPlot
  ; hgtRes@cnLinesOn = True   
  ; ;hgtRes@cnMonoLineColor = False
  ; hgtRes@cnLineColor="blue"
  ; hgtRes@cnLineLabelsOn = True  
  ; hgtRes@cnLineThicknessF=4
  ; hgtRes@cnLineLabelFontHeightF=0.015
  ; hgtRes@cnInfoLabelOn  = False ; set false will not set an information label
  ; hgtRes@cnLevelSelectionMode="ExplicitLevels"
  ; hgtRes@cnLevels       =(/576/)


;----------------绘图----------------------
map       = gsn_csm_map(wks, mpres)
plot_eqt  = gsn_csm_contour(wks,eqt({50000},:,:),eqtres)
plot_pw  = gsn_csm_contour(wks,pw_mask,pwres)
plot_wind = gsn_csm_vector(wks, u({50000},:,:) , v({50000},:,:) , windRes)
; plot_hgt  = gsn_csm_contour(wks,hgt,hgtRes)



overlay(map, plot_pw)
overlay(map, plot_wind)
overlay(map, plot_eqt)
; overlay(map, plot_hgt)


; ;某时刻风暴位置标记
; txres               = True                     ; text mods desired
; txres@txFontHeightF = 0.025                     ; font smaller. default big
; txres@txFontThicknessF=3
; txres@txFontColor="red"
; dum_tc = gsn_add_text(wks,map,"~F37~p",stdlon,stdlat,txres) 

;----------------------绘制波密站点位置-------------------------------
lat_sta=new(5,float)
lon_sta=new(5,float)
lat_sta!0="lat"
lat_sta@units="degrees_north"
lon_sta!0="lon"
lon_sta@units="degrees_east"

lat_sta(0)=27.98 ;错那
lon_sta(0)=91.95
lat_sta(1)=29.86;波密
lon_sta(1)=95.76
lat_sta(2)=28.65 ;察隅
lon_sta(2)=97.46
lat_sta(3)=27.75 ;贡山
lon_sta(3)=98.67
lat_sta(4)=28.45 ;德钦
lon_sta(4)=98.88

;站点标记
mkres               = True                     ; text mods desired
mkres@txFontColor   ="red"
mkres@txFontHeightF = 0.03

dum_sta = gsn_add_text(wks,map,"~F35~z",lon_sta(1),lat_sta(1),mkres) 

;-------高原轮廓线-----------------
shp="/mnt/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,map,shp,lnres) 

;贴上强降水方框
  ypts = (/32.0,  32.0,  28.0, 28.0, 32.0/)
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