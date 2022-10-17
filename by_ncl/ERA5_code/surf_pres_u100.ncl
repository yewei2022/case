;用ERA5资料，时间改两处，读数据，画风暴位置两处
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;代码路径  ncl /mnt/h/d/ncl_related/ERA5_code/surf_pres_u100.ncl

undef("mark_tc")
  procedure mark_tc(wks,YYYYMMDDHH,need_time,map,nrow,lon_tc,lat_tc)
  begin
  tctime_ind = str_match_ind(YYYYMMDDHH, need_time) ; 返回YYYYMMDDHH中和need_time（文件时间）匹配的位置ind
  stdlon=0.
  stdlat=0.
  do i = 0, nrow-1
   if (i.eq.tctime_ind) then ;最方便
     stdlon=lon_tc(i)
     stdlat=lat_tc(i)
   end if
  end do
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.03                     ; font smaller. default big
  txres@txFontThicknessF=6
  txres@txFontColor="red"
  dum_tc = gsn_add_text(wks,map,"~F37~p",stdlon,stdlat,txres) 
  map@$unique_string("dum_tc")$ = dum_tc
  end

undef("mark_topoline")
  procedure mark_topoline(wks,map)
  begin
  shp="/mnt/h/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,map,shp,lnres) 
  map@$unique_string("shp_plot")$ = shp_plot
  end

undef("mark_area")
  procedure mark_area(wks,map)
  begin
  ypts = (/33.0,  33.0,  28.0, 28.0, 33.0/)
  xpts = (/89.0,  98.0,  98.0, 89.0, 89.0/)
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 7.0                       ; thickness of lines
  dum_box = new(4,graphic)  
  ; draw each line separately. Each line must contain two points.
  do i = 0 , 3
    dum_box(i)=gsn_add_polyline(wks,map,xpts(i:i+1),ypts(i:i+1),resp)      
  end do 
  map@$unique_string("dum_box")$ = dum_box
  end

begin

data_dir="/mnt/h/d/data/ERA5/"   
pic_dir ="/mnt/h/d/ncl_related/pictures/surface/"
  data1   = addfile(data_dir+"ERA5_single_level.nc", "r")
  time    =data1->time
  YYYYMMDDHH_f=cd_calendar(time,-3)

  day=  "26"
  hh =  "18"
  need_time="200810"+day+hh
  time_locat = ind(YYYYMMDDHH_f.eq.need_time) ;找该时间所处位置
  sp      =short2flt(data1->sp(time_locat,{10:35},{70:110}))
  sp      =sp*0.01 ;把Pa变成hPa
  sp@units="hPa"

  data2 = addfile(data_dir+"ERA5_wind.nc", "r")
  u     = short2flt(data2 ->u100(time_locat,{10:35},{70:110}))
  v     = short2flt(data2 ->v100(time_locat,{10:35},{70:110}))


;风暴数据,用于标记风暴位置或路径
fiTY = "/mnt/h/d/data/TC/bio042008.txt"
; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
nrow = numAsciiRow(fiTY)
YYYYMMDDHH = new(nrow, "string")
lat = new(nrow, "float")
lon = new(nrow, "float")
cdata = asciiread(fiTY, -1, "string")
YYYYMMDDHH = str_get_field(cdata, 3, ",")
lat_tc = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
lon_tc = stringtofloat(str_get_field(cdata, 8, ",")) *0.1
DateChar = stringtochar(YYYYMMDDHH)
MM = chartostring(DateChar(:,5:6))
DD = chartostring(DateChar(:,7:8))
HH = chartostring(DateChar(:,9:10))

wks = gsn_open_wks("eps", pic_dir+need_time)

res = True
res@gsnDraw     = False
res@gsnFrame    = False
; res@gsnMaximize =True
res@gsnLeftString  = ""
res@gsnRightString = ""
res@tmYROn         =False
res@tmXTOn         =False
;mpres@pmTickMarkDisplayMode  = "Always"
res@tmBorderThicknessF  =2.5    ; border 6x as thick,图框粗细
res@tmXBMajorThicknessF =2.5
res@tmXBMinorThicknessF =2.5
res@tmYLMajorThicknessF =2.5 
res@tmYLMinorThicknessF =2.5   
; ; ;设置坐标轴标签粗细
; res@tmEqualizeXYSizes   =True 
; ; 现在一个坐标轴设置好后，
; ;(tmXBLabelFontHeightF,(tmXBMajorLengthF,(tmXBMinorLengthF
; ; 其他坐标轴就不用重复设置了
; res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号


;设置地图
mpres = res
mpres@mpFillOn                = False
mpres@mpMinLatF               = 25.
mpres@mpMaxLatF               = 35.
mpres@mpMinLonF               = 85.
mpres@mpMaxLonF               = 100.
mpres@gsnMajorLatSpacing    =5. 
mpres@gsnMinorLatSpacing    =1
mpres@gsnMajorLonSpacing    =5. 
mpres@gsnMinorLonSpacing    =1.
mpres@mpDataSetName              = "/mnt/h/d/ChinaMap/Earth..4"
mpres@mpDataBaseVersion          = "MediumRes"
mpres@mpAreaMaskingOn            = True
mpres@mpMaskAreaSpecifiers       = (/"China"/) ;这个和上个两个命令好像是没啥用
mpres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)
; mpres@tmXBLabels=(/"75~S~o~N~E","80~S~o~N~E","85~S~o~N~E","90~S~o~N~E",\
;                   "95~S~o~N~E","100~S~o~N~E","105~S~o~N~E"/)
; mpres@tmYLLabels=(/"20~S~o~N~N","25~S~o~N~N","30~S~o~N~N","35~S~o~N~N"/)
; mpres@pmTickMarkDisplayMode ="Always" ;坐标轴标签上添加符号度,但是次刻度线也变了


;--------------矢量图--------------

windRes = res
;windRes@vcGlyphStyle            = "FillArrow"   ; turn on curly vectors
windRes@vcGlyphStyle    = "LineArrow"
windRes@vcLineArrowThicknessF = 3
windRes@vcMinDistanceF = 0.04
windRes@vcRefMagnitudeF = 2.
windRes@vcRefLengthF    = 0.02
windRes@vcRefAnnoOn = True
windRes@vcRefAnnoFontHeightF =0.015
windRes@vcRefAnnoZone = 1
windRes@vcRefAnnoSide = "Top"
windRes@vcRefAnnoJust = "TopRight"
windRes@vcRefAnnoString1 = "2 m/s"
windRes@vcRefAnnoString2On = False
windRes@vcLineArrowHeadMaxSizeF = 0.01
windRes@vcLineArrowHeadMinSizeF = 0.01
windRes@vcRefAnnoOrthogonalPosF=0
windRes@vcRefAnnoPerimSpaceF=0.08
windRes@vcRefAnnoArrowSpaceF=1

windRes@gsnLeftStringFontHeightF=0.025
windRes@gsnLeftString=day+"/"+hh
windRes@gsnLeftStringOrthogonalPosF=-0.11
windRes@gsnLeftStringParallelPosF =0.01
windRes@txPerimOn=True
windRes@txPerimThicknessF=2.5
windRes@txBackgroundFillColor ="white"


;-------------填色--------------
cnRes                    = res
  cnRes@cnFillOn           = True  ;cn means ContourPlot
  cnRes@cnLinesOn          = False   
  cnRes@cnInfoLabelOn        = False
  cnRes@cnFillPalette       ="precip3_16lev"
  cnRes@cnLevelSelectionMode = "ManualLevels" 
  cnRes@cnLevelSpacingF      = 50 
  cnRes@lbLabelFontHeightF  = 0.015
  cnRes@lbBoxLinesOn        = True ; Label bar
  ; cnRes@lbLabelStride =2 ;每两个一个标签
  ; cnRes@lbLabelAutoStride    = True               ; let NCL determine label spacing
  ; cnRes@lbOrientation        = "Vertical"         ; vertical label bar
  ; cnRes@pmLabelBarWidthF     = 0.05    
  cnRes@pmLabelBarHeightF    = 0.1
  cnRes@pmLabelBarOrthogonalPosF  = 0.1 ;色标离图的距离
  cnRes@lbBoxEndCapStyle = "TriangleBothEnds"


;--------------叠加高度场--------------
hgtRes = res
hgtRes@cnFillOn = False  ;cn means ContourPlot
hgtRes@cnLinesOn = True   
hgtRes@cnLineColor="blue"
hgtRes@cnLineLabelsOn = True  
hgtRes@cnInfoLabelOn = False ; set false will not set an information label
hgtRes@cnLineLabelFontHeightF=0.015
; hgtRes@cnLineLabelBackgroundColor=-1 ;背景变透明
; hgtRes@cnLabelMasking = True ；不穿过线

hgtRes@cnLevelSelectionMode = "ManualLevels" 
; hgtRes@cnMinLevelValF = 1100 
; hgtRes@cnMaxLevelValF = 1250
hgtRes@cnLevelSpacingF = 50
hgtRes@cnLineThicknessF = 6


;绘图
  ; sp = smth9_Wrap(sp,0.5,0.25,False)
  ; sp = smth9_Wrap(sp,0.5,0.25,False)
  map = gsn_csm_map(wks, mpres)
  windPlot = gsn_csm_vector(wks, u , v , windRes)
  contour_height = gsn_csm_contour(wks,sp, cnRes)
  overlay(map, contour_height)
  overlay(map, windPlot)
  mark_tc(wks,YYYYMMDDHH,need_time,map,nrow,lon_tc,lat_tc)
  mark_topoline(wks,map)
  mark_area(wks,map) 

  ;----------------------绘制站点位置-------------------------------
  lat_sta=(/27.98,29.87,28.65/)
  lon_sta=(/91.95,95.77,97.47/)
  ;站点标记
  mkres               = True                     ; text mods desired
  mkres@txFontColor   ="red"
  mkres@txFontHeightF = 0.03
  mkres@txFontThicknessF=3
  dum_sta = gsn_add_text(wks,map,"~F35~z",lon_sta,lat_sta,mkres) 

draw(map)
frame(wks)
end