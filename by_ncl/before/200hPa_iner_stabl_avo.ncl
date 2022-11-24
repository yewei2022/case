;用fnl资料,风矢是箭头，高度场，绝对涡度
;改1个时次，16
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"
;代码路径  ncl /mnt/d/ncl_related/code/200hPa_iner_stabl_avo.ncl

begin
  ;手动设置一个时间，带unit的 timeset_ymdh 
  ;转换为一个字符串 方便文件中读入
   year   = 2008
   month  = 10
   day    = 28
   hour   = 12
   minute = 00
   second = 0
  ;设置时间格式读取数据####################################
   units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
   opt = 0
   timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
   file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)
   path_format = "%Y%N%D%H"   ; 台风路径文件中的时间形式
   path_time  = cd_string(timeset_ymdh,path_format)

  ;数据读取和处理
  data_dir  ="/mnt/d/data/fnl/"
  pic_dir   ="/mnt/d/ncl_related/pictures/200/"
  data      = addfile(data_dir+"fnl_"+file_time+".nc", "r")
  u_in      = data ->UGRD_P0_L100_GLL0({20000},:,:)  ;单位：m/s
  v_in      = data->VGRD_P0_L100_GLL0({20000},:,:)
  lon_in    = data->lon_0
  lat_in    = data->lat_0
  hgt_in    = data->HGT_P0_L100_GLL0({20000},:,:)  ;13={50000}  hgt 的单位：gpm
  avo_in    = data->ABSV_P0_L100_GLL0 ({20000},:,:);26层


  ;读取风暴中心位置   
  fiTY = "/mnt/d/data/TC/bio042008.txt"
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

  ; tctime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
  ; stdlon=0.
  ; stdlat=0.
  ; do i = 0, nrow-1
  ;   if (i.eq.tctime_ind) then ;最方便
  ;     stdlon=lon(i)
  ;     stdlat=lat(i)
  ;   end if
  ; end do

;---------------------------网格点插值--------------------
  lon_out=fspan(55,125,701)
  lat_out=fspan(0, 60, 601)
  u      =linint2_Wrap(lon_in, lat_in(::-1), u_in(::-1,:), True, lon_out, lat_out, 0)
  delete(u_in)
  v      =linint2_Wrap(lon_in, lat_in(::-1), v_in(::-1,:), True, lon_out, lat_out, 0)
  delete(v_in)
  hgt_out=linint2_Wrap(lon_in, lat_in(::-1), hgt_in(::-1,:), True, lon_out, lat_out, 0)
  delete(hgt_in)
  avo_out=linint2_Wrap(lon_in, lat_in(::-1), avo_in(::-1,:), True, lon_out, lat_out, 0)
  delete(avo_in)

hgt  = hgt_out/10 ;位势高度处理
copy_VarCoords(hgt_out, hgt) 
avo =avo_out*1.e+5
copy_VarCoords(avo_out, avo) 

wks = gsn_open_wks("png", pic_dir+file_time+"_avo")
; wks = gsn_open_wks("png", pic_dir+"U_"+file_time)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnLeftString = ""
res@gsnRightString = ""
;res@gsnMaximize=True
res@tmYROn=False
res@tmXTOn=False
;res@pmTickMarkDisplayMode  = "Always" 
;去掉左上角刻度，并且把所有的次刻度都关掉了，过于鸡肋

; ;设置坐标轴标签粗细
; res@tmEqualizeXYSizes =True
; res@tmXBLabelFont=22 ;21是不加粗的helvetica字体
; res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
; res@tmXBMajorLengthF = 0.005

;设置地图
mpres = res
mpres@mpFillOn                = False
mpres@mpMinLatF               = 10.
mpres@mpMaxLatF               = 50.
mpres@mpMinLonF               = 60.
mpres@mpMaxLonF               = 120.
mpres@gsnMajorLonSpacing      =10. 
mpres@gsnMinorLonSpacing      =2.    
mpres@gsnMajorLatSpacing      =10. 
mpres@gsnMinorLatSpacing      =2
mpres@mpDataBaseVersion       = "MediumRes"
mpres@mpDataSetName           = "Earth..4"
mpres@mpOutlineSpecifiers     = (/"China:states"/)


;--------------矢量图--------------

windRes = res
;windRes@vcGlyphStyle            = "FillArrow"   ; turn on curly vectors
windRes@vcGlyphStyle    = "LineArrow"
windRes@vcLineArrowThicknessF = 3
windRes@vcMinDistanceF = 0.04
windRes@vcRefMagnitudeF = 15.
windRes@vcRefLengthF    = 0.02
windRes@vcRefAnnoOn = True
windRes@vcRefAnnoFontHeightF =0.015
windRes@vcRefAnnoZone = 1
windRes@vcRefAnnoSide = "Top"
windRes@vcRefAnnoJust = "TopRight"
windRes@vcRefAnnoString1 = "15 m/s"
windRes@vcRefAnnoString2On = False
windRes@vcLineArrowHeadMaxSizeF = 0.01
windRes@vcLineArrowHeadMinSizeF = 0.01



;--------------涡度--------------

cRes = res
cRes@cnFillOn = True ;cn means ContourPlot
cRes@cnLinesOn = False   
cRes@cnSpanFillPalette = True
cRes@cnFillPalette ="BlueYellowRed"
cRes@cnLineLabelsOn = False  
cRes@cnInfoLabelOn = False ; set false will not set an information label 
; cRes@cnLevelSelectionMode = "ManualLevels" 
; cRes@cnMinLevelValF       = -6             ; set min contour level
; cRes@cnMaxLevelValF       = 6               ; set max contour level
; cRes@cnLevelSpacingF      = 2              ; set contour spacing
cRes@cnLevelSelectionMode   = "ExplicitLevels" 
cRes@cnLevels      =(/-6,-5,-4,-3,-2,-1,0/)
cRes@cnFillColors  =(/5,15,35,45,65,85,115,-1/)

cRes@lbBoxLinesOn = True ; Label bar
cRes@lbLabelAutoStride    = True            ; let NCL determine label spacing
cRes@lbLabelFontHeightF   = 0.015
cRes@lbOrientation        = "Vertical"         ; vertical label bar
cRes@pmLabelBarWidthF = 0.03     
;cRes@pmLabelBarHeightF = 0.05    
;cRes@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
cRes@lbBoxEndCapStyle = "TriangleBothEnds"


;--------------叠加高度场--------------

hgtRes = res
hgtRes@cnFillOn = False  ;cn means ContourPlot
hgtRes@cnLinesOn = True   
;hgtRes@cnMonoLineColor = False
hgtRes@cnLineColor="blue"
;hgtRes@cnSpanLinePalette = True
;hgtRes@cnLinePalette = "rainbow"
hgtRes@cnLineLabelsOn = True  
hgtRes@cnInfoLabelOn = False ; set false will not set an information label

; hgtRes@cnLineLabelBackgroundColor=-1 ;背景变透明
; hgtRes@cnLabelMasking = True ；不穿过线

hgtRes@cnLevelSelectionMode = "ManualLevels" 
; hgtRes@cnMinLevelValF = 1100 
hgtRes@cnMaxLevelValF = 1250
hgtRes@cnLevelSpacingF = 5 
hgtRes@cnLineThicknessF = 6
;hgtRes@gsnSpreadColorStart = 1
;hgtRes@gsnSpreadColorStart =17
;hgtRes@gsnSpreadColors = True ; 6.1.0
; hgtRes@lbBoxLinesOn = False ; Label bar
;hgtRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
;hgtRes@lbOrientation        = "Vertical"         ; vertical label bar



;绘图
  ;1 高度，散度，风矢量
  hgt = smth9_Wrap(hgt,0.5,0.25,False)
  ; hgt = smth9_Wrap(hgt,0.5,0.25,False)
  mapplot = gsn_csm_map(wks, mpres)
  contour_div = gsn_csm_contour(wks,avo,cRes)
  windPlot = gsn_csm_vector(wks, u , v , windRes)
  contour_height = gsn_csm_contour(wks,hgt, hgtRes)
  ; contour_U =gsn_csm_contour(wks,U, uRes)
  ;叠加风矢量和位势高度
  overlay(mapplot, contour_div)
  overlay(mapplot, windPlot)
  overlay(mapplot, contour_height)
  ; overlay(mapplot, contour_U)

  ; ;2 某时刻风暴位置标记---------------------------------
  ; txres               = True                     ; text mods desired
  ; txres@txFontHeightF = 0.03                     ; font smaller. default big
  ; txres@txFontThicknessF=6
  ; txres@txFontColor="red"
  ; dum = gsn_add_text(wks,mapplot,"~F37~p",stdlon,stdlat,txres) 

  ;3 高原轮廓线
  shp="/mnt/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres) 

draw(mapplot)
frame(wks)
end