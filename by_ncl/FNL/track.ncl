; ?籩·??-??????
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;********************************************************

; ncl /mnt/d/case/code/by_ncl/FNL/track.ncl

begin


; ???籩·??
  fiTY = "/mnt/d/case/data/TC/bio042008.txt"
  ; fiTY = "/mnt/f/snow_sts_data/BOB/bio201305.txt" 
  pic_dir="/mnt/d/case/pic/storm/" 
  ; ????????????????Ӧ?Ļ??numAsciiCol???????????
  nrow = numAsciiRow(fiTY)
  ; print(nrow)
  ; exit
 
  YYYYMMDDHH = new(nrow, "string")
  lat = new(nrow, "float")
  lon = new(nrow, "float")
  vmax = new(nrow, "integer")
  mslp = new(nrow, "integer")
 
  data = asciiread(fiTY, -1, "string")
  YYYYMMDDHH = str_get_field(data, 3, ",")
  lat = stringtofloat(str_get_field(data, 7, ",")) *0.1
  lon = stringtofloat(str_get_field(data, 8, ",")) *0.1
  vmax = stringtoint(str_get_field(data, 9, ","))
  mslp = stringtoint(str_get_field(data, 10, ","))
 
  DateChar = stringtochar(YYYYMMDDHH)
  MM = chartostring(DateChar(:,5:6))
  DD = chartostring(DateChar(:,7:8))
  HH = chartostring(DateChar(:,9:10))
 
  ;HHi = mod(stringtoint(YYYYMMDDHH), 100)
  ;DDi = mod(stringtoint(YYYYMMDDHH)/100, 100)
  ;MMi = mod(stringtoint(YYYYMMDDHH)/10000, 100)
 
; plot
wks = gsn_open_wks("eps", pic_dir+"track_paper")
 
res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnAddCyclic=False

  ; res@gsnLeftStringFontHeightF=0.025
  ; res@gsnLeftString="(a)"
  ; res@gsnLeftStringOrthogonalPosF=-0.11
  ; res@gsnLeftStringParallelPosF =0.01
  ; res@txPerimOn=True
  ; res@txPerimThicknessF=2.5
  ; res@txBackgroundFillColor ="white"

  ;?ؿ̶?
  res@tmXTOn = False 
  res@tmYROn = False
  ; ??????߽?
  res@mpMinLatF = 15
  res@mpMaxLatF = 30
  res@mpMinLonF = 80.
  res@mpMaxLonF = 100.
  ; res@gsnMajorLonSpacing    =5. 
  ; res@gsnMinorLonSpacing    =1.    
  ; res@gsnMajorLatSpacing    =5. 
  ; res@gsnMinorLatSpacing    =1
  ; ????½???????????ɫѡȡrainbowɫ??еĵ?60????
  ; ???Ĭ?????
  ; res@mpLandFillColor = 155
  ; res@mpOceanFillColor = "white"
  ; ????澭γ???
  res@mpGridAndLimbOn = "True"
  res@mpGridMaskMode = "MaskNotOcean"
  res@mpGridLineDashPattern = 15
  res@mpGridSpacingF = 5.0
  ; ??ƹ???
  res@mpOutlineOn = True
  res@mpOutlineBoundarySets = "National"
  ; ??????
  res@mpDataBaseVersion = "MediumRes"
  res@mpDataSetName = "/mnt/d/ChinaMap/Earth..4"
    ; res@mpOutlineSpecifiers = "China:States"
  res@tmEqualizeXYSizes   =True  ;XY轴字体都一样大小
  res@tmYLLabelFontHeightF =0.02

  res@pmTickMarkDisplayMode ="Always" ;坐标轴标签上添加符号度


  ; res@tmXBLabels=(/"85~S~o~N~E","90~S~o~N~E","95~S~o~N~E","100~S~o~N~E"/)
  ; res@tmYLLabels=(/"15~S~o~N~N","20~S~o~N~N","25~S~o~N~N","30~S~o~N~N"/)

; ;---????????
;   topores                    = res
;   topores@cnFillOn           = True             ; turn off contour fill
;   topores@cnFillPalette      = "OceanLakeLandSnow"
;   topores@cnLinesOn          = False
;   topores@cnLevelSelectionMode="ManualLevels" 
;   topores@cnLevelSpacingF    = 1000
;   topores@cnMinLevelValF     =0            
;   topores@cnFillMode         = "RasterFill"     ; much faster than AreaFill
;   topores@lbBoxLinesOn       = True ; Label bar
;   topores@lbLabelAutoStride  = True            ; let NCL determine label spacing
;   topores@cnLineLabelsOn     = False            ; turn off line labels
;   topores@cnInfoLabelOn      = False            ; turn off info label
;   topores@pmLabelBarHeightF  =0.08 
;   topores@lbBoxEndCapStyle   ="TriangleBothEnds"
;   topores@pmLabelBarOrthogonalPosF       = 0.15 ;ɫ?????ľ??
;   topores@lbTitleString                  ="elevation/m"
;   topores@lbTitleFontHeightF             =0.02
;   topores@lbTitleJust                    ="BottomRight"
;   topores@lbTitlePosition                ="Bottom"

;--------------------------------??Ƶ??--------------------- 
  ; plot = gsn_csm_contour_map(wks, elev, topores) 

  plot = gsn_csm_map(wks, res)

;-------------------------------------------??????---------------------------------------
  ; ???vmax(??λ????knot??????Сʱ)???????ķ?ٴ?????????Ƚ????֣?
  ; ??????????ͬ?ɫ??
  ; 0~33.5 ??????
  ; 34~64 ????籩
  ;  >64 쫷?
 
gsn_define_colormap(wks, "wgne15")

colors=(/7,13,14/) 
 
; ????
  resLine = True
  resLine@gsLineThicknessF = 3
  dumLine = gsn_add_polyline(wks, plot, lon, lat, resLine)

; ??Ƹ???ʱ?̵ĵ?
  dumDot= new(nrow, graphic)

  resDot = True
  resDot@gsMarkerIndex = 1
  resDot@gsMarkerSizeF = 0.05
  levels  = (/34,64/)
  do i = 0, nrow-1 
    if (vmax(i).lt.levels(0)) then
      resDot@gsMarkerColor = colors(0)
    end if
    if (vmax(i).ge.levels(0).and.vmax(i).lt.levels(1)) then
      resDot@gsMarkerColor = colors(1)
    end if
    if (vmax(i).ge.levels(1)) then
      resDot@gsMarkerColor = colors(2)
    end if
    dumDot(i) = gsn_add_polymarker(wks, plot, lon(i), lat(i), resDot)
  end do
 
 
; ???00ʱ????
  dumDate = new(nrow,graphic)
  resTx = True
  resTx@txFontHeightF = 0.02
  resTx@txFontColor = "black"
  resTx@txJust = "CenterLeft"
 
  do i = 1, nrow-1
    if (HH(i) .ne. "00" ) then
      continue
    end if
    dumDate(i) = gsn_add_text(wks,plot, "2008"+MM(i)+DD(i), lon(i)+0.5, lat(i), resTx)
  end do

;**********************************??Ӹ??????****************
shp="/mnt/f/snow_sts_data/TPshp/DBATP/DBATP_Line.shp"    
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 6.0           ; 2x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,plot,shp,lnres) 


; ??????
  resLg = True 
  resLg@lgItemType = "Markers"
  resLg@lgMonoMarkerIndex = True
  resLg@lgMarkerColors = colors
  resLg@lgMarkerIndex = 1
  resLg@lgMarkerSizeF = 0.05
  resLg@lgMonoDashIndex = True
  resLg@vpWidthF = 0.2
  resLg@vpHeightF = 0.1
  resLg@lgPerimFill = 0
  resLg@lgPerimFillColor = "Background"
  resLg@lgLabelFontHeightF = 0.05
  resLg@lgBoxMinorExtentF=0.2 ;圆点和文字之间的距离

  ; resLg@lgTitleFontHeightF = 0.02
  ; resLg@lgTitleString = "Rashmi Track"
 
  ; lbid = gsn_create_legend(wks, 3, (/" < 17m/s",\
  ; " 17m/s~F34~*~F~31.5m/s","> 31.5m/s"/), resLg) ; ~F34~#~F~64 ;>= ~F34~3~F~

  lbid = gsn_create_legend(wks, 3, (/" TD ",\
  " TS ","H1 or >H1"/), resLg) ; ~F34~#~F~64 ;>= ~F34~3~F~ 
  ;图例位置
  amres = True
  amres@amParallelPosF = 0.375
  amres@amOrthogonalPosF = -0.403 
  dumLg = gsn_add_annotation(plot, lbid, amres)

; ; ?map???ʡ?...
;   letters=(/"Myanmar","Bangladesh","India","Bhutan","Nepal","China"/)
;   nlettes=dimsizes(letters)
;   dumletter = new(nlettes,graphic)
;   lat_letter=(/22,24.5,23,27.3,27.2,29/)
;   lon_letter=(/96.5,90.2,86.5,90.5,86.5,93/)
;   reslet = True
;   reslet@txFontHeightF = 0.02
;   reslet@txFontThicknessF=0.24
;   reslet@txFontColor = "blue"
;   ; reslet@txJust = "CenterLeft"
;   do i=0,nlettes-1
;     dumletter(i)= gsn_add_text(wks,plot,letters(i),\
;                                lon_letter(i), lat_letter(i), reslet)
;   end do

  draw(plot)
  frame(wks)
 
end
