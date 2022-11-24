; 1998年以后的GMS TBB绘图
;  ncl /mnt/d/case/code/by_ncl/FNL/GMS_TBB_1998.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

undef("mark_tc")
  procedure mark_tc(wks,YYYYMMDDHH,path_time,map,nrow,lon_tc,lat_tc)
  begin
  tctime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
  stdlon=0.
  stdlat=0.
  do i = 0, nrow-1
   if (i.eq.tctime_ind) then ;最方便
     stdlon=lon_tc(i)
     stdlat=lat_tc(i)
   end if
  end do
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.02                    ; font smaller. default big
  txres@txFontThicknessF=3
  txres@txFontColor="red"
  dum_tc = gsn_add_text(wks,map,"~F37~p",stdlon,stdlat,txres) 
  map@$unique_string("dum_tc")$ = dum_tc
  end


begin
tbb_dir = "/mnt/d/case/data/TBB/GMS-fxt/TBB_"
suffix_time="08102712"
path_time="20"+suffix_time
pic_dir="/mnt/d/case/pic/TBB/GMS/"  
setfileoption("bin","ReadByteOrder","LittleEndian")
lat = new(1201,"float")
lon = new(1701,"float")
tbb = fbindirread(tbb_dir+suffix_time+".grd",-1,(/1201,1701/),"float")

lat            =fspan(-10, 50, 1201)
lon            =fspan(75, 160, 1701)

tbb!0 ="lat"
tbb&lat = lat
tbb&lat@units ="degrees_north"
tbb!1 ="lon"
tbb&lon = lon
tbb&lon@units ="degrees_east"

;风暴路径数据读取 
    fiTY = "/mnt/d/case/data/TC/bio042008.txt"
    ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
    nrow = numAsciiRow(fiTY) 
    YYYYMMDDHH = new(nrow, "string")
    lat_tc = new(nrow, "float")
    lon_tc = new(nrow, "float")
    vmax = new(nrow, "float")
    mslp = new(nrow, "integer") 
    cdata = asciiread(fiTY, -1, "string")
    YYYYMMDDHH = str_get_field(cdata, 3, ",")
    lat_tc = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
    lon_tc = stringtofloat(str_get_field(cdata, 8, ",")) *0.1
    DateChar = stringtochar(YYYYMMDDHH)
    MM = chartostring(DateChar(:,5:6))
    DD = chartostring(DateChar(:,7:8))
    HH = chartostring(DateChar(:,9:10))


wks = gsn_open_wks("eps",pic_dir+path_time)

res = True

res@gsnDraw = False
res@gsnFrame = False
res@tmYROn=False
res@tmXTOn=False
; res@vpWidthF  = 0.6
; res@vpHeightF = 0.6

res@mpMinLatF = 15
res@mpMaxLatF = 40
res@mpMinLonF = 75
res@mpMaxLonF = 105

res@mpGridAndLimbOn = "True"
; res@mpGridMaskMode = "MaskNotOcean" ;屏蔽海洋/陆地网格
res@mpGridLineDashPattern = 15
res@mpGridSpacingF = 2.5
res@mpOutlineOn = True
; res@mpOutlineBoundarySets = "National"
res@mpFillOn                       = False
res@mpDataSetName                  = "/mnt/d/ChinaMap/Earth..4"
res@mpDataBaseVersion              = "MediumRes"
; res@mpOutlineSpecifiers            = (/"China"/)


res@cnFillOn     = True
res@cnLinesOn    =False
res@gsnAddCyclic  =False         
res@cnLineLabelsOn = False
; res@cnFillMode = "RasterFill"
; res@cnRasterSmoothingOn = True

res@cnLevelSelectionMode           = "ExplicitLevels"
res@cnLevels                       = (/-20,-32,-40,-52,-60,-70/)


res@cnFillPalette                  = "perc2_9lev"   
; res@cnLevelSelectionMode           ="ManualLevels"
res@cnFillColors                   = (/8,7,6,5,4,3,-1/) 

res@lbBoxEndCapStyle               = "TriangleBothEnds"
res@lbTitleString                  = "TBB(~S~o~N~C)"
res@lbTitleJust                    = "TopRight"
res@lbTitlePosition                = "Bottom"
res@lbTitleFontHeightF             = 0.016
res@pmLabelBarOrthogonalPosF       = 0.1 ;色标条与图之间的距离
; res@pmLabelBarWidthF               = 0.08     ;色标变窄
res@pmLabelBarHeightF              = 0.06


res@tiMainString                  =path_time 

plot = gsn_csm_contour_map(wks,tbb,res)
; mark_tc(wks,YYYYMMDDHH,path_time,plot,nrow,lon_tc,lat_tc)


;高原轮廓线
  shp_line = "/mnt/f/snow_sts_data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 4.0           ; 6x thickness          ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,plot,shp_line,lnres)


   
draw(plot)
frame(wks)
end