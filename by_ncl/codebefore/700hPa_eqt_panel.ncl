;700hPa_eqt,将高原直接灰色填充(要比mask掉高原范围内的变量更快捷)
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "/mnt/d/ncl_related/code/shapefile_utils.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"
;  ncl /mnt/d/ncl_related/code/700hPa_eqt_panel.ncl


;某时刻风暴位置标记
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
  txres@txFontHeightF = 0.022                     ; font smaller. default big
  txres@txFontThicknessF=3
  txres@txFontColor="red"
  dum_tc = gsn_add_text(wks,map,"~F37~p",stdlon,stdlat,txres) 
  map@$unique_string("dum_tc")$ = dum_tc
  end


begin
;1 读取数据
  ;1.1 读取再分析资料
  data_dir="/mnt/d/data/fnl/eqt/"   
  pic_dir="/mnt/d/ncl_related/pictures/eqt/700_eqt/"
  nc_files = systemfunc("ls "+data_dir+"fnl_*_*_*.nc");生成文件路径列表
  f = addfiles(nc_files,"r") 
  ListSetType(f,"join")
  f_time = f[:]->RH_P0_L100_GLL0@initial_time; 读取时间，注意原始变量有时间属性

  u_in   = f[:]->UGRD_P0_L100_GLL0 (:,5:25,:,:)  ;单位：m/s lv_ISBL0=26层
  v_in   = f[:]->VGRD_P0_L100_GLL0 (:,5:25,:,:) ;26
  tmp_in = f[:]->TMP_P0_L100_GLL0 (:,5:25,:,: ) ;单位：K 26层
  rh_in  = f[:]->RH_P0_L100_GLL0    ;lv_ISBL4/lat_0/lon_0，21层
  pres11 = f[:]->lv_ISBL4  ;Pa
  pres1  =pres11(0,:)
  lon_in1 = f[:]->lon_0
  lat_in1 = f[:]->lat_0
  lon_in  = lon_in1(0,:)
  lat_in  = lat_in1(0,:)

  ;时间格式的修改,便于匹配风暴位置
  time_set = grib_stime2itime(f_time) 
  ; 将grib格式时间由字符串转化为整型，输出格式为YYYYMMDDHH
  ; 时间变量赋给rh
  rh_in!0="time"
  rh_in&time =time_set
  ;[time | 3] x [lv_ISBL4 | 21] x [lat_0 | 181] x [lon_0 | 360]
  copy_VarCoords(rh_in, u_in)
  copy_VarCoords(rh_in, v_in)
  copy_VarCoords(rh_in, tmp_in)

  ; ;1.2 读取风暴数据,用于标记风暴位置或路径
  ; TC = "/mnt/d/data/TC/bio042008.txt"
  ; ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
  ; nrow = numAsciiRow(TC) 
  ; YYYYMMDDHH = new(nrow, "string")
  ; lat_tc = new(nrow, "float")
  ; lon_tc = new(nrow, "float") 
  ; tc_data = asciiread(TC, -1, "string")
  ; YYYYMMDDHH = str_get_field(tc_data, 3, ",")
  ; lat_tc = stringtofloat(str_get_field(tc_data, 7, ",")) *0.1
  ; lon_tc = stringtofloat(str_get_field(tc_data, 8, ",")) *0.1
  ; DateChar = stringtochar(YYYYMMDDHH)
  ; MM = chartostring(DateChar(:,5:6))
  ; DD = chartostring(DateChar(:,7:8))
  ; HH = chartostring(DateChar(:,9:10))

;2 网格点插值
  lon_out=fspan(55,125,701)
  lat_out=fspan(0, 60, 601)
  u      =linint2_Wrap(lon_in, lat_in(::-1), u_in(:,:,::-1,:), True, lon_out, lat_out, 0)
  delete(u_in)
  v      =linint2_Wrap(lon_in, lat_in(::-1), v_in(:,:,::-1,:), True, lon_out, lat_out, 0)
  delete(v_in)
  tmp    =linint2_Wrap(lon_in, lat_in(::-1), tmp_in(:,:,::-1,:), True, lon_out, lat_out, 0)
  delete(tmp_in)
  rh     =linint2_Wrap(lon_in, lat_in(::-1), rh_in(:,:,::-1,:), True, lon_out, lat_out, 0)
  delete(rh_in)

;3 计算相当位温-----------------------------------
 eqt= pot_temp_equiv(pres1,tmp, rh, 1, "rh")  ;可直接使用这个函数计算,用Pa来算
 copy_VarCoords(tmp, eqt)


;4 打开绘图窗口
;----------------------------------------------------------------------------
wks = gsn_open_wks("ps", pic_dir+"panel")

;5 属性设置
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@tmYROn=False
  res@tmXTOn=False
  ;mpres@pmTickMarkDisplayMode  = "Always"

  ;设置地图
  mpres = res
  mpres@mpFillOn                = False
  mpres@mpMinLatF               = 20.
  mpres@mpMaxLatF               = 30.
  mpres@mpMinLonF               = 85.
  mpres@mpMaxLonF               = 100.
  mpres@gsnMajorLatSpacing    =5. 
  mpres@gsnMinorLatSpacing    =1.
  mpres@gsnMajorLonSpacing    =5. 
  mpres@gsnMinorLonSpacing    =1.    
  mpres@mpDataSetName              = "/mnt/d/ChinaMap/Earth..4"
  mpres@mpDataBaseVersion          = "MediumRes"
  mpres@mpAreaMaskingOn            = True
  mpres@mpMaskAreaSpecifiers       = (/"China"/) ;这和上两个命令好像是没啥用
  mpres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)
  mpres@mpOutlineDrawOrder         ="PreDraw"

  windRes = res
  ;windRes@vcGlyphStyle            = "FillArrow"   ; turn on curly vectors
  windRes@vcGlyphStyle    = "LineArrow"
  windRes@vcLineArrowThicknessF = 1
  windRes@vcMinDistanceF = 0.04
  windRes@vcRefMagnitudeF = 10.
  windRes@vcRefLengthF    = 0.02
  windRes@vcRefAnnoOn = True
  windRes@vcRefAnnoFontHeightF =0.015
  windRes@vcRefAnnoZone = 1
  windRes@vcRefAnnoSide = "Top"
  windRes@vcRefAnnoJust = "TopRight"
  windRes@vcRefAnnoString1 = "10 m/s"
  windRes@vcRefAnnoString2On = False
  windRes@vcLineArrowHeadMaxSizeF = 0.01
  windRes@vcLineArrowHeadMinSizeF = 0.01


  eqtres=res
  eqtres@cnFillOn = False  ;cn means ContourPlot
  eqtres@cnLinesOn = True
  eqtres@cnInfoLabelOn = False ; set false will not set an information label 
  eqtres@cnLineLabelsOn=True
  eqtres@cnLineLabelFontHeightF=0.016
  eqtres@cnLineLabelBackgroundColor="white" ;-1是透明
  ; ;res@cnMonoLineColor = False
  ; ; eqtres@cnLabelMasking = True ;???打开就看不到label了
  ; ; eqtres@cnLineLabelPlacementMode = "Constant"
  eqtres@cnLineDashSegLenF    = 0.3
  eqtres@cnLineThicknessF = 2
  eqtres@cnLevelSelectionMode = "ManualLevels" 
  eqtres@cnLevelSpacingF = 4

  ;高原阴影
  lnres                  =True
  lnres@gsFillColor      = "gray60" 
  shp_fill = "/mnt/d/data/TPshp/DBATP/DBATP_Polygon.shp"


;绘图
map=new(4,graphic)
plot_eqt=new(4,graphic)
plot_wind=new(4,graphic)
plot_topofill=new(4,graphic)

do iplot=0,3
  map(iplot)        = gsn_csm_map(wks, mpres)
  plot_eqt(iplot)   = gsn_csm_contour(wks,eqt(iplot,{70000},:,:),eqtres)
  plot_wind(iplot)  = gsn_csm_vector(wks, u(iplot,{70000},:,:), v(iplot,{70000},:,:),windRes)
  overlay(map(iplot), plot_wind(iplot))
  overlay(map(iplot), plot_eqt(iplot))
  plot_topofill(iplot)     =gsn_add_shapefile_polygons(wks, map(iplot), shp_fill, lnres)
end do
  ; ; ;单独给第1/2张图标记TC位置
  ; path_time=tostring(time_set)
  ; do iplot=0,1
  ;   mark_tc(wks,YYYYMMDDHH,path_time(iplot),map(iplot),nrow,lon_tc,lat_tc)
  ; end do

;7 组图
  pres=True
  ; pres@gsnFrame        = False  ; don't advance frame yet，以便后面使用gsn_legend_ndc   
  pres@gsnMaximize           = True                 ; maximize in frame
  pres@gsnPaperOrientation   = "Portrait"            ;portra和maximize必须同时设置才会生效
  gsn_panel(wks, map, (/2,2/), pres)
  ; draw(mapplot) ；因为panel自动调用draw and frame        
  ; frame(wks)
end