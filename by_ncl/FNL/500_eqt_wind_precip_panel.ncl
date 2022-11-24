;500等压面上相当位温，风矢量，整层可降水量，画组图
;时间改1处,16
;代码路径  ncl /mnt/d/ncl_related/code/500_eqt_wind_precip_panel.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "/mnt/d/ncl_related/code/shapefile_utils.ncl"


 
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

;降水站点标记
undef("mark_station")
  procedure mark_station(wks,map,lon_sta,lat_sta)
  begin
  mkres               = True                     ; text mods desired
  mkres@txFontColor   ="red"
  mkres@txFontHeightF = 0.03
  dum_sta = gsn_add_text(wks,map,"~F35~z",lon_sta,lat_sta,mkres) 
  map@$unique_string("dum_sta")$ = dum_sta
  end

;贴上强降水方框
undef("mark_area")
  procedure mark_area(wks,map)
  begin
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
  map@$unique_string("dum_box")$ = dum_box
  end

begin               

;1 读取数据
  ;1.1 读取再分析资料
  data_dir="/mnt/d/data/fnl/eqt/"   
  pic_dir="/mnt/d/ncl_related/pictures/eqt/500_eqt/"
  nc_files = systemfunc("ls "+data_dir+"fnl_*_*_*.nc");生成文件路径列表
  f = addfiles(nc_files,"r") 
  ListSetType(f,"join")
  f_time = f[:]->RH_P0_L100_GLL0@initial_time; 读取时间，注意原始变量有时间属性

  u_in   = f[:]->UGRD_P0_L100_GLL0 (:,5:25,:,:)  ;单位：m/s lv_ISBL0=26层
  v_in   = f[:]->VGRD_P0_L100_GLL0 (:,5:25,:,:) ;26
  tmp_in = f[:]->TMP_P0_L100_GLL0 (:,5:25,:,: ) ;单位：K 26层
  rh_in  = f[:]->RH_P0_L100_GLL0    ;lv_ISBL4/lat_0/lon_0，21层
  pw_in  = f[:]->PWAT_P0_L200_GLL0  ; (lat_0, lon_0) units:kgm-2  Precipitable water
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
  copy_VarCoords_not_n(rh_in, pw_in,1)

  ;1.2 读取风暴数据,用于标记风暴位置或路径
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

  ;1.3 站点位置数据，便于绘制站点位置
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
  ;1.4 高原轮廓线数据
  shp_line="/mnt/d/data/TPshp/DBATP/DBATP_Line.shp"  


;2网格点插值
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
  pw     =linint2_Wrap(lon_in, lat_in(::-1), pw_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(pw_in)


;3计算相当位温-----------------------------------
 eqt= pot_temp_equiv(pres1,tmp, rh, 1, "rh")  ;可直接使用这个函数计算,用Pa来算
 copy_VarCoords(tmp, eqt)


;4 打开绘图窗口
;-------------------------------------------------------------------------------------------
wks = gsn_open_wks("ps", pic_dir+"panel")

;5 属性设置--------------------------------
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  ; res@gsnMaximize=True
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@tmYROn=False
  res@tmXTOn=False
  ;res@pmTickMarkDisplayMode   = "Always" ;不显示次刻度线
  ; res@tmBorderThicknessF  =3    ; 默认是2，图框粗细,
  ; res@tmXBMajorThicknessF =3
  ; res@tmXBMinorThicknessF =3
  ; res@tmYLMajorThicknessF =3 
  ; res@tmYLMinorThicknessF =3   

  mpres=res
  mpres@mpFillOn                = False
  mpres@mpMinLatF               = 24.
  mpres@mpMaxLatF               = 36.
  mpres@mpMinLonF               = 85.
  mpres@mpMaxLonF               = 100.  
  mpres@gsnMajorLatSpacing    =4. 
  mpres@gsnMinorLatSpacing    =1.
  mpres@gsnMajorLonSpacing    =5. 
  mpres@gsnMinorLonSpacing    =1.  
  mpres@mpDataSetName              = "/mnt/d/ChinaMap/Earth..4"
  mpres@mpDataBaseVersion          = "MediumRes"
  mpres@mpAreaMaskingOn            = True
  mpres@mpMaskAreaSpecifiers       = (/"China"/) ;这和上两个命令好像是没啥用
  mpres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)


  pwres=res
  pwres@cnFillOn = True  ;cn means ContourPlot
  pwres@cnLinesOn = False   
  ; pwres@cnLineColor="skyblue"
  pwres@cnFillPalette ="precip3_16lev" ;"precip3_16lev" 
  ;pwres@cnLineLabelsOn = False 
  ;pwres@cnLineLabelFormat="0@;*.1f" 
  ;pwres@cnLineLabelInterval=0.05 
  ; pwres@cnLevelSelectionMode = "ExplicitLevels" 
  ; pwvalues=fspan(2.5, 25, 10)
  ; pwcolors=ispan(-1, 239, 24)
  ; pwres@cnLevels=pwvalues
  ; pwres@cnFillColors=pwcolors
  pwres@cnLevelSelectionMode = "ManualLevels" 
  pwres@cnMinLevelValF      = 2.5
  pwres@cnMaxLevelValF      = 25
  pwres@cnLevelSpacingF     = 2.5
  pwres@cnInfoLabelOn       = False ; set false will not set an information label
  pwres@lbLabelBarOn        = False

  eqtres=res
  eqtres@cnFillOn = False  ;cn means ContourPlot
  eqtres@cnLinesOn = True 
  eqtres@cnLineColor="black"
  eqtres@cnLineLabelsOn=True
  eqtres@cnLineLabelFontHeightF=0.016
  eqtres@cnLineLabelBackgroundColor="white" ;-1是透明
  ; eqtres@cnLabelMasking = True ;???打开就看不到label了
  ; eqtres@cnLineLabelPlacementMode = "Constant"
  eqtres@cnLineDashSegLenF    = 0.3
  eqtres@cnInfoLabelOn = False ; set false will not set an information label
  eqtres@cnLineThicknessF = 2
  eqtres@cnLevelSelectionMode = "ManualLevels" 
  eqtres@cnLevelSpacingF = 2   

  windRes = res
  ;windRes@vcGlyphStyle            = "FillArrow"   ; turn on curly vectors
  windRes@vcGlyphStyle    = "LineArrow"
  windRes@vcLineArrowThicknessF = 1
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

  ;高原轮廓线
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 4.0           ; 6x thickness 
 

;6 绘图----------------------
  map=new(4,graphic)
  plot_eqt=new(4,graphic)
  plot_pw=new(4,graphic)
  plot_wind=new(4,graphic)
  plot_topoline=new(4,graphic)

  ;mask掉高原以外的pw
  polyshp  ="/mnt/d/data/TPshp/DBATP/DBATP_Polygon.shp"
  pw_mask1 = new(dimsizes(pw),typeof(pw),pw@_FillValue) ;用于后面mask高原以外
  do i=0,3
    pw_mask1(i,:,:) =shapefile_mask_data(pw(i,:,:),polyshp,True) 
  end do
  ;mask掉我定义的框之外的数据
  pw_mask = pw_mask1(:,{28.:32.},{89.:98.}) 

  do iplot=0,3
    map(iplot) = gsn_csm_map(wks, mpres)
    plot_eqt(iplot)   = gsn_csm_contour(wks,eqt(iplot,{50000},:,:),eqtres)
    plot_pw(iplot)    = gsn_csm_contour(wks,pw_mask(iplot,:,:),pwres)
    plot_wind(iplot)  = gsn_csm_vector(wks, u(iplot,{50000},:,:), v(iplot,{50000},:,:),windRes)
    overlay(map(iplot), plot_pw(iplot))
    overlay(map(iplot), plot_wind(iplot))
    overlay(map(iplot), plot_eqt(iplot))
    mark_station(wks,map(iplot),lon_sta(1),lat_sta(1)) ;标记站点位置
    plot_topoline(iplot) = gsn_add_shapefile_polylines(wks,map(iplot),shp_line,lnres)  ;高原轮廓线
    mark_area(wks,map(iplot)) ;强降水方框
  end do
  ; ; ;单独给第1/2张图标记TC位置,没必要
  ; path_time=tostring(time_set)
  ; do iplot=0,1
  ;   mark_tc(wks,YYYYMMDDHH,path_time(iplot),map(iplot),nrow,lon_tc,lat_tc)
  ; end do


;7 组图
  pres=True
  ; pres@gsnFrame        = False  ; don't advance frame yet，以便后面使用gsn_legend_ndc   
  pres@gsnMaximize           = True                 ; maximize in frame
  pres@gsnPaperOrientation   = "Portrait"            ;portra和maximize必须同时设置才会生效
  pres@gsnPanelLabelBar=True    ;共用一个色标
  ; pres@pmLabelBarOrthogonalPosF       = -0.008 ;色标条与图之间的距离
  ; pres@pmLabelBarParallelPosF         = 0.0 
  pres@pmLabelBarWidthF               = 0.5    
  pres@pmLabelBarHeightF              = 0.08
  pres@lbBoxLinesOn                   = True ; 色标之间分界线
  pres@lbBoxEndCapStyle               = "TriangleBothEnds"
  pres@lbLabelFontHeightF             = 0.014
  ; pres@lbTitleOn                      = True
  ; pres@lbTitleFontHeightF             = 0.01
  ; pres@lbTitleString                  = "TBB(~S~o~N~C)"
  ; pres@lbTitleJust                    = "TopRight"
  ; pres@lbTitlePosition                = "Bottom"
  ; pres@gsnPanelBottom = 0.02  ;底部增加空白,用于gsn__ndc
  gsn_panel(wks, map, (/2,2/), pres)
  ; draw(map) ；因为panel自动调用draw and frame        
  ; frame(wks)

end