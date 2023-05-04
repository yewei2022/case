;组图但是！！！和一张一张出来的图不一样 plot_uvgeo.ncl也出现过这个问题
;解决方法是join或者cat数据后，还是需要循环读入数据
;但不是每次都会出错，500hPa_anomaly.ncl就正常,好吧正常性也存疑
;但这里恐怕不是那个原因 因为这里总共就一个数据文件
;使用ERA5数据 500等压面上相当位温，风矢量，整层可降水量，画组图
;  ncl /mnt/e/d/ncl_related/ERA5_code/500_eqt_wind_precip_panel.ncl

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
  map@$unique_string("dum_box")$ = dum_box
  end

begin               

;1 读取数据
  ;1.1 读取再分析资料
  pic_dir  = "/mnt/e/d/ncl_related/pictures/eqt/ERA5/"
  data_dir = "/mnt/e/d/data/ERA5/" 

  data1    = addfile(data_dir+"ERA5.nc", "r") 
  plot_time=(/2008102612,2008102700,2008102706,2008102718/)
  time     = data1->time(0:48:6)
  YYYYMMDDHH_f= cd_calendar(time,-3)
  times=dimsizes(plot_time)
  plot_time_ind = new(dimsizes(plot_time), integer)
  do i =0,times-1
    plot_time_ind(i) = str_match_ind(YYYYMMDDHH_f, plot_time(i))
  end do
  ; print(YYYYMMDDHH_f)
  ; print(plot_time_ind)
  ; exit


  u_short  = data1 ->u(plot_time_ind,{400:600},{20:40},{80:105}) 
  v_short  = data1 ->v(plot_time_ind,{400:600},{20:40},{80:105})
  tmp_short= data1 ->t(plot_time_ind,{400:600},{20:40},{80:105})
  rh_short = data1 ->r(plot_time_ind,{400:600},{20:40},{80:105})
  u     =short2flt(u_short)
  v     =short2flt(v_short)
  tmp   =short2flt(tmp_short)
  rh    =short2flt(rh_short)

  data2    = addfile(data_dir+"ERA5_single_level.nc", "r") 
  pw_short = data2->tcwv(plot_time_ind,{20:40},{80:105});units:kgm-2 Total column water vapour
  pw_in       =short2flt(pw_short)

  level  = data1->level({400:600})  ;hPa
  level1 = level*100
  level1@units="Pa"
  delete([/data1,data2,u_short,v_short,tmp_short,rh_short,pw_short/])

  ;2 pw插值成0.1分辨率
  lon_in =u&longitude
  lat_in =u&latitude
  lon_out=fspan(80,105,251)
  lat_out=fspan(20,40,201)
  lon_out@units="degrees_east"
  lat_out@units="degrees_north"   
  pw     =linint2_Wrap(lon_in, lat_in, pw_in, False, lon_out, lat_out, 0)
  delete(pw_in)


  ;1.2 读取风暴数据,用于标记风暴位置或路径
  TC = "/mnt/e/d/data/TC/bio042008.txt"
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
  lon_sta(1)=95.77
  lat_sta(2)=28.65 ;察隅
  lon_sta(2)=97.47

  ;1.4 高原轮廓线数据
  shp_line="/mnt/d/data/TPshp/DBATP/DBATP_Line.shp"  


;3计算相当位温-----------------------------------
 eqt= pot_temp_equiv(level1,tmp, rh, 1, "rh")  ;可直接使用这个函数计算,用Pa来算
 copy_VarCoords(tmp, eqt)


;4 打开绘图窗口
;-------------------------------------------------------------------------------------------
wks = gsn_open_wks("png", pic_dir+"panel_0317")

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
  pwres@cnMaxLevelValF      = 20
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
  eqtres@cnInfoLabelOn    = False ; set false will not set an information label
  eqtres@cnLineThicknessF = 2
  eqtres@cnLevelSelectionMode = "ManualLevels" 
  eqtres@cnLevelSpacingF = 2   
  eqtres@cnLineLabelInterval=2 
  eqtres@cnLineLabelFormat ="3.3f"

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
  eqt = smth9_Wrap(eqt,0.5,0.25,False)
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
  
  ; ;这里或许需要copy一下坐标
  ; copy_VarCoords(pw, pw_mask1)

  ;mask掉我定义的框之外的数据
  pw_mask = pw_mask1(:,{28.:33.},{89.:98.}) 

  do iplot=0,3
    map(iplot) = gsn_csm_map(wks, mpres)
    plot_eqt(iplot)   = gsn_csm_contour(wks,eqt(iplot,{500},:,:),eqtres)
    plot_pw(iplot)    = gsn_csm_contour(wks,pw_mask(iplot,:,:),pwres)
    plot_wind(iplot)  = gsn_csm_vector(wks, u(iplot,{500},:,:), v(iplot,{500},:,:),windRes)
    overlay(map(iplot), plot_pw(iplot))
    overlay(map(iplot), plot_wind(iplot))
    overlay(map(iplot), plot_eqt(iplot))
    mark_station(wks,map(iplot),lon_sta,lat_sta) ;标记站点位置
    plot_topoline(iplot) = gsn_add_shapefile_polylines(wks,map(iplot),shp_line,lnres)  ;高原轮廓线
    mark_area(wks,map(iplot)) ;强降水方框
  end do
  ; ; ;单独给第1/2张图标记TC位置,没必要
  do iplot=0,2
    mark_tc(wks,YYYYMMDDHH,plot_time(iplot),map(iplot),nrow,lon_tc,lat_tc)
  end do


;7 组图
  pres=True
  ; pres@gsnFrame        = False  ; don't advance frame yet，以便后面使用gsn_legend_ndc   
  pres@gsnMaximize           = True                 ; maximize in frame
  pres@gsnPaperOrientation   = "Portrait"            ;portra和maximize必须同时设置才会生效
  pres@gsnPanelLabelBar = True    ;共用一个色标
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