;读取ERA5整层水汽通量。created by 
;比湿q:kg/kg 速度u:m/s 气压:hPa，g:m/s2; vibeta(q*u/g)→量级100kg/(m*s)
;参考函数：垂直积分 http://www.ncl.ucar.edu/Document/Functions/Built-in/vibeta.shtml
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
; ncl /mnt/d/ncl_related/ERA5_code/vertical_integral_vf_panel.ncl


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

;高原轮廓线
undef("mark_topoline")
  procedure mark_topoline(wks,map)
  begin
  shp="/mnt/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 4.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,map,shp,lnres) 
  map@$unique_string("shp_plot")$ = shp_plot
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


data_dir="/mnt/d/data/ERA5/"   
pic_dir ="/mnt/d/ncl_related/pictures/vapour/ERA5/"
data1   = addfile(data_dir+"ERA5_single_levels.nc", "r")
vint_qu0 =data1->$"p71.162"$
vint_qv0 =data1->$"p72.162"$

;1.找时间,读取数据
  time    =data1->time
  YYYYMMDDHH_f=cd_calendar(time,-3)
  path_time=tostring(YYYYMMDDHH_f(0:24:24))
  ; print(path_time)
  ; exit
  ; test_time=2008102600
  ; time_locat = ind(YYYYMMDDHH_f.eq.test_time) ;找该时间所处位置
  ; print(time_locat)
  ; exit
  vint_qu      = short2flt(vint_qu0(0:24:24,:,:))  ;单位：m/s,time, lat, lon
  vint_qv      = short2flt(vint_qv0(0:24:24,:,:))  ;单位：m/s,time, lat, lon
  mm=sqrt(vint_qu^2+vint_qv^2) ;水汽通量大小
  copy_VarMeta(vint_qu,mm)


  ;风暴数据,用于标记风暴位置或路径 
  fiTY = "/mnt/d/data/TC/bio042008.txt"
  ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
  nrow = numAsciiRow(fiTY)
  YYYYMMDDHH = new(nrow, "string")
  lat = new(nrow, "float")
  lon = new(nrow, "float")
  ;以字符串形式读入
  data3 = asciiread(fiTY, -1, "string")
  YYYYMMDDHH = str_get_field(data3, 3, ",")
  lat_tc = stringtofloat(str_get_field(data3, 7, ",")) *0.1
  lon_tc = stringtofloat(str_get_field(data3, 8, ",")) *0.1
  DateChar = stringtochar(YYYYMMDDHH)
  MM = chartostring(DateChar(:,5:6))
  DD = chartostring(DateChar(:,7:8))
  HH = chartostring(DateChar(:,9:10))


wks = gsn_open_wks("x11", pic_dir+"vap_panel")

;2 属性设置
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  ; res@gsnMaximize=True
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@tmYROn=False
  res@tmXTOn=False
  ;res@pmTickMarkDisplayMode   = "Always" ;不显示次刻度线
  ; res@tmBorderThicknessF  =2.5    ; border 6x as thick,图框粗细
  ; res@tmXBMajorThicknessF =2.5
  ; res@tmXBMinorThicknessF =2.5
  ; res@tmYLMajorThicknessF =2.5 
  ; res@tmYLMinorThicknessF =2.5 


  ; ;设置坐标轴标签粗细，图例在左时候设置的
  ; res@tmEqualizeXYSizes   =True 
  ; ; 上个命令，(tmXBLabelFontHeightF,(tmXBMajorLengthF,(tmXBMinorLengthF设置好后
  ; ; 其他坐标轴就不用重复设置了
  ; res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号

  ;设置地图
  mpres                        =res
  mpres@mpFillOn               = False
  mpres@mpMinLatF              = 16.
  mpres@mpMaxLatF              = 34.
  mpres@mpMinLonF              = 84.
  mpres@mpMaxLonF              = 104.  
  mpres@gsnMajorLatSpacing     =4. 
  mpres@gsnMinorLatSpacing     =1.
  mpres@gsnMajorLonSpacing     =4. 
  mpres@gsnMinorLonSpacing     =1.  
  mpres@mpDataSetName          = "/mnt/d/ChinaMap/Earth..4"
  mpres@mpDataBaseVersion      = "MediumRes"
  mpres@mpOutlineSpecifiers    = (/"China","China:Provinces"/)
  ;边界线
  mpres@mpOutlineOn                = True
  mpres@mpOutlineBoundarySets      = "National"
  mpres@mpNationalLineColor        = "black"
  mpres@mpProvincialLineColor      = "black"
  mpres@mpGeophysicalLineColor     = "black"
  mpres@mpNationalLineThicknessF   = 1
  mpres@mpProvincialLineThicknessF = 1
  mpres@mpGeophysicalLineThicknessF= 1

  ;填色-----------------------------------------
  cnres=res
  cnres@lbLabelBarOn    = False
  cnres@cnFillOn        = True                    ; turn on color
  cnres@cnLinesOn       = False                   ; turn off contour lines
  cnres@cnLineLabelsOn  = False                   ; tuen off line labels
  cnres@cnInfoLabelOn   = False 
  ;水汽通量散度填色
  ; cnres@cnFillPalette       = "BlueWhiteOrangeRed"  ;水汽通量散度填色
  ; cnres@cnLevelSelectionMode= "ManualLevels" 
  ; cnres@cnMinLevelValF      =-10
  ; cnres@cnMaxLevelValF      = 10
  ; cnres@cnLevelSpacingF     = 1
  ;水汽通量填色
  cnres@cnFillPalette         = "MPL_PuBuGn"
  cnres@cnLevelSelectionMode  ="ExplicitLevels" 
  cnres@cnLevels              =(/100,200,300,400,500,600,700,800,900,1000/)
  cnres@cnFillColors          =(/20,30,40,50,60,70,80,90,100,110,120/)



  ;--------------水汽通量矢量--------------
  vRes = res
  vRes@vcGlyphStyle            = "CurlyVector"   ; turn on curly vectors
  ;vRes@vcGlyphStyle           = "LineArrow"
  vRes@vcLineArrowThicknessF   = 1.5
  vRes@vcMinDistanceF          = 0.015
  vRes@vcLineArrowHeadMaxSizeF = 0.005
  vRes@vcLineArrowHeadMinSizeF = 0.005
  vRes@vcRefMagnitudeF      = 400.
  vRes@vcRefLengthF         = 0.02
  vRes@vcRefAnnoOn          = True
  vRes@vcRefAnnoFontHeightF = 0.015
  vRes@vcRefAnnoZone = 1.5
  vRes@vcRefAnnoSide = "Top"
  vRes@vcRefAnnoJust = "TopRight" ;默认值
  ;vRes@vcRefAnnoString1 = "10 gcm~S~-1~N~hPa~S~-1~N~s~S~-1~N~"
  vRes@vcRefAnnoString1 = "400"
  vRes@vcRefAnnoString2On = False
  ;vRes@vcRefAnnoParallelPosF=0.33
  ;vRes@vcRefAnnoOrthogonalPosF=-0.1 ;调整上下位置
  ;vRes@vcRefAnnoBackgroundColor=-1 ;将矢量注释的填充色设为透明
  ;vRes@vcRefAnnoPerimOn=False ;将矢量注释的边框去掉


;3 绘图
  map=new(2,graphic)
  plot_vc=new(2,graphic)
  plot_cn=new(2,graphic)
  do iplot=0,1
  ; iplot=1
    map(iplot)=gsn_csm_map(wks, mpres)
    plot_vc = gsn_csm_vector(wks, vint_qu(iplot,:,:) , vint_qv(iplot,:,:) , vRes)
    plot_cn = gsn_csm_contour(wks, mm(iplot,:,:), cnres)
    overlay(map(iplot),plot_cn(iplot))
    overlay(map(iplot),plot_vc(iplot))
    mark_tc(wks,YYYYMMDDHH,path_time(iplot),map(iplot),nrow,lon_tc,lat_tc)
    mark_topoline(wks,map(iplot))
    mark_area(wks,map(iplot)) ;贴上强降水方框
  end do

; 4 组图
  plres=True
  ; plres@gsnMaximize   = True
  ; plres@gsnFrame      = False  ; don't advance frame yet，以便后面使用gsn_legend_ndc   
  plres@gsnPanelLabelBar=True    ;共用一个色标
  ; ; plres@pmLabelBarOrthogonalPosF       = -0.008 ;色标条与图之间的距离
  ; ; plres@pmLabelBarParallelPosF         = 0.0 
  ; plres@pmLabelBarWidthF               = 0.25     
  plres@pmLabelBarHeightF              = 0.06
  plres@lbLabelStride         = 2         ; 两个level只标一个标签
  ; plres@lbLabelAutoStride   = True            ; let NCL determine label spacing
  plres@lbBoxLinesOn          = True ; 色标之间分界线
  plres@lbBoxEndCapStyle      = "TriangleBothEnds"
  plres@lbLabelFontHeightF    = 0.014 
  ; plres@lbTitleOn                      = True
  ; plres@lbTitleString                  = "TBB(~S~o~N~C)"
  ; plres@lbTitleJust                    = "TopRight"
  ; plres@lbTitlePosition                = "Bottom"
  ; plres@lbTitleFontHeightF             = 0.01
  ; plres@gsnPanelBottom = 0.02   ;底部增加空白
  gsn_panel(wks, map, (/1,2/), plres)  ;panel叠加gsn_csm_map的底图才不会出问题

;   ; ;绘制≥25mm的站点位置--------------
; ; gsres               = True                     ; text mods desired
; ; gsres@txFontColor   ="red"
; ; gsres@txFontHeightF = 0.03
; pmid = gsn_add_text(wks,plot,"~F35~z",lon_new,lat_new,gsres)

; draw(plot_cn)
; frame(wks)
end