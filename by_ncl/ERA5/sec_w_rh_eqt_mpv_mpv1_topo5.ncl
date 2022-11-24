;波密95.8E的垂直速度矢量，相对湿度图案填充，相当位温等值线，mpv、mpv1填色
;created by yewei ，2020
;地形剖面填充灰色，法二：用叠图的方式，推荐
; 参考;http://www.ncl.ucar.edu/Applications/Scripts/mask_14.ncl
;湿层填充图案参考 http://www.ncl.ucar.edu/Document/Graphics/Interfaces/gsn_contour_shade.shtml
;气压层单位是Pa和hPa的时候画的图不一致，Pa的时候明显不对，上升运动区的矢量箭头向下
;！！！遇到的问题： y轴标签始终不能设置成等距的等压坐标；800hPa和1000hPa的距离始终比别的200hPa间隔大
;解决上述问题参考 https://www.ncl.ucar.edu/Applications/vertterp.shtml vert_2.ncl
;  ncl /mnt/e/d/ncl_related/ERA5_code/sec_w_rh_eqt_mpv_mpv1_topo5.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;绘图属性设置------------------------------------------
undef("set_all_res")
  procedure set_all_res(res,etores,rhres,optres,eqtres,mpvres,mpv1res,wres)
  local res,etores,rhres,optres,eqtres,mpvres,mpv1res,wres
  begin
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize   =True
  res@tmXTOn = False ;关刻度
  res@tmYROn = False
  res@gsnLeftString = "" ;左上角字符
  res@gsnRightString=""
  res@tiYAxisString ="P(hPa)"
  ; res@tmYRMode     = "Automatic"          ; turn off special labels on right axis 
  ; res@tmBorderThicknessF       = 6.0    ; border 6x as thick,图框粗细
  res@trXMinF              = 25    ; Set min/max of X axis
  res@trXMaxF              = 35
  res@trYReverse           = True   ; Reverse Y axis
  ; res@trYLog               = True  ;1.对数气压坐标,默认False,这里用不了
  ; res@trYAxisType          = "LogAxis"        ;2.对数气压坐标,这里用不了
  res@gsnYAxisIrregular2Linear = True    ;3.Log Convert Y axis to logarithmic;Linear 
  ; res@tmYLMode             = "Manual"
  ; res@tmYLTickSpacingF     = 50
  ; res@tmYLMinorPerMajor    = 1 ;次要刻度线
  ; res@tmYLMode              = "Explicit"
  ; res@tmYLValues            = (/200,400,600,800,900,1000/)
  ; res@tmYLLabels            = "" + res@tmYLValues
  ; res@tmYLStyle               ="Irregular" ;无效
  ; res@tmYLDataBottomF         =150
  ; res@tmYLDataTopF            =1000
  res@tmYLMinorPerMajor=1
  res@tmBorderThicknessF  =2.5    ; border 6x as thick,图框粗细
  res@tmXBMajorThicknessF =2.5
  res@tmXBMinorThicknessF =2.5
  res@tmYLMajorThicknessF =2.5 
  res@tmYLMinorThicknessF =2.5   
  res@gsnPaperOrientation   = "Portrait"           
   ;portra和maximize必须同时设置才会生效

  ;地形填充属性
  etores=True
  etores@gsnMaximize      = True
  etores@trYMaxF          = 1000
  etores@gsnYRefLine      = 1000.          ; create a reference line 
  etores@gsnBelowYRefLineColor = "gray60"              ; above ref line fill red
  etores@gsnDraw        = False             ; don't draw the plots yet
  etores@gsnFrame       = False        ; don't advance the frame

  ;相当位温
  eqtres=res
  eqtres@cnFillOn             = False  ;cn means ContourPlot
  eqtres@cnLinesOn            = True 
  eqtres@cnLineColor          ="black"
  ; eqtres@cnLabelMasking       = True
  eqtres@cnLineLabelFontHeightF   =0.018
  eqtres@cnLineLabelPlacementMode = "Constant" ;使标签是线的一部分
  eqtres@cnLineDashSegLenF    = 0.3
  eqtres@cnLineLabelBackgroundColor = "white" ;设置标签背景透明-1
  eqtres@cnInfoLabelOn        = False ; set false will not set an information label
  eqtres@cnLineThicknessF     = 4
  eqtres@cnLevelSelectionMode = "ManualLevels" 
  eqtres@cnLevelSpacingF      = 2
  eqtres@cnLineDrawOrder   ="PreDraw"

  ;相对湿度
  rhres=res
  rhres@cnFillOn          = False                  ; turn on color
  rhres@cnLinesOn         = False
  rhres@cnLineLabelsOn    = False
  ; rhres@cnFillPattern     = 17
  ; rhres@cnFillDotSizeF    = 0.003            ; increase dot size of contours
  ; rhres@cnFillScaleF      = 0.3 ; densities
  rhres@cnInfoLabelOn     = False
  rhres@cnLevelSelectionMode = "ExplicitLevels" 
  rhres@cnLevels          =(/90/)
  rhres@cnFillColor       ="royalblue"
  rhres@cnLineColor       ="royalblue"
  rhres@cnLineThicknessF  = 4
  ; rhres@cnLineLabelsOn    = False
  rhres@cnLineDrawOrder   ="PreDraw"
  rhres@cnFillDrawOrder   ="PreDraw"

  optres                      = True
  optres@gsnShadeFillType     = "pattern"
  optres@gsnShadeHigh         = 17          ; stipple,gsnShadelow代表填充大小值之间的范围
  optres@gsnShadeFillScaleF   = 2         ; add extra density
  optres@gsnShadeFillDotSizeF = 0.006      ; make dots larger
  ;如果要改变图案颜色在cnfillon=false，cnfillcolor可以设置

  ;MPV
  mpvres=res
  mpvres@cnFillOn   = True  ;cn means ContourPlot
  mpvres@cnLinesOn  = False   
  mpvres@cnFillPalette = "MPL_PiYG"
  mpvres@cnInfoLabelOn = False ; set false will not set an information label

  mpvres@cnLevelSelectionMode = "ExplicitLevels" 
  mpvres@cnLevels     =(/-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0/)
  ; mpvres@cnLevels     =(/-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0/)
  mpvres@cnFillColors =(/120,115,110,105,100,95,90,85,80,-1/)

  mpvres@lbBoxLinesOn = True ; Label bar
  ; mpvres@lbLabelAutoStride  = True            ; let NCL determine label spacing
  mpvres@lbLabelStride =2
  ; mpvres@pmLabelBarWidthF     = 0.05     ;色标变窄
  mpvres@pmLabelBarHeightF  =0.07 
  mpvres@lbBoxEndCapStyle   ="TriangleBothEnds"
  mpvres@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
  mpvres@cnFillDrawOrder    = "Predraw"  ;调整画图的顺序 

  ; MPV1
  mpv1res=res
  mpv1res@cnFillOn  = False  ;cn means ContourPlot
  mpv1res@cnLinesOn = True   
  mpv1res@cnLineColor         ="red"
  mpv1res@cnLineLabelFontColor="red" 
  mpv1res@cnLineLabelFormat   ="f" ;"0@;*.1f" 为0.0
  mpv1res@cnLineLabelFontHeightF  =0.018
  ; mpv1res@cnLevelSelectionMode  = "ManualLevels"   
  ; ;mpv1res@cnLevelSpacingF      =0.5
  ; ;mpv1res@cnLineDashPattern    =14
  ; mpv1res@gsnContourZeroLineThicknessF=8
  mpv1res@cnLineThicknessF     = 4
  mpv1res@cnLevelSelectionMode = "ExplicitLevels" 
  mpv1res@cnLevels             =(/0/)
  ; mpv1res@gsnContourPosLineDashPattern =14 ; sets positive contours to dash pattern 14
  mpv1res@cnLabelMasking       = True ;线条不穿过标签
  mpv1res@cnLineLabelBackgroundColor   = -1 ;contour label背景透明
  mpv1res@cnLineLabelPlacementMode     = "Constant"
  mpv1res@cnLineDashSegLenF    = 0.2 
  mpv1res@cnInfoLabelOn        = False ; set false will not set an information label
  mpv1res@cnLineDrawOrder      = "Predraw"  ;调整画图的顺序

  ;风矢量
  wres=res
  wres@vcMapDirection        = False
  wres@vcGlyphStyle          = "LineArrow"
  wres@vcLineArrowThicknessF = 3
  wres@vcMinDistanceF        = 0.025
  wres@vcRefMagnitudeF       = 15 ;设置单位量级
  wres@vcRefLengthF          = 0.02
  wres@vcLineArrowHeadMaxSizeF = 0.01
  wres@vcLineArrowHeadMinSizeF = 0.01
  wres@vcRefAnnoOn           = True ;anno代表矢量注释
  wres@vcRefAnnoFontHeightF  = 0.018
  wres@vcRefAnnoZone         = 1
  wres@vcRefAnnoSide         = "Top"
  wres@vcRefAnnoJust         = "TopRight"
  ;res@vcRefAnnoString1 = "5 g cm~S~-1~N~hPa~S~-1~Ns~S~-1~N)"
  wres@vcRefAnnoString1      = "15"
  wres@vcRefAnnoString2On    = False
  wres@vcVectorDrawOrder     ="PreDraw"

  end

;读地形文件----------------------------------------------
undef("read_height_data")
  function read_height_data(topo_file)
  local nlat,nlon,topo_file,lat,lon
  begin
  nlat = 2160
  nlon = 4320
  setfileoption("bin","ReadByteOrder","BigEndian")
  elev = tofloat(cbinread(topo_file,(/nlat,nlon/),"short"))
  
  lat = fspan(90,-90,nlat)
  lon = fspan(0,360,nlon)
  lat!0 = "lat"
  lon!0 = "lon"
  lat@units = "degrees_north"
  lon@units = "degrees_east"
  lat&lat = lat
  lon&lon = lon

  elev!0 = "lat"
  elev!1 = "lon"
  elev&lat = lat
  elev&lon = lon
  
  return(elev)
  end


;--------------------main code-------------
begin

;1 读变量-----------------------------------------------------
  pic_dir ="/mnt/d/ncl_related/pictures/eqt/ERA5/solo/"
  data_dir= "/mnt/d/data/ERA5/"   
  data    = addfile(data_dir+"ERA5.nc", "r")  
  time    =data->time
  YYYYMMDDHH_f=cd_calendar(time,-3)

  need_time=2008102706
  time_locat= ind(YYYYMMDDHH_f.eq.need_time) ;找该时间所处位置
  path_time=tostring(need_time)

  u  = short2flt(data ->u(time_locat,{100:1000},{20:40},{85:105})) 
  v  = short2flt(data ->v(time_locat,{100:1000},{20:40},{85:105}))
  tmp= short2flt(data ->t(time_locat,{100:1000},{20:40},{85:105}))
  vvel  =short2flt(data->w(time_locat,{100:1000},{20:40},{85:105}))  
  rh    =short2flt(data->r(time_locat,{100:1000},{20:40},{85:105}))  
  vo    =short2flt(data->vo(time_locat,{100:1000},{20:40},{85:105})) 
  lv    = data->level({100:1000})  ;hPa
  lon   = data->longitude({85:105})
  lat   = data->latitude({20:40})
  lv1      =lv*100
  lv1@units="Pa"

;3 读地形文件并处理-----------------------------------------
  stdlon =95.8;选择固定某一角度 波密95.77,29.87    
  datadir  ="/mnt/d/data/other/"   ;地形数据路径
  elev     = read_height_data(datadir+"ETOPO5.DAT")
  lat1     = fspan(90,-90,2160)  ;用于gsn_xy 地形图
  ;; Convert terrain data from units "m" to "hPa", it is described as a high pressure formula
  elev = 1013.25*(1-elev*0.0065/288.15)^5.25145
  ;变量的气压维度单位是Pa，这里乘100,这里改了记得etores改
  elev1=elev(:,{stdlon})       ;固定某一经度的地形数据，用于叠加高分辨率的灰色地形

;4 rh设定 --------------------------------------------------
  ;; Determine the terrain
  rh_lon =rh(:,:,{stdlon})

;5 速度设定;将Pa/s转换成m/s,垂直速度转换 -------------------
  w1 = omega_to_w(vvel, conform(vvel,lv1, 0), tmp)  ; w[*][*][*] (m/s)，要求温度K，气压Pa
  w =w1*100.
  copy_VarCoords(w1, w)
  w_lon =w(:,:,{stdlon})
  v_lon =v(:,:,{stdlon})


;6 计算相当位温---------------------------------------------
  eqt           = pot_temp_equiv(lv1,tmp, rh, 0, "rh")  ;可直接使用这个函数计算,用K,Pa来算
  eqt_lon =eqt(:,:,{stdlon})


;7 计算湿位涡-----------------------------------------------
  g      =9.8
  dTdp   = center_finite_diff_n(eqt,lv1,False,0,0)
  dudp   = center_finite_diff_n(u,lv1,False,0,0)
  dvdp   = center_finite_diff_n(v,lv1,False,0,0)

  nlat=dimsizes(lat)
  mlon=dimsizes(lon)
  dlon   = (lon(2)-lon(1))*0.0174533 ; convert to radians 0.017=3.14/180
  dTdX   = new ( dimsizes(eqt),typeof(eqt)); pre-allocate space
  do nl=0,nlat-1                      ; loop over each latitude
    dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
    dTdX(:,nl:nl,:) = center_finite_diff_n (eqt(:,nl:nl,:), dX , False,0,2) ;False代表not cyclic
  end do
  dY   =6378388.*0.0174533*(lat(2)-lat(1))
  dTdY = new ( dimsizes(eqt),typeof(eqt)); pre-allocate space
  do nl=0,mlon-1                      ; loop over each latitude
    dTdY(:,:,nl:nl) = center_finite_diff_n (eqt(:,:,nl:nl), dY , False,0,1)
  end do

  f1d =new(dimsizes(lat), float)
  f1d =coriolis_param(lat)
  f2d =conform_dims(dimsizes(vo),f1d, 1) 
  copy_VarCoords(vo, f2d)
  avo =vo+f2d
  copy_VarCoords(vo, avo)

  mpv1=(-g*avo*dTdp)*1.e+6 
  mpv2=(g*(dvdp*dTdX-dudp*dTdY))*1.e+6
  copy_VarCoords(rh, mpv1)
  copy_VarCoords(rh, mpv2)
  mpv =mpv1+mpv2
  copy_VarCoords(rh, mpv) 

  mpv_lon =mpv(:,:,{stdlon})
  mpv1_lon =mpv1(:,:,{stdlon})


;8 绘图---------------------------------
wks = gsn_open_wks("ps",pic_dir+stdlon+"E_w_rh_eqt_mpv_"+path_time) 

res    = True
etores = True
rhres  = True
optres  = True
eqtres  = True
mpvres  = True 
mpv1res = True
wres    =True
set_all_res(res,etores,rhres,optres,eqtres,mpvres,mpv1res,wres)


plot_mpv = gsn_csm_contour(wks,mpv_lon({150:1000},:),mpvres)
plot_hgt = gsn_csm_xy(wks,lat1,elev1,etores) 
; plot_mpv1 = gsn_csm_contour(wks,mpv1_lon({150:1000},:),mpv1res)
plot_eqt = gsn_csm_contour(wks,eqt_lon({150:1000},:),eqtres)
plot_rh = gsn_csm_contour(wks,rh_lon({150:1000},:),rhres)
plot_rh = gsn_contour_shade(plot_rh,0,90, optres)
plot_w=gsn_csm_vector(wks,v_lon({150:1000},:),w_lon({150:1000},:),wres)


overlay(plot_mpv,plot_rh)
overlay(plot_mpv,plot_w)
overlay(plot_mpv,plot_eqt)
; overlay(plot_mpv,plot_mpv1)
overlay(plot_mpv,plot_hgt)

draw(plot_mpv)
frame(wks)
end