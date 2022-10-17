;波密95.8E的垂直速度矢量，相对湿度填色，相当位温等值线填色
;created by yewei ，2020
;地形剖面填充灰色，法二：用叠图的方式，推荐
; 参考;http://www.ncl.ucar.edu/Applications/Scripts/mask_14.ncl
;湿层填充图案参考 http://www.ncl.ucar.edu/Document/Graphics/Interfaces/gsn_contour_shade.shtml
;气压层单位是Pa和hPa的时候画的图不一致，Pa的时候明显不对，上升运动区的矢量箭头向下
;！！！遇到的问题： y轴标签始终不能设置成等距的等压坐标；
; 800hPa和1000hPa的距离始终比别的200hPa间隔大
;解决上述问题参考 https://www.ncl.ucar.edu/Applications/vertterp.shtml vert_2.ncl
;  ncl /mnt/h/d/ncl_related/ERA5_code/sec_w_rh_eqt_mpv_mpv1_topo2.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"




;--------------------main code-------------
begin

;1 读变量-----------------------------------------------------
  pic_dir ="/mnt/h/d/ncl_related/pictures/eqt/ERA5/solo/topo2/"
  data_dir= "/mnt/h/d/data/ERA5/"   
  data    = addfile(data_dir+"ERA5.nc", "r")  
  time    =data->time
  YYYYMMDDHH_f=cd_calendar(time,-3)

count=3
  need_time=(/2008102612,2008102612,2008102706,2008102706/)
  plot_str=(/"(a)26/12","(b)26/12","(c)27/06","(d)27/06"/)
  info_str=(/"","Zayu(2331m)","Cona(4281m)","Bomi(2737m)"/)
  plot_lon =(/94.0,97.5,92.0,95.8/)
  plot_lat =(/15.0,28.7,28.0,29.9/)
  stdlon     = plot_lon(count)  ;选择固定某一角度 波密95.8,29.9
  stdlat     =plot_lat(count)    

  time_locat= ind(YYYYMMDDHH_f.eq.need_time(count)) ;找该时间所处位置
  path_time=tostring(need_time(count))

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
  topo        = addfile("/mnt/h/d/data/other/ETOPO2v2c_f4.nc", "r")
  elev0       = topo->z 
  y           = topo->y
  lat1        = y  ;用于gsn_xy 地形图
  ;; Convert terrain data from units "m" to "hPa", it is described as a high pressure formula
  ; elev        = 1013.25*(1-elev*0.0065/288.15)^5.25145 ;手动转换 高度转气压 hPa
  elev_transf=stdatmus_z2tdp(elev0) ;函数转换 两种得到的高度画出来是相同的
  elev=elev_transf(2,:,:)
  copy_VarCoords(elev0, elev)
  ;如果上面的转换乘100变成Pa之后,记得etores也要改

  ;==========固定某一经度的地形数据，用于叠加高分辨率的灰色地形==================
  elev1=elev(:,{stdlon})     

  ; ;==========固定取三站的地形平均值，用于叠加高分辨率的灰色地形=================
  ; elev_rec=elev(:,{plot_lon(1:3)})
  ; elev1=dim_avg_n_Wrap(elev_rec, 1)
  ; ; printVarSummary(elev1)
  ; ; exit  



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


;8 ---------------------------------绘图---------------------------------
wks = gsn_open_wks("x11",pic_dir+stdlon+"E_w_rh_eqt_"+path_time+"_topo2") 

res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize   =True
  res@tmXTOn = False ;关刻度
  res@tmYROn = False
  res@gsnLeftString = "" ;左上角字符
  res@gsnRightString=""
  res@tiYAxisString ="P(hPa)"   
  ; ; 下面这个是默认的
  res@trXMinF              = 24    ; Set min/max of X axis
  res@trXMaxF              = 33
  res@tmXBTickSpacingF     = 2 ;主刻度线间隔
  res@tmXBLabels=(/"24~S~o~N~N","26~S~o~N~N","28~S~o~N~N","30~S~o~N~N",\
                    "32~S~o~N~N"/)


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
  res@gsnPaperOrientation   = "Portrait"            ;portra和maximize必须同时设置才会生效

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

  eqtres@gsnLeftStringFontHeightF=0.02
  eqtres@gsnLeftString=plot_str(count)
  eqtres@gsnLeftStringOrthogonalPosF=-0.085
  eqtres@gsnLeftStringParallelPosF =0.01
  eqtres@txPerimOn=True
  eqtres@txPerimThicknessF=2.5
  eqtres@txBackgroundFillColor ="white"


  ;相对湿度
  rhres=res
  rhres@cnFillOn          = True                  ; turn on color
  rhres@cnLinesOn         = False
  ; rhres@cnFillPattern     = 17
  ; rhres@cnFillDotSizeF    = 0.003            ; increase dot size of contours
  ; rhres@cnFillScaleF      = 0.3 ; densities
  rhres@cnInfoLabelOn     = False
  rhres@cnLevelSelectionMode = "ExplicitLevels" 
  rhres@cnLevels          =(/90,98/)
  rhres@cnFillColors       =(/"white","skyblue","khaki"/)
  ; rhres@cnLineColor       ="royalblue"
  ; rhres@cnLineThicknessF  = 4
  ; rhres@cnLineLabelsOn    = False
  ; rhres@cnLineDrawOrder   ="PreDraw"
  rhres@cnFillDrawOrder   ="PreDraw"
  rhres@lbLabelStrings=(/"90%","98%"/)
  rhres@lbBoxEndCapStyle     = "TriangleBothEnds"
  rhres@pmLabelBarHeightF =0.1
  rhres@pmLabelBarWidthF =0.3



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
  wres@vcRefAnnoOrthogonalPosF=0

print(eqt_lon({700},{24.5}))
plot_rh = gsn_csm_contour(wks,rh_lon({150:1000},:),rhres)
; plot_mpv = gsn_csm_contour(wks,mpv_lon({150:1000},:),mpvres)
plot_hgt = gsn_csm_xy(wks,lat1,elev1,etores) 
; plot_mpv1 = gsn_csm_contour(wks,mpv1_lon({150:1000},:),mpv1res)
plot_eqt = gsn_csm_contour(wks,eqt_lon({150:1000},:),eqtres)
; plot_rh = gsn_contour_shade(plot_rh,0,90, optres)
plot_w=gsn_csm_vector(wks,v_lon({150:1000},:),w_lon({150:1000},:),wres)


overlay(plot_rh,plot_w)
overlay(plot_rh,plot_eqt)
overlay(plot_rh,plot_hgt)

;某时刻风暴位置标记
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.04                     ; font smaller. default big
  txres@txFontColor   ="red"
  dum1 = gsn_add_text(wks,plot_rh,"~F35~z",stdlat,1000,txres)
  txres@txFontHeightF = 0.025                     ; font smaller. default big
  dum2 = gsn_add_text(wks,plot_rh,info_str(count),stdlat+0.1,960,txres)  
  ;站点名称右挪一点点

draw(plot_rh)
frame(wks)
end