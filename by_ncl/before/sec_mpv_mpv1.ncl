;计算相当位温，垂直环流，MPV1,MPV,画区域平均沿996E的剖面图
;These files are loaded by default in NCL V6.2.0 and newer
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"


;;-------------------绘图属性----------------------------
undef("set_all_res")
procedure set_all_res(res,etores,MPVres,MPV1res,MPV2res,eqtres,vores)
 local res,etores,MPVres,MPV1res,MPV2res,eqtres,vores
 begin
 res = True
    res@gsnDraw  = False
    res@gsnFrame = False
    res@tmXTOn   = False ;关刻度
    res@tmYROn   = False
    res@gsnMaximize   =True
    res@gsnLeftString = "" ;左上角字符
    res@gsnRightString= "" ;左上角字符
    ;res@gsnMajorLatSpacing   = 2 ;主刻度纬度间隔为2
    res@gsnMajorLonSpacing    =5. 
    res@gsnMinorLonSpacing    =1.    
    res@gsnMajorLatSpacing    =5. 
    res@gsnMinorLatSpacing    =1
    res@tiYAxisString ="P/hPa"
    res@tiYAxisFont   =22 ;坐标轴标题的字体
    ;res@tiYAxisFontHeightF=0.5
    res@tmYRMode      = "Automatic" 
    ; 画高度剖面图时 turn off special labels on right axis


    ;设置坐标轴刻度线边上标签粗细
    res@tmXBLabelFont=22
    res@tmYLLabelFont=22
    res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
    res@tmYLLabelFontHeightF =0.02



 ;地形填充属性
 etores=True
    ;etores@gsnBelowYRefLineColor = "black" 
    ;etores@gsnXYBelowFillColors="black"
    etores@gsnMaximize      = True
    etores@trYMaxF=1000
    etores@gsnYRefLine      = 1000.          ; create a reference line 
    etores@gsnBelowYRefLineColor = "black"              ; above ref line fill red
    etores@gsnDraw        = False             ; don't draw the plots yet
    etores@gsnFrame       = False        ; don't advance the frame

  MPVres=res
    MPVres@cnFillOn   = True  ;cn means ContourPlot
    MPVres@cnLinesOn  = True   
    MPVres@cnLineColor="black"
    MPVres@cnFillPalette = "WhiteBlue"
    ;res@cnLineLabelsOn = False 
    ;res@cnLineLabelFormat="0@;*.1f" 
    ;res@cnLineLabelInterval=0.05
    MPVres@cnInfoLabelOn = False ; set false will not set an information label
    MPVres@cnLevelSelectionMode = "ExplicitLevels" 
    MPVres@cnLevels     =(/-1,-0.8,-0.4,-0.2,0/)
    MPVres@cnFillColors =(/240,200,160,120,80,0/)
    MPVres@lbBoxLinesOn = True ; Label bar
    MPVres@lbLabelAutoStride    = True            ; let NCL determine label spacing
    ;res@lbOrientation        = "Vertical"         ; vertical label bar
    ;res@pmLabelBarWidthF     = 0.05     ;色标变窄
    MPVres@pmLabelBarHeightF  =0.05 
    MPVres@lbBoxEndCapStyle   ="TriangleBothEnds"
    MPVres@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离

    MPVres@cnFillDrawOrder    = "Predraw"  ;调整画图的顺序 

    ;矢量叠加在contour上，属性设置合并
    MPVres@vcGlyphStyle          = "LineArrow"
    MPVres@vcLineArrowThicknessF = 4
    MPVres@vcMinDistanceF        = 0.025
    MPVres@vcRefMagnitudeF       = 15 ;设置单位量级
    MPVres@vcRefLengthF          = 0.02
    MPVres@vcLineArrowHeadMaxSizeF = 0.01
    MPVres@vcLineArrowHeadMinSizeF = 0.01
    MPVres@vcRefAnnoOn           = True ;anno代表矢量注释
    MPVres@vcRefAnnoFontHeightF  = 0.015
    MPVres@vcRefAnnoZone         = 1
    MPVres@vcRefAnnoSide         = "Top"
    MPVres@vcRefAnnoJust         = "TopRight"
    ;res@vcRefAnnoString1 = "5 g cm~S~-1~N~hPa~S~-1~Ns~S~-1~N)"
    MPVres@vcRefAnnoString1      = "15"
    MPVres@vcRefAnnoString2On    = False

  MPV1res=res
    MPV1res@cnFillOn  = False  ;cn means ContourPlot
    MPV1res@cnLinesOn = True   
    MPV1res@cnLineColor         ="red"
    MPV1res@cnLineLabelFontColor="red"
    MPV1res@cnLevelSelectionMode= "ManualLevels"   
    ;MPV1res@cnLevelSpacingF     =0.5
    ;MPV1res@cnLineDashPattern  =14
    MPV1res@gsnContourZeroLineThicknessF=8
    MPV1res@cnLineThicknessF    = 5
    ;MPV1res@cnLevelSelectionMode = "ExplicitLevels" 
    ;MPV1res@cnLevels=(/0,0.5,1,1.5,2/)
    MPV1res@gsnContourPosLineDashPattern =14 ; sets positive contours to dash pattern 14
    MPV1res@cnLabelMasking      = True ;线条不穿过标签
    MPV1res@cnLineLabelBackgroundColor   = -1 ;contour label背景透明
    MPV1res@cnLineLabelPlacementMode     = "Constant"
    MPV1res@cnLineDashSegLenF   = 0.2 
    MPV1res@cnInfoLabelOn       = False ; set false will not set an information label
    MPV1res@cnLineDrawOrder     = "Predraw"  ;调整画图的顺序

  MPV2res=res
    MPV2res@cnFillOn  = False  ;cn means ContourPlot
    MPV2res@cnLinesOn = True   
    MPV2res@cnLineColor         ="purple"
    MPV2res@cnLineLabelFontColor="purple"
    MPV2res@cnLevelSelectionMode= "ManualLevels"  
    ;MPV2res@cnLevelSpacingF     =0.5
    ;MPV2res@cnLineDashPattern  =14
    MPV2res@gsnContourZeroLineThicknessF=8
    MPV2res@cnLineThicknessF    = 5 
    ;MPV2res@cnLevelSelectionMode = "ExplicitLevels" 
    ;MPV2res@cnLevels=(/0,0.5,1,1.5,2/)
    MPV2res@gsnContourPosLineDashPattern =14 ; sets positive contours to dash pattern 14
    MPV2res@cnLabelMasking      = True ;线条不穿过标签
    MPV2res@cnLineLabelBackgroundColor   = -1 ;contour label背景透明
    MPV2res@cnLineLabelPlacementMode     = "Constant"
    MPV2res@cnLineDashSegLenF   = 0.2 
    MPV2res@cnInfoLabelOn       = False ; set false will not set an information label
    MPV2res@cnLineDrawOrder     = "Predraw"  ;调整画图的顺序

 ;相当位温 
 eqtres=res
    eqtres@cnFillOn = False  ;cn means ContourPlot
    eqtres@cnLinesOn = True 
    ;res@cnMonoLineColor = False
    eqtres@cnLineColor="black"
    eqtres@cnLabelMasking = True
    eqtres@cnLineLabelPlacementMode = "Constant" ;使标签是线的一部分
    eqtres@cnLineDashSegLenF    = 0.3
    eqtres@cnLineLabelBackgroundColor = -1 ;设置标签背景透明
    eqtres@cnInfoLabelOn = False ; set false will not set an information label
    eqtres@cnLineThicknessF = 6
    eqtres@cnLevelSelectionMode = "ManualLevels" 
    eqtres@cnLevelSpacingF = 2   

 vores=res
    vores@cnFillOn = True  ;cn means ContourPlot
    vores@cnLinesOn = False  
    vores@cnInfoLabelOn = False ; set false will not set an information label 
    vores@cnFillPalette = "BlueWhiteOrangeRed"
    ;vores@cnLevelSelectionMode = "ExplicitLevels" 
    ;vores@cnLevels=(/-0.3,-0.2,-0.1,-0.05,0/)
    ;vores@cnFillColors=(/240,200,160,120,80,0/)
    vores@cnLevelSelectionMode = "ManualLevels" 
    ;res@cnLevelSpacingF = 4 
    ;res@cnLineThicknessF = 4
    vores@lbBoxLinesOn = True ; Label bar
    vores@lbLabelAutoStride     = True            ; let NCL determine label spacing
    ;res@lbOrientation        = "Vertical"         ; vertical label bar
    ;res@pmLabelBarWidthF = 0.05     ;色标变窄
    vores@pmLabelBarHeightF =0.05 
    vores@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离

    vores@lbBoxEndCapStyle = "TriangleBothEnds"
    vores@cnFillDrawOrder      = "Predraw"  ;调整画图的顺序


end

undef("read_elev_data")
function read_elev_data(topo_file)
local nlat, nlon, topo_file, lat, lon
begin
 ;---Read data as a straight binary file
  nlat = 2160
  nlon = 4320
  setfileoption("bin","ReadByteOrder","BigEndian")
  elev= cbinread(topo_file,(/nlat,nlon/),"short")

 ;---Create 1D coordinate arrays
  lat       = fspan(90,-90,nlat)
  lon       = fspan(0,360,nlon)
  lat!0     = "lat"
  lon!0     = "lon"
  lat@units = "degrees_north"
  lon@units = "degrees_east"
  lat&lat   = lat
  lon&lon   = lon

 ;---Attach the coordinate arrays
  elev!0    = "lat"
  elev!1    = "lon"
  elev&lat  = lat
  elev&lon  = lon
 
  return(elev)
end


;-------------main code----------------------------------

begin

 pic_dir ="/cygdrive/D/ncl_related/pictures/MPV/solo/"

 ;读地形文件----------------------------
   ;lon1       fspan(0,360,4320)  ;gsn_xy    
   lat1     =fspan(90,-90,2160)  ;用于gsn_xy 地形图
   data1_dir  ="/cygdrive/D/ncl_related/data/other/"       ;地形数据路径
   setfileoption("nc","MissingToFillValue",False)    ;设置缺省值
   elev= read_elev_data(data1_dir+"ETOPO5.DAT")
   ;printVarSummary(elev)

 ;处理地形数据，将地形数据的单位：m转化为hPa
   elev0 = 1013.25*(1-elev*0.0065/288.15)^5.25145  
   copy_VarCoords(elev,elev0) 
   ;printVarSummary(elev0)

   stdlon=96;选择固定某一角度
    elev1=elev0(:,{stdlon})       ;固定某一经度的地形数据
    elev1@_FillValue= -9.96921e+36


;读变量-----------------------------------------------------
  ;;;>>>手动设置一个时间 用各个数据生成一个带unit的 timeset_ymdh 并且转换为一个字符串 方便文件中读入
   year   = 2008
   month  = 10
   day    = 26
   hour   = 18
   minute = 00
   second = 00
  ;设置时间格式读取数据####################################
   units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
   opt = 0
   timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
   file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)
   data2_dir= "/cygdrive/D/ncl_related/data/fnl/"  ; data1_dir是地形，这里略去了  
   data   = addfile(data2_dir+"fnl_"+file_time+".nc", "r")                        

    uData  = data ->UGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s lv_ISBL0=26层
    vData  = data->VGRD_P0_L100_GLL0(5:25,:,:) ;26
    vvel   = data->VVEL_P0_L100_GLL0( :, :, :);float VVEL_P0_L100_GLL0 ( lv_ISBL4, lat_0, lon_0 ),21层
    tmpdata=data->TMP_P0_L100_GLL0 (5:25,:,: ) ;单位：K 26层
    avo    =data->ABSV_P0_L100_GLL0 (5:25,:,:);26层
    rh     =data->RH_P0_L100_GLL0 ( :, :, :) ;lv_ISBL4=21层 %
    pres1   =data->lv_ISBL4  ;Pa
    lon    =data->lon_0
    lat    =data->lat_0
    pres   =pres1/100
    pres@units="hPa"
    tmp    =tmpdata-273.15
    copy_VarCoords(rh, tmp)

;printVarSummary(uData)
;print(pres)
;exit 

;计算垂直速度矢量
  w1 = omega_to_w(vvel, conform(vvel,pres1, 0), tmpdata)  ; w[*][*][*] (m/s)，要求温度K，气压Pa
  w =w1*100
  copy_VarCoords(w1, w)
  h_lw = w(:,:,{stdlon})
  h_lw!0="level"
  h_lw&level=pres

;设定经向风速v经度
  h_lv = vData(:,:,{stdlon})
  h_lv!0="level"
  h_lv&level=pres

;----------------------------计算相当位温-----------------------------------
 eqt= pot_temp_equiv(pres1,tmpdata, rh, 0, "rh")  ;可直接使用这个函数计算,用Pa来算
 copy_VarCoords(tmp, eqt)
  ;printVarSummary(eqt)
  ;print(pres)
  ;exit 
  h_leqt  =eqt(:,:,{stdlon})
  ;h_leqt = dim_avg_n_Wrap(eqt(:,:,{92:98}), 2)
  h_leqt!0="level"
  h_leqt&level=pres


  
;-----------------------------计算湿位涡----------------------------    
 g      =9.8
 dTdp   = center_finite_diff_n(eqt,pres1,False,0,0)
 dudp   = center_finite_diff_n(uData,pres1,False,0,0)
 dvdp   = center_finite_diff_n(vData,pres1,False,0,0)

 nlat=dimsizes(lat)
 mlon=dimsizes(lon)

 dlon   = (lon(2)-lon(1))*0.0174533 ; convert to radians 0.017=3.14/180
 dTdX   = new ( dimsizes(eqt),typeof(eqt)); pre-allocate space
 do nl=0,nlat-1                      ; loop over each latitude
  dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
  dTdX(:,nl:nl,:) = center_finite_diff_n (eqt(:,nl:nl,:), dX , True,0,2)
 end do

 dY   =6378388.*0.0174533
 dTdY = new ( dimsizes(eqt),typeof(eqt)); pre-allocate space
 do nl=0,mlon-1                      ; loop over each latitude
  dTdY(:,:,nl:nl) = center_finite_diff_n (eqt(:,:,nl:nl), dY , False,0,1)
 end do

 mpv1=-g*avo*dTdp 
 mpv2=g*(dvdp*dTdX-dudp*dTdY)
 mpv =mpv1+mpv2
 copy_VarCoords(tmp, mpv1)
 copy_VarCoords(tmp, mpv2)
 copy_VarCoords(tmp, mpv)
 copy_VarCoords(tmp, mpv)

 ;h_ldTdp   = dim_avg_n_Wrap(dTdp(:,:,{92:98}), 2)


 ;h_lmpv       =dim_avg_n_Wrap(mpv(:,:,{92:98}), 2)*1.e+6
 h_lmpv      =mpv(:,:,{stdlon})*1.e+6
 copy_VarCoords(tmp, h_lmpv)
 h_lmpv!0     ="level"
 h_lmpv&level =pres

 ;h_lmpv1      =dim_avg_n_Wrap(mpv1(:,:,{92:98}), 2)*1.e+6
 h_lmpv1      =mpv1(:,:,{stdlon})*1.e+6
 copy_VarCoords(tmp, h_lmpv1)
 h_lmpv1!0    ="level"
 h_lmpv1&level=pres

 ;h_lmpv2      =dim_avg_n_Wrap(mpv2(:,:,{92:98}), 2)*1.e+6
 h_lmpv2      =mpv2(:,:,{stdlon})*1.e+6
 copy_VarCoords(tmp, h_lmpv2)
 h_lmpv2!0    ="level"
 h_lmpv2&level=pres

;----------------------------计算垂直涡度------------------
  vr = uv2vr_cfd (uData,vData,lat,lon, 2) 
  copy_VarCoords(rh, vr)

  h_lvr  =vr(:,:,{stdlon})*1.e+5
  ;h_lvr = 1.e+5*dim_avg_n_Wrap(vr(:,:,{92:98}), 2)
  copy_VarCoords(rh, h_lvr)
  h_lvr!0="level"
  h_lvr&level=pres
  h_lvr!1="lat"
  h_lvr&lat=lat
  ;printVarSummary(h_lvr)
  ;exit

;------------绘图---------------------

wks = gsn_open_wks("png",pic_dir+ stdlon+"E_sec_mpv_mpv1_"+file_time) 
;wks = gsn_open_wks("png",pic_dir+"96E_sec_vr_mpv2_"+file_time) 

;gsn_define_colormap(wks,"BlueWhiteOrangeRed") 


res    = True
 etores=True
 MPVres = True
 MPV1res = True
 MPV2res = True
 eqtres = True
 vores  = True
 set_all_res(res,etores,MPVres,MPV1res,MPV2res,eqtres,vores)


strlat=15
endlat=35
plot_eqt   = gsn_csm_pres_hgt(wks, h_leqt(:,{strlat:endlat}), eqtres ) 
;plot_vr= gsn_csm_pres_hgt(wks, h_lvr(:,{25:35}), vores ) 
;plot_vr= gsn_csm_pres_hgt_vector(wks,h_lvr(:,{25:35}),h_lv(:,{25:35}),h_lw(:,{25:35}),vores)
;plot_MPV  = gsn_csm_pres_hgt(wks, h_lmpv(:,{25:35}), MPVres )
plot_MPV  = gsn_csm_pres_hgt_vector(wks,h_lmpv(:,{strlat:endlat}),h_lv(:,{strlat:endlat}),h_lw(:,{strlat:endlat}),MPVres)  
plot_MPV1 = gsn_csm_pres_hgt(wks, h_lmpv1(:,{strlat:endlat}), MPV1res ) 
;plot_MPV2  = gsn_csm_pres_hgt(wks, h_lmpv2(:,{strlat:endlat}), MPV2res ) 
plot_hgt = gsn_csm_xy(wks,lat1,elev1,etores) 



overlay(plot_MPV,plot_eqt)
overlay(plot_MPV,plot_MPV1)
overlay(plot_MPV,plot_hgt)
;overlay(plot_vr,plot_MPV2)
;overlay(plot_vr,plot_eqt)
;overlay(plot_vr, plot_hgt)



draw(plot_MPV)
;draw(plot_vr)
frame(wks)
end