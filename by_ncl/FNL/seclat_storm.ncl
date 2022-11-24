;将1*1插值成0.1*0.1
;过各时刻台风中心纬向剖面，目前，环流，温度距平，高度距平
;垂直速度矢量，相当位温等值线，相对湿度>90%填色,温度距平，高度距平
;换时次时改两处，需改动61和最后风暴标记的位置
;由于纬向移动较小，横坐标可以统一用80-100的范围
;代码路径  ncl /mnt/d/ncl_related/code/seclat_storm.ncl


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

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


begin

;1. 读地形文件
  lon1       =fspan(0,360,4320)  ;gsn_xy    
  ;lat1     =fspan(90,-90,2160)  ;用于gsn_xy 地形图
  ;datadir="/mnt/d/data/"
  datadir  ="/mnt/d/data/other/"   ;地形数据路径
  setfileoption("nc","MissingToFillValue",False)    ;设置缺省值
  elev= read_elev_data(datadir+"ETOPO5.DAT")
  ;printVarSummary(elev)
  elev0 = 1013.25*(1-elev*0.0065/288.15)^5.25145 ;将地形数据的单位：m转化为hPa
  copy_VarCoords(elev,elev0) 
  ;printVarSummary(elev0)

;2 读变量
  ;手动设置一个时间 用各个数据生成一个带unit的 timeset_ymdh 
  ;并且转换为一个字符串 方便文件中读入
  year   = 2008
  month  = 10
  day    = 26
  hour   = 12
  minute = 00
  second = 0
  ;设置时间格式读取数据####################################
  units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
  opt = 0
  timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
  file_format  = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
  file_time    = cd_string(timeset_ymdh,file_format)
  path_format  = "%Y%N%D%H"   ; 台风路径文件中的时间形式     
  path_time    = cd_string(timeset_ymdh,path_format)

  path_rawdata= "/mnt/d/data/fnl/"  ; 变量数据路径
  data       = addfile(path_rawdata+"fnl_"+file_time+".nc", "r") 
  rh_in      = data->RH_P0_L100_GLL0   ;lv_ISBL4/lat_0/lon_0,21层
  tmp_in = data->TMP_P0_L100_GLL0(5:25,:,: ) ;单位：K,lv_ISBL0/lat_0/lon_0,26层，选择对应的21层
  u_in       = data->UGRD_P0_L100_GLL0(5:25,:,:) 
  v_in       = data->VGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
  vvel_in    = data->VVEL_P0_L100_GLL0(:,:,: ) ;Pa/s,lv_ISBL4, lat_0, lon_0 
  hgtData    =data->HGT_P0_L100_GLL0(5:25,:,:)  ;hgt 的单位：gpm
  hgt_in     = hgtData/10.
  copy_VarCoords(hgtData, hgt_in)
  pres1      =data->lv_ISBL4  ;Pa
  lon_in     =data->lon_0
  lat_in     =data->lat_0
  pres       =pres1/100
  pres@units ="hPa"

;3 读取风暴中心位置  
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

  ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
  tctime_ind = str_match_ind(YYYYMMDDHH, path_time) 
  stdlon=0.
  stdlat=0.
  do i = 0, nrow-1
    ;if ((DD(i).eq.day).and.(HH(i).eq.hour)) then
    ; if ((DD(i).eq."27").and.(HH(i).eq."00")) then
    if (i.eq.tctime_ind) then ;最方便
      stdlon=lon(i)
      stdlat=lat(i)
    end if
  end do
; lon0=floattoint(stdlon);用于绘图范围
; lat0=floattoint(stdlat)
;4 选择固定某一角度
  elev1=elev0({stdlat},:)       ;固定某一经度的地形数据    
  elev1@_FillValue= -9.96921e+36


;5 网格点插值
  lon_out=fspan(70,110,401)
  lat_out=fspan(10,50, 401)
  u      =linint2_Wrap(lon_in, lat_in(::-1), u_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(u_in)
  v      =linint2_Wrap(lon_in, lat_in(::-1), v_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(v_in)
  vvel   =linint2_Wrap(lon_in, lat_in(::-1), vvel_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(vvel_in)
  rh     =linint2_Wrap(lon_in, lat_in(::-1), rh_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(rh_in)
  tmp    =linint2_Wrap(lon_in, lat_in(::-1), tmp_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(tmp_in)
  hgt      =linint2_Wrap(lon_in, lat_in(::-1), hgt_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(hgt_in)

;6 设定垂直速度矢量
  ;设定纬向风速，将Pa/s转换成m/s   
  w1 = omega_to_w(vvel, conform(vvel,pres1, 0), tmp)  ; w[*][*][*] (m/s)，要求温度K，气压Pa
  w =w1*100.
  copy_VarCoords(w1, w)
  h_lw = w(:,{stdlat},:)
  h_lw!0="level" ;有这步没这步画出来的图都不一样
  h_lw&level=pres
  ;设定经向风速v经度
  h_lu = u(:,{stdlat},:)
  h_lu!0="level"
  h_lu&level=pres

; ;7 计算相当位温
;   eqt= pot_temp_equiv(pres1,tmp, rh, 0, "rh")  ;可直接使用这个函数计算,用Pa来算
;   copy_VarCoords(tmp, eqt)
;   h_leqt  =eqt(:,{stdlat},:)
;   h_leqt!0="level"
;   h_leqt&level=pres

;8 设定相对湿度
;   h_lrh  =rh(:,{stdlat},:)
;   h_lrh!0="level"
;   h_lrh&level=pres

; ;9 计算垂直涡度
;   vo = uv2vr_cfd (u,v,lat_out,lon_out, 2) 
;   copy_VarCoords(u, vo)
;   vo=vo*1.e5
;   h_lvo  =vo(:,{stdlat},:)
;   h_lvo!0="level"
;   h_lvo&level=pres

;10 计算距平
  min_lon=stdlon-10
  max_lon=stdlon+10
  min_lat=stdlat-10
  max_lat=stdlat+10
  ;温度距平
  tmp_area= tmp(:,{min_lat:max_lat},{min_lon:max_lon}) 
  tmp_avg=dim_avg_n_Wrap(tmp_area, (/1,2/)) 
  tmp_dev=tmp-conform(tmp,tmp_avg,0)
  copy_VarCoords(tmp, tmp_dev)
  h_ltmp_dev  =tmp_dev(:,{stdlat},:)
  h_ltmp_dev!0="level"
  h_ltmp_dev&level=pres
  ;高度距平
  hgt_area= hgt(:,{min_lat:max_lat},{min_lon:max_lon}) 
  hgt_avg=dim_avg_n_Wrap(hgt_area, (/1,2/)) 
  hgt_dev=hgt-conform(hgt,hgt_avg,0)
  copy_VarCoords(hgt, hgt_dev)
  h_lhgt_dev  =hgt_dev(:,{stdlat},:)
  h_lhgt_dev!0="level"
  h_lhgt_dev&level=pres


;11 打开绘图窗口
pic_dir ="/mnt/d/ncl_related/pictures/storm/sec_eqt/lat/"
wks = gsn_open_wks("png",pic_dir+stdlat+"N"+"_eqt_hgtdev_wind_"+file_time) 


;12 绘图属性设置
  res = True
  res@gsnDraw       = False
  res@gsnFrame      = False
  res@gsnMaximize   = True
  res@tmXTOn        = False ;关刻度
  res@tmYROn        = False
  res@gsnLeftString = "" ;左上角字符
  res@gsnRightString=""
  res@tiYAxisString ="P(hPa)"
  res@tmYRMode      = "Automatic" ;左右Y轴均显示气压层
  res@tmBorderThicknessF  =4    ; border 6x as thick,图框粗细
  res@tmXBMajorThicknessF =4
  res@tmXBMinorThicknessF =4
  res@tmYLMajorThicknessF =4 
  res@tmYLMinorThicknessF =4   

  ;设置坐标轴
  res@tmXBMode       ="Explicit"
  res@tmXBValues     =(/80,85,90,95,100/)
  res@tmXBLabels     =(/"80E","85E","90E","95E","100E"/)
  res@tmXBMinorValues=(/81,82,83,84,86,87,88,89,91,92,93,94,96,97,98,99/)

  ; ;相当位温 
  ; eqtres=res
  ; eqtres@cnFillOn      = False  ;cn means ContourPlot
  ; eqtres@cnLinesOn     = True 
  ; ;res@cnMonoLineColor = False
  ; eqtres@cnLineColor   ="black"
  ; eqtres@cnLabelMasking= True
  ; eqtres@cnLineLabelPlacementMode   = "Constant" ;使标签是线的一部分
  ; eqtres@cnLineDashSegLenF          = 0.3
  ; eqtres@cnLineLabelBackgroundColor = -1 ;设置标签背景透明
  ; eqtres@cnInfoLabelOn              = False ; set false will not set an information label
  ; eqtres@cnLineThicknessF           = 4
  ; eqtres@cnLevelSelectionMode       = "ManualLevels" 
  ; eqtres@cnLevelSpacingF            = 2  

  ; ;垂直速度矢量叠加在contour上，属性设置合并
  ; eqtres@vcMapDirection        = False ;绘制垂直剖面图要设置为False
  ; eqtres@vcGlyphStyle          = "LineArrow"
  ; eqtres@vcLineArrowThicknessF = 4
  ; eqtres@vcMinDistanceF        = 0.02
  ; eqtres@vcRefMagnitudeF       = 15 ;设置单位量级
  ; eqtres@vcRefLengthF          = 0.02
  ; eqtres@vcLineArrowHeadMaxSizeF = 0.01
  ; eqtres@vcLineArrowHeadMinSizeF = 0.01
  ; eqtres@vcRefAnnoOn           = True ;anno代表矢量注释
  ; eqtres@vcRefAnnoFontHeightF  = 0.015
  ; eqtres@vcRefAnnoZone         = 1
  ; eqtres@vcRefAnnoSide         = "Top"
  ; eqtres@vcRefAnnoJust         = "TopRight"
  ; ;eqtres@vcRefAnnoString1 = "5 g cm~S~-1~N~hPa~S~-1~Ns~S~-1~N)"
  ; eqtres@vcRefAnnoString1      = "15"
  ; eqtres@vcRefAnnoString2On    = False

  ; ;---------相对湿度-------------
  ; rhres=res
  ; rhres@cnFillOn          = True                  ; turn on color
  ; rhres@cnLinesOn         = False
  ; rhres@cnFillPalette     = "WhiteBlue"             ; set color map
  ; rhres@cnInfoLabelOn     =False
  ; rhres@cnLevelSelectionMode = "ExplicitLevels" 
  ; rhres@cnLevels          =(/90/)
  ; rhres@cnFillColors      =(/-1,80/)
  ; rhres@cnLineLabelsOn    =False
  ; rhres@lbLabelBarOn      =False
  ; ;rhres@cnLevelSelectionMode ="ManualLevels"
  ; rhres@cnFillDrawOrder   ="PreDraw"

  ; vores=res
  ; vores@cnFillOn      = True  ;cn means ContourPlot
  ; vores@cnLinesOn     = False  
  ; vores@cnInfoLabelOn = False ; set false will not set an information label 
  ; vores@cnFillPalette = "BlueWhiteOrangeRed"
  ; vores@cnLevelSelectionMode   = "ExplicitLevels" 
  ; vores@cnLevels      =(/-4,-2,0,2,4,6,8,10,12,14,16,18/)
  ; vores@cnFillColors  =(/45,65,75,148,168,178,188,198,208,218,228,235,240/)
  ; ;vores@cnLevelSelectionMode = "ManualLevels"
  ; ;vores@cnMinLevelValF  =-4
  ; ;vores@cnMaxLevelValF  =20
  ; ;vores@cnLevelSpacingF = 2
  ; ;res@cnLineThicknessF  = 4
  ; vores@lbBoxLinesOn           = True ; Label bar
  ; vores@lbLabelAutoStride      = True            ; let NCL determine label spacing
  ; ;res@lbOrientation           = "Vertical"         ; vertical label bar
  ; ;res@pmLabelBarWidthF        = 0.05     ;色标变窄
  ; vores@pmLabelBarHeightF      =0.05 
  ; vores@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
  ; vores@lbBoxEndCapStyle          = "TriangleBothEnds"
  ; vores@cnFillDrawOrder           = "Predraw"  ;调整画图的顺序

  tmpres=res
  tmpres@cnFillOn      = True  ;cn means ContourPlot
  tmpres@cnLinesOn     = False  
  tmpres@cnInfoLabelOn = False ; set false will not set an information label 
  tmpres@cnFillPalette = "BlueWhiteOrangeRed"
  tmpres@cnLevelSelectionMode ="ExplicitLevels" 
  tmpres@cnLevels      =(/-4,-3,-2,-1,0,1,2,3,4/)
  tmpres@cnFillColors  =(/25,45,65,75,110,130,178,188,200,230/)
  ; tmpres@cnMinLevelValF  =-4
  ; tmpres@cnMaxLevelValF  =4
  ; tmpres@cnLevelSpacingF = 1
  ;res@cnLineThicknessF  = 4
  tmpres@lbBoxLinesOn           = True ; Label bar
  tmpres@lbLabelAutoStride      = True            ; let NCL determine label spacing
  ;res@lbOrientation           = "Vertical"         ; vertical label bar
  tmpres@pmLabelBarWidthF       = 0.5     ;色标长度
  tmpres@pmLabelBarHeightF      = 0.07 
  ; tmpres@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
  tmpres@lbBoxEndCapStyle          = "TriangleBothEnds"
  tmpres@cnFillDrawOrder           = "Predraw"  ;调整画图的顺序

  hgtres=res
  hgtres@cnFillOn      = False  ;cn means ContourPlot
  hgtres@cnLinesOn     = True 
  ;res@cnMonoLineColor = False
  ; hgtres@cnLineColor   ="royalblue"
  ; hgtres@cnLineLabelFontColor="deepskyblue"
  hgtres@cnLineLabelFontHeightF=0.015
  ; hgtres@cnLineLabelPlacementMode   = "Constant" ;将标签绘制为虚线的一部分；默认Randomized
  hgtres@cnLineDashSegLenF          = 0.1 ;调整每根虚线中短线密集程度
  hgtres@cnLineLabelDensityF         =0.8
  ; hgtres@cnLabelMasking             = True
  hgtres@cnLineLabelBackgroundColor = "white" ;设置标签背景透明-1
  hgtres@cnInfoLabelOn              = False ; set false will not set an information label
  hgtres@cnLineThicknessF           = 6
  hgtres@gsnContourNegLineDashPattern=16
  hgtres@cnLevelSelectionMode       = "ManualLevels" 
  hgtres@cnLevelSpacingF            = 1 
  hgtres@cnLineDrawOrder            ="PreDraw"

  ;垂直速度矢量叠加在contour上，属性设置合并
  hgtres@vcMapDirection        = False ;绘制垂直剖面图要设置为False
  hgtres@vcGlyphStyle          = "LineArrow"
  hgtres@vcLineArrowThicknessF = 3
  hgtres@vcMinDistanceF        = 0.02
  hgtres@vcRefMagnitudeF       = 5 ;设置单位量级
  hgtres@vcRefLengthF          = 0.02
  hgtres@vcLineArrowHeadMaxSizeF = 0.01
  hgtres@vcLineArrowHeadMinSizeF = 0.01
  hgtres@vcRefAnnoOn           = True ;anno代表矢量注释
  hgtres@vcRefAnnoFontHeightF  = 0.015
  hgtres@vcRefAnnoZone         = 1
  hgtres@vcRefAnnoSide         = "Top"
  hgtres@vcRefAnnoJust         = "TopRight"
  ;hgtres@vcRefAnnoString1 = "5 g cm~S~-1~N~hPa~S~-1~Ns~S~-1~N)"
  hgtres@vcRefAnnoString1      = "5"
  hgtres@vcRefAnnoString2On    = False
  hgtres@vcVectorDrawOrder     ="PreDraw"


  etores=True
  etores@gsnMaximize      = True
  etores@gsnDraw          = False           ; don't draw the plots yet
  etores@gsnFrame         = False           ; don't advance the frame
  etores@trYMaxF          =1000
  etores@gsnYRefLine      =1000.          ; create a reference line 
  etores@gsnBelowYRefLineColor = "gray60"  ; above ref line fill red

;13 绘图
strlon=80
endlon=100
plot_base= gsn_csm_pres_hgt_vector(wks,h_lhgt_dev({200:1000},{strlon:endlon}),\
                                    h_lu({200:1000},{strlon:endlon}),\
                                    h_lw({200:1000}, {strlon:endlon}),hgtres)
; plot_rh   = gsn_csm_pres_hgt(wks, h_lrh({150:1000},{strlon:endlon}), rhres ) 
; plot_vo   = gsn_csm_pres_hgt(wks, h_lvo({150:1000},{strlon:endlon}), vores ) 
; plot_hgt   = gsn_csm_pres_hgt(wks, h_lhgt_dev({200:1000},{strlon:endlon}), hgtres ) 
plot_tmp   = gsn_csm_pres_hgt(wks, h_ltmp_dev({200:1000},{strlon:endlon}), tmpres ) 
plot_topo  = gsn_csm_xy(wks,lon1,elev1,etores) 

; overlay(plot_vector,plot_rh)
; overlay(plot_vector,plot_vo)
overlay(plot_tmp,plot_base)
overlay(plot_tmp,plot_topo)
draw(plot_tmp)

;14 某时刻风暴位置标记
  txres                 = True                     ; text mods desired
  txres@txFontHeightF   = 0.03                     ; font smaller. default big
  txres@txFontThicknessF= 4
  txres@txFontColor     ="red"
  ;dum = gsn_add_text(wks,plot_vector,"~F37~p",1000,stdlat,txres) 
  ;page coordinates are normalized. go from 0->1 not 0->8" etc. 第一个为左右，左小右大，第二个为上下，下小上大
  gsn_text_ndc(wks,"~F37~p",0.525,0.22,txres)  
  ;2706-0.59;2700-0.58;2600-0.43 2612-0.48 2618-54
frame(wks)
end
;2008102600,182N,881E;2008102606,187N,884E;
;2008102612,198N,888E;2008102618,210N,896E;
;2008102700,234N,904E;2008102706,250N,906E,