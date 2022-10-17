;过26日00时、27日00时的台风中心纬度剖面：垂直速度矢量，相当位温等值线，散度填色
;换时次时，需改动148,186，最后标记的位置
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;---------------------------------------绘图属性设置------------------------------------------
undef("set_all_res")
procedure set_all_res(res,etores,cnres,eqtres)
 local res,etores,cnres,eqtres
 begin
 
  ;--------------------------------------------------
   res = True
   res@gsnDraw       = False
   res@gsnFrame      = False
   res@gsnMaximize   = True
   res@tmXTOn        = False ;关刻度
   res@tmYROn        = False
   res@gsnLeftString = "" ;左上角字符
   res@gsnRightString=""
   res@tiYAxisString ="P/hPa"
   res@tiYAxisFont   =22
   res@tmYRMode      = "Automatic" ;关掉右边的height（km）字样
   res@gsnMajorLonSpacing    =5. 
   res@gsnMinorLonSpacing    =1.    
   res@gsnMajorLatSpacing    =5. 
   res@gsnMinorLatSpacing    =1



  ;设置坐标轴标签粗细
    res@tmXBLabelFont=22
    res@tmYLLabelFont=22
    res@tmXBLabelFontHeightF = 0.025 ; x轴标签字号
    res@tmYLLabelFontHeightF =0.025



  ;--------------------散度---------------------------
   cnres=res
   cnres@cnFillOn          = True                  ; turn on color
   cnres@cnInfoLabelOn     = False ; set false will not set an information label
   cnres@cnFillPalette     = "NEO_div_vegetation_c"            ; set color map
   cnres@cnLinesOn         = False                 ; no contour lines
   cnres@cnLevelSelectionMode = "ManualLevels" 
   cnres@cnMinLevelValF       = -8             ; set min contour level
   cnres@cnMaxLevelValF       = 8               ; set max contour level
   cnres@cnLevelSpacingF      = 2              ; set contour spacing
   cnres@lbBoxLinesOn         = True ; Label bar
   cnres@lbLabelAutoStride    = True            ; let NCL determine label spacing
   ;cnres@lbOrientation       = "vertical"            ; 垂直坐标
   cnres@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
   cnres@lbBoxEndCapStyle     = "TriangleBothEnds"
   ;cnres@pmLabelBarWidthF    = 0.05     ;色标变窄
   cnres@pmLabelBarHeightF    =0.05
   cnres@cnFillDrawOrder      = "PreDraw"

  ;---------vector-----------------
    ;垂直速度矢量叠加在contour上，属性设置合并
    cnres@vcGlyphStyle          = "LineArrow"
    cnres@vcLineArrowThicknessF = 4
    cnres@vcMinDistanceF        = 0.02
    cnres@vcRefMagnitudeF       = 15 ;设置单位量级
    cnres@vcRefLengthF          = 0.02
    cnres@vcLineArrowHeadMaxSizeF = 0.01
    cnres@vcLineArrowHeadMinSizeF = 0.01
    cnres@vcRefAnnoOn           = True ;anno代表矢量注释
    cnres@vcRefAnnoFontHeightF  = 0.015
    cnres@vcRefAnnoZone         = 1
    cnres@vcRefAnnoSide         = "Top"
    cnres@vcRefAnnoJust         = "TopRight"
    ;cnres@vcRefAnnoString1 = "5 g cm~S~-1~N~hPa~S~-1~Ns~S~-1~N)"
    cnres@vcRefAnnoString1      = "15"
    cnres@vcRefAnnoString2On    = False


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


  ;-----------地形---------------------------
   etores=True
    etores@gsnMaximize      = True
    etores@gsnDraw          = False           ; don't draw the plots yet
    etores@gsnFrame         = False           ; don't advance the frame
    etores@trYMaxF          =1000
    etores@gsnYRefLine      =1000.          ; create a reference line 
    etores@gsnBelowYRefLineColor = "black"  ; above ref line fill red
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


begin

;----------------------------------------------读地形文件----------------------------
   lon1       =fspan(0,360,4320)  ;gsn_xy    
   ;lat1     =fspan(90,-90,2160)  ;用于gsn_xy 地形图
   ;datadir="/cygdrive/D/ncl_related/data/"
   datadir  ="/cygdrive/D/ncl_related/data/other/"   ;地形数据路径
   setfileoption("nc","MissingToFillValue",False)    ;设置缺省值
   elev= read_elev_data(datadir+"ETOPO5.DAT")
   ;printVarSummary(elev)
   elev0 = 1013.25*(1-elev*0.0065/288.15)^5.25145        ;将地形数据的单位：m转化为hPa
   copy_VarCoords(elev,elev0) 
   ;printVarSummary(elev0)

;----------------------------------------------选择固定某一角度--------------------------------
  ;stdlon=88  
  stdlat=23   
  ;2008102600,182N,881E;2008102606,187N,884E;
  ;2008102612,198N,888E;2008102618,210N,896E;
  ;2008102700,234N,904E;2008102706,250N,906E,
    elev1=elev0({stdlat},:)       ;固定某一经度的地形数据    
    elev1@_FillValue= -9.96921e+36

;---------------------------------------------读变量-----------------------------------------------------
  ;手动设置一个时间 用各个数据生成一个带unit的 timeset_ymdh 并且转换为一个字符串 方便文件中读入
   year   = 2008
   month  = 10
   day    = 27
   hour   = 00
   minute = 00
   second = 0
  ;设置时间格式读取数据####################################
   units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
   opt = 0
   timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
   file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)


   path_rawdata= "/cygdrive/D/ncl_related/data/fnl/"  ; 变量数据路径
   data   = addfile(path_rawdata+"fnl_"+file_time+".nc", "r") 
   rh     = data->$"RH_P0_L100_GLL0"$   ;lv_ISBL4/lat_0/lon_0,21层
   tmpdata= data->TMP_P0_L100_GLL0(5:25,:,: ) ;单位：K,lv_ISBL0/lat_0/lon_0,26层，选择对应的21层
   u      = data->UGRD_P0_L100_GLL0(5:25,:,:) 
   v      = data->VGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
   vvel   = data->VVEL_P0_L100_GLL0(:,:,: ) ;Pa/s,lv_ISBL4, lat_0, lon_0 
   pres1  =data->lv_ISBL4  ;Pa
   lon    =data->lon_0
   lat    =data->lat_0
   pres   =pres1/100
   pres@units="hPa"

  ;-------------------------计算各层散度----------------------------------------------------------
  div=uv2dvF_Wrap(u(:,::-1,:),v(:,::-1,:)) ;用这个球面间谐波函数计算散度更精确
   ;di = uv2dv_cfd (u,v,lat,lon, 3) 
   div=div*1.e5
   h_ldiv=div(:,{stdlat},:)

  ;----------------------------------设定垂直速度矢量-------------------------------------------------
  ;设定纬向风速，将Pa/s转换成m/s   
   w1 = omega_to_w(vvel, conform(vvel,pres1, 0), tmpdata)  ; w[*][*][*] (m/s)，要求温度K，气压Pa
   w =w1*100.
   copy_VarCoords(w1, w)
   h_lw = w(:,{stdlat},:)
   h_lw!0="level" ;有这步没这步画出来的图都不一样
   h_lw&level=pres
  ;设定经向风速v经度
   h_lu = u(:,{stdlat},:)
   h_lu!0="level"
   h_lu&level=pres

 ;----------------------------计算相当位温---------------------------------------------
 eqt= pot_temp_equiv(pres1,tmpdata, rh, 0, "rh")  ;可直接使用这个函数计算,用Pa来算
 copy_VarCoords(tmpdata, eqt)
  ;printVarSummary(eqt)
  ;print(pres)
  ;exit 
  h_leqt  =eqt(:,{stdlat},:)
  ;h_leqt = dim_avg_n_Wrap(eqt(:,:,{92:98}), 2)
  h_leqt!0="level"
  h_leqt&level=pres



pic_dir ="/cygdrive/D/ncl_related/pictures/storm/sec_eqt_div/lat/"
wks = gsn_open_wks("png",pic_dir+stdlat+"N"+file_time) 


;------------绘图---------------------

res     = True
 etores = True
 cnres  = True
 eqtres = True
 set_all_res(res,etores,cnres,eqtres)


;print(min(h_lvvel))
;print(max(h_lrh))
;exit

strlon=80
endlon=100
plot_vector= gsn_csm_pres_hgt_vector(wks,h_ldiv(:,{strlon:endlon}),\
                                     h_lu(:,{strlon:endlon}),h_lw(:,{strlon:endlon}),cnres)
plot_eqt  = gsn_csm_pres_hgt(wks, h_leqt(:,{strlon:endlon}), eqtres ) 
plot_hgt = gsn_csm_xy(wks,lon1,elev1,etores) 
overlay(plot_vector,plot_eqt)
overlay(plot_vector,plot_hgt)
draw(plot_vector)

;某时刻风暴位置标记
 txres               = True                     ; text mods desired
 txres@txFontHeightF = 0.03                     ; font smaller. default big
 txres@txFontThicknessF=6
 txres@txFontColor   ="red"
 ;dum = gsn_add_text(wks,plot_vector,"~F37~p",1000,stdlat,txres) 
;page coordinates are normalized. go from 0->1 not 0->8" etc. 第一个为左右，左小右大，第二个为上下，下小上大
gsn_text_ndc(wks,"~F37~p",0.57,0.20,txres)  
;2700-0.57;2600-0.49
frame(wks)
end