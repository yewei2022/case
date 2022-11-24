;0.1*0.1插值
;26日00时、27日00时的台风中心剖面：经向垂直速度矢量，相当位温等值线
;换时次时，需改动151，最后风暴标记位置两处
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;---------------------------------------绘图属性设置------------------------------------------
undef("set_all_res")
procedure set_all_res(res,etores,eqtres,rhres)
 local res,etores,eqtres,rhres
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
   ;res@gsnMajorLonSpacing    =5. 
   ;res@gsnMinorLonSpacing    =1.    
   ;res@gsnMajorLatSpacing    =5. 
   ;res@gsnMinorLatSpacing    =1



  ;设置坐标轴
    res@tmXBMode       ="Explicit"
    res@tmXBValues     =(/16,20,24,28,32/)
    res@tmXBLabels     =(/"16N","20N","24N","28N","32N"/)
    res@tmXBMinorValues=(/17,18,19,21,22,23,25,26,27,29,30,31/)
    res@tmXBLabelFont=22
    res@tmYLLabelFont=22
    res@tmXBLabelFontHeightF =0.025 ; x轴标签字号
    res@tmYLLabelFontHeightF =0.025

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
      ;---------vector-----------------
    ;垂直速度矢量叠加在contour上，属性设置合并
    eqtres@vcMapDirection        = False ;绘制垂直剖面图要设置为False
    eqtres@vcGlyphStyle          = "LineArrow"
    eqtres@vcLineArrowThicknessF = 4
    eqtres@vcMinDistanceF        = 0.02
    eqtres@vcRefMagnitudeF       = 15 ;设置单位量级
    eqtres@vcRefLengthF          = 0.02
    eqtres@vcLineArrowHeadMaxSizeF = 0.01
    eqtres@vcLineArrowHeadMinSizeF = 0.01
    eqtres@vcRefAnnoOn           = True ;anno代表矢量注释
    eqtres@vcRefAnnoFontHeightF  = 0.015
    eqtres@vcRefAnnoZone         = 1
    eqtres@vcRefAnnoSide         = "Top"
    eqtres@vcRefAnnoJust         = "TopRight"
    ;eqtres@vcRefAnnoString1 = "5 g cm~S~-1~N~hPa~S~-1~Ns~S~-1~N)"
    eqtres@vcRefAnnoString1      = "15"
    eqtres@vcRefAnnoString2On    = False

  ;---------相对湿度-------------
   rhres=res
   rhres@cnFillOn          = True                  ; turn on color
   rhres@cnLinesOn         = False
   rhres@cnFillPalette     = "WhiteBlue"             ; set color map
   rhres@cnInfoLabelOn     =False
   rhres@cnLevelSelectionMode = "ExplicitLevels" 
   rhres@cnLevels          =(/90/)
   rhres@cnFillColors      =(/-1,80/)
   rhres@cnLineLabelsOn    =False
   rhres@lbLabelBarOn      =False
   ;rhres@cnLevelSelectionMode ="ManualLevels"
   rhres@cnFillDrawOrder   ="PreDraw"


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
   ;lon1       fspan(0,360,4320)  ;gsn_xy    
   lat1     =fspan(90,-90,2160)  ;用于gsn_xy 地形图
   datadir  ="/cygdrive/D/ncl_related/data/other/"   ;地形数据路径
   setfileoption("nc","MissingToFillValue",False)    ;设置缺省值
   elev= read_elev_data(datadir+"ETOPO5.DAT")
   ;printVarSummary(elev)
   elev0 = 1013.25*(1-elev*0.0065/288.15)^5.25145        ;将地形数据的单位：m转化为hPa
   copy_VarCoords(elev,elev0) 
   ;printVarSummary(elev0)


;---------------------------------------------读变量-----------------------------------------------------
  ;手动设置一个时间 用各个数据生成一个带unit的 timeset_ymdh 并且转换为一个字符串 方便文件中读入
   year   = 2008
   month  = 10
   day    = 26
   hour   = 00
   minute = 00
   second = 0
  ;设置时间格式读取数据####################################
   units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
   opt = 0
   timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
   file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)
   path_format  = "%Y%N%D%H"   ; 台风路径文件中的时间形式     
   path_time    = cd_string(timeset_ymdh,path_format)


   path_rawdata= "/cygdrive/D/ncl_related/data/fnl/"  ; 变量数据路径
   data       = addfile(path_rawdata+"fnl_"+file_time+".nc", "r") 
   rh_in      = data->RH_P0_L100_GLL0   ;lv_ISBL4/lat_0/lon_0,21层
   tmpdata_in = data->TMP_P0_L100_GLL0(5:25,:,: ) ;单位：K,lv_ISBL0/lat_0/lon_0,26层，选择对应的21层
   u_in       = data->UGRD_P0_L100_GLL0(5:25,:,:) 
   v_in       = data->VGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
   vvel_in    = data->VVEL_P0_L100_GLL0(:,:,: ) ;Pa/s,lv_ISBL4, lat_0, lon_0 
   pres1      =data->lv_ISBL4  ;Pa
   lon_in     =data->lon_0
   lat_in     =data->lat_0
   pres       =pres1/100
   pres@units ="hPa"

;--------------------读取风暴中心位置-------------------------------------------------   
fiTY = "/cygdrive/D/ncl_related/data/TC/bio042008.txt"
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

tctime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
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
;----------------------------------------------选择固定某一角度-------------------------------- 
;2008102600,182N,881E;2008102606,187N,884E;
;2008102612,198N,888E;2008102618,210N,896E;
;2008102700,234N,904E;2008102706,250N,906E,
elev1=elev0(:,{stdlon})       ;固定某一经度的地形数据    
elev1@_FillValue= -9.96921e+36


;---------------------------网格点插值--------------------
lon_out=fspan(70,110,401)
lat_out=fspan(0, 50, 501)
   u      =linint2_Wrap(lon_in, lat_in(::-1), u_in(:,::-1,:), True, lon_out, lat_out, 0)
   delete(u_in)
   v      =linint2_Wrap(lon_in, lat_in(::-1), v_in(:,::-1,:), True, lon_out, lat_out, 0)
   delete(v_in)
   vvel   =linint2_Wrap(lon_in, lat_in(::-1), vvel_in(:,::-1,:), True, lon_out, lat_out, 0)
   delete(vvel_in)
   rh     =linint2_Wrap(lon_in, lat_in(::-1), rh_in(:,::-1,:), True, lon_out, lat_out, 0)
   delete(rh_in)
   tmpdata=linint2_Wrap(lon_in, lat_in(::-1), tmpdata_in(:,::-1,:), True, lon_out, lat_out, 0)
   delete(tmpdata_in)
   

  ;----------------------------------设定垂直速度矢量-------------------------------------------------
  ;设定纬向风速，将Pa/s转换成m/s   
   w1     = omega_to_w(vvel, conform(vvel,pres1, 0), tmpdata)  ; w[*][*][*] (m/s)，要求温度K，气压Pa
   w      =w1*100.
   copy_VarCoords(w1, w)
   h_lw   = w(:,:,{stdlon})
   h_lw!0 ="level" ;有这步没这步画出来的图都不一样
   h_lw&level=pres
  ;设定经向风速v经度
   h_lv   = v(:,:,{stdlon})
   h_lv!0 ="level"
   h_lv&level=pres

 ;----------------------------计算相当位温---------------------------------------------
 eqt         = pot_temp_equiv(pres1,tmpdata, rh, 0, "rh")  ;可直接使用这个函数计算,用Pa来算
 copy_VarCoords(tmpdata, eqt)
 h_leqt      =eqt(:,:,{stdlon})
 ;h_leqt     =dim_avg_n_Wrap(eqt(:,:,{92:98}), 2)
 h_leqt!0    ="level"
 h_leqt&level=pres

 ;--------------------设定相对湿度--------------------------------
 h_lrh  =rh(:,:,{stdlon})
 h_lrh!0="level"
 h_lrh&level=pres



pic_dir ="/cygdrive/D/ncl_related/pictures/storm/sec_eqt/lon/"
wks = gsn_open_wks("png",pic_dir+stdlon+"E"+"_eqt_"+file_time) 

;------------绘图---------------------

res     = True
 etores = True
 eqtres = True
 rhres  = True
 set_all_res(res,etores,eqtres,rhres)


;print(min(h_lvvel))
;print(max(h_lrh))
;exit

strlat=16
endlat=32
plot_vector= gsn_csm_pres_hgt_vector(wks,h_leqt({150:1000},{strlat:endlat}),\
                                     h_lv({150:1000},{strlat:endlat}),h_lw({150:1000},{strlat:endlat}),eqtres)
; plot_rh   = gsn_csm_pres_hgt(wks, h_lrh({150:1000},{strlat:endlat}), rhres ) 
plot_hgt = gsn_csm_xy(wks,lat1,elev1,etores) 
; overlay(plot_vector,plot_rh)
overlay(plot_vector,plot_hgt)
draw(plot_vector)

;某时刻风暴位置标记
 txres               = True                     ; text mods desired
 txres@txFontHeightF = 0.03                     ; font smaller. default big
 txres@txFontThicknessF=6
 txres@txFontColor   ="red"
 ;dum = gsn_add_text(wks,plot_vector,"~F37~p",1000,stdlat,txres) 
;page coordinates are normalized. go from 0->1 not 0->8" etc. 第一个为左右，左小右大，第二个为上下，下小上大
gsn_text_ndc(wks,"~F37~p",0.31,0.185,txres)  
;2700-0.545;2600-0.31 2612-0.38
frame(wks)
end