;0305 绘制96E垂直速度，水平水汽通量散度垂直剖面
;计算等压面上的比湿 mixhum_ptrh，
;网址 http://www.ncl.ucar.edu/Document/Functions/Built-in/mixhum_ptrh.shtml
;将Pa/s转换成m/s，用于与v合成
;网址 https://www.ncl.ucar.edu/Document/Functions/Contributed/omega_to_w.shtml
;根据气压计算高度http://www.ncl.ucar.edu/Document/Functions/Built-in/stdatmus_p2tdz.shtml
;气压纬度矢量图gsn_csm_pres_hgt_vector 
;http://www.ncl.ucar.edu/Document/Graphics/Interfaces/gsn_csm_pres_hgt_vector.shtml

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;---------------------------------------绘图属性设置------------------------------------------
undef("set_all_res")
procedure set_all_res(res,etores)
 local res,etores
 begin
 
  ;--------------------------------------------------
   res = True
   res@gsnDraw = False
   res@gsnFrame = False
   res@tmXTOn = False ;关刻度
   res@tmYROn = False
   res@tmYRMode = "Automatic"  ; turn off special labels on right axis   
   res@gsnLeftString = "" ;左上角字符
   res@gsnRightString=""
   res@tiYAxisString="P/hPa"
   res@tiYAxisFont=22

  ;设置坐标轴标签粗细
    ;res@gsnMaximize=True
    res@tmXBLabelFont=22
    res@tmYLLabelFont=22
    res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
    res@tmYLLabelFontHeightF =0.02


  ;---------contour-------------

   ;cres                    = res
   res@cnFillOn           = True ;cn means ContourPlot
   res@cnLinesOn          = False   
   ;res@cnSpanFillPalette  = True
   res@cnFillPalette      = "BlueWhiteOrangeRed"
   res@cnInfoLabelOn      = False ; set false will not set an information label 
   res@cnLevelSelectionMode = "ManualLevels" 
   res@cnMinLevelValF     =-4
   res@cnMaxLevelValF     =4
   res@cnLevelSpacingF    =0.5
   ;res@cnLevelSelectionMode = "ExplicitLevels" 
   ;res@cnLevels=(/-8,-6,-4,-2,0,2,4,6,8/)
   ;res@cnFillColors=(/2,24,36,48,90,140,180,200,220,250/)
   res@lbLabelBarOn       =True
   res@lbBoxLinesOn       = True ; Label bar
   res@lbLabelAutoStride = True            ; let NCL determine label spacing
   ;res@lbOrientation     = "Vertical"         ; vertical label bar
   ;res@pmLabelBarWidthF  = 0.05     ;色标变窄
   res@pmLabelBarHeightF  = 0.05    ;色标变短
   res@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
   res@lbBoxEndCapStyle   = "TriangleBothEnds"
   res@cnFillDrawOrder    = "PreDraw"


  ;---------vector-----------------
    ;矢量叠加在contour上，属性设置合并
    res@vcGlyphStyle          = "LineArrow"
    res@vcLineArrowThicknessF = 4
    res@vcMinDistanceF        = 0.01
    res@vcRefMagnitudeF       = 15 ;设置单位量级
    res@vcRefLengthF          = 0.02
    res@vcLineArrowHeadMaxSizeF = 0.01
    res@vcLineArrowHeadMinSizeF = 0.01
    res@vcRefAnnoOn           = True ;anno代表矢量注释
    res@vcRefAnnoFontHeightF  = 0.015
    res@vcRefAnnoZone         = 1
    res@vcRefAnnoSide         = "Top"
    res@vcRefAnnoJust         = "TopRight"
    ;res@vcRefAnnoString1 = "5 g cm~S~-1~N~hPa~S~-1~Ns~S~-1~N)"
    res@vcRefAnnoString1      = "15"
    res@vcRefAnnoString2On    = False

  ;-----------地形---------------------------
   etores=True
    etores@gsnMaximize      = True
    etores@trYMaxF=1000
    etores@gsnYRefLine      = 1000.          ; create a reference line 
    etores@gsnBelowYRefLineColor = "black"              ; above ref line fill red
    etores@gsnDraw        = False             ; don't draw the plots yet
    etores@gsnFrame       = False        ; don't advance the frame
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
 
 data1_dir  ="/cygdrive/D/ncl_related/data/fnl/" ;fnl数据路径
 data2_dir  ="/cygdrive/D/ncl_related/data/other/"    ;地形数据路径
 pic_dir  ="/cygdrive/D/ncl_related/pictures/vapour/sec/" ;图片存储路径

;读地形文件----------------------------
   ;lon1       fspan(0,360,4320)  ;gsn_xy    
   lat1     =fspan(90,-90,2160)  ;用于gsn_xy 地形图
   setfileoption("nc","MissingToFillValue",False)    ;设置缺省值
   elev= read_elev_data(data2_dir+"ETOPO5.DAT")
   ;printVarSummary(elev)
   elev0 = 1013.25*(1-elev*0.0065/288.15)^5.25145        ;将地形数据的单位：m转化为hPa
   copy_VarCoords(elev,elev0) 
   ;printVarSummary(elev0)

  stdlon=96;选择固定某一角度
    elev1=elev0(:,{stdlon})       ;固定某一经度的地形数据
    elev1@_FillValue= -9.96921e+36

; 数据处理-----------------------------------------------------
  ;;;>>>手动设置一个时间 用各个数据生成一个带unit的 timeset_ymdh 并且转换为一个字符串 方便文件中读入
   year   = 2008
   month  = 10
   day    = 26
   hour   = 06
   minute = 00
   second = 0
  ;设置时间格式读取数据####################################
   units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
   opt = 0
   timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
   file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)

  ;计算等压面上的比湿
   data   = addfile(data1_dir+"fnl_"+file_time+".nc", "r") 
   rh     = data->RH_P0_L100_GLL0   ;lv_ISBL4/lat_0/lon_0，21层
   tmpdata= data->TMP_P0_L100_GLL0(5:25,:,: ) ;单位：K,lv_ISBL0/lat_0/lon_0,26层，选择对应的21层
   u      = data->UGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
   v      = data->VGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
   vvel   = data->VVEL_P0_L100_GLL0(:,:,: ) ;Pa/s,lv_ISBL4, lat_0, lon_0 
   pres1  =data->lv_ISBL4  ;Pa,21层
   lon    =data->lon_0
   lat    =data->lat_0
   pres   =pres1*0.01
   pres@units="hPa"
   q      = mixhum_ptrh(conform(tmpdata,pres,0),tmpdata,rh,-2)
   ; specific humidity (g/kg)，要求温度K,气压hPa
   copy_VarCoords(rh,q)
   q@units="g/kg"

  ;将Pa/s转换成m/s   
   w1 = omega_to_w(vvel, conform(vvel,pres1, 0), tmpdata)  ; w[*][*][*] (m/s)，要求温度K，气压Pa
   w =w1*100.
   copy_VarCoords(w1, w)
   h_lw = w(:,:,{stdlon})
   h_lw!0="level"
   h_lw&level=pres

  ;设定经向风速v经度
   h_lv = v(:,:,{stdlon})
   h_lv!0="level"
   h_lv&level=pres

  ;计算垂直水汽通量和各层水平~~~散度 
   qu = new((/21,181,360/),float)
   qv = new((/21,181,360/),float)
   qu = q*u/9.8
   qv = q*v/9.8
   ;copy_VarCoords_not_n(q,qv,(/2/)) ; 拷贝除了lon以外的坐标信息
   copy_VarCoords(q,qu)
   copy_VarCoords(q,qv)
   vapdiv = new((/21,181,360/),float)
   copy_VarCoords(q,vapdiv)
   vapdiv=uv2dv_cfd(qu,qv,lat,lon,2)             ; u,v ==> divergence
   vapdiv=vapdiv*1.e5
   copy_VarCoords(rh,vapdiv) 
   ;10^(-5) "g∙" 〖"cm" 〗^"-2"  〖"hPa" 〗^"-1"  "s" ^"-1
    h_lvapdiv=vapdiv(:,:,{stdlon})
    copy_VarCoords(tmpdata,h_lvapdiv)
    h_lvapdiv!0     ="level"
    h_lvapdiv&level =pres
 





;dir="/cygdrive/d/ncl_related/data/"
wks = gsn_open_wks("png",pic_dir+"vapdiv_w"+file_time)


;------------绘图---------------------

res    = True
 ;cres  = True
 etores=True
 set_all_res(res,etores)


;print(min(h_lvvel))
;print(max(h_lrh))
;exit

plot_vector= gsn_csm_pres_hgt_vector(wks,h_lvapdiv(:,{25:35}),h_lv(:,{25:35}),h_lw(:,{25:35}),res)
plot_hgt = gsn_csm_xy(wks,lat1,elev1,etores) 
overlay(plot_vector,plot_hgt)


draw(plot_vector)
frame(wks)
end