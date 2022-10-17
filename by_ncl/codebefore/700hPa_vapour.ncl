;1种数据，用fnl的风速和相对湿度算水汽通量和水汽通量散度
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"

begin

;手动设置一个时间,需要改两处
;用各个数据生成一个带unit的 timeset_ymdh 
;转换为一个字符串 方便文件中读入
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
   file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)

;数据读取
data_dir="/cygdrive/D/ncl_related/data/fnl/"
  pic_dir="/cygdrive/D/ncl_related/pictures/700/"
  data   = addfile(data_dir+"fnl_"+file_time+".nc", "r") 
  rh     = data->$"RH_P0_L100_GLL0"$   ;lv_ISBL4/lat_0/lon_0，21层
  tmpdata= data->TMP_P0_L100_GLL0(5:25,:,: ) ;单位：K,lv_ISBL0/lat_0/lon_0,26层，选择对应的21层
  u      = data->UGRD_P0_L100_GLL0({70000},:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
  v      = data->VGRD_P0_L100_GLL0({70000},:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
  hgtData= data->HGT_P0_L100_GLL0({70000},:,:)  ;13={50000}  hgt 的单位：gpm
  pres1  =data->lv_ISBL4  ;Pa,21层
  lon1    =data->lon_0
  lat1    =data->lat_0

;数据处理
  pres   =pres1*0.01
  pres@units="hPa"
  q      = mixhum_ptrh(conform(tmpdata,pres,0),tmpdata,rh,-2) ;每层都算了
   ; specific humidity (g/kg)，要求温度K,气压hPa
  copy_VarCoords(rh,q)
  q@units="g/kg"
  q1=q({70000},:,:)
  copy_VarCoords(u, q1) 
  qu = new((/181,360/),float)
  qv = new((/181,360/),float)
  qu= q1*u/9.8
  qv= q1*v/9.8
  copy_VarCoords(q1,qu)
  copy_VarCoords(q1,qv)
  vapdiv = new((/181,360/),float)
  copy_VarCoords(q1,vapdiv)
  ;vapdiv=uv2dv_cfd(qu,qv,lat1,lon1,2)             ; u,v ==> divergence
  vapdiv=uv2dvF_Wrap(qu(::-1,:), qv(::-1,:)) ;为什么这两个函数算出来完全相反
  ;因为用uv2dvF函数时要求qu，qv纬度从小到大排列
  vapdiv=vapdiv*1.e5
  hgt = hgtData/10
  copy_VarCoords(hgtData, hgt)

;max1=min(vapdiv({25:30},{90:100}))
;print(max1)
;exit


;风暴数据,用于标记风暴位置或路径 
fiTY = "/cygdrive/D/ncl_related/data/TC/bio042008.txt"
; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
nrow = numAsciiRow(fiTY)
YYYYMMDDHH = new(nrow, "string")
lat = new(nrow, "float")
lon = new(nrow, "float")
;以字符串形式读入
data3 = asciiread(fiTY, -1, "string")
YYYYMMDDHH = str_get_field(data3, 3, ",")
lat = stringtofloat(str_get_field(data3, 7, ",")) *0.1
lon = stringtofloat(str_get_field(data3, 8, ",")) *0.1
DateChar = stringtochar(YYYYMMDDHH)
MM = chartostring(DateChar(:,5:6))
DD = chartostring(DateChar(:,7:8))
HH = chartostring(DateChar(:,9:10))


wks = gsn_open_wks("png", pic_dir+day+hour+"test")
gsn_define_colormap(wks,"BlueWhiteOrangeRed") 


res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize=True
res@gsnLeftString = ""
res@gsnRightString = ""
res@tmYROn=False
res@tmXTOn=False
;res@pmTickMarkDisplayMode   = "Always" ;不显示次刻度线
;res@tmXBMajorLengthF = 0.005  ;主刻度线长度
;res@vpWidthF=0.6
;res@vpHeightF=0.4


;设置坐标轴标签粗细
res@tmXBLabelFont=22
res@tmYLLabelFont=22
res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
res@tmYLLabelFontHeightF =0.02




;设置地图
cRes = res
cRes@mpFillOn                = False
cRes@mpMinLatF               = 10.
cRes@mpMaxLatF               = 40.
cRes@mpMinLonF               = 70.
cRes@mpMaxLonF               = 110.
cRes@mpDataBaseVersion       = "MediumRes"
cRes@mpDataSetName           = "Earth..4"
cRes@mpOutlineSpecifiers     = (/"China:states"/)




;------water vapor flux divergence---------
cRes                    = res
cRes@cnFillOn           = True ;cn means ContourPlot
cRes@cnLinesOn          = False   
cRes@cnSpanFillPalette  = True
cRes@cnFillPalette      = "BlueWhiteOrangeRed"
cRes@cnLineLabelsOn     = False  
cRes@cnInfoLabelOn      = False ; set false will not set an information label 
cRes@cnLevelSelectionMode = "ManualLevels" 
cRes@cnMinLevelValF     =-5
cRes@cnMaxLevelValF     = 5
cRes@cnLevelSpacingF    = 1
;cRes@cnLineThicknessF  = 4
;res@cnLineLabelInterval  = 1 ；all labeled 
cRes@lbBoxLinesOn       = True; Label bar
cRes@lbLabelAutoStride  = True            ; let NCL determine label spacing
cRes@lbOrientation      = "Vertical"         ; vertical label bar
cRes@pmLabelBarWidthF   = 0.05     ;色标变窄
;cRes@pmLabelBarHeightF = 0.5    ;色标变短
cRes@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
cRes@lbBoxEndCapStyle   = "TriangleBothEnds"

;--------------水汽通量矢量--------------


vRes = res
;vRes@vcGlyphStyle            = "CurlyVector"   ; turn on curly vectors
vRes@vcGlyphStyle    = "LineArrow"
vRes@vcLineArrowThicknessF = 3
vRes@vcMinDistanceF = 0.05
vRes@vcLineArrowHeadMaxSizeF = 0.01
vRes@vcLineArrowHeadMinSizeF = 0.01
vRes@vcRefMagnitudeF = 10.
vRes@vcRefLengthF    = 0.02
vRes@vcRefAnnoOn = True
vRes@vcRefAnnoFontHeightF =0.015
vRes@vcRefAnnoZone = 1.5
vRes@vcRefAnnoSide = "Top"
vRes@vcRefAnnoJust = "TopRight" ;默认值
;vRes@vcRefAnnoString1 = "10 gcm~S~-1~N~hPa~S~-1~N~s~S~-1~N~"
vRes@vcRefAnnoString1 = "10"
vRes@vcRefAnnoString2On = False
;vRes@vcRefAnnoParallelPosF=0.33
;vRes@vcRefAnnoOrthogonalPosF=-0.1 ;调整上下位置
;vRes@vcRefAnnoBackgroundColor=-1 ;将矢量注释的填充色设为透明
;vRes@vcRefAnnoPerimOn=False ;将矢量注释的边框去掉


;--------------叠加高度场--------------

heightRes = res
heightRes@cnFillOn = False  ;cn means ContourPlot
heightRes@cnLinesOn = True   
heightRes@cnLineColor="blue"
heightRes@cnLineLabelsOn = True  
heightRes@cnInfoLabelOn = False ; set false will not set an information label
heightRes@cnLineLabelBackgroundColor=-1
heightRes@cnLabelMasking = True
heightRes@cnLevelSelectionMode = "ManualLevels" 
heightRes@cnLevelSpacingF = 1
heightRes@cnLineThicknessF = 4


;平滑
vapdiv = smth9_Wrap(vapdiv,0.5,0.25,False)
vapdiv = smth9_Wrap(vapdiv,0.5,0.25,False)
hgt = smth9_Wrap(hgt,0.5,0.25,False)
hgt = smth9_Wrap(hgt,0.5,0.25,False)

;contour on contour
;plot = gsn_csm_contour_map_overlay(wks,vf,hgt,cRes, heightRes)

plot=gsn_csm_contour_map(wks, vapdiv, cRes)
vectorplot= gsn_csm_vector(wks, qu , qv , vRes)
heightplot=gsn_csm_contour(wks,hgt, heightRes)
;mapplot=gsn_csm_map(wks, mpres)


;叠加矢量和等值线
overlay(plot,vectorplot)
overlay(plot,heightplot)
;overlay(mapplot,vectorplot)

;某时刻风暴位置标记
txres               = True                     ; text mods desired
txres@txFontHeightF = 0.03                     ; font smaller. default big
txres@txFontColor="red"
do i = 0, nrow-1
   if ((HH(i) .eq."12").and.(DD(i) .eq."26")) then
      dum = gsn_add_text(wks,plot,"~F37~p",lon(i), lat(i),txres) 
    end if
end do
draw(plot)
frame(wks)
end