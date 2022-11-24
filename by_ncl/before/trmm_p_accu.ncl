load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
; 降水数据
data_dir     = "/cygdrive/D/ncl_related/data/trmm_rain/"  ;参数文件路径
picture_dir ="/cygdrive/D/ncl_related/pictures/rain/" ;图片存储路径


a = addfile(data_dir+"p_25.nc","r")
p1 = a ->HQPRECIPITATION (0,:,:)  

b = addfile(data_dir+"p_26.nc","r")
p2 = b ->HQPRECIPITATION (0,:,:)  

c = addfile(data_dir+"p_27.nc","r")
p3 = c ->HQPRECIPITATION (0,:,:)

d = addfile(data_dir+"p_28.nc","r")
p4 = d ->HQPRECIPITATION (0,:,:)  

p=p1(:,:)+p2(:,:)+p3(:,:)+p4(:,:)
copy_VarCoords(p1,p)

;print(max(p({27:32},{91:97})))
;exit


;风暴数据,用于标记风暴位置或路径
 
fiTY ="/cygdrive/D/ncl_related/data/TC/bio042008.txt"
 
; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
nrow = numAsciiRow(fiTY)
 
YYYYMMDDHH = new(nrow, "string")
lat = new(nrow, "float")
lon = new(nrow, "float")
 
data = asciiread(fiTY, -1, "string")
YYYYMMDDHH = str_get_field(data, 3, ",")
lat = stringtofloat(str_get_field(data, 7, ",")) *0.1
lon = stringtofloat(str_get_field(data, 8, ",")) *0.1

DateChar = stringtochar(YYYYMMDDHH)
MM = chartostring(DateChar(:,5:6))
DD = chartostring(DateChar(:,7:8))
HH = chartostring(DateChar(:,9:10))


wks = gsn_open_wks("png", picture_dir+"p_accu")
;gsn_define_colormap(wks,"prcp_1") 

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnLeftString = "" ;左方字符串
res@gsnRightString = ""
res@gsnMaximize=True

res@vpWidthF=0.6
res@vpHeightF=0.4

;设置坐标轴标签粗细
res@gsnMaximize=True
res@tmXBLabelFont=22
res@tmYLLabelFont=22
res@tmXBLabelFontHeightF = 0.015 ; x轴标签字号
res@tmYLLabelFontHeightF =0.015

res@pmTickMarkDisplayMode  = "Always"


;设置地图
preRes = res
preRes@mpFillOn                = False
;res@mpFillColor               ="white"
preRes@mpMinLatF               = 15.
preRes@mpMaxLatF               = 35.
preRes@mpMinLonF               = 80.
preRes@mpMaxLonF               = 105.
preRes@mpDataBaseVersion       = "MediumRes"
preRes@mpDataSetName           = "Earth..4"
preRes@mpOutlineSpecifiers     = (/"China:states"/)

preRes@mpOutlineBoundarySets       = "National"
preRes@mpNationalLineCol           = "black"
preRes@mpProvincialLineColor       = "black"
preRes@mpGeophysicalLineColor      = "black"
preRes@mpNationalLineThicknessF    = 3
preRes@mpProvincialLineThicknessF  = 3
preRes@mpGeophysicalLineThicknessF = 3

;preRes@tmXBMajorLengthF = 0.005 下方标签条刻度长度
;map = gsn_csm_map(wks,preRes)
;preRes@tiMainString ="p_2008_10_26" ;标题
;res@tiMainFontColor="red"           ;主题字体颜色
;res@gsnCenterString="center title"  ;图片正中央字符串
 
;res@tiXAxisString="X axis" ;x轴坐标名称
;res@tiYAxisString="Y axis"
; res@tiXAxisFontColor="black" ；字体颜色

;青藏高原轮廓线
;shptp="/home/YEWEI/draw/add_TP/DBATP_Line.shp"     ;修改为你的shp文件存放路径
;tpres                  = True
;tpres@gsLineColor      = "forestgreen"     ;设置廓线颜色
;tpres@gsLineThicknessF = 3              ;设置廓线宽度
;map                    = gsn_add_shapefile_polylines(wks,map,shptp,tpres)  ;map为底图


;画等值线图
preRes = res
preRes@cnFillOn = True  ;cn means ContourPlot
preRes@cnLinesOn = True   
preRes@cnMonoLineColor = False ;线条不使用同一个颜色
preRes@cnSpanLinePalette = True
preRes@cnFillPalette = "wh-bl-gr-ye-re"
;preRes@cnLinePalette = "prcp_1" ;线条的颜色
preRes@cnLineLabelsOn = False   
preRes@cnInfoLabelOn = False ; set false will not set an information label
preRes@cnLevelSelectionMode = "ManualLevels" 
;preRes@cnMinLevelValF = 5500
;preRes@cnMaxLevelValF = 5800
;preRes@cnLevelSelectionMode = "ExplicitLevels"
;preRes@cnLevels = (/10, 25, 30, 40, 50/)
;preRes@cnFillColors=(/2,3,4,5,6/)
;preRes@cnMonoLineThickness = False 所有等值线不设置成同样粗细
;res@cnLineThicknesses = (/1, 1, 3, 1, 1, 1, 1, 1, 1/)

preRes@cnLineThicknessF = 2
;preRes@cnLevelSpacingF = 10 
;preRes@gsnSpreadColors = True ; 
;preRes@gsnSpreadColorStart = 1
;preRes@gsnSpreadColorStart =17 

;设置色标条
preRes@lbBoxLinesOn = True ; Label bar
preRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
preRes@lbOrientation        = "Vertical"         ; vertical label bar 色标垂直
preRes@pmLabelBarWidthF = 0.035     ;色标变窄
;preRes@pmLabelBarHeightF = 0.02    ;色标变细
preRes@lbBoxEndCapStyle = "TriangleBothEnds"
preRes@pmLabelBarOrthogonalPosF=0.02





;某时刻风暴位置标记
;txres               = True                     ; text mods desired
;txres@txFontHeightF = 0.03                     ; font smaller. default big
;txres@txFontColor="red"
;do i = 0, nrow-1
;   if (HH(i) .eq. "00".and.DD(i) .eq. "27") then
;      dum = gsn_add_text(wks,plot,"~F37~p",lon(i), lat(i),txres) 
;    end if
;end do


;绘分布图
p = smth9_Wrap(p,0.5,0.25,False)
;p = smth9_Wrap(p,0.5,0.25,False)
plot = gsn_csm_contour_map(wks,p,preRes)   

;风暴路径绘制
resLine =True
resLine@gsLineThicknessF = 4
resLine@gsLineColor = "black"
txres = True                     ; text mods desired
txres@txFontHeightF   = 0.02                     ; font smaller. default big
txres@txFontThicknessF=4
txres@txFontColor     ="red"


;画线

dumLine = new(nrow, graphic)
do i = 0, nrow-2
    xx = (/ lon(i), lon(i+1)/)
    yy = (/ lat(i), lat(i+1)/)
    dumLine(i) = gsn_add_polyline(wks, plot, xx, yy, resLine)
end do 

;画点
do i = 0, nrow-1
    dum = gsn_add_text(wks,plot,"~F37~p",lon(i), lat(i),txres) 
end do


; 标注00时的日期

  dumDate = new(nrow,graphic)
  resTx = True
  resTx@txFontHeightF = 0.03
  resTx@txFontColor = "black"
  ;resTx@txBackgroundFillColor="white"
  resTx@txJust = "BottomLeft"
 
  do i = 1, nrow-1
    if (HH(i) .ne. "00" ) then
      continue
    end if
    dumDate(i) = gsn_add_text(wks,plot, MM(i)+DD(i), lon(i)+0.7, lat(i)-0.3, resTx)
  end do
draw(plot)
frame(wks)
end