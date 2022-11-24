load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"
begin

;手动设置一个时间 
;用各个数据生成一个带unit的 timeset_ymdh 
;转换为一个字符串 方便文件中读入
   year   = 2008
   month  = 10
   day    = 28
   hour   = 12
   minute = 00
   second = 0
  ;设置时间格式读取数据####################################
   units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
   opt = 0
   timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
   file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)

data_dir="/cygdrive/D/ncl_related/data/fnl/"   
pic_dir ="/cygdrive/D/ncl_related/pictures/500/"
data    = addfile(data_dir+"fnl_"+file_time+".nc", "r")

uData = data ->UGRD_P0_L100_GLL0(13,:,:)  ;单位：m/s
vData = data->VGRD_P0_L100_GLL0(13,:,:)
hgtData = data->HGT_P0_L100_GLL0(13,:,:)  ;13={50000}  hgt 的单位：gpm
tmpData=data->TMP_P0_L100_GLL0 (13,:,: )  ;单位：K

hgt = hgtData/10
copy_VarCoords(hgtData, hgt)
tmp=tmpData-273.15
copy_VarCoords(tmpData, tmp)


;风暴数据,用于标记风暴位置或路径
 
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


wks = gsn_open_wks("png", pic_dir+file_time)
gsn_define_colormap(wks,"rainbow") 


res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize=True
res@gsnLeftString = ""
res@gsnRightString = ""
res@tmYROn=False
res@tmXTOn=False
;res@vpHeightF = 1
;res@vpwidthF = 1
;mpres@pmTickMarkDisplayMode  = "Always"
;mpres@tmXBMajorLengthF = 0.005
;mpres@tiMainString = "500hPa_2900"
;mpres@gsnDraw = False

;设置坐标轴标签粗细
res@tmXBLabelFont=22
res@tmYLLabelFont=22
res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
res@tmYLLabelFontHeightF =0.02


;设置地图
mpres = res
mpres@mpFillOn                = False
mpres@mpMinLatF               = 0.
mpres@mpMaxLatF               = 60.
mpres@mpMinLonF               = 50.
mpres@mpMaxLonF               = 140.
mpres@gsnMajorLonSpacing    =10. 
mpres@gsnMinorLonSpacing    =2.    
mpres@gsnMajorLatSpacing    =10. 
mpres@gsnMinorLatSpacing    =2
mpres@mpDataBaseVersion       = "MediumRes"
mpres@mpDataSetName           = "Earth..4"
mpres@mpOutlineSpecifiers     = (/"China:states"/)

mapplot = gsn_csm_map(wks, mpres)

;--------------风矢量图--------------

windRes = res
;windRes@vcGlyphStyle            = "FillArrow"   ; turn on curly vectors
windRes@vcGlyphStyle    = "LineArrow"
windRes@vcLineArrowThicknessF = 3
windRes@vcMinDistanceF = 0.03
windRes@vcRefMagnitudeF = 15.
windRes@vcRefLengthF    = 0.02
windRes@vcRefAnnoOn = True
windRes@vcRefAnnoFontHeightF =0.015
windRes@vcRefAnnoZone = 1
windRes@vcRefAnnoSide = "Top"
windRes@vcRefAnnoJust = "TopRight"
windRes@vcRefAnnoString1 = "15 m/s"
windRes@vcRefAnnoString2On = False
windRes@vcLineArrowHeadMaxSizeF = 0.01
windRes@vcLineArrowHeadMinSizeF = 0.01



;--------------叠加高度场--------------

heightRes = res
heightRes@cnFillOn = False  ;cn means ContourPlot
heightRes@cnLinesOn = True   
;heightRes@cnMonoLineColor = False
heightRes@cnLineColor="blue"
;heightRes@cnSpanLinePalette = True
;heightRes@cnLinePalette = "rainbow"



heightRes@cnLineLabelsOn = True  
heightRes@cnInfoLabelOn = False ; set false will not set an information label
heightRes@cnLevelSelectionMode = "ManualLevels" 
;heightRes@cnMinLevelValF = 552
;heightRes@cnMaxLevelValF = 588
heightRes@cnLevelSpacingF = 1 
heightRes@cnLineThicknessF = 4
;heightRes@gsnSpreadColorStart = 1
;heightRes@gsnSpreadColorStart =17
;heightRes@gsnSpreadColors = True ; 6.1.0
heightRes@lbBoxLinesOn = False ; Label bar
heightRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
heightRes@lbOrientation        = "Vertical"         ; vertical label bar


;--------------叠加温度场--------------

tmRes = res
tmRes@cnFillOn = False  ;cn means ContourPlot
tmRes@cnLinesOn = True   
tmRes@cnLineColor="red"
tmRes@cnLineDashPattern="1"
;tmRes@cnMonoLineColor = False
;tmRes@cnSpanLinePalette = True
;tmRes@cnLinePalette = "rainbow"



tmRes@cnLineLabelsOn = False  
tmRes@cnInfoLabelOn = False ; set false will not set an information label
tmRes@cnLevelSelectionMode = "ManualLevels" 
tmRes@cnLevelSpacingF = 4 
tmRes@cnLineThicknessF = 3
tmRes@lbBoxLinesOn = False ; Label bar
tmRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
tmRes@lbOrientation        = "Vertical"         ; vertical label bar


;绘图
hgt = smth9_Wrap(hgt,0.5,0.25,False)
hgt = smth9_Wrap(hgt,0.5,0.25,False)
tmp = smth9_Wrap(tmp,0.5,0.25,False)
tmp = smth9_Wrap(tmp,0.5,0.25,False)



contour_height = gsn_csm_contour(wks,hgt, heightRes)
contour_tmp = gsn_csm_contour(wks,tmp, tmRes)
windPlot = gsn_csm_vector(wks, uData , vData , windRes)


;叠加风矢量和位势高度
overlay(mapplot, contour_height)
overlay(mapplot, contour_tmp)
overlay(mapplot, windPlot)

;某时刻风暴位置标记
txres               = True                     ; text mods desired
txres@txFontHeightF = 0.03                     ; font smaller. default big
txres@txFontColor="red"
do i = 0, nrow-1
   if (HH(i) .eq. "00".and.DD(i) .eq. "28") then
      dum = gsn_add_text(wks,mapplot,"~F37~p",lon(i), lat(i),txres) 
    end if
end do

;高原轮廓线
shp="/cygdrive/D/ncl_related/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres) 

draw(mapplot)
frame(wks)
end