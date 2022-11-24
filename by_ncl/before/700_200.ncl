load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"

;改时间，第10,21，63,189行

begin

data = addfile("/cygdrive/D/ncl_related/data/fnl/fnl_20081027_06_00.nc", "r")
data1= addfile("/home/YEWEI/listen/ERA5_pl.nc", "r")
pic_dir="/home/YEWEI/listen/"

time = data1->time ; time是一个变量，units :hours since 1900-01-01 00:00:00.0
time1= cd_calendar(time, -3)
;print(time1)
;exit
;输出(0)-2008102600...(11)-2008102818

u = data ->UGRD_P0_L100_GLL0({70000},:,:)  ;单位：m/s
v = data->VGRD_P0_L100_GLL0({70000},:,:)
hgtData = data->HGT_P0_L100_GLL0({70000},:,:)  ;13={50000}  hgt 的单位：gpm

;time=hgtData@initial_time ;变量的属性initial_time :   10/27/2008 (00:00)
;time1=cd_calendar(time, -3)
;print(time1)
;exit


t=5
di=short2flt(data1->d(t,1,:,: )) ; time(0-11)/level:0-200hPa,1-700hPa/latitude/longitude

di=di*1.e5
hgt = hgtData/10
copy_VarCoords(hgtData, hgt)

;max1=min(di({20:40},{80:100}))
;print(max1)
;exit


;m=sqrt(u^2+v^2)
;copy_VarCoords(u, m)
;max1=max(m({20:40},{80:100}))
;print(max1)
;exit



;风暴数据,用于标记风暴位置或路径
 
fiTY ="/cygdrive/D/ncl_related/data/TC/bio042008.txt"

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


wks = gsn_open_wks("png", pic_dir+"700_200_2706")
gsn_define_colormap(wks,"NEO_div_vegetation_c") 


res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnLeftString = ""
res@gsnRightString = ""


;设置坐标轴标签粗细
;res@gsnMaximize=True
res@tmXBLabelFont=22
res@tmYLLabelFont=22
res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
res@tmYLLabelFontHeightF =0.02



;设置地图
mpres = res
mpres@mpFillOn                = False
mpres@mpMinLatF               = 0.
mpres@mpMaxLatF               = 40.
mpres@mpMinLonF               = 60.
mpres@mpMaxLonF               = 120.
mpres@mpDataBaseVersion       = "MediumRes"
mpres@mpDataSetName           = "Earth..4"
mpres@mpOutlineSpecifiers     = (/"China:states"/)
mpres@pmTickMarkDisplayMode   = "Always"
mpres@tmXBMajorLengthF = 0.005
;mpres@tiMainString = "500hPa_2900"
mapplot = gsn_csm_map(wks, mpres)

;--------------矢量图--------------

windRes = res
;windRes@vcGlyphStyle            = "FillArrow"   ; turn on curly vectors
windRes@vcGlyphStyle    = "LineArrow"
windRes@vcLineArrowThicknessF = 3
windRes@vcMinDistanceF = 0.05
windRes@vcRefMagnitudeF = 10
windRes@vcRefLengthF    = 0.02
windRes@vcRefAnnoOn = True
windRes@vcRefAnnoFontHeightF =0.015
windRes@vcRefAnnoZone = 1
windRes@vcRefAnnoSide = "Top"
windRes@vcRefAnnoJust = "TopRight"
windRes@vcRefAnnoString1 = "10 m/s"
windRes@vcRefAnnoString2On = False
windRes@vcLineArrowHeadMaxSizeF = 0.01
windRes@vcLineArrowHeadMinSizeF = 0.01



;--------------散度--------------

cRes = res
cRes@cnFillOn = True ;cn means ContourPlot
cRes@cnLinesOn = False   
;cRes@cnMonoLineColor = False
;cRes@cnLineColor="blue"
cRes@cnSpanFillPalette = True
cRes@cnFillPalette = "NEO_div_vegetation_c"
cRes@cnLineLabelsOn = False  
cRes@cnInfoLabelOn = False ; set false will not set an information label 
cRes@cnLevelSelectionMode = "ManualLevels" 
;cRes@cnMaxLevelValF=0
;cRes@cnLevelSpacingF = 4 
cRes@cnLineThicknessF = 4
cRes@lbBoxLinesOn = False ; Label bar
cRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
cRes@lbOrientation        = "Vertical"         ; vertical label bar
cRes@pmLabelBarWidthF = 0.05     ;色标变窄
;cRes@pmLabelBarHeightF = 0.02    ;色标变短
cRes@lbBoxEndCapStyle = "TriangleBothEnds"


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
;heightRes@cnMinLevelValF = 552
;heightRes@cnMaxLevelValF = 588
heightRes@cnLevelSpacingF = 1
heightRes@cnLineThicknessF = 4
heightRes@lbBoxLinesOn = False ; Label bar
heightRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
heightRes@lbOrientation        = "Vertical"         ; vertical label bar


;绘图
di = smth9_Wrap(di,0.5,0.25,False)
di = smth9_Wrap(di,0.5,0.25,False)


hgt = smth9_Wrap(hgt,0.5,0.25,False)
hgt = smth9_Wrap(hgt,0.5,0.25,False)



contour_div = gsn_csm_contour(wks,di,cRes)
windPlot = gsn_csm_vector(wks, u , v , windRes)
contour_height = gsn_csm_contour(wks,hgt, heightRes)

;叠加风矢量和位势高度
overlay(mapplot, contour_div)
overlay(mapplot, windPlot)
overlay(mapplot, contour_height)

;某时刻风暴位置标记
txres               = True                     ; text mods desired
txres@txFontHeightF = 0.03                     ; font smaller. default big
txres@txFontColor="red"
do i = 0, nrow-1
   if (HH(i) .eq. "06".and.DD(i) .eq. "27") then
      dum = gsn_add_text(wks,mapplot,"~F37~p",lon(i), lat(i),txres) 
    end if
end do
draw(mapplot)
frame(wks)
end