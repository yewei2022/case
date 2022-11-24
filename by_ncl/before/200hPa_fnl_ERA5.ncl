;;fnl+ERA5资料，效果更好
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"
begin
;改时间要改3处，10行，17行，197行
data_dir="/cygdrive/D/ncl_related/data/fnl/"
pic_dir="/home/YEWEI/listen/"
fil_time="fnl_20081027_00_00"
data = addfile(data_dir+fil_time+".nc", "r")
data1= addfile("/home/YEWEI/listen/ERA5_pl.nc", "r")

time = data1->time
time1= cd_calendar(time, -3);(0)-2008102600...(11)-2008102818

t=4
u = data ->UGRD_P0_L100_GLL0({20000},:,:)  ;单位：m/s
v = data->VGRD_P0_L100_GLL0({20000},:,:)
hgtData = data->HGT_P0_L100_GLL0({20000},:,:)  ;13={50000}  hgt 的单位：gpm
di=short2flt(data1->d(t,0,:,: )) ; time(0-11)/level:0-200hPa,1-700hPa/latitude/longitude
di=di*1.e5

hgt = hgtData/10
copy_VarCoords(hgtData, hgt)

;max1=max(di({0:40},{80:100}))
;print(max1)
;exit


;m=sqrt(u^2+v^2)
;copy_VarCoords(u, m)
;max1=max(m({20:40},{80:100}))
;copy_VarCoords(u, max1)
;print(max1)
;exit



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


wks = gsn_open_wks("png", pic_dir+"200hPa_2700")
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
mpres@mpMaxLatF               = 60.
mpres@mpMinLonF               = 50.
mpres@mpMaxLonF               = 140.
mpres@mpDataBaseVersion       = "MediumRes"
mpres@mpDataSetName           = "Earth..4"
mpres@mpOutlineSpecifiers     = (/"China:states"/)
mpres@pmTickMarkDisplayMode  = "Always"
mpres@tmXBMajorLengthF = 0.005
;mpres@tiMainString = "500hPa_2900"
mapplot = gsn_csm_map(wks, mpres)

;--------------矢量图--------------

windRes = res
;windRes@vcGlyphStyle            = "FillArrow"   ; turn on curly vectors
windRes@vcGlyphStyle    = "LineArrow"
windRes@vcLineArrowThicknessF = 3
windRes@vcMinDistanceF = 0.05
windRes@vcRefMagnitudeF = 20.
windRes@vcRefLengthF    = 0.02
windRes@vcRefAnnoOn = True
windRes@vcRefAnnoFontHeightF =0.015
windRes@vcRefAnnoZone = 1
windRes@vcRefAnnoSide = "Top"
windRes@vcRefAnnoJust = "TopRight"
windRes@vcRefAnnoString1 = "20 m/s"
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
;cRes@cnMinLevelValF       = 0               ; set min contour level
;res@cnMaxLevelValF      = 20               ; set max contour level
;cRes@cnLevelSpacingF      = 5              ; set contour spacing
;cRes@cnLineThicknessF = 4
;cRes@cnLevelSelectionMode = "ExplicitLevels"
;cRes@cnLevels = (/0, 5, 10, 15,20/)
;cRes@cnFillColors=(/120,150,180,210,240,250/)

cRes@lbBoxLinesOn = True ; Label bar
cRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
;cRes@lbLabelFontHeightF   =0.02
;cRes@lbOrientation        = "Vertical"         ; vertical label bar
;cRes@pmLabelBarWidthF = 0.05     
cRes@pmLabelBarHeightF = 0.05    
cRes@pmLabelBarOrthogonalPosF  = 0.03 ;色标离图的距离
cRes@lbBoxEndCapStyle = "TriangleBothEnds"


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

heightRes@cnLineLabelBackgroundColor=-1
heightRes@cnLabelMasking = True

heightRes@cnLevelSelectionMode = "ManualLevels" 
heightRes@cnMinLevelValF =1100 
heightRes@cnMaxLevelValF = 1250
heightRes@cnLevelSpacingF = 4 
heightRes@cnLineThicknessF = 4
;heightRes@gsnSpreadColorStart = 1
;heightRes@gsnSpreadColorStart =17
;heightRes@gsnSpreadColors = True ; 6.1.0
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
   if (HH(i) .eq. "00".and.DD(i) .eq. "27") then
      dum = gsn_add_text(wks,mapplot,"~F37~p",lon(i), lat(i),txres) 
    end if
end do
draw(mapplot)
frame(wks)
end