;用fnl资料,区域放大，风矢换成了全风羽
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"

begin
;手动设置一个时间 
;用各个数据生成一个带unit的 timeset_ymdh 转换为一个字符串 方便文件中读入
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

;数据读取
  data_dir="/cygdrive/D/ncl_related/data/fnl/"
  pic_dir ="/cygdrive/D/ncl_related/pictures/200/"
  data    = addfile(data_dir+"fnl_"+file_time+".nc", "r")
  u       = data ->UGRD_P0_L100_GLL0({20000},:,:)  ;单位：m/s
  v       = data->VGRD_P0_L100_GLL0({20000},:,:)
  lon1    =data->lon_0
  lat1    =data->lat_0
  hgtData = data->HGT_P0_L100_GLL0({20000},:,:)  ;13={50000}  hgt 的单位：gpm
  hgt     = hgtData/10
  copy_VarCoords(hgtData, hgt)
  di=uv2dvF_Wrap(u(::-1,:),v(::-1,:)) ;用这个球面间谐波函数计算散度更精确
  ;di = uv2dv_cfd (u,v,lat1,lon1, 3) 
  di=di*1.e5
;max1=max(di({30:40},{90:100}))
;print(max1)
;exit


;m=sqrt(u^2+v^2)
;copy_VarCoords(u, m)
;max1=max(m({20:40},{80:100}))
;copy_VarCoords(u, max1)
;print(max1)
;exit



wks = gsn_open_wks("png", pic_dir+"new_"+file_time)
gsn_define_colormap(wks,"NEO_div_vegetation_c") 


res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnLeftString = ""
res@gsnRightString = ""
;res@gsnMaximize=True
res@tmYROn=False
res@tmXTOn=False
;res@pmTickMarkDisplayMode  = "Always" ;去掉左上角刻度，并且把所有的次刻度都关掉了，过于鸡肋
;res@tiMainString = "500hPa_2900"
;res@tmXBMajorLengthF = 0.005
;res@vpWidthF=0.8  ;画布大小
;res@vpHeightF=0.7

;设置坐标轴标签粗细
res@tmXBLabelFont=22
res@tmYLLabelFont=22
res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
res@tmYLLabelFontHeightF =0.02


;设置地图
mpres = res
mpres@mpFillOn                = False
mpres@mpMinLatF               = 10.
mpres@mpMaxLatF               = 50.
mpres@mpMinLonF               = 60.
mpres@mpMaxLonF               = 120.
mpres@gsnMajorLonSpacing    =10. 
mpres@gsnMinorLonSpacing    =2.    
mpres@gsnMajorLatSpacing    =10. 
mpres@gsnMinorLatSpacing    =2
mpres@mpDataBaseVersion       = "MediumRes"
mpres@mpDataSetName           = "Earth..4"
mpres@mpOutlineSpecifiers     = (/"China:states"/)


;--------------矢量图--------------

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
cRes@cnMinLevelValF       = -8             ; set min contour level
cRes@cnMaxLevelValF       = 8               ; set max contour level
cRes@cnLevelSpacingF      = 2              ; set contour spacing
;cRes@cnLineThicknessF = 4
;cRes@cnLevelSelectionMode = "ExplicitLevels"
;cRes@cnLevels = (/0, 5, 10, 15,20/)
;cRes@cnFillColors=(/120,150,180,210,240,250/)

cRes@lbBoxLinesOn = True ; Label bar
cRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
cRes@lbLabelFontHeightF   =0.02
cRes@lbOrientation        = "Vertical"         ; vertical label bar
cRes@pmLabelBarWidthF = 0.03     
;cRes@pmLabelBarHeightF = 0.05    
;cRes@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
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
heightRes@cnLineThicknessF = 6
;heightRes@gsnSpreadColorStart = 1
;heightRes@gsnSpreadColorStart =17
;heightRes@gsnSpreadColors = True ; 6.1.0
heightRes@lbBoxLinesOn = False ; Label bar
;heightRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
;heightRes@lbOrientation        = "Vertical"         ; vertical label bar


;绘图
di = smth9_Wrap(di,0.5,0.25,False)
di = smth9_Wrap(di,0.5,0.25,False)
hgt = smth9_Wrap(hgt,0.5,0.25,False)
hgt = smth9_Wrap(hgt,0.5,0.25,False)


mapplot = gsn_csm_map(wks, mpres)
contour_div = gsn_csm_contour(wks,di,cRes)
windPlot = gsn_csm_vector(wks, u , v , windRes)
contour_height = gsn_csm_contour(wks,hgt, heightRes)

;叠加风矢量和位势高度
overlay(mapplot, contour_div)
overlay(mapplot, windPlot)
overlay(mapplot, contour_height)

;----------------------------------某时刻风暴位置标记-----------------------------

;读风暴数据,用于标记风暴位置或路径
 
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

txres               = True                     ; text mods desired
txres@txFontHeightF = 0.03                     ; font smaller. default big
txres@txFontThicknessF=6
txres@txFontColor="red"
do i = 0, nrow-1
   if (HH(i) .eq. "00".and.DD(i) .eq. "27") then
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