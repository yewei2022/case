load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;************************************************
;数据中时间转换的参考函数：http://www.ncl.ucar.edu/Document/Functions/date.shtml 
;用到的参考函数：https://www.ncl.ucar.edu/Document/Functions/Built-in/cd_inv_calendar.shtml
begin
; ;手动设置一个时间 
;用各个数据生成一个带unit的 timeset_ymdh 
;转换为一个字符串 方便文件中读入
   year   = 2008
   month  = 10
   day    = 25
   hour   = 00
   minute = 00
   second = 0
  ;设置时间格式读取数据####################################
   units  = "days since 0001-01-01 00:00:00"  ; 时间设置，注意单位
   opt = 0
   timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
   file_format = "%D" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)

;读取数据
file_dir="/cygdrive/D/ncl_related/data/trmm_rain/"
pic_dir="/cygdrive/D/ncl_related/pictures/rain/"
a = addfile(file_dir+"p_"+file_time+".nc","r")
p = a ->HQPRECIPITATION (0,:,:)  
;latt=a->LAT
;lonn=a->LON

;latt!0="latitude"
;latt@units="degrees_north"
;lonn!0="longitude"
;lonn@units="degrees_east"




;风暴数据,用于标记风暴位置或路径
 
fiTY =  "/cygdrive/D/ncl_related/data/TC/bio042008.txt"
 
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


wks = gsn_open_wks("png", pic_dir+"trmm"+file_time)
;gsn_define_colormap(wks,"prcp_1") 


;---Set up some map resources.
  mpres                       = True
  mpres@mpFillOn              = False
  ;mpres@gsnMaximize           = True             ; Maximize plot in frame.
  mpres@gsnDraw               = False            ; Will draw later
  mpres@gsnFrame              = False            ; Don't advance the frame
  mpres@gsnLeftString = ""
  mpres@gsnRightString = ""
  mpres@tmYROn=False
  mpres@tmXTOn=False
  ;mpres@pmTickMarkDisplayMode = "Always"

  mpres@mpMinLatF             = 16.
  mpres@mpMaxLatF             = 36.
  mpres@mpMinLonF             = 80
  mpres@mpMaxLonF             = 105
  mpres@gsnMajorLonSpacing    =5. 
  mpres@gsnMinorLonSpacing    =1.    
  mpres@gsnMajorLatSpacing    =4. 
  mpres@gsnMinorLatSpacing    =1


  ;设置坐标轴标签粗细
  mpres@tmXBLabelFont=22
  mpres@tmYLLabelFont=22
  mpres@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
  mpres@tmYLLabelFontHeightF =0.02

  mpres@vpWidthF=1.0
  mpres@vpHeightF=0.6

  mpres@mpDataBaseVersion       = "MediumRes"
  mpres@mpDataSetName           = "Earth..4"
  mpres@mpOutlineSpecifiers     = (/"China:states"/)
  ;边界线
  mpres@mpOutlineBoundarySets           = "National"
  mpres@mpNationalLineColor                        = "black"
  mpres@mpProvincialLineColor                = "black"
  mpres@mpGeophysicalLineColor          = "black"
  mpres@mpNationalLineThicknessF         = 3
  mpres@mpProvincialLineThicknessF  = 3
  mpres@mpGeophysicalLineThicknessF = 3


;画等值线图
preRes = mpres
preRes@cnFillOn = True  ;cn means ContourPlot
preRes@cnLinesOn = True   
preRes@cnMonoLineColor = False ;线条不使用同一个颜色
preRes@cnSpanLinePalette = True
preRes@cnFillPalette = "prcp_1"
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
preRes@cnLevelSpacingF = 10 

;preRes@gsnSpreadColors = True ; 
;preRes@gsnSpreadColorStart = 1
;preRes@gsnSpreadColorStart =17 
preRes@lbBoxLinesOn = True ; Label bar
preRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
preRes@lbOrientation        = "Vertical"         ; vertical label bar 色标垂直
preRes@pmLabelBarWidthF = 0.05     ;色标变窄
;preRes@pmLabelBarHeightF = 0.02    ;色标变细
preRes@lbBoxEndCapStyle = "TriangleBothEnds"


;某时刻风暴位置标记
;txres               = True                     ; text mods desired
;txres@txFontHeightF = 0.03                     ; font smaller. default big
;txres@txFontColor="red"
;do i = 0, nrow-1
;   if (HH(i) .eq. "00".and.DD(i) .eq. "27") then
;      dum = gsn_add_text(wks,plot,"~F37~p",lon(i), lat(i),txres) 
;    end if
;end do


;绘图
p = smth9_Wrap(p,0.5,0.25,False)
;p = smth9_Wrap(p,0.5,0.25,False)
plot = gsn_csm_contour_map(wks,p,preRes)   



;某日风暴路径绘制
resDot = True
resLine = True
resDot@gsMarkerIndex = 1
resDot@gsMarkerSizeF = 0.05
resDot@gsMarkerColor = "red"
resLine@gsLineThicknessF = 4
resLine@gsLineColor = "black"

nnrow = 4
dumDot= new(nnrow, graphic)
dumLine =new(nnrow, graphic)
lon1=new(nnrow,float)
lat1=new(nnrow,float)

ii = 0
   do i=0,nrow-1
     if (DD(i) .eq. "25") then
       lon1(ii)=lon(i)
      lat1(ii)=lat(i)
     ii = ii+1
 end if 
    end do

;画线
do i = 0, nnrow-2
    xx = (/ lon1(i), lon1(i+1)/)
    yy = (/ lat1(i), lat1(i+1)/)
    dumLine(i) = gsn_add_polyline(wks, plot, xx, yy, resLine)
end do 

;画点
do i = 0, nnrow-1
  dumDot(i) = gsn_add_polymarker(wks, plot, lon1(i), lat1(i), resDot)
end do

draw(plot)
frame(wks)
end