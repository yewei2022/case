; 气压梯度力
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

data_dir="/cygdrive/D/ncl_related/data/fnl/"   
pic_dir ="/cygdrive/D/ncl_related/pictures/500/big/"
data    = addfile(data_dir+"fnl_"+file_time+".nc", "r")

uData = data ->UGRD_P0_L100_GLL0({50000} ,:,:)  ;单位：m/s
vData = data->VGRD_P0_L100_GLL0({50000} ,:,:)
hgtData = data->HGT_P0_L100_GLL0({50000} ,:,:)  ;13={50000}  hgt 的单位：gpm
lat= data->lat_0
lon= data->lon_0

;计算气压气度力
 nlat=dimsizes(lat)
 mlon=dimsizes(lon)

 dlon   = (lon(2)-lon(1))*0.0174533 ; convert to radians 0.017=3.14/180
 dfdX   = new ( dimsizes(hgtData),float); pre-allocate space
 do nl=0,nlat-1                      ; loop over each latitude
  dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
  dfdX(nl:nl,:) = center_finite_diff_n (hgtData(nl:nl,:), dX , True,0,1)
 end do

 dY   =6378388.*0.0174533
 dfdY = new ( dimsizes(hgtData),float); pre-allocate space
 do nl=0,mlon-1                      ; loop over each latitude
  dfdY(:,nl:nl) = center_finite_diff_n (hgtData(:,nl:nl), dY , False,0,0)
 end do

copy_VarCoords(hgtData, dfdX)
copy_VarCoords(hgtData, dfdY)

grad_f=sqrt(dfdX^2+dfdY^2)
copy_VarCoords(hgtData, grad_f)

grad_p=grad_f*1.e5
copy_VarCoords(hgtData, grad_p)
;print(grad_p({30},{100}))
;exit

wks = gsn_open_wks("png", pic_dir+file_time)
;gsn_define_colormap(wks,"WhiteBlueGreenYellowRed") 


res = True
res@gsnDraw = False
res@gsnFrame = False
;res@gsnMaximize=True
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
mpres@mpMinLatF               = 15.
mpres@mpMaxLatF               = 45.
mpres@mpMinLonF               = 70.
mpres@mpMaxLonF               = 110.
mpres@gsnMajorLonSpacing    =5. 
mpres@gsnMinorLonSpacing    =1.    
mpres@gsnMajorLatSpacing    =5. 
mpres@gsnMinorLatSpacing    =1
mpres@mpDataBaseVersion       = "MediumRes"
mpres@mpDataSetName           = "Earth..4"
mpres@mpOutlineSpecifiers     = (/"China:states"/)

mapplot = gsn_csm_map(wks, mpres)

;--------------风矢量图--------------

windRes = res
windRes@vcGlyphStyle    = "WindBarb"
windRes@vcWindBarbScaleFactorF=2.5 ;转换为2m/s一短杆，4m/s一长杆
windRes@vcWindBarbTickLengthF=0.6 ;指示风速的杆子长短，默认是0.3
windRes@vcWindBarbTickSpacingF=0.125 ;默认分杆之间的间隔为0.125
windRes@vcWindBarbLineThicknessF=3
windRes@vcMinDistanceF = 0.03
windRes@vcRefMagnitudeF = 4.;量级
windRes@vcRefLengthF    = 0.03 ;指示方向的杆子的长短
windRes@vcRefAnnoOn = True
windRes@vcRefAnnoFontHeightF =0.015
windRes@vcRefAnnoZone = 1
windRes@vcRefAnnoSide = "Top"
windRes@vcRefAnnoJust = "TopRight"
windRes@vcRefAnnoString1 = "4 m/s"
windRes@vcRefAnnoString2On = False



;--------------叠加气压梯度力--------------------------------

heightRes = res
heightRes@cnFillOn = True  ;cn means ContourPlot
heightRes@cnLinesOn = False   
heightRes@cnFillPalette = "MPL_BuPu"
heightRes@cnInfoLabelOn = False ; set false will not set an information label
heightRes@cnLineLabelBackgroundColor=-1
heightRes@cnLabelMasking = True
heightRes@cnLevelSelectionMode = "ManualLevels" 
heightRes@cnMinLevelValF =0
heightRes@cnMaxLevelValF = 20
heightRes@cnLevelSpacingF = 2 
;heightRes@cnLineThicknessF = 4
;heightRes@gsnSpreadColorStart = 1
;heightRes@gsnSpreadColorStart =17
;heightRes@gsnSpreadColors = True ; 6.1.0
heightRes@lbLabelBarOn=True
heightRes@lbBoxLinesOn = True ; Label bar
heightRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
heightRes@lbLabelFontHeightF   =0.02
heightRes@lbOrientation        = "Vertical"         ; vertical label bar
heightRes@pmLabelBarWidthF = 0.05     
;cRes@pmLabelBarHeightF = 0.05    
;cRes@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
heightRes@lbBoxEndCapStyle = "TriangleBothEnds"





;绘图



contour_height = gsn_csm_contour(wks,grad_p, heightRes)
windPlot = gsn_csm_vector(wks, uData , vData , windRes)


;叠加风矢量和位势高度
overlay(mapplot, contour_height)
overlay(mapplot, windPlot)

;-----------------某时刻风暴位置标记------------------------------------------------
;读取风暴数据,用于标记风暴位置或路径
 
fiTY = "/cygdrive/D/ncl_related/data/TC/bio042008.txt"

; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
nrow = numAsciiRow(fiTY)
 
YYYYMMDDHH = new(nrow, "string")
latt = new(nrow, "float")
lonn = new(nrow, "float")
 
cdata = asciiread(fiTY, -1, "string")
YYYYMMDDHH = str_get_field(cdata, 3, ",")
latt = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
lonn = stringtofloat(str_get_field(cdata, 8, ",")) *0.1

DateChar = stringtochar(YYYYMMDDHH)
MM = chartostring(DateChar(:,5:6))
DD = chartostring(DateChar(:,7:8))
HH = chartostring(DateChar(:,9:10))
txres               = True                     ; text mods desired
txres@txFontHeightF = 0.03                     ; font smaller. default big
txres@txFontColor="red"
do i = 0, nrow-1
   if (HH(i) .eq. "00".and.DD(i) .eq. "27") then
      dum = gsn_add_text(wks,mapplot,"~F37~p",lonn(i), latt(i),txres) 
    end if
end do


;----------------------绘制三个站点位置-------------------------------

lat2=new(5,float)
lon2=new(5,float)
lat2!0="lat"
lat2@units="degrees_north"
lon2!0="lon"
lon2@units="degrees_east"

lat2(0)=27.98 ;错那
lon2(0)=91.95
lat2(1)=29.86;波密
lon2(1)=95.76
lat2(2)=28.65 ;察隅
lon2(2)=97.46
lat2(3)=27.75 ;贡山
lon2(3)=98.67
lat2(4)=28.45 ;德钦
lon2(4)=98.88

;站点标记
mkres               = True                     ; text mods desired
mkres@txFontColor="red"
mkres@txFontHeightF = 0.025

do i =0,4
dum1 = gsn_add_text(wks,mapplot,"~F35~z",lon2(i),lat2(i),mkres) 
end do
;--------------------------高原轮廓线------------------------------------------------
shp="/cygdrive/D/ncl_related/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres) 

draw(mapplot)
frame(wks)
end