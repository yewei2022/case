;用ERA5资料,风矢是箭头，高度场，散度，惯性稳定度，不插值,没用
;时间改三处，14 data，18 need_time，52 找中心经纬度
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"

begin

data_dir="/cygdrive/D/ncl_related/data/ERA5/"   
pic_dir ="/cygdrive/D/ncl_related/pictures/200/ERA5/"
; data   = addfile(data_dir+"ERA5_24_26.nc", "r")
data    = addfile(data_dir+"ERA5_27_29.nc", "r")

;找时间
  need_time=2008102700
  time    =data->time
  YYYYMMDDHH_f=cd_calendar(time,-3)
  time_locat = ind(YYYYMMDDHH_f.eq.need_time) ;找该时间所处位置

;数据读取
 u = short2flt(data ->u(time_locat,{200},:,:))  ;单位：m/s
 v = short2flt(data->v(time_locat,{200},:,:))
 hgtData = short2flt(data->z(time_locat,{200},:,:))  ; hgt 的单位：m**2 s**-2
 hgt     = hgtData/10/9.8
 copy_VarCoords(hgtData, hgt) ;现在dgpm
 lon_out=data->longitude    ;1440
 lat_out=data->latitude ;721

;读取风暴中心位置   
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

stdlon=0.
stdlat=0.
do i = 0, nrow-1
  ;if ((DD(i).eq.day).and.(HH(i).eq.hour)) then
  if ((DD(i).eq."27").and.(HH(i).eq."00")) then
  stdlon=lon(i)
  stdlat=lat(i)
  end if
end do

 lon0=floattoint(stdlon) ;将台风应用于固定出图的经纬度范围
 lat0=floattoint(stdlat)
  

  ;----------------------------计算惯性稳定度---------------------------------------------
 ;1.设置三维的距离变量,地转参数f
  r= new ( dimsizes(u),typeof(u)); pre-allocate space
  copy_VarCoords(u, r)
  f= new(dimsizes(u), typeof(u))
  copy_VarCoords(u, f)
  x= new ( dimsizes(u),typeof(u))
  copy_VarCoords(u, x) 
  y= new(dimsizes(u), typeof(u))
  copy_VarCoords(u, y)
 ;2.给f、x，y，r 变量赋值
  nlat=dimsizes(lat_out)
  mlon=dimsizes(lon_out)
  dY   =6378388.*0.0174533 ;一个纬度的单位弧长
  dX= new(dimsizes(u), float)
  copy_VarCoords(u, dX)
  ;dlon = (lon_out(2)-lon_out(1))*0.0174533 ; convert to radians 0.017=3.14/180
  ;不对，上面这个式子代表0.1个经度的弧长
  ;应该改成1个经度的
  dlon=0.0174533 ;R*1°/180°
  do nl=0,nlat-1                      ; loop over each latitude
    f(nl,:) = 2*7.292*sin(lat_out(nl)*3.14159/180.0)*0.00001  ;f=2wsin0
    y(nl,:) = dY*lat_out(nl)
    dX(nl,:)= 6378388.*cos(0.0174533*lat_out(nl))*dlon  ; constant at this latitude，R*COS()*单位经度弧长
  end do
  do nl=0,mlon-1                      ; loop over each latitude
    x(:,nl) =dX(:,nl)*lon_out(nl)
  end do
 ;3.计算每个格点相对风暴中心的x,y,r
  X=x
  Y=y
  X0= x({stdlat},{stdlon})
  Y0= y({stdlat},{stdlon})

  X=x-conform(X, X0, 0) ;存疑。这里应该是X=x-conform(X, X0, 1)？其实没啥区别
  Y=y-conform(Y, Y0, 0)

  r=sqrt(X^2+Y^2)
  r@units="m"

  
 ;4计算垂直涡度
  vr = uv2vr_cfd (u,v,lat_out,lon_out, 2) 
  copy_VarCoords(u, vr)

 ;5.计算切向风
  rinv = 1. / where(r.ne.0,r,r@_FillValue)
  ; Vt=v*abs(X)*rinv-u*abs(Y)*rinv
  Vt=v*X*rinv-u*Y*rinv
  copy_VarCoords(u, Vt) 

 ;6. 计算惯性稳定度
  iner_stabl= (vr+f)*(2.*Vt*rinv+f)
  copy_VarCoords(u, iner_stabl)
  iner_stabl=iner_stabl*1.e9    ;
  
  levs  =fspan(-4, 4, 9)
  colors=ispan(20, 220, 25)

  ; print(max(iner_stabl({30},{92})))
  ; exit



wks = gsn_open_wks("png", pic_dir+"_iner_stabl_"+need_time)
;gsn_define_colormap(wks, color_map) 
;wkColorMap contains special elements (0 and 1) for 
;the Background and Foreground colors,
; whereas the palette-type resources do not.



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
mpres@mpMinLatF               =16
mpres@mpMaxLatF               =36
mpres@mpMinLonF               =75
mpres@mpMaxLonF               =105
; mpres@mpMinLatF               =lat0-5
; mpres@mpMaxLatF               =lat0+5
; mpres@mpMinLonF               =lon0-5
; mpres@mpMaxLonF               =lon0+5
mpres@gsnMajorLonSpacing    =4. 
mpres@gsnMinorLonSpacing    =2.    
mpres@gsnMajorLatSpacing    =4. 
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


;--------------叠加高度场--------------

hgtRes = res
hgtRes@cnFillOn = False  ;cn means ContourPlot
hgtRes@cnLinesOn = True   
;hgtRes@cnMonoLineColor = False
hgtRes@cnLineColor="blue"
;hgtRes@cnSpanLinePalette = True
;hgtRes@cnLinePalette = "rainbow"
hgtRes@cnLineLabelsOn = True  
hgtRes@cnInfoLabelOn = False ; set false will not set an information label

hgtRes@cnLineLabelBackgroundColor=-1
hgtRes@cnLabelMasking = True

hgtRes@cnLevelSelectionMode = "ManualLevels" 
hgtRes@cnMinLevelValF = 1100 
hgtRes@cnMaxLevelValF = 1250
hgtRes@cnLevelSpacingF = 4 
hgtRes@cnLineThicknessF = 6
;hgtRes@gsnSpreadColorStart = 1
;hgtRes@gsnSpreadColorStart =17
;hgtRes@gsnSpreadColors = True ; 6.1.0
hgtRes@lbBoxLinesOn = False ; Label bar
;hgtRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
;hgtRes@lbOrientation        = "Vertical"         ; vertical label bar

;--------------------惯性稳定度---------------------------
iner_stablRes=res
iner_stablRes@cnLinesOn         =False
iner_stablRes@cnFillOn          =True                  ; turn on color
iner_stablRes@cnInfoLabelOn         = False ; set false will not set an information label
iner_stablRes@cnFillPalette         = "MPL_PiYG"
iner_stablRes@cnLevelSelectionMode  = "ManualLevels" 
iner_stablRes@cnLevelSpacingF       = 1             ; set contour spacing,2600-0.1,2700-0.2
iner_stablRes@cnMaxLevelValF        = 4               ; set max contour level
iner_stablRes@cnMinLevelValF        =-4
; iner_stablRes@cnLevelSelectionMode   = "ExplicitLevels"
; iner_stablRes@cnLevels               = levs
; iner_stablRes@cnFillColors           = colors ;(/20,45,70,95,120,145,170,195,220/)
iner_stablRes@lbBoxLinesOn         = True ; Label bar
iner_stablRes@lbLabelAutoStride    = True            ; let NCL determine label spacing
;iner_stablRes@lbLabelFontHeightF   = 0.02
iner_stablRes@lbOrientation        = "Vertical"         ; vertical label bar
iner_stablRes@pmLabelBarWidthF     = 0.03     
;iner_stablRes@pmLabelBarHeightF   = 0.05    
;iner_stablRes@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
iner_stablRes@lbBoxEndCapStyle     = "TriangleBothEnds"
iner_stablRes@cnFillDrawOrder      ="PreDraw"


;绘图
hgt = smth9_Wrap(hgt,0.5,0.25,False)
hgt = smth9_Wrap(hgt,0.5,0.25,False)

mapplot = gsn_csm_map(wks, mpres)
windPlot = gsn_csm_vector(wks, u , v , windRes)
; contour_height = gsn_csm_contour(wks,hgt, hgtRes)
contour_iner_stabl=gsn_csm_contour(wks,iner_stabl, iner_stablRes)

;叠加风矢量和位势高度
overlay(mapplot, contour_iner_stabl)
overlay(mapplot, windPlot)
; overlay(mapplot, contour_height)


;------------------------某时刻风暴位置标记---------------------------------


txres               = True                     ; text mods desired
txres@txFontHeightF = 0.03                     ; font smaller. default big
txres@txFontThicknessF=6
txres@txFontColor="red"

dum = gsn_add_text(wks,mapplot,"~F37~p",stdlon,stdlat,txres) 

; do i = 0, nrow-1
;    if (HH(i) .eq. "00".and.DD(i) .eq. "26") then
;       dum = gsn_add_text(wks,mapplot,"~F37~p",lon(i), lat(i),txres) 
;     end if
; end do

;高原轮廓线
shp="/cygdrive/D/ncl_related/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres) 

draw(mapplot)
frame(wks)
end