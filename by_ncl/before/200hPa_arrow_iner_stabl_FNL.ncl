;用fnl资料,风矢是箭头，惯性稳定度，时间改一处，读取变量那儿,cmap的颜色可颠倒
;!!!画指定风圈
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"
;代码路径  ncl /mnt/d/ncl_related/code/200hPa_arrow_iner_stabl.ncl

begin
;手动设置一个时间 
;用各个数据生成一个带unit的 timeset_ymdh 
;转换为一个字符串 方便文件中读入
   year   = 2008
   month  = 10
   day    = 26
   hour   = 00
   minute = 00
   second = 0
  ;设置时间格式读取数据####################################
   units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
   opt = 0
   timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
   file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)
   path_format = "%Y%N%D%H"   ; 台风路径文件中的时间形式
   path_time  = cd_string(timeset_ymdh,path_format)
   ;dd=tostring(day)
   ;hh=tostring(hour)

;数据读取和处理
  data_dir  ="/mnt/d/ncl_related/data/fnl/"
  pic_dir   ="/mnt/d/ncl_related/pictures/200/"
  data      = addfile(data_dir+"fnl_"+file_time+".nc", "r")
  u_in      = data ->UGRD_P0_L100_GLL0({20000},:,:)  ;单位：m/s
  v_in      = data->VGRD_P0_L100_GLL0({20000},:,:)
  lon_in    = data->lon_0
  lat_in    = data->lat_0
  hgt_in    = data->HGT_P0_L100_GLL0({20000},:,:)  ;13={50000}  hgt 的单位：gpm

;--------------------读取风暴中心位置-------------------------------------------------   
fiTY = "/mnt/d/ncl_related/data/TC/bio042008.txt"
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

tctime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
stdlon=0.
stdlat=0.
do i = 0, nrow-1
  if (i.eq.tctime_ind) then 
    stdlon=lon(i)
    stdlat=lat(i)
  end if
end do

lon0=floattoint(stdlon)
lat0=floattoint(stdlat)

;-----------------------------------------网格点插值--------------------
  lon_out=fspan(55,125,701)
  lat_out=fspan(0, 60, 601)
   u      =linint2_Wrap(lon_in, lat_in(::-1), u_in(::-1,:), True, lon_out, lat_out, 0)
   delete(u_in)
   v      =linint2_Wrap(lon_in, lat_in(::-1), v_in(::-1,:), True, lon_out, lat_out, 0)
   delete(v_in)
   hgt_out=linint2_Wrap(lon_in, lat_in(::-1), hgt_in(::-1,:), True, lon_out, lat_out, 0)
   delete(hgt_in)

  hgt     = hgt_out/10 ;位势高度处理
  copy_VarCoords(hgt_out, hgt)
  

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

  X=x-conform(X, X0, 0) ;存疑。这里应该是X=x-conform(X, X0, 1)？
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
  iner_stabl0= (vr+f)*(2.*Vt*rinv+f)
  copy_VarCoords(u, iner_stabl0)
  iner_stabl0=iner_stabl0*1.e9    ;[LAT_0 | 601] x [LON_0 | 701]


;-------------------------显示风暴圆圈图------------------------------
  clat = new(100,float)    ; Create arrays to hold great circle.
  clon = new(100,float)    ; As of NCL 6.5.0, can be float or double

  nggcog(stdlat,stdlon,10.,clat,clon) ;计算10个经纬度范围内的圆的经纬度

;---For each center and radius: extract extent of each circle [ reduce gc_inout time ]
;---Use NCL's 'reassignment syntax [ := ] to accomodate possiblr changing array sizes

  min_lat  = min(clat)            ; min of current latitude polygon
  max_lat  = max(clat)            ; max       "
  min_lon  = min(clon)            ; min       "    longitude polygon
  max_lon  = max(clon)

;---Extract the desired rectangle of data 

  iner_stabl = iner_stabl0({min_lat:max_lat},{min_lon:max_lon}) ;[LAT_0 | 200] x [LON_0 | 219]

; ;---Set points that are outside of the circle of data to missing
  lat2d := conform(iner_stabl, iner_stabl&LAT_0, 0)
  lon2d := conform(iner_stabl, iner_stabl&LON_0, 1)
  latlon_circle := gc_inout(lat2d,lon2d, clat,clon)
  iner_stabl = where(latlon_circle,iner_stabl, iner_stabl@_FillValue)

  ; printVarSummary(iner_stabl)
  ; exit

wks = gsn_open_wks("png", pic_dir+"_iner_stabl_"+file_time)
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
res@tmXBLabelFontHeightF = 0.016 ; x轴标签字号
res@tmYLLabelFontHeightF =0.016


;设置地图
mpres = res
mpres@mpFillOn                = False
mpres@mpMinLatF               =lat0-11
mpres@mpMaxLatF               =lat0+12
mpres@mpMinLonF               =lon0-15
mpres@mpMaxLonF               =lon0+15
mpres@gsnMajorLonSpacing    =4. 
mpres@gsnMinorLonSpacing    =1.    
mpres@gsnMajorLatSpacing    =4. 
mpres@gsnMinorLatSpacing    =1
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

;--------------------惯性稳定度---------------------------
iner_stablRes=res
iner_stablRes@cnLinesOn         =False
iner_stablRes@cnFillOn          =True                  ; turn on color
iner_stablRes@cnInfoLabelOn         = False ; set false will not set an information label
; iner_stablRes@cnFillPalette         = "MPL_PiYG"
;--------------------颜色设置----------------------
cmap  = read_colormap_file("MPL_PiYG")
cmap = cmap(::-1,:)
iner_stablRes@cnFillPalette = cmap
iner_stablRes@cnLevelSelectionMode  = "ManualLevels" 
iner_stablRes@cnLevelSpacingF       = 1             ; set contour spacing,2600-0.1,2700-0.2
iner_stablRes@cnMaxLevelValF        = 5               ; set max contour level
iner_stablRes@cnMinLevelValF        =-5
; iner_stablRes@cnLevelSelectionMode   = "ExplicitLevels"
; iner_stablRes@cnLevels               = levs
; iner_stablRes@cnFillColors           = colors ;(/20,45,70,95,120,145,170,195,220/)
iner_stablRes@lbBoxLinesOn         = True ; Label bar
iner_stablRes@lbLabelAutoStride    = True            ; let NCL determine label spacing
;iner_stablRes@lbLabelFontHeightF   = 0.02
iner_stablRes@lbOrientation        = "Vertical"         ; vertical label bar
iner_stablRes@pmLabelBarWidthF     = 0.04     
;iner_stablRes@pmLabelBarHeightF   = 0.05    
;iner_stablRes@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
iner_stablRes@lbBoxEndCapStyle     = "TriangleBothEnds"
iner_stablRes@cnFillDrawOrder      ="PreDraw"



mapplot = gsn_csm_map(wks, mpres)
windPlot = gsn_csm_vector(wks, u , v , windRes)
contour_iner_stabl=gsn_csm_contour(wks,iner_stabl, iner_stablRes)

;叠加
overlay(mapplot, contour_iner_stabl)
overlay(mapplot, windPlot)


;------------------------某时刻风暴位置标记---------------------------------
txres               = True                     ; text mods desired
txres@txFontHeightF = 0.03                     ; font smaller. default big
txres@txFontThicknessF=6
txres@txFontColor="red"
dum = gsn_add_text(wks,mapplot,"~F37~p",stdlon,stdlat,txres) 

;高原轮廓线
shp="/mnt/d/ncl_related/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres) 
;画圆圈
  plres=True
  plres@gsLineColor      ="black"
  plres@gsLineThicknessF = 3
  circle=gsn_add_polyline(wks, mapplot, clon, clat, plres)

draw(mapplot)
frame(wks)
end