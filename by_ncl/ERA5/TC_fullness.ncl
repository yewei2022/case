;用ERA5资料,风矢是箭头，高度场，散度，惯性稳定度，不插值,没用
;时间改三处，14 data，18 need_time，52 找中心经纬度
;  ncl /mnt/h/d/ncl_related/ERA5_code/TC_fullness.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin

data_dir="/mnt/h/d/data/ERA5/"   
pic_dir ="/mnt/h/d/ncl_related/pictures/storm/ERA5/fullness/new/"
data    = addfile(data_dir+"ERA5.nc", "r")

;找时间
  need_time=2008102618
  path_time=tostring(need_time)
  time    =data->time
  YYYYMMDDHH_f=cd_calendar(time,-3)
  time_locat = ind(YYYYMMDDHH_f.eq.need_time) ;找该时间所处位置

;数据读取
  u_in = short2flt(data ->u(time_locat,{1000},{10:30},{80:100}))  ;单位：m/s
  v_in = short2flt(data ->v(time_locat,{1000},{10:30},{80:100}))
  lon_in=data->longitude({80:100})   
  lat_in=data->latitude({10:30}) 
;网格点插值成0.1*0.1
  lon_out=fspan(80,100,201)
  lat_out=fspan(10, 30, 201)
  lon_out@units="degrees_east"
  lat_out@units="degrees_north"
  u      =linint2_Wrap(lon_in, lat_in, u_in, False, lon_out, lat_out, 0)
  delete(u_in)
  v      =linint2_Wrap(lon_in, lat_in, v_in, False, lon_out, lat_out, 0)
  delete(v_in)
;计算风速
  U = sqrt(u^2+v^2) 
  copy_VarCoords(u, U)

;读取风暴中心位置   
  fiTY = "/mnt/h/d/data/TC/bio042008.txt"
  ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
  nrow = numAsciiRow(fiTY)
  YYYYMMDDHH = new(nrow, "string")
  lat_tc = new(nrow, "float")
  lon_tc = new(nrow, "float")
  cdata = asciiread(fiTY, -1, "string")
  YYYYMMDDHH = str_get_field(cdata, 3, ",")
  lat_tc = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
  lon_tc = stringtofloat(str_get_field(cdata, 8, ",")) *0.1
  rmw = stringtofloat(str_get_field(cdata, 20, ","))*1852.

  stdlon=0.
  stdlat=0.
  rmw_need=0.
  tctime_ind = str_match_ind(YYYYMMDDHH, path_time) 
  ; str can not be judge by .eq.
  ; return YYYYMMDDHH and path_time's matched position 
  do i = 0, nrow-1
    if (i.eq.tctime_ind) then
      stdlon=lon_tc(i)
      stdlat=lat_tc(i)
      rmw_need=rmw(i)
    end if
  end do

;计算丰满度 r17
  ;1. 函数计算r
  ;挑出U17的格点的经纬度
  lat2d= conform_dims(dimsizes(u), lat_out, 0)
  copy_VarCoords(u, lat2d)
  lon2d= conform_dims(dimsizes(u), lon_out, 1)
  copy_VarCoords(u, lon2d)
  lat1d =ndtooned(lat2d)
  lon1d =ndtooned(lon2d)
  U1D=ndtooned(U)
  U1D = where(U1D.ge.16.99.and.U1D.le.17.01,U1D, U1D@_FillValue)
  indexes = ind(.not.ismissing(U1D))
  ; printVarSummary(indexes)
  ; exit
  npts_range = dimsizes(indexes)   ; # of points in this range.
  lat17=new(npts_range, float)
  lon17=new(npts_range, float)
  lat17(0:npts_range-1) = lat1d(indexes)
  lon17(0:npts_range-1) = lon1d(indexes)

  dist = gc_latlon(stdlat,stdlon,lat17,lon17,2,-3)
  r=new(npts_range, float)
  r=dist@spacing 
  r17=avg(r)

  ;2 计算tcf
  tcf=1-rmw_need/r17
  print(r17)
  print(rmw_need)
  print(tcf)
  exit


  ;3 用于绘图，看U17大概位置
  u = where(U.ge.16.99.and.U.le.17.01,u, u@_FillValue) ;U17's u
  v = where(U.ge.16.99.and.U.le.17.01,v, v@_FillValue) ;V17's v


; ;计算r10
;   ;挑出U10的格点的经纬度
;   lat2d= conform_dims(dimsizes(u), lat_out, 0)
;   copy_VarCoords(u, lat2d)
;   lon2d= conform_dims(dimsizes(u), lon_out, 1)
;   copy_VarCoords(u, lon2d)
;   lat1d =ndtooned(lat2d)
;   lon1d =ndtooned(lon2d)
;   U1D=ndtooned(U)
;   U1D = where(U1D.ge.9.99.and.U1D.le.10.01,U1D, U1D@_FillValue)
;   indexes = ind(.not.ismissing(U1D))
;   ; printVarSummary(indexes)
;   ; exit
;   npts_range = dimsizes(indexes)   ; # of points in this range.
;   lat10=new(npts_range, float)
;   lon10=new(npts_range, float)
;   lat10(0:npts_range-1) = lat1d(indexes)
;   lon10(0:npts_range-1) = lon1d(indexes)

;   dist = gc_latlon(stdlat,stdlon,lat10,lon10,2,-3)
;   r=new(npts_range, float)
;   r=dist@spacing 
;   r10=avg(r)
;   print(r10)


;   ;3 用于绘图，看U10大概位置
;   u = where(U.ge.9.99.and.U.le.10.01,u, u@_FillValue) ;U10's u
;   v = where(U.ge.9.99.and.U.le.10.01,v, v@_FillValue) ;V10's v


wks = gsn_open_wks("x11", pic_dir+"fullness_"+path_time)

  res = True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@gsnMaximize = True
  res@gsnLeftString  = ""
  res@gsnRightString = ""
  res@tiMainString   =path_time
  res@tmYROn         =False
  res@tmXTOn         =False
  ; res@pmTickMarkDisplayMode  = "Always" ;turn ° on or off

  ;设置地图
  mpres = res
  mpres@mpFillOn                = False
  mpres@mpMinLatF               = 15.
  mpres@mpMaxLatF               = 25.
  mpres@mpMinLonF               = 85.
  mpres@mpMaxLonF               = 95.
  ; mpres@gsnMajorLatSpacing    =5. 
  ; mpres@gsnMinorLatSpacing    =1
  ; mpres@gsnMajorLonSpacing    =5. 
  ; mpres@gsnMinorLonSpacing    =1. 
  mpres@mpDataSetName              = "/mnt/h/d/ChinaMap/Earth..4"
  mpres@mpDataBaseVersion          = "MediumRes"
  mpres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)
  ; mpres@mpOutlineOn = True  ;国界
  ; mpres@mpOutlineBoundarySets = "National" 

  windRes = res
  ;windRes@vcGlyphStyle            = "FillArrow"   ; turn on curly vectors
  windRes@vcGlyphStyle    = "LineArrow"
  windRes@vcLineArrowThicknessF = 3
  ; windRes@vcMinDistanceF = 0.04  ;default value 0.0
  windRes@vcRefMagnitudeF = 10.
  windRes@vcRefLengthF    = 0.03
  windRes@vcRefAnnoOn = True
  windRes@vcRefAnnoFontHeightF =0.015
  windRes@vcRefAnnoZone = 1
  windRes@vcRefAnnoSide = "Top"
  windRes@vcRefAnnoJust = "TopRight"
  windRes@vcRefAnnoString1 = "10 m/s"
  windRes@vcRefAnnoString2On = False
  windRes@vcLineArrowHeadMaxSizeF = 0.01
  windRes@vcLineArrowHeadMinSizeF = 0.01

;绘图
  mapplot = gsn_csm_map(wks, mpres)
  windplot = gsn_csm_vector(wks, u , v , windRes)
  overlay(mapplot, windplot)

;TC position
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.03                     ; font smaller. default big
  txres@txFontThicknessF=6
  txres@txFontColor="red"
  dum_tc = gsn_add_text(wks,mapplot,"~F37~p",stdlon,stdlat,txres) 
 
draw(mapplot)
frame(wks) 
end