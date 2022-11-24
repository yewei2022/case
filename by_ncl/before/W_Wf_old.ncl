; 等压面速度与地形强迫垂直速度做对比，sigmal不行，用近地面的呢
; 4.10 读取 所有 时次文件
;    目的：绘制time变化图， lat 26~32 lon 92~100 区域平均
; 4.15 将时次文件处理平均之后的时间剖面图
;    通过 restick属性设置XB的time label
;      参考函数：https://www.ncl.ucar.edu/Document/Functions/User_contributed/time_axis_labels.shtml
;    通过 grib_stime2time 将字符串时间转化为YYYYCCDDHHMMSS形式 用整除分段 cd_inv_cdalendar转化为Julian历
;      参考函数: https://www.ncl.ucar.edu/Document/Functions/Built-in/cd_inv_calendar.shtml
;    通过 gsn_csm_pres_hgt 绘制不同高度的剖面图
;      参考绘图集: https://www.ncl.ucar.edu/Applications/height_time.shtml 
; 6.15 插值 linint2要求经纬度单增
;参考https://www.ncl.ucar.edu/Document/Functions/Built-in/linint2.shtml
;load "~/offical_sample_script/download/pvs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("read_elev_data") ;读地形文件的函数
function read_elev_data(topo_file)
local nlat, nlon, topo_file, lat, lon
begin
  ;---Read data as a straight binary file
  nlat = 2160
  nlon = 4320
  setfileoption("bin","ReadByteOrder","BigEndian")
  elev= cbinread(topo_file,(/nlat,nlon/),"short")

  ;---Create 1D coordinate arrays
  lat       = fspan(90,-90,nlat)
  lon       = fspan(0,360,nlon)
  lat!0     = "lat"
  lon!0     = "lon"
  lat@units = "degrees_north"
  lon@units = "degrees_east"
  lat&lat   = lat
  lon&lon   = lon

  ;---Attach the coordinate arrays
  elev!0    = "lat"
  elev!1    = "lon"
  elev&lat  = lat
  elev&lon  = lon
 
  return(elev)
end


begin

;-------------1.变量数据读入--------------------

  path_rawdata = "/cygdrive/D/ncl_related/data/fnl/"
  path_picture ="/cygdrive/D/ncl_related/pictures/terrain/"  
  nc_files     = systemfunc("ls "+path_rawdata+"fnl_*_*_*.nc")
  f            = addfiles(nc_files,"r") ; f[0-12]
  ListSetType(f,"join")    ;打开多个文件并且返回文件句柄(是列表),指定如何返回文件变量
  f_time       =  f[:]->LAND_P0_L1_GLL0@initial_time; 注意原始数据中有些变量如 LAND_P0_L1_GLL0 有时间属性

  vvel1_in     =f[:]->$"VVEL_P0_L100_GLL0"$ ;float VVEL_P0_L100_GLL0 ( lv_ISBL4, lat_0, lon_0 )
  u_sig_in     =f[:]->$"UGRD_P0_L104_GLL0"$ ;10m风速UGRD_P0_L103_GLL0，sigmalfloat UGRD_P0_L104_GLL0 ( lat_0, lon_0 )
  v_sig_in     =f[:]->$"VGRD_P0_L104_GLL0"$ ;10m风速VGRD_P0_L103_GLL0，sigmal层 VGRD_P0_L104_GLL0 ( lat_0, lon_0 )

  tmp_in       =f[:]->$"TMP_P0_L100_GLL0"$
  pres1        =f[:]->$"lv_ISBL4"$
  pres         =pres1(0,:) ;Pa
  pres_hpa     =pres*0.01
  pres_hpa@units="hPa"
  lon_in       =u_sig_in&lon_0
  lat_in       =u_sig_in&lat_0

  ; 修改维度名
  tmp_in!0 = "time"
  vvel1_in!0="time"
  u_sig_in!0="time" 
  v_sig_in!0="time"

;---------------------2.时间格式的修改-----------------------------
  time_set = grib_stime2itime(f_time) ; 将grib格式时间由字符串转化为整型
  time_set_yy = time_set/1000000
  time_set_mm = time_set/10000%100
  time_set_dd = time_set/100%100
  time_set_hh = time_set%100
  time_set_min = time_set_hh
  time_set_sec = time_set_hh
  time_set_min = 0
  time_set_sec = 0
  ; print(time_set_hh)
  ; print(time_set)
  ; 时间格式重新写为 ju格式
  u_sig_in&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)

  v_sig_in&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)

  vvel1_in&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)

  tmp_in&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)



;---------------------------------------3.读地形文件------------------------------------------
  ;其实fnl资料有地形高度数据，float HGT_P0_L1_GLL0 ( lat_0, lon_0 )
  ;lon1       fspan(0,360,4320)  ;gsn_xy    
  ;lat1     =fspan(90,-90,2160)  ;用于gsn_xy 地形图
  datadir  ="/cygdrive/D/ncl_related/data/other/"      ;数据、图片所在路径D:\ncl_related\data\other
  setfileoption("nc","MissingToFillValue",False)    ;设置缺省值 
  ;参考网址 http://www.ncl.ucar.edu/Document/Functions/Built-in/setfileoption.shtml
  elev_old= read_elev_data(datadir+"ETOPO5.DAT")
  elev_float=short2flt(elev_old)
  copy_VarCoords(elev_old, elev_float)
  elev_float=where(elev_float.ge.0., elev_float, elev_float@_FillValue)


;----------------------------------------------网格点插值----------------------------------------
lon_out=fspan(80,110,301)
lat_out=fspan(0, 40, 401)
lon_out@units="degrees_east"
lat_out@units="degrees_north"
   
  u_sig     =linint2_Wrap(lon_in, lat_in(::-1), u_sig_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(u_sig_in)
  v_sig     =linint2_Wrap(lon_in, lat_in(::-1), v_sig_in(:,::-1,:), True, lon_out, lat_out, 0)
  delete(v_sig_in)
  vvel1     =linint2_Wrap(lon_in, lat_in(::-1), vvel1_in(:,:,::-1,:), True, lon_out, lat_out, 0)
  delete(vvel1_in)
  tmp       =linint2_Wrap(lon_in, lat_in(::-1), tmp_in(:,5:25,::-1,:), True, lon_out, lat_out, 0)
  delete(tmp_in)
  elev      = linint2_Wrap(elev_float&lon,elev_float&lat(::-1),elev_float(::-1,:),True,lon_out, lat_out,0)
  ;[LAT | 401] x [LON | 301]
  delete(elev_old)
  delete(elev_float)


; ;------------------------数据处理 ：等压面垂直速度------------------------------------
;   vvel2 = omega_to_w(vvel1,conform(vvel1,pres, 1), tmp)  ; w[*][*][*] (m/s)，要求温度K，气压Pa
;   copy_VarCoords(vvel1, vvel2)  ;[time | 12] x [lv_ISBL4 | 21] x [LAT_0 | 401] x [LON_0 | 301]
;   LEVEL=60000
;   vvel_areas= vvel2(time|:,lv_ISBL4|:,{LAT_0|28:30},{LON_0 |90:98}) ;确定范围
;   vvel_area_ave= dim_avg_n_Wrap(vvel_areas, (/2,3/)) ; 对lat 和 lon 做区域平均
;   vvel=vvel_area_ave(:,{LEVEL})
;   ; printMinMax(vvel2, 0); vertical velocity (m/s) : min=-0.352249   max=0.583592
;   ; exit
;   printMinMax(u_sig(:,{28:30},{90:98}), 0);  U-component of wind (m s-1) : min=-12.5   max=15.5
;   exit ;   V-component of wind (m s-1) : min=-8.2   max=9.8

  

;--------------------------数据处理: 地形强迫垂直速度------------------

;算地形梯度
 nlat=dimsizes(elev&LAT)
 mlon=dimsizes(elev&LON)
 lat    =elev&LAT
 lon    =elev&LON

 dlon   = (lon(2)-lon(1))*0.0174533 ; convert to radians 0.017=3.14/180
 dHdX   = new ( dimsizes(elev),typeof(elev)); pre-allocate space
 do nl=0,nlat-1                      ; loop over each latitude
  dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
  dHdX(nl:nl,:) = center_finite_diff_n (elev(nl:nl,:), dX , True,0,1) ; 1指的经度那一维
 end do

 dY   =(lat(2)-lat(1))*6378388.*0.0174533 
 ;一个纬度的长度，如果格点是1*1做差分时可以直接用这个dY，但如果不是1*1需要变一下
 dHdY = new ( dimsizes(elev),typeof(elev)); pre-allocate space
 do nl=0,mlon-1                      ; loop over each latitude
  dHdY(:,nl:nl) = center_finite_diff_n (elev(:,nl:nl), dY , False,0,0)
 end do

;算速度
 w1=u_sig*conform(u_sig,dHdX,(/1,2/))+v_sig*conform(v_sig,dHdY,(/1,2/))
 copy_VarCoords(u_sig, w1)
 ;[time | 12] x [LAT_0 | 401] x [LON_0 | 301]
 ; printMinMax(w1, 0); min=-0.999717   max=1.40928
 ; exit
;算区域平均
  w_areas = w1(time|:,{LAT_0|28:30},{LON_0|90:98}) ;确定范围
  w_area_ave = dim_avg_n_Wrap(w_areas, (/1,2/)) ; 对lat 和 lon 做区域平均
  ; w_area_ave=  w_area_ave*(-1.0)  
  copy_VarCoords(w_areas,w_area_ave)
  w_area_ave!0="time"  

;将变量整合到一个数组
 data  = new((/2,12/),float)

 data(0,:) = vvel
 data(1,:) = w_area_ave
 data!1="time"
 data&time=u_sig&time
 ;print(data&time)
 ;exit 



wks  = gsn_open_wks ("png",path_picture+"W_Wf_"+LEVEL)          ; send graphics to PNG file



;---------------------------------------绘图属性设置------------------------------------------
res = True
   res@gsnDraw = False
   res@gsnFrame = False

 
  ;res@gsnMaximize=True
  res@tmXTOn = False ;关刻度
  res@tmYROn = False
  ;res@gsnLeftString = "hPa" ;左上角字符
  res@tiYAxisString="W / W~B~f~N~  (m/s)"
  res@tiYAxisFont=22
  ;设置坐标轴
  res@tmXBLabelFont=22
  res@tmXBLabelAngleF = 48.   ;角度
  res@tmYLLabelFont=22
  res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
  res@tmYLLabelFontHeightF =0.02
  ;设置画图区域大小
  res@vpWidthF=0.7
  res@vpHeightF=0.4

  ;>>>时间轴设置
  trestick = True
  trestick@ttmFormat = "%D/%H" ; 时间轴的格式
  trestick@ttmAxis   = "XB"
  trestick@ttmNumTicks = 12    ; 显示所有12个点


vvelres=res
  vvelres@xyDashPattern  = 0                         ; Make curves all solid
  vvelres@xyLineThicknessF= 8
  vvelres@xyMarkLineMode = "MarkLines"               ; Markers *and* lines
  vvelres@xyMarker      = 16               ; 2 different markers
  vvelres@xyMarkerColors = (/"blue","red"/) ; 2 different colors
  vvelres@xyLineColors =(/"blue","red"/)
  vvelres@trYMinF=min(data)-0.002



;------------------------------------------绘图-----------------------------------------------------


time_axis_labels(data&time,vvelres,trestick) ; call the formatting procedure 
plot  = gsn_csm_xy (wks,data&time,data,vvelres) ; create plot

;-------------------------------------------------------------------------------------
; 绘制图例
  resLg = True
  colours=(/"blue","red"/)
 
  resLg@lgItemType = "MarkLines"

  resLg@lgMonoMarkerIndex = True
  resLg@lgMarkerColors = colours
  resLg@lgMarkerIndex = 1
  resLg@lgMarkerSizeF = 0.04
 
 
  resLg@lgMonoDashIndex = True
  resLg@lgDashIndex = 0
  resLg@lgLineColors = colours
  resLg@lgLineThicknessF = 3
 
  resLg@vpWidthF = 0.14
  resLg@vpHeightF = 0.1
 
  resLg@lgPerimFill = 0
  resLg@lgPerimFillColor = "Background"
 
  resLg@lgLabelFontHeightF = 0.08
 
  ;resLg@lgTitleFontHeightF = 0.015
  ;resLg@lgTitleString = ""
 
  lbid = gsn_create_legend(wks, 2, (/"W","W~B~f~N~"/), resLg)
 
; 将图例放置在图中
  amres = True
  amres@amParallelPosF = 0.398
  amres@amOrthogonalPosF = -0.37
  dumLg = gsn_add_annotation(plot, lbid, amres)

draw(plot)
frame(wks)
end