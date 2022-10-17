; 计算某时刻地形强迫垂直速度，近地面风速取sigma
; 6.15 插值了物理量 linint2要求经纬度单增 但地形高度为保持高分辨率，没有插值
;时间改1处
; 参考https://www.ncl.ucar.edu/Document/Functions/Built-in/linint2.shtml
; load "~/offical_sample_script/download/pvs.ncl"
;代码路径  ncl /mnt/d/ncl_related/code/Wf.ncl

; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" ;不知道这个是干嘛用的
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"


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
;手动设置一个时间 
;用各个数据生成一个带unit的 timeset_ymdh 
;转换为一个字符串 方便文件中读入
  year   = 2008
  month  = 10
  day    = 27
  hour   = 06
  minute = 00
  second = 0
  ;设置时间格式读取数据####################################
  units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
  opt = 0
  timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
  file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
  file_time  = cd_string(timeset_ymdh,file_format)

  data_dir="/mnt/d/data/fnl/"   
  pic_dir ="/mnt/d/ncl_related/pictures/terrain/"
  data    = addfile(data_dir+"fnl_"+file_time+".nc", "r")
  ; u_in  = data ->UGRD_P0_L100_GLL0 ;单位：m/s 等压面速度lv_ISBL0
  ; v_in  = data->VGRD_P0_L100_GLL0
  u_in  = data->UGRD_P0_L104_GLL0  ;单位：m/s sigma层 104,10m风速103
  v_in  = data->VGRD_P0_L104_GLL0 
  vvel_in = data->VVEL_P0_L104_GLL0   ;sigmal 层垂直速度
  pres    = data->lv_ISBL0
  pres_hpa     =pres*0.01
  pres_hpa@units="hPa"
  pres_x       =pres_hpa(11:17)  ;7个，400-700
  pres_x@units ="hPa"
  lon_in       =u_in&lon_0
  lat_in       =u_in&lat_0



;---------------------------------------3.读地形文件------------------------------------------
  ;lon1       fspan(0,360,4320)  ;gsn_xy    
  ;lat1     =fspan(90,-90,2160)  ;用于gsn_xy 地形图
  topo_dir  ="/mnt/d/data/other/"      ;数据、图片所在路径D:\ncl_related\data\other
  setfileoption("nc","MissingToFillValue",False)    ;设置缺省值 
  ;参考网址 http://www.ncl.ucar.edu/Document/Functions/Built-in/setfileoption.shtml
  elev_old= read_elev_data(topo_dir+"ETOPO5.DAT")
  elev_float=short2flt(elev_old)
  copy_VarCoords(elev_old, elev_float)
  elev_float=where(elev_float.ge.0., elev_float, elev_float@_FillValue)


;----------------------------------------------网格点插值----------------------------------------
  lon_out=fspan(85,105,241)
  lat_out=fspan(10, 35, 241)
  lon_out@units="degrees_east"
  lat_out@units="degrees_north"
   
  u     =linint2_Wrap(lon_in, lat_in(::-1), u_in(::-1,:), True, lon_out, lat_out, 0)
  delete(u_in)
  v     =linint2_Wrap(lon_in, lat_in(::-1), v_in(::-1,:), True, lon_out, lat_out, 0)
  delete(v_in)
  vvel     =linint2_Wrap(lon_in, lat_in(::-1), vvel_in(::-1,:), True, lon_out, lat_out, 0)
  delete(vvel_in)
  ;地形高度网格插值
  ; elev      = linint2_Wrap(elev_float&lon,elev_float&lat(::-1),elev_float(::-1,:),True,lon_out, lat_out,0)
  ; ;[LAT | 401] x [LON | 301]
  ; elev_lon=elev_float&lon
  ; elev_lat=elev_float&lat
  ; 高分辨率插值到低分辨率
  ; elev = area_hi2lores_Wrap (elev_lon,elev_lat(::-1),elev_float(::-1,:), True, 1, lon_out,lat_out, False)
  ; delete(elev_old)
  ; delete(elev_float)

  ;printVarSummary(u) 
  ;Dimensions and sizes:   [lv_ISBL0 | 26] x [LAT_0 | 401] x [LON_0 | 301]
  ;Coordinates:lv_ISBL0: [1000..100000] LAT_0: [ 0..40] LON_0: [80..110]
  ;printVarSummary(elev) 
  ;Dimensions and sizes:   [LAT | 401] x [LON | 301]
  ;Coordinates:LAT: [ 0..40] LON: [80..110]

;--------------------------数据处理: 地形强迫垂直速度------------------
  ;算地形梯度
  ;不对地形高度进行网格插值
  elev=elev_float({10:35},{85:105}) 
  nlat=dimsizes(elev&lat)
  mlon=dimsizes(elev&lon)
  lat    =elev&lat
  lon    =elev&lon

  ; nlat=dimsizes(elev&LAT)
  ; mlon=dimsizes(elev&LON)
  ; lat    =elev&LAT
  ; lon    =elev&LON

  dlon   = (lon(2)-lon(1))*0.0174533 ; convert to radians 0.017=3.14/180
  dHdX   = new ( dimsizes(elev),typeof(elev)); pre-allocate space
  copy_VarCoords(elev, dHdX)
  do nl=0,nlat-1                      ; loop over each latitude
  dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
  dHdX(nl:nl,:) = center_finite_diff_n (elev(nl:nl,:), dX , False,0,1) ; 1指的经度那一维
  end do

  dY   =(lat(2)-lat(1))*6378388.*0.0174533 ;标量
  ;一个纬度的长度，如果格点是1*1做差分时可以直接用这个dY，但如果不是1*1需要变一下
  dHdY = new ( dimsizes(elev),typeof(elev)); pre-allocate space
  copy_VarCoords(elev, dHdY)
  ; do nl=0,mlon-1                      ; loop over each latitude，其实可以不用循环
  ; dHdY(:,nl:nl) = center_finite_diff_n (elev(:,nl:nl), dY , False,0,0)
  ; end do
  dHdY= center_finite_diff_n (elev, dY , False,0,0)


;算速度
;不对地形高度进行网格插值
dHdY_cuo = dHdY({28.0},{92.0})
dHdX_cuo = dHdX({28.0},{92.0})
wf_cuo    = v({28.0},{92.0})*dHdY_cuo;+u({28.0},{92.0})*dHdX_cuo 
w_cuo    = vvel({28.0},{92.0})


dHdY_bomi= dHdY({29.9},{95.8})
dHdX_bomi= dHdX({29.9},{95.8})
wf_bomi   = v({29.9},{95.8})*dHdY_bomi;+u({29.9},{95.8})*dHdX_cuo 
w_bomi    = vvel({29.9},{95.8})


;试一试用10m风速和sigma的风速

;   ; wf=u*conform(u,dHdX,(/1,2/))+v*conform(v,dHdY,(/1,2/))
;   wf=v*conform(v,dHdY,(/1,2/))
;   copy_VarCoords(u, wf)
;   ;[lv_ISBL0 | 26] x [LAT_0 | 401] x [LON_0 | 301]
;   ; printMinMax(w1, 0);   min=-6.73754   max=21.4344
;   ; exit
; ;定点，错那，92.0E,28.0N;波密，95.8E,29.9N
;   w_cuo = wf({lv_ISBL0|70000},{LAT_0|26.5},{LON_0|92.0}) 
;   w_bomi = wf({lv_ISBL0|70000},{LAT_0|28.4},{LON_0|95.8}) 
;   dHdY_cuo=dHdY({26.5},{92.0})
;   dHdY_bomi=dHdY({28.4},{95.8})

  print(wf_cuo)
  print(w_cuo) 
  print(wf_bomi) 
  print(w_bomi) 
  ; print(dHdY_cuo)
  ; print(dHdY_bomi) 
  ; exit

end