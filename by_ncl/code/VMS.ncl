;用fnl资料,计算垂直风切变VMS
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"
;代码路径  ncl /mnt/d/ncl_related/code/VMS.ncl


undef("calculate_vms")
  function calculate_vms(min_lon,max_lon,min_lat,max_lat,vms_all)
  local vms_area,vms_1d,vms
  begin
  vms_area= vms_all({min_lat:max_lat},{min_lon:max_lon}) 
  vms_1d=ndtooned(vms_area) ;不用dim_avg_n_Wrap的原因：
  ;输入的维数要大于求平均的维数
  vms = avg(vms_1d)

  return(vms)
  end


begin

  fil_dir = "/mnt/d/data/fnl/"
  pic_dir ="/mnt/d/ncl_related/pictures/500/" 
  nc_files = systemfunc("ls "+fil_dir+"fnl_*_*_*.nc"); 
  f        = addfiles(nc_files,"r") ; f[0-11],12个文件
  ListSetType(f,"join"); 打开多个文件并且返回文件句柄(是列表),指定如何返回文件变量
  f_time   =  f[:]->LAND_P0_L1_GLL0@initial_time
  ; 注意原始数据中有些变量如 LAND_P0_L1_GLL0 有时间属性
  ; print(f_time) ; 是每个文件的具体时间,具有单位, initial_time : 10/27/2008 (00:00)
  u_in = f[:]->$"UGRD_P0_L100_GLL0"$ 
  v_in = f[:]->$"VGRD_P0_L100_GLL0"$ 
  ;第一个维度是 join 格式（列表）
  u_in!0 = "time" ;修改维度名
  v_in!1 = "lv"
  pres_in  = f[:]->$"lv_ISBL0"$  ;Pa,26层 9*26
  pres     = pres_in(0,:)  ;输出是一维[lv_ISBL0 | 26]

  lon_in = f[:]->$"lon_0"$
  lat_in = f[:]->$"lat_0"$
  lon = lon_in(0,:)  ;输出是一维[lv_ISBL0 | 26]
  lat = lat_in(0,:)  ;输出是一维[lv_ISBL0 | 26]

  ;修改第一维的时间格式
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
  u_in&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec,\
   "hour since 2008-10-26 00:00:00",0)
  ; print(tmp&time) ; [time | 9] x [lv_ISBL0 | 26] x [lat_0 | 181] x [lon_0 | 360]
  copy_VarCoords(u_in, v_in)
  path_format = "%Y%N%D%H"   ; 台风路径文件中的时间形式
  path_time   = cd_string(u_in&time,path_format)
  ; print(path_time)
  ; exit

;--------------------读取风暴中心位置---------------------------------------   
fiTY = "/mnt/d/data/TC/bio042008.txt"
; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
nrow = numAsciiRow(fiTY)
YYYYMMDDHH = new(nrow, "string")
lat_tc = new(nrow, "float")
lon_tc = new(nrow, "float")
cdata  = asciiread(fiTY, -1, "string")
YYYYMMDDHH = str_get_field(cdata, 3, ",")
lat_tc = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
lon_tc = stringtofloat(str_get_field(cdata, 8, ",")) *0.1

stdlon=new(12,float)
stdlat=new(12,float)
stdlon(0:5)=lon_tc(7:12)
stdlat(0:5)=lat_tc(7:12)
stdlon(6:)=lon_tc(12)
stdlat(6:)=lat_tc(12)

;------------------------计算风切变----------------------------------
u200 =u_in(:,{20000},:,:)
v200 =v_in(:,{20000},:,:)
u850 =u_in(:,{85000},:,:)
v850 =v_in(:,{85000},:,:)

vms_all=sqrt((u200-u850)^2+(v200-v850)^2)
copy_VarCoords_not_n(u_in, vms_all,1)
;[time | 12] x [lat_0 | 181] x [lon_0 | 360]

min_lon=new(12,float)
min_lat=new(12,float)
max_lon=new(12,float)
max_lat=new(12,float)
min_lon=stdlon-10
max_lon=stdlon+10
min_lat=stdlat-10
max_lat=stdlat+10

vms=new(12, float)
do i=0,11
  vms(i)=calculate_vms(min_lon(i),max_lon(i),min_lat(i),max_lat(i),vms_all(i,:,:))
end do
print(vms)
exit

; min_lat_tem=0.
; max_lat_tem=0.
; min_lon_tem=0.
; max_lon_tem=0.
; vms_area=new(12, float)
; do i =0,11
;   min_lon_tem=min_lon(i)
;   max_lon_tem=max_lon(i)
;   min_lat_tem=min_lat(i)
;   max_lat_tem=max_lat(i)  
;   vms_area(i) = vms_all(i, {min_lat_tem:max_lat_tem},{min_lon_tem:max_lon_tem}) 
; end do
; vms = dim_avg_n_Wrap(vms_area, (/1,2/))
; print(vms)
; exit

end