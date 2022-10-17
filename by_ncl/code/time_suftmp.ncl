;绘制micaps地面气温随时间变化图
;文件读取方式 http://www.ncl.ucar.edu/Applications/list_io.shtml
;插值函数http://www.ncl.ucar.edu/Document/Functions/interp.shtml

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"


begin

;1.读取数据
  filepath = "/cygdrive/D/ncl_related/data/tmp/"  ;参数文件路径
  picture_path="/cygdrive/D/ncl_related/pictures/tmp/"
  data_cuo = asciiread(filepath+"surtmp55690.txt",-1,"string") ;以字符串形式读取参数文件入数组
  tmp_cuo  = stringtointeger(str_get_field(data_cuo,2," ")) ;从数组中获取温度值
  data_bo = asciiread(filepath+"surtmp56227.txt",-1,"string") ;以字符串形式读取参数文件入数组
  tmp_bo  = stringtointeger(str_get_field(data_bo,2," ")) ;从数组中获取温度值
  data_cha = asciiread(filepath+"surtmp56434.txt",-1,"string") ;以字符串形式读取参数文件入数组
  tmp_cha  = stringtointeger(str_get_field(data_cha,2," ")) ;从数组中获取温度值

  ;print(tmp_cuo)
  ;exit

;2.将变量整合到一个数组,并设置时间属性
  data      = new((/3,9/),integer)
  data(0,:) = tmp_cuo
  data(1,:) = tmp_cha
  data(2,:) = tmp_bo


  ;设置时间属性，从nc文件中获取
  path_rawdata = "/cygdrive/D/ncl_related/data/tmp/"
  nc_files = systemfunc("ls "+path_rawdata+"fnl_*_*_*.nc")
  f = addfiles(nc_files,"r") ; f[0-8],9个文件
  ListSetType(f,"join"); 打开多个文件并且返回文件句柄(是列表),指定如何返回文件变量
  f_time =  f[:]->TMP_P0_L100_GLL0@initial_time; 注意原始数据中有些变量如 LAND_P0_L1_GLL0 有时间属性
  ;print(f_time) ; 是每个文件的具体时间,具有单位, initial_time : 10/27/2008 (00:00)
  ;exit
  tmpdata= f[:]->$"TMP_P0_L100_GLL0"$ 
  tmpdata!0 = "time"

  ;;;>>>时间格式的修改
  ; tmpdata&time = f_time
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
  ; tmp的时间格式重新写为 ju格式
  tmpdata&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)


  data!1="time"
  data&time=tmpdata&time
  ;print(data&time)
  ;exit 




wks = gsn_open_wks("png", picture_path+"micap_suf_tmp")


 ;------------------------------绘图属性设置------------------------------------------
res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize=True
  res@tmXTOn = False ;关刻度
  res@tmYROn = False
  res@tiYAxisString="T (~S~o~N~C)"
  res@tiYAxisFont=22 ;y轴标题字体
  res@tmXBLabelFont=22
  res@tmXBLabelAngleF = 48.   ;角度
  res@tmYLLabelFont=22
  res@tmXBLabelFontHeightF = 0.015 ; x轴标签字号
  res@tmYLLabelFontHeightF =0.015

  res@trYMaxF=max(data)+12 ;扩大y轴便签显示的最大值


  ;设置画图区域大小
  res@vpWidthF=0.5
  res@vpHeightF=0.4

  ;时间轴设置
  trestick = True
  trestick@ttmFormat = "%D/%H" ; 时间轴的格式
  trestick@ttmAxis   = "XB"
  trestick@ttmNumTicks = 9 ; 显示所有9个点


  Tres=res
  ;Tres@xyDashPatterns  = (/1,2,3/)                        ; Make curves all solid
  Tres@xyDashPattern =0
  Tres@xyLineThicknessF= 8
  Tres@xyMarkLineMode = "MarkLines"               ; Markers *and* lines
  Tres@xyMarker =16
  Tres@xyMarkerThicknessF=8
  Tres@xyLineColors =(/"blue","red","purple"/) ;cuo,cha,bo

  ;图例
  Tres@pmLegendDisplayMode    = "Always"                              ;add lengend
  Tres@xyExplicitLegendLabels = (/"Cuona","Chayu","Bomi"/)          ;lengend title by english
  Tres@pmLegendWidthF         = 0.11                                   ;lengend width
  Tres@pmLegendHeightF        = 0.1                                   ;lengend height

  Tres@lgLabelFontHeightF = 0.015

  ;放在图中
  Tres@pmLegendOrthogonalPosF  = -1.21                              ;
  Tres@pmLegendParallelPosF    = 0.865


;绘图
time_axis_labels(data&time,Tres,trestick) ; call the formatting procedure 
plot  = gsn_csm_xy (wks,data&time,data,Tres) ; create plot


draw(plot)
frame(wks)
end
