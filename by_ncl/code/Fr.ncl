; 某时次高度-湿弗劳德数图，v用850hpa-500hPa高度的每层风速
; 6.15 插值 linint2要求经纬度单增
;时间改1处
;参考https://www.ncl.ucar.edu/Document/Functions/Built-in/linint2.shtml
;代码路径  ncl /mnt/d/ncl_related/code/Fr.ncl

; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"



begin

;1.变量数据读入--------------------

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

	data_dir = "/mnt/d/data/fnl/"
	pic_dir ="/mnt/d/ncl_related/pictures/terrain/"  
	data    = addfile(data_dir+"fnl_"+file_time+".nc", "r")

	u_in  = data ->UGRD_P0_L100_GLL0 ;单位：m/s 等压面速度lv_ISBL0
	v_in  = data->VGRD_P0_L100_GLL0
	tmp_in= data->TMP_P0_L100_GLL0 
	rh_in = data->RH_P0_L100_GLL0   ;21层
	pres  = data->lv_ISBL4({40000:100000})  ;Pa
	pres_hpa     =pres*0.01
	pres_hpa@units="hPa"
	lon_in       =u_in&lon_0
	lat_in       =u_in&lat_0

;2 网格点插值------------
	lon_out=fspan(80,110,301)
	lat_out=fspan(0, 40, 401)
	lon_out@units="degrees_east"
	lat_out@units="degrees_north"   
	u_ipl_grid   =linint2_Wrap(lon_in, lat_in(::-1), u_in({40000:100000},::-1,:), True, lon_out, lat_out, 0)
	delete(u_in)
	v_ipl_grid   =linint2_Wrap(lon_in, lat_in(::-1), v_in({40000:100000},::-1,:), True, lon_out, lat_out, 0)
	delete(v_in) 
	tmp_ipl_grid =linint2_Wrap(lon_in, lat_in(::-1), tmp_in({40000:100000},::-1,:), True, lon_out, lat_out, 0)
	delete(tmp_in)
	rh_ipl_grid  =linint2_Wrap(lon_in, lat_in(::-1), rh_in({40000:100000},::-1,:), True, lon_out, lat_out, 0)
	delete(rh_in)

;3 气压层插值，间隔25hPa
	;int2p_n_Wrap 等压面插值； stdatmus_p2tdz 气压换算成高度
	pres_ipl=fspan(40000,100000,25) ;插值25个格点
	pres_ipl@units="Pa"
	linlog=2
	u = int2p_n_Wrap(pres,u_ipl_grid,pres_ipl,linlog,0)
	v = int2p_n_Wrap(pres,v_ipl_grid,pres_ipl,linlog,0)
	tmp = int2p_n_Wrap(pres,tmp_ipl_grid,pres_ipl,linlog,0)
	rh = int2p_n_Wrap(pres,rh_ipl_grid,pres_ipl,linlog,0)


;4 数据处理：湿Froude数----------------------------------------
  ;4.1 计算虚位温
 	pres_ipl_hpa      = pres_ipl*0.01
	pres_ipl_hpa@units= "hPa"
	pres_ipl_hpa!0="pres_ipl_hpa"
	pres_ipl_hpa&pres_ipl_hpa=pres_ipl_hpa
	w = mixhum_ptrh((conform(tmp,pres_ipl_hpa,0)),tmp, rh, 1) 
	; 1代表混合比，mixing ratio (kg/kg),要求气压mb
	tvw= tmp*(1.0+0.61*w)  ;虚温
	copy_VarCoords(tmp, tvw)
	theta_v=pot_temp(pres_ipl, tvw, 0, False)  ; 虚位温，Pa,K

  ;4.2 湿静力稳定度
	tdz1 = stdatmus_p2tdz(pres_ipl_hpa) ; 将气压(要求hPa)转换成高度 tdz(3,9) ,tdz(2,:)-height
	tdz  =tdz1(2,:)
	copy_VarCoords(pres_ipl, tdz) ;让变量和高度的坐标一致，均为Pa
	dtdz = center_finite_diff_n (theta_v,tdz,False,0,0)
	g=9.8
	Nw2=g*dtdz/theta_v
	Nw=sqrt(Nw2) 
	copy_VarCoords(tmp, Nw)

  ;4.3 湿弗劳德数
	copy_VarCoords(tmp, u)
	copy_VarCoords(tmp, v)
	hmax=conform_dims(dimsizes(tdz), 4000., 0) ;转换成和tdz一样的坐标变量
	h=hmax-tdz;(9)从上往下第9才为正值，对应625hPa
	Nw@_FillValue=9.e+10
	h@_FillValue=9.e+10
	h      = where(h.ge.0,h,h@_FillValue)	
	; print(pres_ipl_hpa(ind(ismissing(h))))
	; exit
	Nw_inv = 1. / where(Nw.ne.0,Nw,Nw@_FillValue)
	h_inv  = 1. / where(h.ne.0,h,h@_FillValue)
	Fr=v*Nw_inv*conform(v, h_inv, 0) 
	copy_VarCoords(tmp, Fr)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    

; 	; Fr_area= Fr(time|:,{LAT_0|20:25},{LON_0|90:98}) ;确定范围
; 	; Fr_ave= dim_avg_n_Wrap(Fr_area, (/1,2/)) ; 对lv_ISBL0做平均
	Fr_A=Fr({65000:100000},{26.5},{92.0}) 
	Fr_B=Fr({65000:100000},{28.4},{95.8}) ;16层，

	;4.4 将变量整合到一个数组
	Fr_plot  = new((/2,15/),float)
	Fr_plot(0,:) = Fr_A
	Fr_plot(1,:) = Fr_B
	pres_x=pres_ipl_hpa({650:1000})  ;7个，400-700
	pres_x@units="hPa"
	Fr_plot!1="pres_x" ;横坐标维
	Fr_plot&pres_x=pres_x
; print(Fr_plot(0,:))
; exit

wks  = gsn_open_wks ("x11",pic_dir+"H_Fr_"+file_time)          ; send graphics to PNG file

;5 绘图属性设置
  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  ;res@gsnMaximize=True
  ;设置坐标轴
  res@tmXTOn = False ;关刻度
  res@tmYROn = False
  ; res@tmYLMinorOn=False 打开次要刻度线
  res@tmYLMinorPerMajor=3 ;次要刻度线数目
  res@tiYAxisString="P (hPa)"
  res@tiXAxisString="Fr~B~m~N~" ; m s~S~-1~N~
  res@tmBorderThicknessF  =4    ; border 6x as thick,图框粗细
  res@tmXBMajorThicknessF =4
  res@tmXBMinorThicknessF =4
  res@tmYLMajorThicknessF =4 
  res@tmYLMinorThicknessF =4   
  ; res@tmEqualizeXYSizes   =True 
  ; res@tmXBLabelFontHeightF = 0.02 
  ;设置画图区域大小
  res@vpWidthF=0.4
  res@vpHeightF=0.6

  xyres=res
  xyres@xyDashPattern   = 0                         ; Make curves all solid
  xyres@xyLineThicknessF= 8
  xyres@xyMarkLineMode  = "Lines"               ; "MarkLines" Markers *and* lines
  ; xyres@xyMarker        = 16               ; 2 different markers
  ; xyres@xyMarkerColors  = (/"blue","purple"/) ; 2 different colors
  xyres@xyLineColors    = (/"blue","purple"/)
  xyres@trYReverse      = True
  ; xyres@trXMinF         = -0.2
  xyres@trYMinF         = 625
  xyres@trXMinF         =-0.5
  ;添加辅助线
  xyres@gsnXRefLine=1
  xyres@gsnXRefLineDashPattern=4
  xyres@gsnXRefLineThicknessF=6
  xyres@gsnXRefLineColor  ="black"




;6 绘图

plot  = gsn_csm_xy (wks,Fr_plot,Fr_plot&pres_x,xyres) ; create plot

;7 绘制图例
  resLg = True
  colours=(/"blue","purple"/)
  resLg@lgItemType = "Lines"
  ; resLg@lgMonoMarkerIndex = True
  ; resLg@lgMarkerColors = colours
  ; resLg@lgMarkerIndex = 1
  ; resLg@lgMarkerSizeF = 0.04
  resLg@lgMonoDashIndex = True
  resLg@lgDashIndex = 0
  resLg@lgLineColors = colours
  resLg@lgLineThicknessF = 6
  resLg@vpWidthF = 0.16
  resLg@vpHeightF = 0.1
  resLg@lgPerimOn   =False ;关掉图例边框
  ; resLg@lgPerimFill = 0
  ; resLg@lgPerimFillColor = "Background"
  resLg@lgLabelFontHeightF = 0.07
  ;resLg@lgTitleFontHeightF = 0.015
  ;resLg@lgTitleString = ""
  lbid = gsn_create_legend(wks, 2, (/"A","B"/), resLg)
 
  ; 将图例放置在图中
  amres = True
  amres@amParallelPosF   = 0.30
  amres@amOrthogonalPosF = 0.37
  dumLg = gsn_add_annotation(plot, lbid, amres)

draw(plot)
frame(wks)
end