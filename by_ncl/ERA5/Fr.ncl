; 一审前的图
; 利用ERA5数据 某时次高度-湿弗劳德数图，v用850hpa-500hPa高度的每层风速
;时间改1处
;参考https://www.ncl.ucar.edu/Document/Functions/Built-in/linint2.shtml
;  ncl /mnt/h/d/ncl_related/ERA5_code/Fr.ncl

; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"



begin

;1.变量数据读入--------------------

	data_dir = "/mnt/h/d/data/ERA5/"   
	pic_dir ="/mnt/h/d/ncl_related/pictures/terrain/ERA5/"  
	data   = addfile(data_dir+"ERA5.nc", "r")  
	time    =data->time
  YYYYMMDDHH_f=cd_calendar(time,-3)

  need_time=2008102706
  time_locat = ind(YYYYMMDDHH_f.eq.need_time) ;找该时间所处位置
  path_time =tostring(need_time)
	pres_hpa  = data->level({400:1000}) ;hPa 18个点
	pres     =pres_hpa*100
	pres@units="Pa"
	u_short  = data ->u(time_locat,{400:1000},{10:40},{85:100}) 
	v_short  = data ->v(time_locat,{400:1000},{10:40},{85:100})
	tmp_short= data ->t(time_locat,{400:1000},{10:40},{85:100})
	rh_short = data ->r(time_locat,{400:1000},{10:40},{85:100})
	u_in     =short2flt(u_short)
	v_in     =short2flt(v_short)
	tmp_in   =short2flt(tmp_short)
	rh_in    =short2flt(rh_short)
	delete([/u_short,v_short,tmp_short,rh_short/])

;2 气压层插值，间隔25hPa
	pres_ipl=fspan(400,1000,25) ;插值25个格点
	pres_ipl@units="hPa"
	pres_ipl!0="pres_ipl"
	pres_ipl&pres_ipl=pres_ipl
	linlog=2
	u = int2p_n_Wrap(pres_hpa,u_in,pres_ipl,linlog,0)
	v = int2p_n_Wrap(pres_hpa,v_in,pres_ipl,linlog,0)
	tmp = int2p_n_Wrap(pres_hpa,tmp_in,pres_ipl,linlog,0)
	rh = int2p_n_Wrap(pres_hpa,rh_in,pres_ipl,linlog,0)


 ; tdz_cha = stdatmus_p2tdz(pres_ipl)
 ; print(pres_ipl)
 ; print(tdz_cha(2,:))
 ; exit

;3 数据处理：湿Froude数----------------------------------------
  ;3.1 计算虚位温
	w = mixhum_ptrh((conform(tmp,pres_ipl,0)),tmp, rh, 1) 
	; 1代表混合比，mixing ratio (kg/kg),要求气压mb
	tvw= tmp*(1.0+0.61*w)  ;虚温
	copy_VarCoords(tmp, tvw)
	theta_v=pot_temp(pres_ipl, tvw, 0, False)  ; 虚位温，Pa,K

  ;3.2 湿静力稳定度
	tdz1 = stdatmus_p2tdz(pres_ipl) ; 将气压(要求hPa)转换成高度 tdz(3,9) ,tdz(2,:)-height
	tdz  =tdz1(2,:)
	copy_VarCoords(pres_ipl, tdz) ;让变量和高度的坐标一致，均为Pa
	dtdz = center_finite_diff_n (theta_v,tdz,False,0,0)
	g=9.8
	Nw2=g*dtdz/theta_v
	Nw=sqrt(Nw2) 
	copy_VarCoords(tmp, Nw)

  ;3.3 湿弗劳德数
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

  lat_A=26.5 ;27.98 ;错那
  lon_A=92.0 ;91.95
  lat_B=28.4 ;29.87;波密
  lon_B=95.8 ;95.77
  lat_C=27.2 ;28.65 ;察隅
  lon_C=97.5 ;97.47

  ; lat_A=27.0 ;27.98 ;错那
  ; lon_A=92.0 ;91.95
  ; lat_B=28.9 ;29.87;波密
  ; lon_B=95.8 ;95.77
  ; lat_C=27.7 ;28.65 ;察隅
  ; lon_C=97.5 ;97.47

	Fr_A=Fr({625:1000},{lat_A},{lon_A}) 
	Fr_B=Fr({625:1000},{lat_B},{lon_B}) 
  Fr_C=Fr({625:1000},{lat_C},{lon_C}) 

	;3.4 将变量整合到一个数组
	Fr_plot  = new((/3,16/),float)
	Fr_plot(0,:) = Fr_A
	Fr_plot(1,:) = Fr_B
  Fr_plot(2,:) = Fr_C
	pres_x=pres_ipl({625:1000})  ;7个，400-700
	pres_x@units="hPa"
	Fr_plot!1="pres_x" ;横坐标维
	Fr_plot&pres_x=pres_x
; print(Fr_plot(0,:))
; exit

wks  = gsn_open_wks ("x11",pic_dir+"Fr_h_"+path_time)          ; send graphics to PNG file
cmap = read_colormap_file("wgne15")
colors=new((/3,4/), float)
; printVarSummary(cmap)
; exit
colors(0,:)=cmap(12,:)
colors(1,:)=cmap(6,:)
colors(2,:)=cmap(3,:)

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
  ; res@tmEqualizeXYSizes   =True 
  ; res@tmXBLabelFontHeightF = 0.02 
  ;设置画图区域大小
  res@vpWidthF=0.4
  res@vpHeightF=0.6

  xyres=res
  xyres@xyDashPattern   = 0                         ; Make curves all solid
  xyres@xyLineThicknessF= 4
  xyres@xyMarkLineMode  = "Lines"               ; "MarkLines" Markers *and* lines
  ; xyres@xyMarker        = 16               ; 2 different markers
  ; xyres@xyMarkerColors  = (/"blue","purple"/) ; 2 different colors
  xyres@xyLineColors    = colors
  xyres@trYReverse      = True
  ; xyres@trXMinF         = -0.2
  xyres@trYMinF         = 625
  xyres@trXMinF         =-0.5
  xyres@trXMaxF         =7
  ;添加辅助线
  xyres@gsnXRefLine=1
  xyres@gsnXRefLineDashPattern=4
  xyres@gsnXRefLineThicknessF=4
  xyres@gsnXRefLineColor  ="black"




;6 绘图

plot  = gsn_csm_xy (wks,Fr_plot,Fr_plot&pres_x,xyres) ; create plot

;7 绘制图例
  resLg = True
  colours=colors
  resLg@lgItemType = "Lines"
  ; resLg@lgMonoMarkerIndex = True
  ; resLg@lgMarkerColors = colours
  ; resLg@lgMarkerIndex = 1
  ; resLg@lgMarkerSizeF = 0.04
  resLg@lgMonoDashIndex = True
  resLg@lgDashIndex = 0
  resLg@lgLineColors = colours
  resLg@lgLineThicknessF = 4
  resLg@vpWidthF = 0.16
  resLg@vpHeightF = 0.1
  resLg@lgPerimOn   =False ;关掉图例边框
  ; resLg@lgPerimFill = 0
  ; resLg@lgPerimFillColor = "Background"
  resLg@lgLabelFontHeightF = 0.07
  ;resLg@lgTitleFontHeightF = 0.015
  ;resLg@lgTitleString = ""
  lbid = gsn_create_legend(wks, 3, (/"A","B","C"/), resLg)
 
  ; 将图例放置在图中
  amres = True
  amres@amParallelPosF   = 0.30
  amres@amOrthogonalPosF = 0.37
  dumLg = gsn_add_annotation(plot, lbid, amres)

draw(plot)
frame(wks)
end