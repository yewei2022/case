;用fnl资料,风矢是箭头，高度场，散度
;改1个时次，16
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"
;代码路径  ncl /mnt/d/ncl_related/code/200hPa_div_arrow.ncl

begin
	;手动设置一个时间，带unit的 timeset_ymdh 
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
	 path_format = "%Y%N%D%H"   ; 台风路径文件中的时间形式
	 path_time  = cd_string(timeset_ymdh,path_format)

	;数据读取和处理
	data_dir  ="/mnt/d/data/fnl/"
	pic_dir   ="/mnt/d/ncl_related/pictures/200/"
	data      = addfile(data_dir+"fnl_"+file_time+".nc", "r")
	u_in      = data ->UGRD_P0_L100_GLL0({20000},:,:)  ;单位：m/s
	v_in      = data->VGRD_P0_L100_GLL0({20000},:,:)
	lon_in    = data->lon_0
	lat_in    = data->lat_0
	hgt_in    = data->HGT_P0_L100_GLL0({20000},:,:)  ;13={50000}  hgt 的单位：gpm

;读取风暴中心位置   
	fiTY = "/mnt/d/data/TC/bio042008.txt"
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
		if (i.eq.tctime_ind) then ;最方便
			stdlon=lon(i)
			stdlat=lat(i)
		end if
	end do

;---------------------------网格点插值--------------------
	lon_out=fspan(55,125,701)
	lat_out=fspan(0, 60, 601)
	u      =linint2_Wrap(lon_in, lat_in(::-1), u_in(::-1,:), True, lon_out, lat_out, 0)
	delete(u_in)
	v      =linint2_Wrap(lon_in, lat_in(::-1), v_in(::-1,:), True, lon_out, lat_out, 0)
	delete(v_in)
	hgt_out=linint2_Wrap(lon_in, lat_in(::-1), hgt_in(::-1,:), True, lon_out, lat_out, 0)
	delete(hgt_in)

hgt  = hgt_out/10 ;位势高度处理
copy_VarCoords(hgt_out, hgt) 
;di=uv2dvF_Wrap(u,v) ;散度计算，球面简谐波函数计算散度
di = uv2dv_cfd(u,v,lat_out,lon_out, 2)  ;中央差分计算散度
copy_VarCoords(u, di)
di = di*1.e5
; ;风速大小
; U=sqrt(u^2+v^2)
; copy_VarCoords(u, U)



wks = gsn_open_wks("png", pic_dir+file_time)
; wks = gsn_open_wks("png", pic_dir+"U_"+file_time)

res = True
res@gsnDraw     = False
res@gsnFrame    = False
res@gsnMaximize =True
res@gsnLeftString  = ""
res@gsnRightString = ""
res@tmYROn         =False
res@tmXTOn         =False
;mpres@pmTickMarkDisplayMode  = "Always"
res@tmBorderThicknessF  =4    ; border 6x as thick,图框粗细
res@tmXBMajorThicknessF =4
res@tmXBMinorThicknessF =4
res@tmYLMajorThicknessF =4 
res@tmYLMinorThicknessF =4   
; ;设置坐标轴标签粗细
res@tmEqualizeXYSizes   =True 
; 现在一个坐标轴设置好后，
;(tmXBLabelFontHeightF,(tmXBMajorLengthF,(tmXBMinorLengthF
; 其他坐标轴就不用重复设置了
res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号


;设置地图
mpres = res
mpres@mpFillOn                = False
mpres@mpMinLatF               = 20.
mpres@mpMaxLatF               = 40.
mpres@mpMinLonF               = 70.
mpres@mpMaxLonF               = 105.
mpres@gsnMajorLatSpacing    =5. 
mpres@gsnMinorLatSpacing    =1
mpres@gsnMajorLonSpacing    =5. 
mpres@gsnMinorLonSpacing    =1. 
mpres@mpDataSetName              = "/mnt/d/ChinaMap/Earth..4"
mpres@mpDataBaseVersion          = "MediumRes"
mpres@mpAreaMaskingOn            = True
mpres@mpMaskAreaSpecifiers       = (/"China"/) ;这个和上个两个命令好像是没啥用
mpres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)


;--------------矢量图--------------

windRes = res
;windRes@vcGlyphStyle            = "FillArrow"   ; turn on curly vectors
windRes@vcGlyphStyle    = "LineArrow"
windRes@vcLineArrowThicknessF = 3
windRes@vcMinDistanceF = 0.04
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



;--------------散度--------------

cRes = res
cRes@cnFillOn = True ;cn means ContourPlot
cRes@cnLinesOn = False   
cRes@cnSpanFillPalette = True
cRes@cnFillPalette = "NEO_div_vegetation_c"
cRes@cnLineLabelsOn = False  
cRes@cnInfoLabelOn = False ; set false will not set an information label 
cRes@cnLevelSelectionMode = "ManualLevels" 
cRes@cnMinLevelValF       = -8             ; set min contour level
cRes@cnMaxLevelValF       = 8               ; set max contour level
cRes@cnLevelSpacingF      = 2              ; set contour spacing
;cRes@cnLineThicknessF = 4
;cRes@cnLevelSelectionMode = "ExplicitLevels"
;cRes@cnLevels = (/0, 5, 10, 15,20/)
;cRes@cnFillColors=(/120,150,180,210,240,250/)

cRes@lbBoxLinesOn = True ; Label bar
cRes@lbLabelAutoStride    = True            ; let NCL determine label spacing
; cRes@lbLabelFontHeightF   = 0.015
; cRes@lbOrientation        = "Vertical"         ; vertical label bar
; cRes@pmLabelBarWidthF = 0.03     
cRes@pmLabelBarHeightF = 0.1    
cRes@pmLabelBarOrthogonalPosF  = 0.15 ;色标离图的距离
cRes@lbBoxEndCapStyle = "TriangleBothEnds"


;--------------叠加高度场--------------

hgtRes = res
hgtRes@cnFillOn = False  ;cn means ContourPlot
hgtRes@cnLinesOn = True   
hgtRes@cnLineColor="blue"
hgtRes@cnLineLabelsOn = True  
hgtRes@cnInfoLabelOn = False ; set false will not set an information label
hgtRes@cnLineLabelFontHeightF=0.015
; hgtRes@cnLineLabelBackgroundColor=-1 ;背景变透明
; hgtRes@cnLabelMasking = True ；不穿过线

hgtRes@cnLevelSelectionMode = "ManualLevels" 
; hgtRes@cnMinLevelValF = 1100 
hgtRes@cnMaxLevelValF = 1250
hgtRes@cnLevelSpacingF = 5 
hgtRes@cnLineThicknessF = 6

; ;--------------叠加风速等值线，看急流位置--------------
; uRes = res
; uRes@cnFillOn    = False  ;cn means ContourPlot
; uRes@cnLinesOn   = True   
; uRes@cnLineColor = "red"
; uRes@cnLineLabelsOn       = True  
; uRes@cnInfoLabelOn        = False ; set false will not set an information label
; uRes@cnLevelSelectionMode = "ManualLevels" 
; uRes@cnMinLevelValF   =30 
; uRes@cnLevelSpacingF  =5 
; uRes@cnLineThicknessF = 3



;绘图
	;1 高度，散度，风矢量
	; hgt = smth9_Wrap(hgt,0.5,0.25,False)
	; hgt = smth9_Wrap(hgt,0.5,0.25,False)

	mapplot = gsn_csm_map(wks, mpres)
	contour_div = gsn_csm_contour(wks,di,cRes)
	windPlot = gsn_csm_vector(wks, u , v , windRes)
	contour_height = gsn_csm_contour(wks,hgt, hgtRes)
	; contour_U =gsn_csm_contour(wks,U, uRes)
	;叠加风矢量和位势高度
	overlay(mapplot, contour_div)
	overlay(mapplot, windPlot)
	overlay(mapplot, contour_height)
	; overlay(mapplot, contour_U)

	;2 某时刻风暴位置标记---------------------------------
	txres               = True                     ; text mods desired
	txres@txFontHeightF = 0.03                     ; font smaller. default big
	txres@txFontThicknessF=6
	txres@txFontColor="red"
	dum = gsn_add_text(wks,mapplot,"~F37~p",stdlon,stdlat,txres) 

	;3 高原轮廓线
	shp="/mnt/d/data/TPshp/DBATP/DBATP_Line.shp"  
	lnres                  = True     
	lnres@gsLineColor      ="orange"
	lnres@gsLineThicknessF = 8.0           ; 6x thickness  
	shp_plot               = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres) 

	;4 贴上强降水区域的方框
	ypts = (/32.0,  32.0,  28.0, 28.0, 32.0/)
	xpts = (/89.0,  98.0,  98.0, 89.0, 89.0/)

	resp                  = True                      ; polyline mods desired
	resp@gsLineColor      = "black"                     ; color of lines
	resp@gsLineThicknessF = 7.0                       ; thickness of lines
	dum_box = new(4,graphic)  
	; draw each line separately. Each line must contain two points.
	do i = 0 , 3
		dum_box(i)=gsn_add_polyline(wks,mapplot,xpts(i:i+1),ypts(i:i+1),resp)      
	end do 
draw(mapplot)
frame(wks)
end
