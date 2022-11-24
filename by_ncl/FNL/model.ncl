;某时刻风暴位置，30min后TBB填色，读取其后高原6h降水(剔除了≤0.1mm)
; 不用挑选高原站点 即不用mask
;高原轮廓线，3张组图
;写函数，若需要多次调用的，最好封装一下
;参考网址http://www.ncl.ucar.edu/Support/talk_archives/2012/0899.html
;labelbar的标标题设置参考 http://www.ncl.ucar.edu/Document/Graphics/Resources/lb.shtml#lbTitlePosition
;代码路径  ncl /mnt/h/d/ncl_related/code/model.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


undef("plot_track")
	procedure plot_track(wks,mapplot,lat_tc,lon_tc,nrow,MM,DD,HH)
	local i,plot_nrow,dumline,dumdot,dumdate,xx,yy
	begin

	;画线
	resLine = True
	resLine@gsLineThicknessF = 6

  resLine@gsLineColor="yellow"
	dumLine = gsn_add_polyline(wks, mapplot, lon_tc, lat_tc, resLine)

	; ;画点
	; resdot               = True                     ; text mods desired
	; resdot@txFontHeightF = 0.024                     ; font smaller. default big
	; resdot@txFontThicknessF=4
	; resdot@txFontColor="red"
	; dumdot = gsn_add_text(wks,mapplot,"~F37~p",lon_tc, lat_tc,resdot) 

	; ;标注日期
	; dumDate = new(nrow,graphic)
	; resTx = True
	; resTx@txFontHeightF = 0.02
	; resTx@txFontColor = "black"
	; resTx@txJust = "CenterLeft"

	; do i = 1, nrow-1
	; 	dumDate(i) = gsn_add_text(wks,mapplot, "2008"+MM(i)+DD(i)+HH(i), \
	; 	                          lon_tc(i)+0.5,lat_tc(i), resTx)
	; end do

	mapplot@$unique_string("dumLine")$ = dumLine
	; mapplot@$unique_string("dumdot")$ = dumdot
	; mapplot@$unique_string("dumDate")$ = dumDate

	;以上两行是为了解决报错
	;warning:TransformPostDraw: tfPolyDrawList element 0 is invalid
	;当函数需要重复调用时需要这样，为了保证dumline，dumdot作为局部变量在返回id后不会被删除
 
	end 

; Give a contour/map plot, retrieve the contour levels associated with it.
undef("get_contour_levels")
  function get_contour_levels(plot)
  begin
    getvalues plot@contour
      "cnLevels" : levels
  end getvalues
  return(levels)
  end

; Give a contour/map plot, retrieve the contour colors associated with the levels
undef("get_contour_colors")
  function get_contour_colors(plot)
  begin
    getvalues plot@contour
      "cnLineColors" : colors
  end getvalues
  return(colors)
  end

;----------------------------------------------------------------------
; Loop through each contour level, get the isolines for that level, and
; add them to the plot only if they fit in an area larger than the
; given X/Y range.
;----------------------------------------------------------------------
undef("add_isolines_by_range")
  function add_isolines_by_range(wks,plot_orig,plot_map,x_range_min,y_range_min)
  local plres, i, j, b, e, x, y, i, nlevels, iso
  begin
    levels = get_contour_levels(plot_orig)
    colors = get_contour_colors(plot_orig)
    nlevels = dimsizes(levels)

  ;---Loop through each level, get isolines, and add them to the plot.
  plres = True
  plres@gsLineThicknessF = 3
  do i = 0, nlevels-1
    iso := get_isolines(plot_orig@contour,levels(i))
    plres@gsLineColor = colors(i)
    do j = 0, iso@segment_count -1
      b = iso@start_point(j)
      e = b + iso@n_points(j) - 1
      y := iso(0,b:e)
      x := iso(1,b:e)
      if((max(y)-min(y)).lt.y_range_min.and.\
         (max(x)-min(x)).lt.x_range_min) then
        continue
      end if
      plot_map@$unique_string("isolines")$ = gsn_add_polyline(wks,plot_map,x,y,plres)
    end do
  end do

  return(plot_map)
  end



undef("mark_tc")
  procedure mark_tc(wks,YYYYMMDDHH,path_time,map,nrow,lon_tc,lat_tc)
  begin
  tctime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
  stdlon=0.
  stdlat=0.
  do i = 0, nrow-1
   if (i.eq.tctime_ind) then ;最方便
     stdlon=lon_tc(i)
     stdlat=lat_tc(i)
   end if
  end do
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.03                    ; font smaller. default big
  txres@txFontThicknessF=3
  txres@txFontColor="red"
  dum_tc = gsn_add_text(wks,map,"~F37~p",stdlon,stdlat,txres) 
  map@$unique_string("dum_tc")$ = dum_tc
  end




;高原三条等高线图例
  ;每张子图中都加图例，用函数方便
undef("add_topoline_label")
  procedure add_topoline_label(wks,plot)
  begin
  resLg = True
  resLg@lgItemType = "Lines" 
  resLg@lgMonoDashIndex = True
  resLg@lgPerimOn    =False
  resLg@lgDashIndex  = 0
  resLg@lgLineColors = (/"gray0","gray40","gray70"/)
  resLg@lgLineThicknessF = 3
  resLg@vpWidthF  = 0.18
  resLg@vpHeightF = 0.15
  resLg@lgPerimFill = 0
  resLg@lgPerimFillColor = "Background"
  resLg@lgLabelFontHeightF = 0.08
  resLg@lgTitleString = ""
  lbid = gsn_create_legend(wks, 3, (/"2 km",\
  "3 km","4 km"/), resLg)
 
 ; 将图例放置在图中
  amres = True
  amres@amParallelPosF = -0.36 ;左右，正右左负
  amres@amOrthogonalPosF = 0.34 ;上下，正下负上
  dumLg = gsn_add_annotation(plot, lbid, amres)

  tmpid1 = "anno"+unique_string("id")
  tmpid2 = "lege"+unique_string("id")
  plot@$tmpid1$ = dumLg
  plot@$tmpid2$ = lbid

  end




; ============================main code====================================

begin

pic_dir   ="/mnt/h/d/ncl_related/pictures/"  
path_time="20162706"

;读地形文件
  topo        = addfile("/mnt/h/d/data/other/ETOPO2v2c_f4.nc", "r")
  elev        = topo->z({15:50},{70:110}) 
  x           = topo->x({70:110}) 
  y           = topo->y({15:50}) 

  lat_topo    = y
  lon_topo    = x
  lat_topo!0     = "lat"
  lon_topo!0     = "lon"
  lat_topo@units = "degrees_north"
  lon_topo@units = "degrees_east"
  lat_topo&lat   = lat_topo
  lon_topo&lon   = lon_topo

  ;---Attach the coordinate arrays
  elev!0    = "lat"
  elev!1    = "lon"
  elev&lat  = lat_topo
  elev&lon  = lon_topo

;风暴路径数据读取 
	fiTY = "/mnt/h/d/data/TC/bio042008.txt"
	; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
	nrow = numAsciiRow(fiTY) 
	YYYYMMDDHH = new(nrow, "string")
	lat_tc = new(nrow, "float")
	lon_tc = new(nrow, "float")
	vmax = new(nrow, "float")
	mslp = new(nrow, "integer") 
	cdata = asciiread(fiTY, -1, "string")
	YYYYMMDDHH = str_get_field(cdata, 3, ",")
	lat_tc = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
	lon_tc = stringtofloat(str_get_field(cdata, 8, ",")) *0.1
	DateChar = stringtochar(YYYYMMDDHH)
	MM = chartostring(DateChar(:,5:6))
	DD = chartostring(DateChar(:,7:8))
	HH = chartostring(DateChar(:,9:10))

;---- Begin to plot.
wks = gsn_open_wks("eps",pic_dir+"model")

res = True
	res@gsnDraw                         = False
	res@gsnFrame                        = False
	;res@gsnMaximize                    = True
	;res@gsnPaperOrientation            = "portrait"
	res@gsnAddCyclic                   = False  ; regional data, don't add longitude cyclic point
	;res@pmTickMarkDisplayMode          = "Always" ;关掉不需要的刻度,次刻度线和右上刻度线
	res@gsnLeftString=""
	res@gsnRightString=""
	res@tmYROn=False 
	res@tmXTOn=False 

	; ;设置坐标轴标签粗细，图例在左时候设置的
	; res@tmEqualizeXYSizes   =True 
	; ; 上个命令，(tmXBLabelFontHeightF,(tmXBMajorLengthF,(tmXBMinorLengthF设置好后
	; ; 其他坐标轴就不用重复设置了
	; res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
	; res@tmXBLabels=(/"85~S~o~N~E","90~S~o~N~E","95~S~o~N~E","100~S~o~N~E"/)
	; res@tmYLLabels=(/"20~S~o~N~N","25~S~o~N~N","30~S~o~N~N","35~S~o~N~N"/)


	;绘制省界
	mpres                                = res
	mpres@mpFillOn                       = False
	mpres@mpOutlineOn                    = False
	mpres@mpDataSetName              = "/mnt/h/d/ChinaMap/Earth..4"
	mpres@mpDataBaseVersion          = "MediumRes"
	; mpres@mpAreaMaskingOn            = True
	; mpres@mpMaskAreaSpecifiers       = (/"China"/) ;这和上两个命令好像是没啥用
	; mpres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)
	;绘制国界
	; mpres@mpOutlineBoundarySets = "National"

	;地图边界
	mpres@tmXBLabelsOn=False
	mpres@tmYLLabelsOn=False
	mpres@mpMinLonF                    = 80
	mpres@mpMaxLonF                    = 100
	mpres@mpMinLatF                    = 10
	mpres@mpMaxLatF                    = 30
	mpres@pmTickMarkDisplayMode ="Always" ;坐标轴标签上添加符号度
	; mpres@gsnMajorLonSpacing    =5. 
	; mpres@gsnMinorLonSpacing    =1.    
	; mpres@gsnMajorLatSpacing    =5. 
	; mpres@gsnMinorLatSpacing    =1


;---高原3km,4km，5km等值线
  mpres1                    = mpres
  mpres1@cnFillOn           = False             ; turn off contour fill
  mpres1@cnLinesOn          = True
  mpres1@cnMonoLineColor    = False ;这个一定要打开才能设置不同的颜色
  mpres1@cnLevelSelectionMode="ExplicitLevels"
  mpres1@cnLevels           =(/2000,3000,4000/)
  mpres1@cnLineColors       =(/"gray0","gray40","gray70"/)
  ; mpres1@cnLineThicknessF   = 2. ;无效，在add lines那儿有效
  mpres1@cnLineLabelsOn     = False            ; turn off line labels
  mpres1@cnInfoLabelOn      = False            ; turn off info label

  ;对地形高度数据做平滑
  elev = smth9_Wrap(elev,0.5,0.25,False)
  elev = smth9_Wrap(elev,0.5,0.25,False)
  elev = smth9_Wrap(elev,0.5,0.25,False)
  elev = smth9_Wrap(elev,0.5,0.25,False)  


;------------------------绘图-----------------------
    ; map1 = gsn_csm_contour_map(wks,elev,mpres1)

    mapplot = gsn_csm_map(wks,mpres)

    ; lon_range = 4. ;切除了过小的闭合等值线，使not busy
    ; lat_range = 4.
    ; mapplot= add_isolines_by_range(wks,map1,\
    ;                                     mapplot,lon_range,lat_range)
    ; add_topoline_label(wks,mapplot)

	;3 某时刻风暴位置
	; mark_tc(wks,YYYYMMDDHH,path_time,mapplot,nrow,lon_tc,lat_tc)

	; ;4 高原轮廓线
	;高原轮廓属性设置
	; shp="/mnt/h/d/data/TPshp/DBATP/DBATP_Line.shp"  
	; lnres                  = True     
	; lnres@gsLineColor      ="orange"
	; lnres@gsLineThicknessF = 4.0           ; 6x thickness
	; shp_plot = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres)

  ; lnres                  =True
  ; lnres@gsFillColor      = "gray70" 
  ; shp_fill = "/mnt/h/d/data/TPshp/DBATP/DBATP_Polygon.shp"
  ; fill_plot              =gsn_add_shapefile_polygons(wks, mapplot, shp_fill, lnres)


; ;5 贴上强降水区域的方框
; 	ypts = (/33.0,  33.0,  28.0, 28.0, 33.0/)
; 	xpts = (/89.0,  98.0,  98.0, 89.0, 89.0/)

; 	resp                  = True                      ; polyline mods desired
; 	resp@gsLineColor      = "black"                     ; color of lines
; 	resp@gsLineThicknessF = 3.0                       ; thickness of lines
; 	dum_box3 = new(4,graphic)  
; 	; draw each line separately. Each line must contain two points.
; 	do i = 0 , 3
; 	  dum_box3(i)=gsn_add_polyline(wks,mapplot(2),xpts(i:i+1),ypts(i:i+1),resp)      
; 	end do  

	; 	单独绘制TC路径
	plot_track(wks,mapplot,lat_tc,lon_tc,nrow,MM,DD,HH)


  ; ;站点标记
  ; lon_sta=(/91.95, 95.77, 97.47/)
  ; lat_sta=(/27.98, 29.87, 28.65/)
  ; mkres               = True                     ; text mods desired
  ; mkres@txFontColor   ="black"
  ; mkres@txFontHeightF = 0.01
  ; mkres@txFontThicknessF=1.5
  ; dum_sta = gsn_add_text(wks,mapplot,"~F35~8",lon_sta(0:1),lat_sta(0:1),mkres) 

	; ; 每张图标上BoB
	; reslet = True
	; reslet@txFontHeightF = 0.024
	; reslet@txFontThicknessF=0.24
	; reslet@txFontColor = "black"
	; ; reslet@txJust = "CenterLeft"
	; dumcount= gsn_add_text(wks,mapplot,"BoB",\
	;                            91., 18., reslet)
	draw(mapplot)
  frame(wks)                ; now advance the frame

end 
