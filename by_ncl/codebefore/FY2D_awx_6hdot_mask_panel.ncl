;TBB填色，其后6h降水(剔除了≤0.1mm以及高原以外站点降水)，;高原轮廓线，3张组图；
; 某日风暴路径
;写函数，若需要多次调用的，最好封装一下
;参考网址http://www.ncl.ucar.edu/Support/talk_archives/2012/0899.html
;labelbar的标标题设置参考 http://www.ncl.ucar.edu/Document/Graphics/Resources/lb.shtml#lbTitlePosition
;代码路径  ncl /mnt/d/ncl_related/code/FY2D_awx_6hdot_mask_panel.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("read_tbb_data")
	function read_tbb_data(tbb_file)
	local nlat, nlon,tbb, tbb_file, latt,lonn,tbb_1d,ind_notValid
	begin
	;---Read data as a straight binary file
	nlat           =1201
	nlon           =1201
	setfileoption("bin", "ReadByteOrder", "LittleEndian")
	tbb            = fbindirread (tbb_file,0,(/nlat,nlon/), "integer")
	tbb!0          ="lat"
	tbb!1          ="lon"
	latt            =fspan(-60.0, 60.0, 1201)
	lonn            =fspan(26.5, 146.5, 1201)
	latt@units      = "degrees_north"
	lonn@units      = "degrees_east"

	;---- Add meta data for TBB.
	tbb&lat        = latt
	tbb&lon        = lonn
	tbb@units      = "degree Kelvin"
	tbb@long_name  = "Temperature of Bright Blackbody"
	tbb@coordinates= "lat lon"
	tbb@_FillValue = -999
	tbb            = tbb-273        ;(convert degree Kelvin to degree Celsius)
	tbb@units      = "degree Celsius"

	;---- Replace data outside the valid range with missing value.      
	tbb_1d                 = ndtooned(tbb) ;二维变一维
	ind_notValid           = ind(tbb_1d.gt.-32)
	tbb_1d(ind_notValid)   = tbb@_FillValue
	tbb                    = onedtond(tbb_1d,dimsizes(tbb))

	return(tbb)
	end

undef("rain_pick")
	function rain_pick(R_old,lat,lon,mrb_lat,mrb_lon,npts)
	begin
	R = new(dimsizes(R_old),typeof(R_old),R_old@_FillValue)
  do i=0,npts-1
		R(i)=where(gc_inout(lat(i),lon(i),mrb_lat,mrb_lon), R_old(i), R_old@_FillValue)
  end do
  R=where(R.gt.0.1, R, R@_FillValue) ;剔除了小于等于0.1mm的降水

  return(R)
  end

undef("mark_lat_lon")
	function mark_lat_lon(num_distinct_markers,R,levels,lat,lon,npts,labels)
	local num_distinct_markers,R,levels,lat,lon,npts,labels,i,indexes,\
	npts_range,lat_new,lon_new
	; Create X and Y arrays to hold the points for each range and initialize
	; them to missing values.  We want to use num_distinct_markers
	; different colors, so we need num_distinct_markers sets of X and
	; Y points.
	begin
	; Group the points according to which range they fall in. At the
	; same time, create the label that we will use later in the labelbar  
	lat_new = new((/num_distinct_markers,npts/),float,-999);设置成二维，标记*站点数
	lon_new = new((/num_distinct_markers,npts/),float,-999)
	do i = 0, num_distinct_markers-1
		if (i.eq.0) then ;第一个范围
			indexes = ind(R.lt.levels(0)) ;定位，且小于第一个level的点会画出来
			labels(i) = "x < " + levels(0)
		end if
		if (i.eq.num_distinct_markers-1) then ;最后一个范围
			indexes = ind(R.ge.max(levels))
			labels(i) = "x >= " + max(levels)
		end if
		if (i.gt.0.and.i.lt.num_distinct_markers-1) then    ;中间的范围，因为格式不一样，要分别设置
			indexes = ind(R.ge.levels(i-1).and.R.lt.levels(i))
			labels(i) = levels(i-1) + " <= x < " + levels(i)
		end if
	 ; Now that we have the set of indexes whose values fall within 
	 ; the given range, take the corresponding lat/lon values and store
	 ; them, so later we can color this set of markers with the appropriate color.
		if (.not.any(ismissing(indexes))) then 
			npts_range = dimsizes(indexes)   ; # of points in this range.
			lat_new(i,0:npts_range-1) = lat(indexes)
			lon_new(i,0:npts_range-1) = lon(indexes)
		end if
		delete(indexes)            ; Necessary b/c "indexes" may be a different
		; size next time.
	end do

	return([/lat_new,lon_new/])
 end

undef("plot_scatter")
	procedure plot_scatter(wks,size,thick,pmid,num_distinct_markers,\
												lat_new,lon_new,gsres,mapplot)
	local i
	begin
	do i = 0, num_distinct_markers-1
		if (.not.ismissing(lat_new(i,0)))
			gsres@gsMarkerSizeF      = size(i)
			gsres@gsMarkerThicknessF = thick(i)
			pmid(i) = gsn_add_polymarker(wks,mapplot,lon_new(i,:),lat_new(i,:),gsres)      
		end if
	end do
	mapplot@$unique_string("pmid")$ = pmid
	end

undef("plot_track")
	procedure plot_track(wks,mapplot,lat_tc,lon_tc,plot_nrow,date)
	local i,plot_nrow,dumline,dumdot,dumdate,xx,yy
	begin

	;画线
	resline = True
	resline@gsLineThicknessF = 6
	resline@gsLineColor = "red"
	dumline = new(plot_nrow, graphic) 
	i=0
	  do i = 0, plot_nrow-2
		xx = (/ lon_tc(i), lon_tc(i+1)/) ;从第8个点开始，即2600时开始
		yy = (/ lat_tc(i), lat_tc(i+1)/) 
		dumline(i) = gsn_add_polyline(wks, mapplot, xx, yy, resline)
	  end do

	;画点
	resdot               = True                     ; text mods desired
	resdot@txFontHeightF = 0.024                     ; font smaller. default big
	resdot@txFontThicknessF=4
	resdot@txFontColor="red"
	dumdot = gsn_add_text(wks,mapplot,"~F37~p",lon_tc, lat_tc,resdot) 

	;标注日期
  dumdate = new(plot_nrow,graphic)
  resTx = True
  resTx@txFontHeightF = 0.02
  resTx@txFontColor = "black"
  resTx@txJust = "CenterLeft"
  dumdate(0)= gsn_add_text(wks,mapplot, date(0), \
                           lon_tc(0)+0.8, lat_tc(0)-0.4, resTx)
  dumdate(1:)= gsn_add_text(wks,mapplot, date(1:),\
                            lon_tc(1:)+0.8, lat_tc(1:), resTx)
	mapplot@$unique_string("dumline")$ = dumline
	mapplot@$unique_string("dumdot")$ = dumdot
	mapplot@$unique_string("dumdate")$ = dumdate

	;以上两行是为了解决报错
	;warning:TransformPostDraw: tfPolyDrawList element 0 is invalid
	;当函数需要重复调用时需要这样，为了保证dumline，dumdot作为局部变量在返回id后不会被删除
 
	end 


;----------------------------------------main code-------------------------
begin

pictures_dir   ="/mnt/d/ncl_related/pictures/TBB/"  

;风暴路径数据读取 
	fiTY = "/mnt/d/data/TC/bio042008.txt"
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
	vmax = stringtofloat(str_get_field(cdata, 9, ","))*0.5
	vmax_str = tostring_with_format(vmax,"%3.1f")
	mslp = stringtoint(str_get_field(cdata, 10, ","))
	DateChar = stringtochar(YYYYMMDDHH)
	MM = chartostring(DateChar(:,5:6))
	DD = chartostring(DateChar(:,7:8))
	HH = chartostring(DateChar(:,9:10))
	plot_lat_tc=lat_tc(7:)
	plot_lon_tc=lon_tc(7:)
	plot_date=DD(7:)+HH(7:)+","+vmax_str(7:)\
	+"m/s"+","+tostring(mslp(7:))+"hPa"
	plot_nrow=6 ;只画6个点
;读TBB并处理
	filename_tbb_dir="/mnt/d/data/TBB/filename_ncl.txt"
	filename_tbb    = asciiread(filename_tbb_dir,-1,"string") ;以字符串形式读取参数文件入数组argu 
	tbb             = new((/3,1201,1201/), integer)
	do i=0,2
	tbb(i,:,:)    = read_tbb_data(filename_tbb(i)) ;读第3个文件，即file_name的第二行
	end do
 
;读降水并处理
	;高原轮廓，用于mask高原以外区域
	tp_line = addfile("/mnt/d/data/TPshp/DBATP/DBATP_Polygon.shp", "r")
	mrb_lon = tp_line->x
	mrb_lat = tp_line->y

	rain_dir = "/mnt/d/data/precipitation/6h/rain12_py.txt"  ;数据
	f        = asciiread(rain_dir,-1,"string")
	station  = stringtofloat(str_get_field(f,1," ")) ;第1列中获取站号,列数的索引不用减1
	lat      = stringtofloat(str_get_field(f,2," ")) ;获取纬度值
	lon      = stringtofloat(str_get_field(f,3," ")) ;获取经度值
	R1    = stringtofloat(str_get_field(f,6," ")) ;获取2712的过去6h降水
	R2    = stringtofloat(str_get_field(f,8," ")) ;获取2712的过去6h降水
	R3    = stringtofloat(str_get_field(f,9," ")) ;获取2712的过去6h降水
	npts  = numAsciiRow(rain_dir) ;获取站点数(文本行数)，注意这里是文本rain_dir，而不是字符串f
	R_old =new((/3,npts/),float)
	R_old(0,:)=R1 
	R_old(1,:)=R2
	R_old(2,:)=R3
	R_old@_FillValue =  9999.
	R_old=where(R_old.gt.0.1, R_old, R_old@_FillValue) ;剔除了小于等于0.1mm的降水
	R      = new(dimsizes(R_old),typeof(R_old),R_old@_FillValue)
	levels = (/0.,5.,10.,20.,30.,40.,50./)
	nlevels= dimsizes(levels)
	num_distinct_markers = nlevels+1        ; number of distinct markers
	labels   = new(num_distinct_markers,string)  ; Labels for legend. 
	lat_new3 = new((/3,num_distinct_markers,npts/),float,-999);设置成二维，标记*站点数
	lon_new3 = new((/3,num_distinct_markers,npts/),float,-999)
	do i= 0,2
	  R(i,:)=rain_pick(R_old(i,:),lat,lon,mrb_lat,mrb_lon,npts)
		ref   =mark_lat_lon(num_distinct_markers,R(i,:),levels,lat,lon,npts,labels)
		lat_new3(i,:,:)=ref[0]
		lon_new3(i,:,:)=ref[1]
	end do

;---- Begin to plot.
wkstype="pdf"
wks = gsn_open_wks(wkstype,pictures_dir+"panel")

res = True
 res@gsnDraw                        = False
 res@gsnFrame                       = False
 ;res@gsnMaximize                    = True
 ;res@gsnPaperOrientation            = "portrait"
 ;res@gsnAddCyclic                   = False  ; regional data, don't add longitude cyclic point
 ;res@pmTickMarkDisplayMode          = "Always" ;关掉不需要的刻度,次刻度线和右上刻度线
 res@gsnLeftString=""
 res@gsnRightString=""
 res@tmYROn=False
 res@tmXTOn=False
 res@tmBorderThicknessF  =2.5    ; border 6x as thick,图框粗细
 res@tmXBMajorThicknessF =2.5
 res@tmXBMinorThicknessF =2.5
 res@tmYLMajorThicknessF =2.5 
 res@tmYLMinorThicknessF =2.5   

;设置坐标轴标签粗细，图例在左时候设置的
	res@tmEqualizeXYSizes   =True 
	; 上个命令，(tmXBLabelFontHeightF,(tmXBMajorLengthF,(tmXBMinorLengthF设置好后
	; 其他坐标轴就不用重复设置了
	res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号

 ;绘制省界
 mpres                                =res
 mpres@mpFillOn                       = False
 mpres@mpDataSetName                  = "Earth..4"
 mpres@mpDataBaseVersion              = "MediumRes"
 mpres@mpOutlineSpecifiers            = (/"China","China:Provinces"/)
 mpres@mpNationalLineColor            = "black"  ;我国边界
 mpres@mpProvincialLineColor          = "black"  ;省界
 mpres@mpGeophysicalLineColor         = "black" ;海陆边界
 mpres@mpNationalLineThicknessF       = 1  ;边界加粗
 mpres@mpProvincialLineThicknessF     = 1
 mpres@mpGeophysicalLineThicknessF    = 1
 ;绘制国界
 mpres@mpOutlineOn = True
 mpres@mpOutlineBoundarySets = "National"

 ;地图边界
 mpres@mpMinLonF                    = 85
 mpres@mpMaxLonF                    = 104
 mpres@mpMinLatF                    = 17
 mpres@mpMaxLatF                    = 36
 mpres@gsnMajorLonSpacing    =5. 
 mpres@gsnMinorLonSpacing    =1.    
 mpres@gsnMajorLatSpacing    =5. 
 mpres@gsnMinorLatSpacing    =1


;tbb填图
	tbbres                                = res
	tbbres@cnFillOn                       = True
	tbbres@cnFillMode                     = "RasterFill" ;矩形栅格填充
	tbbres@cnLinesOn                      = False
	tbbres@cnLineLabelsOn                 = False
	tbbres@cnInfoLabelOn                  = False
	;----颜色设置
	tbbres@cnFillPalette= "MPL_cool" 
	tbbres@cnLevelSelectionMode = "ExplicitLevels"
	tbbres@cnLevels = (/-32,-40,-52,-60,-70,-80/) 
	;tbbcolors=ispan(0, 125, 18)
	tbbres@cnFillColors=(/0,18,36,54,90,108,-1/)
	;色标设置
	tbbres@lbLabelBarOn                   = False

	;降水散点图的属性设置
	gsres               = True
	gsres@gsMarkerIndex = 16          ; Use filled dots for markers.
	gsres@gsMarkerColor = "black"
	; Loop through each grouping of markers, and draw them one set at
	; a time, assigning the proper color and size with gsn_marker.
	size=new(num_distinct_markers,float) ;将size存起来制作图例
	thick=new(num_distinct_markers,float)
	base_size = 0.01
	do i = 0, num_distinct_markers-1
		size(i)=base_size * (i+1)/8.
		thick(i)=0.5*(i+1)
	end do

	;高原轮廓属性设置
	shp="/mnt/d/data/TPshp/DBATP/DBATP_Line.shp"  
	lnres                  = True     
	lnres@gsLineColor      ="orange"
	lnres@gsLineThicknessF = 8.0           ; 6x thickness
;------------------------绘图-----------------------
mapplot=new(3,graphic)
tbb_plot=new(3,graphic)
shp_plot=new(3,graphic)
do iplot=0,2
	;1 地图，tbb填色
	mapplot(iplot) = gsn_csm_map(wks, mpres)
	tbb_plot(iplot) = gsn_csm_contour(wks,tbb(iplot,:,:),tbbres)
	overlay(mapplot(iplot), tbb_plot(iplot))

	;;2 降水散点图
		pmid = new(num_distinct_markers,graphic)
		plot_scatter(wks,size,thick,pmid,num_distinct_markers,\
													lat_new3(iplot,:,:),lon_new3(iplot,:,:),gsres,mapplot(iplot))
	;3 风暴路径
		plot_track(wks,mapplot(iplot),plot_lat_tc,plot_lon_tc,plot_nrow,plot_date)
	;4 高原轮廓线
		shp_plot(iplot) = gsn_add_shapefile_polylines(wks,mapplot(iplot),shp,lnres)
end do

;5 贴上强降水区域的方框
	ypts = (/32.0,  32.0,  28.0, 28.0, 32.0/)
	xpts = (/89.0,  98.0,  98.0, 89.0, 89.0/)

	resp                  = True                      ; polyline mods desired
	resp@gsLineColor      = "black"                     ; color of lines
	resp@gsLineThicknessF = 6.0                       ; thickness of lines
	dum_box = new(4,graphic)  
	; draw each line separately. Each line must contain two points.
	do i = 0 , 3
	  dum_box(i)=gsn_add_polyline(wks,mapplot(2),xpts(i:i+1),ypts(i:i+1),resp)      
	end do 

;6 组图
	pres=True
	pres@gsnFrame        = False  ; don't advance frame yet，以便后面使用gsn_legend_ndc   
	pres@gsnPanelLabelBar=True  	;共用一个色标
	; pres@gsnPanelCenter  = False   ; don't center ,3张图才能左对齐
	; pres@gsnMaximize           = True                 ; maximize in frame
	; pres@gsnPaperOrientation   = "Portrait"            ;portra和maximize必须同时设置才会生效
	;此处可以不必设置以上两行
	;------------------------------------------------------------
	;色标在右
	; pres@lbOrientation                  = "vertical"   	  
	; pres@pmLabelBarOrthogonalPosF       = 0.02 ;色标条与图之间的距离
	; pres@pmLabelBarParallelPosF         = 0.45 ;把lb上下移动，上正下负
	; pres@pmLabelBarWidthF               = 0.08     ;色标变窄
	; pres@pmLabelBarHeightF              = 0.6
	; pres@lbTitleFontHeightF             = 0.02

	;-----------------------------------------------------------------
	; 色标在下 
	 pres@pmLabelBarOrthogonalPosF       = -0.008 ;色标条与图之间的距离
	 pres@pmLabelBarParallelPosF         = 0.0 
	 pres@pmLabelBarWidthF               = 0.25     
	 pres@pmLabelBarHeightF              = 0.06
	 pres@lbTitleFontHeightF             = 0.01
	 pres@lbLabelFontHeightF             = 0.01
	;-----------------------------------------------------------------

	pres@lbBoxLinesOn                   = True ; 色标之间分界线
	pres@lbBoxEndCapStyle               = "TriangleBothEnds"
	pres@lbTitleOn                      = True
	pres@lbTitleString                  = "TBB(~S~o~N~C)"
	pres@lbTitleJust                    = "TopRight"
	pres@lbTitlePosition                = "Bottom"
	pres@gsnPanelBottom = 0.08   ;底部增加空白
	gsn_panel(wks, mapplot, (/3,1/), pres)
	; draw(mapplot) ；因为panel自动调用draw and frame				
	; frame(wks)

;7 绘制降水散点图图例                
; Manually specify location of legends,先调整内部间距，再整体平移
  xleg = (/0.15,0.29,0.43,0.57,0.15,0.29,0.43/)   ; Location of
  ; xtxt = (/0.22,0.32,0.42,0.52,0.22,0.32,0.42/)  ; legend markers
  xleg=xleg+0.12
  xtxt=xleg+0.06
  yleg = (/0.18,0.18,0.18,0.18,0.15,0.15,0.15/)   ; and text
  yleg=yleg-0.12
  ytxt = yleg   ; strings.

  mkres=True
  mkres@gsMarkerIndex=16

  txres               = True
  txres@txFontHeightF = 0.01

	size_lg  =new(num_distinct_markers-1,float) ;将size存起来制作图例
	thick_lg =new(num_distinct_markers-1,float)
	labels_lg= new(num_distinct_markers-1,string)  ; Labels for legend. 
	size_lg=size(1:)
	thick_lg=thick(1:)
	labels_lg=labels(1:)
	do i = 0,num_distinct_markers-2
		 mkres@gsMarkerSizeF = size_lg(i)
		 mkres@gsMarkerThicknessF =thick_lg(i)
		 gsn_polymarker_ndc(wks,xleg(i),yleg(i),mkres)
		 gsn_text_ndc      (wks,labels_lg(i),xtxt(i),ytxt(i),txres)
	end do
  frame(wks)                ; now advance the frame

end 
