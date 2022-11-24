;某时刻风暴位置，30min后TBB填色，读取其后高原6h降水(剔除了≤0.1mm)
;高原轮廓线，3张组图
;写函数，若需要多次调用的，最好封装一下
;参考网址http://www.ncl.ucar.edu/Support/talk_archives/2012/0899.html
;labelbar的标标题设置参考 http://www.ncl.ucar.edu/Document/Graphics/Resources/lb.shtml#lbTitlePosition
;代码路径  ncl /mnt/h/d/ncl_related/code/FY2D_awx_6hdot_mask_panel.ncl

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

; undef("rain_pick")
; 	function rain_pick(R_old,lat,lon,mrb_lat,mrb_lon,npts)
; 	begin
; 	R = new(dimsizes(R_old),typeof(R_old),R_old@_FillValue)
;   do i=0,npts-1
; 		R(i)=where(gc_inout(lat(i),lon(i),mrb_lat,mrb_lon), R_old(i), R_old@_FillValue)
;   end do

;   return(R)
;   end

undef("mark_lat_lon_phase")
	function mark_lat_lon_phase(num_distinct_markers,R,phase,levels,lat,lon,npts,labels)
	local num_distinct_markers,R,levels,lat,lon,npts,labels,i,indexes,\
	npts_range,lat_new,lon_new,phase_new
	; Create X and Y arrays to hold the points for each range and initialize
	; them to missing values.  We want to use num_distinct_markers
	; different colors, so we need num_distinct_markers sets of X and
	; Y points.
	begin
	; Group the points according to which range they fall in. At the
	; same time, create the label that we will use later in the labelbar  
	lat_new = new((/num_distinct_markers,npts/),float,-999);设置成二维，标记*站点数
	lon_new = new((/num_distinct_markers,npts/),float,-999)
	phase_new = new((/num_distinct_markers,npts/),integer,-999)
	do i = 0, num_distinct_markers-1
		if (i.eq.0) then ;第一个范围
			indexes = ind(R.le.levels(0)) ;定位，且小于等于第一个level的点会画出来,
			;但是其实没有小于等于0的点,因为都设成缺测了
			labels(i) = "x<" + levels(0)
		end if
		if (i.eq.num_distinct_markers-1) then ;最后一个范围
			indexes = ind(R.gt.max(levels))
			labels(i) = "x>" + max(levels) ;x > 
		end if
		if (i.gt.0.and.i.lt.num_distinct_markers-1) then    ;中间的范围，因为格式不一样，要分别设置
			indexes = ind(R.gt.levels(i-1).and.R.le.levels(i))
			labels(i) = levels(i-1) + "<x~F34~#~F~" + levels(i) ;< x < =
		end if
	 ; Now that we have the set of indexes whose values fall within 
	 ; the given range, take the corresponding lat/lon values and store
	 ; them, so later we can color this set of markers with the appropriate color.
		if (.not.any(ismissing(indexes))) then 
			npts_range = dimsizes(indexes)   ; # of points in this range.
			lat_new(i,0:npts_range-1) = lat(indexes)
			lon_new(i,0:npts_range-1) = lon(indexes)
			phase_new(i,0:npts_range-1) = phase(indexes)
		end if
		delete(indexes)            ; Necessary b/c "indexes" may be a different
		; size next time.
	end do

	return([/lat_new,lon_new,phase_new/])
 end

undef("plot_scatter")
	procedure plot_scatter(wks,size,thick,num_distinct_markers,\
												lat_new,lon_new,phase_new,mapplot,phase,npts)
	local i,j,gsres,pmid
	begin
	pmid = new((/num_distinct_markers,npts/),graphic)
	gsres               = True
	gsres@gsMarkerColor = "black"
	do i = 0, num_distinct_markers-1
		; if (.not.ismissing(lat_new(i,0)))
		gsres@gsMarkerSizeF      = size(i)
		gsres@gsMarkerThicknessF = thick(i)
		do j=0,npts-1
			if((.not.ismissing(phase_new(i,j))).and.phase(i,j).eq.2) then
				gsres@gsMarkerIndex=3
			else if((.not.ismissing(phase_new(i,j))).and.phase(i,j).eq.1)  then
				gsres@gsMarkerIndex=16
			else
				continue
			end if
			end if
			pmid(i,j) = gsn_add_polymarker(wks,mapplot,lon_new(i,j),lat_new(i,j),gsres)
			; else if ((.not.ismissing(lat_new(i,j))).and.(phase(j).eq.1)) then 
			; 	gsres@gsMarkerIndex = 16
			; 	pmid(i,j) = gsn_add_polymarker(wks,mapplot,lon_new(i,j),lat_new(i,j),gsres) 
			; else
			; 	continue
			; end if
			; end if
		end do   
		; end if
	end do
	mapplot@$unique_string("pmid")$ = pmid
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

begin

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
;读TBB并处理
	filename_tbb_dir="/mnt/h/d/data/TBB/filename_ncl.txt"
	filename_tbb    = asciiread(filename_tbb_dir,-1,"string") ;以字符串形式读取参数文件入数组argu 
	tbb             = new((/3,1201,1201/), integer)
	do i=0,2
	tbb(i,:,:)    = read_tbb_data(filename_tbb(i)) ;读第3个文件，即file_name的第二行
	end do
 
;读降水并处理
	; ;高原轮廓，用于mask高原以外区域
	; tp_line = addfile("/mnt/h/d/data/TPshp/DBATP/DBATP_Polygon.shp", "r")
	; mrb_lon = tp_line->x
	; mrb_lat = tp_line->y

	rain_dir = "/mnt/h/d/data/precipitation/6h/tp_rain12_py.txt"  ;数据
	f        = asciiread(rain_dir,-1,"string")
	station  = stringtofloat(str_get_field(f,1," ")) ;第1列中获取站号,列数的索引不用减1
	lat      = stringtofloat(str_get_field(f,2," ")) ;获取纬度值
	lon      = stringtofloat(str_get_field(f,3," ")) ;获取经度值
	R1    = stringtofloat(str_get_field(f,6," ")) ;获取2618的过去6h降水
	R2    = stringtofloat(str_get_field(f,8," ")) ;获取2706的过去6h降水
	R3    = stringtofloat(str_get_field(f,9," ")) ;获取2712的过去6h降水
	npts  = numAsciiRow(rain_dir) ;获取站点数(文本行数)，注意这里是文本rain_dir，而不是字符串f
	R_old =new((/3,npts/),float)
	R_old(0,:)=R1 
	R_old(1,:)=R2
	R_old(2,:)=R3
	R_old@_FillValue =  9999.
	R_old=where(R_old.gt.0.1, R_old, R_old@_FillValue) ;剔除了小于等于0.1mm的降水

	R=R_old ;不用挑选高原站点
	; R      = new(dimsizes(R_old),typeof(R_old),R_old@_FillValue)

	phase_dir = "/mnt/h/d/data/precipitation/6h/phase.txt"  ;数据
	phasefile = asciiread(phase_dir,-1,"string")
	phase =new((/3,npts/),integer)
	phase(0,:)=stringtoint(str_get_field(phasefile,4," ")) 
	phase(1,:)=stringtoint(str_get_field(phasefile,6," "))
	phase(2,:)=stringtoint(str_get_field(phasefile,7," "))	

	levels = (/0.,5.,10.,20.,30.,40.,50./)
	nlevels= dimsizes(levels)
	num_distinct_markers = nlevels+1        ; number of distinct markers
	labels   = new(num_distinct_markers,string)  ; Labels for legend. 
	lat_new3 = new((/3,num_distinct_markers,npts/),float,-999);设置成二维，标记*站点数
	lon_new3 = new((/3,num_distinct_markers,npts/),float,-999)
	phase_new3 = new((/3,num_distinct_markers,npts/),integer,-999)

	do i= 0,2
	  ; R(i,:)=rain_pick(R_old(i,:),lat,lon,mrb_lat,mrb_lon,npts)
		ref1   =mark_lat_lon_phase(num_distinct_markers,R(i,:),\
		                           phase(i,:),levels,lat,lon,npts,labels)
		lat_new3(i,:,:)=ref1[0]
		lon_new3(i,:,:)=ref1[1]
		phase_new3(i,:,:)=ref1[2]
	end do

;---- Begin to plot.
pic_dir   ="/mnt/h/d/ncl_related/pictures/TBB/"  
wks = gsn_open_wks("x11",pic_dir+"panel")

res = True
	res@gsnDraw                         = False
	res@gsnFrame                        = False
	;res@gsnMaximize                    = True
	;res@gsnPaperOrientation            = "portrait"
	;res@gsnAddCyclic                   = False  ; regional data, don't add longitude cyclic point
	;res@pmTickMarkDisplayMode          = "Always" ;关掉不需要的刻度,次刻度线和右上刻度线
	res@gsnLeftString=""
	res@gsnRightString=""
	res@tmYROn=False 
	res@tmXTOn=False 

	;设置坐标轴标签粗细，图例在左时候设置的
	res@tmEqualizeXYSizes   =True 
	; 上个命令，(tmXBLabelFontHeightF,(tmXBMajorLengthF,(tmXBMinorLengthF设置好后
	; 其他坐标轴就不用重复设置了
	res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
	res@tmXBLabels=(/"85~S~o~N~E","90~S~o~N~E","95~S~o~N~E","100~S~o~N~E","105~S~o~N~E"/)
	res@tmYLLabels=(/"20~S~o~N~N","25~S~o~N~N","30~S~o~N~N","35~S~o~N~N"/)


	;绘制省界
	mpres                                = res
	mpres@mpFillOn                       = False
	mpres@mpDataSetName              = "/mnt/h/d/ChinaMap/Earth..4"
	mpres@mpDataBaseVersion          = "MediumRes"
	mpres@mpAreaMaskingOn            = True
	mpres@mpMaskAreaSpecifiers       = (/"China"/) ;这和上两个命令好像是没啥用
	mpres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)
	;绘制国界
	mpres@mpOutlineBoundarySets = "National"

	;地图边界
	mpres@mpMinLonF                    = 85
	mpres@mpMaxLonF                    = 100
	mpres@mpMinLatF                    = 20
	mpres@mpMaxLatF                    = 35
	mpres@gsnMajorLonSpacing    =5. 
	mpres@gsnMinorLonSpacing    =1.    
	mpres@gsnMajorLatSpacing    =5. 
	mpres@gsnMinorLatSpacing    =1

	;tbb
	tbbres                                = res
	tbbres@cnFillOn                       = True
	; tbbres@cnFillMode                     = "RasterFill" ;矩形栅格填充
	;RasterFill这个命令画出来的ps图会被黑色填充，得去掉
	tbbres@cnLinesOn                      = False
	tbbres@cnLineLabelsOn                 = False
	tbbres@cnInfoLabelOn                  = False
	;----颜色设置
	tbbres@cnLevelSelectionMode = "ExplicitLevels"
	tbbres@cnLevels = (/-32,-40,-52,-60,-70,-80/) 
	; tbbres@cnFillPalette= "MPL_cool" 
	; tbbres@cnFillColors=(/0,18,36,54,90,108,-1/)
	tbbres@cnFillPalette= "perc2_9lev" 
	tbbres@cnFillColors=(/8,7,6,5,4,2,-1/)
	;色标设置
	tbbres@lbLabelBarOn                   = False

	; Loop through each grouping of markers, and draw them one set at
	; a time, assigning the proper color and size with gsn_marker.
	size=new(num_distinct_markers,float) ;将size存起来制作图例
	thick=new(num_distinct_markers,float)
	do i = 0, num_distinct_markers-1
		size(i)=0.02/8.*(sqrt(i*1.5)+(i+4)*0.08)
		thick(i)=0.8*sqrt(i+1)
	end do

	;高原轮廓属性设置
	shp="/mnt/h/d/data/TPshp/DBATP/DBATP_Line.shp"  
	lnres                  = True     
	lnres@gsLineColor      ="orange"
	lnres@gsLineThicknessF = 4.0           ; 6x thickness

;------------------------绘图-----------------------
mapplot=new(3,graphic)
tbb_plot=new(3,graphic)
shp_plot=new(3,graphic)
path_time=(/"2008102612","2008102700","2008102706"/)
do iplot=0,2
	;1 地图，tbb填色
	mapplot(iplot) = gsn_csm_map(wks, mpres)
	tbb_plot(iplot) = gsn_csm_contour(wks,tbb(iplot,:,:),tbbres)
	overlay(mapplot(iplot), tbb_plot(iplot))

	;2 降水散点图
	plot_scatter(wks,size,thick,num_distinct_markers,lat_new3(iplot,:,:),\
	        lon_new3(iplot,:,:),phase_new3(iplot,:,:),mapplot(iplot),phase_new3(iplot,:,:),npts)
	;3 某时刻风暴位置
	mark_tc(wks,YYYYMMDDHH,path_time(iplot),mapplot(iplot),nrow,lon_tc,lat_tc)
	;4 高原轮廓线
	shp_plot(iplot) = gsn_add_shapefile_polylines(wks,mapplot(iplot),shp,lnres)
end do

;5 贴上强降水区域的方框
	ypts = (/33.0,  33.0,  28.0, 28.0, 33.0/)
	xpts = (/89.0,  98.0,  98.0, 89.0, 89.0/)

	resp                  = True                      ; polyline mods desired
	resp@gsLineColor      = "black"                     ; color of lines
	resp@gsLineThicknessF = 3.0                       ; thickness of lines
	dum_box1 = new(4,graphic)  
	; draw each line separately. Each line must contain two points.
	do i = 0 , 3
	  dum_box1(i)=gsn_add_polyline(wks,mapplot(0),xpts(i:i+1),ypts(i:i+1),resp)      
	end do
	dum_box2 = new(4,graphic) 
	do i = 0 , 3
	  dum_box2(i)=gsn_add_polyline(wks,mapplot(1),xpts(i:i+1),ypts(i:i+1),resp)      
	end do
	dum_box3 = new(4,graphic) 
	do i = 0 , 3
	  dum_box3(i)=gsn_add_polyline(wks,mapplot(2),xpts(i:i+1),ypts(i:i+1),resp)      
	end do  

;6 组图
	pres=True
	pres@gsnFrame        = False  ; don't advance frame yet，以便后面使用gsn_legend_ndc   
	pres@gsnPanelLabelBar= True  	;共用一个色标
	; pres@gsnPanelCenter  = False   ; don't center ,3张图才能左对齐
	; pres@gsnMaximize           = True                 ; maximize in frame
	; pres@gsnPaperOrientation   = "Portrait"            ;portra和maximize必须同时设置才会生效
	;此处可以不必设置以上两行
	; ;------------------------------------------------------------
	; ; 色标在右
	; pres@lbOrientation                  = "vertical"   	  
	; pres@pmLabelBarOrthogonalPosF       = 0.02 ;色标条与图之间的距离
	; pres@pmLabelBarParallelPosF         = -0.025 ;把lb上下移动，上正下负
	; pres@pmLabelBarWidthF               = 0.06     ;色标变窄
	; pres@pmLabelBarHeightF              = 0.45
	; pres@lbTitleFontHeightF             = 0.01
	; pres@lbLabelFontHeightF             = 0.01

	; -----------------------------------------------------------------
	; 色标在下 
	 pres@pmLabelBarOrthogonalPosF       = 0.005 ;色标条与图之间的距离
	 pres@pmLabelBarParallelPosF         = -0.005
	 pres@pmLabelBarWidthF               = 0.23     
	 pres@pmLabelBarHeightF              = 0.05
	 pres@lbTitleFontHeightF             = 0.01
	 pres@lbLabelFontHeightF             = 0.01
	;-----------------------------------------------------------------

	pres@lbBoxLinesOn                   = True ; 色标之间分界线
	pres@lbBoxEndCapStyle               = "TriangleBothEnds"
	pres@lbTitleOn                      = True
	pres@lbTitleString                  = "TBB(~S~o~N~C)"
	pres@lbTitleJust                    = "CenterRight"
	pres@lbTitlePosition                = "Right"
	pres@lbTitleDirection               = "Across"
	pres@gsnPanelBottom = 0.08   ;底部增加空白
	pres@gsnPanelFigureStrings=(/"(a)26/12","(b)27/00","(c)27/06"/)
	pres@gsnPanelFigureStringsFontHeightF=0.01
	pres@amJust="topLeft"
	gsn_panel(wks, mapplot, (/3,1/), pres)
	; draw(mapplot) ；因为panel自动调用draw and frame				
	; frame(wks)

;7 绘制降水散点图图例                
; Manually specify location of legends,先调整内部间距，再整体平移
	xleg = (/0.15,0.23,0.31,0.15,0.23,0.31,0.15/)   ; Location of
	xleg=xleg+0.248
	xtxt=xleg+0.045
	yleg = (/0.18,0.18,0.18,0.16,0.16,0.16,0.14/)   ; and text
	yleg=yleg-0.1
	ytxt = yleg   ; strings.

	rainmkres=True
	rainmkres@gsMarkerIndex=16
	snowmkres=True
	snowmkres@gsMarkerIndex=3

	txres               = True
	txres@txFontHeightF = 0.01

	size_lg  =new(num_distinct_markers-1,float) ;将size存起来制作图例
	thick_lg =new(num_distinct_markers-1,float)
	labels_lg= new(num_distinct_markers-1,string)  ; Labels for legend. 
	size_lg=size(1:)
	thick_lg=thick(1:)
	labels_lg=labels(1:)
	do i = 0,num_distinct_markers-2
		 rainmkres@gsMarkerSizeF = size_lg(i)
		 rainmkres@gsMarkerThicknessF =thick_lg(i)
		 snowmkres@gsMarkerSizeF = size_lg(i)
		 snowmkres@gsMarkerThicknessF =thick_lg(i)
		 gsn_polymarker_ndc(wks,xleg(i),yleg(i),rainmkres)
		 gsn_polymarker_ndc(wks,xleg(i)+0.012,yleg(i),snowmkres)
		 gsn_text_ndc      (wks,labels_lg(i),xtxt(i),ytxt(i),txres)
	end do
  frame(wks)                ; now advance the frame

end 
