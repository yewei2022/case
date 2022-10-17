;读取站点过程累计 降水量 降雨量，降雪量，组图
;画降水彩色散点图(降水数据本身只包含高原站点)，剔除≤0.1的站点
;绘制地形3/4/5km等值线,并将等值线变稀疏（切除+平滑）
;文件读取方式 http://www.ncl.ucar.edu/Applications/list_io.shtml
;散点图参考 http://www.ncl.ucar.edu/Applications/polyg.shtml polyg_8_lbar.ncl
;代码路径  ncl /mnt/h/d/ncl_related/code/rainfall_snowfall_topo2_panel.ncl


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


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

;加彩色散点的图例,
;在panel中加公共图例(gsn_labelbar_ndc)和
;在子图中(gsn_create_labelbar+gsn_add_annotation)
undef("add_labelbar")
  procedure add_labelbar(wks,levels,colors)
  local lbres, labels
  begin
  nlevels = dimsizes(levels)
  ; Draw a labelbar
  lbres                    = True
  lbres@vpWidthF           = 0.36           ; width
  lbres@vpHeightF          = 0.06            ; height
  lbres@lbPerimOn          = False            ; Turn off perimeter.
  lbres@lbOrientation      = "Horizontal"     ; Default is vertical.
  lbres@lbLabelAlignment   = "InteriorEdges"  ; Default is "BoxCenters",标签如何与标签栏对齐
  lbres@lbFillColors       = colors           ; Colors for boxes.
  lbres@lbMonoFillPattern  = True             ; Fill them all solid.
  lbres@lbLabelFontHeightF = 0.012            ; label font height
  lbres@lbBoxLinesOn       =True
  lbres@lbBoxEndCapStyle   = "TriangleBothEnds"
  lbres@lbTitleString      ="(mm)"
  lbres@lbTitleFontHeightF = 0.012
  lbres@lbTitleJust        ="BottomRight"
  lbres@lbTitlePosition    ="Bottom"

  labels = sprintf("%2.0f",levels)

  gsn_labelbar_ndc(wks,nlevels+1,labels,0.342,0.09,lbres)
  frame(wks)                                   ; now advance the frame

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

;根据降水量将经纬度数组分类
undef("mark_lat_lon")
  function mark_lat_lon(num_distinct_markers,R,levels,lat,lon,npts)
  local num_distinct_markers,R,levels,lat,lon,npts,i,indexes,\
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
      indexes = ind(R.le.levels(0)) ;定位，且小于等于第一个level的点会画出来
    end if
    if (i.eq.num_distinct_markers-1) then ;最后一个范围
      indexes = ind(R.ge.max(levels))
    end if
    if (i.gt.0.and.i.lt.num_distinct_markers-1) then    ;中间的范围，因为格式不一样，要分别设置
      indexes = ind(R.ge.levels(i-1).and.R.lt.levels(i))
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

undef("plot_rain_scatter")
  procedure plot_rain_scatter(wks,colors,num_distinct_markers,\
                        lat_new,lon_new,mapplot)
  local i,gsres,pmid
  begin
  pmid = new(num_distinct_markers,graphic)
  gsres               = True
  gsres@gsMarkerIndex = 16  ;雨16 雪3          ; Use filled dots for markers.
  gsres@gsMarkerSizeF      = 0.015
  do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
      gsres@gsMarkerColor      = colors(i,:)
      pmid(i) = gsn_add_polymarker(wks,mapplot,lon_new(i,:),lat_new(i,:),gsres)      
    end if
  end do
  mapplot@$unique_string("pmid")$ = pmid
  end

undef("plot_snow_scatter")
  procedure plot_snow_scatter(wks,colors,num_distinct_markers,\
                        lat_new,lon_new,mapplot)
  local i,gsres,pmid
  begin
  pmid = new(num_distinct_markers,graphic)
  gsres               = True
  gsres@gsMarkerIndex = 3  ;雨16 雪3    
  gsres@gsMarkerSizeF      = 0.015
  gsres@gsMarkerThicknessF = 5 ;若为雨，注释掉,不需要那么粗
  do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
      gsres@gsMarkerColor      = colors(i,:)
      pmid(i) = gsn_add_polymarker(wks,mapplot,lon_new(i,:),lat_new(i,:),gsres)      
    end if
  end do
  mapplot@$unique_string("pmid")$ = pmid
  end

undef("mark_area")
  procedure mark_area(wks,map)
  begin
  ypts = (/33.0,  33.0,  28.0, 28.0, 33.0/)
  xpts = (/89.0,  98.0,  98.0, 89.0, 89.0/)
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 4.0                       ; thickness of lines
  dum_box = new(4,graphic)  
  ; draw each line separately. Each line must contain two points.
  do i = 0 , 3
    dum_box(i)=gsn_add_polyline(wks,map,xpts(i:i+1),ypts(i:i+1),resp)      
  end do 
  map@$unique_string("dum_box")$ = dum_box
  end



;-----------------------------main code------------------------------------
begin
  pic_path="/mnt/h/d/ncl_related/pictures/rain/"
;读地形文件
  topo        = addfile("/mnt/h/d/data/other/ETOPO2v2c_f4.nc", "r")
  elev        = topo->z({15:50},{70:110}) 
  ; printVarSummary(elev)
  ; exit 
  lat_topo    = fspan(15,50,1050)
  lon_topo    = fspan(70,110,1200)
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

;读站点数据
  filepath= "/mnt/h/d/data/precipitation/6h/sta_r_s_total_py2600_2800.txt"  ;参数文件路径
  f1       = asciiread(filepath,-1,"string")
  sta     = stringtointeger(str_get_field(f1(1:),1," "))  ;读取站号
  lat     = stringtofloat(str_get_field(f1(1:),2," ")) 
  lon     = stringtofloat(str_get_field(f1(1:),3," ")) 
  rain_total= stringtofloat(str_get_field(f1(1:),4," "))
  snow_total= stringtofloat(str_get_field(f1(1:),5," "))
  prep_total= stringtofloat(str_get_field(f1(1:),6," "))
  ;rain_total 4;snow_total 5; prep_total 6; rain_oneday 7;snow_oneday 8；prep_oneday 9

  phase   ="prep_rain_snow"  ;用于修改图片名称,改画雨雪的时候注意307,319行也要改
  npts  = numAsciiRow(filepath)-1 ;获取站点数(文本行数)，注意这里是文本rain_dir，而不是字符串f
  R =new((/3,npts/),float)
  R(0,:)=prep_total 
  R(1,:)=rain_total
  R(2,:)=snow_total
  R@_FillValue = 9999.0
  ;注意：这里的txt文件最后一行数据下面需要空白行，才能读出文件的正确行数
  ;大多数txt文件都满足要求，但还是需要核实一下
  R=where(R.gt.0.1, R, R@_FillValue) ;剔除了小于等于0.1mm的降水
  ; print(num(.not.ismissing(R)))
  ; exit
  ; index=ind(prep_total.ge.50.) 
  ; print(sta(index))
  ; exit 
  ;日降雪量≥5mm的站点有30个，日降雨量超过25mm的有6个

  levels = (/0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100./) ; bin settings (bin0 = < 0., ; bin1 = 0.:4.999, etc.)
  nlevels = dimsizes(levels)
  colors  = span_color_rgba("wh-bl-gr-ye-re",nlevels+1)
  ;Given the number of desired color values, 
  ;return an array of RGB triplets or RGBA quadruplets 
  ;that nicely span the given color map.
  num_distinct_markers = nlevels+1        ; number of distinct markers
  lat_new3 = new((/3,num_distinct_markers,npts/),float,-999)
  lon_new3 = new((/3,num_distinct_markers,npts/),float,-999)
  do i= 0,2
    ref   =mark_lat_lon(num_distinct_markers,R(i,:),levels,lat,lon,npts)
    lat_new3(i,:,:)=ref[0]
    lon_new3(i,:,:)=ref[1]
  end do

;===========================================================================
; Begin plotting section.
  wks = gsn_open_wks("x11",pic_path+phase+"_2600_2800_topo2")              ; send graphics to PNG file

;---Set up some map resources.
  res                       = True
  ;res@gsnMaximize           = True             ; Maximize plot in frame.
  res@gsnDraw               = False            ; Will draw later
  res@gsnFrame              = False            ; Don't advance the frame
  res@gsnLeftString         = ""
  res@gsnRightString        = ""
  res@tmYROn                =False
  res@tmXTOn                =False
  ; res@tmBorderThicknessF  =4    ; border 6x as thick,图框粗细
  ; res@tmXBMajorThicknessF =4
  ; res@tmXBMinorThicknessF =4
  ; res@tmYLMajorThicknessF =4 
  ; res@tmYLMinorThicknessF =4 
  res@gsnAddCyclic          =False
  ; res@pmTickMarkDisplayMode = "Always"
  res@tmEqualizeXYSizes   =True 
  res@tmXBLabelFontHeightF = 0.025 ; x轴标签字号
  res@tmXBLabels=(/"85~S~o~N~E","90~S~o~N~E","95~S~o~N~E","100~S~o~N~E","105~S~o~N~E"/)
  res@tmYLLabels=(/"24~S~o~N~N","28~S~o~N~N","32~S~o~N~N"/)

  ; res@vpWidthF=1.0
  ; res@vpHeightF=0.6

  ;-----------地图属性-----------------
  mpres                       =res
  mpres@mpMinLatF             = 24.
  mpres@mpMaxLatF             = 35.
  mpres@mpMinLonF             = 85.
  mpres@mpMaxLonF             = 100.
  mpres@gsnMajorLonSpacing    =5. 
  mpres@gsnMinorLonSpacing    =1.    
  mpres@gsnMajorLatSpacing    =4. 
  mpres@gsnMinorLatSpacing    =1
  mpres@mpFillOn              =False
  mpres@mpDataBaseVersion       = "MediumRes"
  mpres@mpDataSetName           = "/mnt/d/ChinaMap/Earth..4"
  mpres@mpOutlineSpecifiers     = (/"China:states"/)
  ; mpres@mpOutlineOn = True
  ; mpres@mpOutlineBoundarySets           = "National"
  ; mpres@mpNationalLineColor             = "black"
  ; mpres@mpProvincialLineColor           = "black"
  ; mpres@mpGeophysicalLineColor          = "black"
  ; mpres@mpNationalLineThicknessF        = 1
  ; mpres@mpProvincialLineThicknessF      = 1
  ; mpres@mpGeophysicalLineThicknessF     = 1


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

;-----------------------------绘图-------------------------
  ;对地形高度数据做平滑
  elev = smth9_Wrap(elev,0.5,0.25,False)
  elev = smth9_Wrap(elev,0.5,0.25,False)
  elev = smth9_Wrap(elev,0.5,0.25,False)
  elev = smth9_Wrap(elev,0.5,0.25,False)  

  map1=new(3,graphic)
  map2=new(3,graphic)

  ;1 绘制第一张图 过程累计降水量
    mpres@mpOutlineOn = True
    mpres@mpOutlineBoundarySets           = "National"
    mpres@mpNationalLineColor             = "black"
    mpres@mpNationalLineThicknessF        = 3

    map2(0) = gsn_csm_map(wks,mpres)

    ;降雨散点
    plot_rain_scatter(wks,colors,num_distinct_markers,\
                        lat_new3(0,:,:),lon_new3(0,:,:),map2(0))
    ;强降水方框
    mark_area(wks,map2(0))

  ;2 单独绘制第2张日降雨量
    delete(mpres@mpDataBaseVersion)
    delete(mpres@mpDataSetName)
    delete(mpres@mpOutlineOn)
    delete(mpres@mpOutlineBoundarySets)

    map1(1) = gsn_csm_contour_map(wks,elev,mpres1)
    map2(1) = gsn_csm_map(wks,mpres)
    lon_range = 4. ;切除了过小的闭合等值线，使not busy
    lat_range = 4.
    map2(1) = add_isolines_by_range(wks,map1(1),\
                                        map2(1),lon_range,lat_range)
    add_topoline_label(wks,map2(1))
    plot_rain_scatter(wks,colors,num_distinct_markers,\
                        lat_new3(1,:,:),lon_new3(1,:,:),map2(1))
    mark_area(wks,map2(1)) 

  ;2 单独绘制第3张降雪量

    map1(2) = gsn_csm_contour_map(wks,elev,mpres1)
    map2(2) = gsn_csm_map(wks,mpres)
    lon_range = 4. ;切除了过小的闭合等值线，使not busy
    lat_range = 4.
    map2(2) = add_isolines_by_range(wks,map1(2),\
                                        map2(2),lon_range,lat_range)
    add_topoline_label(wks,map2(2))
    plot_snow_scatter(wks,colors,num_distinct_markers,\
                        lat_new3(2,:,:),lon_new3(2,:,:),map2(2))
    mark_area(wks,map2(2)) 

;3 单独在第一张图map2(0)上标记站点的站点位置
  letters=(/"Cona","L","Bomi","Zayu","TB","QH","SC"/)
  nletters=dimsizes(letters)
  dumletter = new(nletters,graphic)
  lat_letter=(/27.6, 29.1, 29.1, 27.8, 34,35.5,33.6/)
  lon_letter=(/91.4, 94.2, 96.8, 96.4, 87, 94, 99.2/)
  reslet = True
  reslet@txFontHeightF = 0.023
  reslet@txFontThicknessF=0.24
  reslet@txFontColor = "black"
  ; reslet@txJust = "CenterLeft"
  ; 站点名称
  do i=0,3
    dumletter(i)= gsn_add_text(wks,map2(0),letters(i),\
                               lon_letter(i), lat_letter(i)+1.2, reslet)
  end do
  ; 省份名称
  do i=4,nletters-1
    dumletter(i)= gsn_add_text(wks,map2(0),letters(i),\
                               lon_letter(i), lat_letter(i)-1.2, reslet)
  end do

  ; 单独在第一张图上标上三角形
  lat_triangle=(/27.98, 29.57, 29.87,  28.65/)
  lon_triangle=(/91.95, 94.47,95.77, 97.47/)
  gsres               = True
  gsres@gsMarkerIndex = 7   
  gsres@gsMarkerSizeF      = 0.02
  gsres@gsMarkerThicknessF = 1.5 ;若为雨，注释掉,不需要那么粗
  gsres@gsMarkerColor      = "black"
  pmid = gsn_add_polymarker(wks,map2(0),lon_triangle,lat_triangle+0.1,gsres)      


; 单独在第一张图标上国家
  country=(/"Myanmar","Bangladesh","India","Bhutan","Nepal","China"/)
  ncount=dimsizes(country)
  dumcount = new(ncount,graphic)
  lat_count=(/26, 24.6, 25.5, 27.3, 27.2, 31/)
  lon_count=(/97, 90.2, 93.5, 90.5, 86.5, 87/)
  reslet = True
  reslet@txFontHeightF = 0.022
  reslet@txFontThicknessF=0.24
  reslet@txFontColor = "red"
  ; reslet@txJust = "CenterLeft"
  do i=0,ncount-1
    dumcount(i)= gsn_add_text(wks,map2(0),country(i),\
                               lon_count(i), lat_count(i), reslet)
  end do

  ;-单独在第一张图上绘制高原轮廓线-----------------
  shp="/mnt/h/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 5.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,map2(0),shp,lnres) 

;4 组图
  pres=True
  pres@gsnFrame        = False  ; don't advance frame yet，以便后面使用gsn_(legend)_ndc 
  pres@gsnPanelBottom = 0.1   ;底部增加空白
  pres@gsnPanelFigureStrings=(/"(a)","(b)","(c)"/)
  pres@gsnPanelFigureStringsFontHeightF=0.015
  pres@amJust="topLeft"
  ; pres@gsnPanelFigureStringsPerimOn=False
  gsn_panel(wks,(/map2(0),map2(1),map2(2)/), (/3,1/), pres) 

;5 绘制降水散点图图例--------------------------
  ;---Draw labelbar and advance frame.
  add_labelbar(wks,levels,colors) 

;   draw(map2)
;   frame(wks)
end     
 

