;读取站点累计降水(剔除了<0.1mm以及mask掉高原以外站点降水)，画降水彩色散点图
;绘制地形3/4/5km等值线,并将等值线变稀疏（切除+平滑）
;文件读取方式 http://www.ncl.ucar.edu/Applications/list_io.shtml
;散点图参考 http://www.ncl.ucar.edu/Applications/polyg.shtml polyg_8_lbar.ncl
;代码路径  ncl /mnt/d/ncl_related/code/prep_3day_mask_topoline.ncl


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("read_elev_data")
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
  plres@gsLineThicknessF = 8
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

;加彩色散点的图例
undef("add_labelbar")
procedure add_labelbar(wks,plot,levels,colors)
local lbres, labels
begin
  nlevels = dimsizes(levels)
  ; Draw a labelbar
  lbres                    = True
  lbres@vpWidthF           = 0.7           ; width
  lbres@vpHeightF          = 0.07            ; height
  lbres@lbPerimOn          = False            ; Turn off perimeter.
  lbres@lbOrientation      = "Horizontal"     ; Default is vertical.
  lbres@lbLabelAlignment   = "InteriorEdges"  ; Default is "BoxCenters",标签如何与标签栏对齐
  lbres@lbFillColors       = colors           ; Colors for boxes.
  lbres@lbMonoFillPattern  = True             ; Fill them all solid.
  lbres@lbLabelFontHeightF = 0.02            ; label font height
  lbres@lbBoxLinesOn       =True
  lbres@lbBoxEndCapStyle   = "TriangleBothEnds"
  lbres@lbTitleString      ="precipitation(mm)"
  lbres@lbTitleFontHeightF = 0.02
  lbres@lbTitleJust        ="BottomRight"
  lbres@lbTitlePosition    ="Bottom"

  labels = sprintf("%2.0f",levels)

  lbid = gsn_create_labelbar(wks,nlevels+1,labels,lbres)
   
  ; Now, create some annotation resources indicating how we want to
  ; attach the labelbar to the plot. Here, we are using the top right
  ; corner of the labelbar as the point which we are going to position
  ; it, and then we use amParallelPosF and amOrthogonalPosF to indicate
  ; where we want to place it.
  ; amParallelPosF/amOrthogonalPosF                                       ;
  ;    0.0/ 0.0  -  annotation in dead center of plot
  ;    0.5/ 0.5  -  annotation at bottom right of plot
  ;    0.5/-0.5  -  annotation at top right of plot
  ;   -0.5/-0.5  -  annotation at top left of plot
  ;   -0.5/ 0.5  -  annotation at bottom left of plot
  amres                  = True
  amres@amJust           = "TopCenter" ;"TopCenter"
  amres@amParallelPosF   =  0.0   ; keep labelbar centered
  amres@amOrthogonalPosF =  0.7  ; move down 

  ; Give both annotation id and labelbar id unique names.
  ; Attaching them to plot with unique names ensures that
  ; labelbar "lives" outside this procedure.
  
  tmpid1 = "anno"+unique_string("id")
  tmpid2 = "lbar"+unique_string("id")

  plot@$tmpid1$ = gsn_add_annotation(plot,lbid,amres)
  plot@$tmpid2$ = lbid
end
;-----------------------------main code------------------------------------
begin
;-----------------读地形文件----------------------------------------------------
  data1_dir  ="/mnt/d/data/other/"     ;地形数据路径
  setfileoption("nc","MissingToFillValue",False)        ;设置缺省值
  elev_old= read_elev_data(data1_dir+"ETOPO5.DAT")
  elev=tofloat(elev_old)
  copy_VarCoords(elev_old, elev)
  elev@_FillValue=9.e+10

;--------------------读站点数据--------
  filepath = "/mnt/d/data/precipitation/6h/"  ;参数文件路径
  name="rain12_py"
  picture_path="/mnt/d/ncl_related/pictures/rain/"
  f       = asciiread(filepath+name+".txt",-1,"string")
  sta     = stringtointeger(str_get_field(f,1," "))  ;读取站号
  lat     = stringtofloat(str_get_field(f,2," ")) ;获取纬度值lon
  lon     = stringtofloat(str_get_field(f,3," ")) ;获取经度值lat
  R_old   = stringtofloat(str_get_field(f,17," ")) ;获取第16列站点累计降水

 ;设置缺测值，mask掉高原以外的数据
  R_old@_FillValue = 9999.0
  ;注意：这里的txt文件最后一行数据下面需要空白行，才能读出文件的正确行数
  ;大多数txt文件都满足要求，但还是需要核实一下
  npts = numAsciiRow(filepath+name+".txt") ;获取行数-站点数
  ;npts=3280 ;2512-2612 3275;  2612-2712是3280；2712-2812 3162

  ;---Open shapefile and read Mississippi River Basin lat/lon values.
  f1      = addfile("/mnt/d/data/TPshp/DBATP/DBATP_Polygon.shp", "r")
  mrb_lon = f1->x
  mrb_lat = f1->y
  ;---Start with data all missing
  R = new(dimsizes(R_old),typeof(R_old),R_old@_FillValue)
  do i=0,npts-1
    R(i)=where(gc_inout(lat(i),lon(i),mrb_lat,mrb_lon), R_old(i), R_old@_FillValue)
  end do

R=where(R.ge.0.1, R, R@_FillValue) ;剔除了小于0.1mm的降水,剩下66
; print(num(.not.ismissing(R)))
; exit

  ; index=ind(R.ge.50.) 
  ; print(sta(index))
  ; exit ;降水量≥50mm的站点/降水量 ,7个
  ; (0)55690 错那* (1)56227 波密* (2)56317 米林 (3)56434 察隅 (4)56444 德钦* (5)56533 贡山* (6)56643 怒江州*
  ; (0)     107  (1)     111    (2)     57.4 (3)     96    (4)     101.7(5)     132.4 (6)     139.6

;------------------------------ 
; Create X and Y arrays to hold the points for each range and initialize
; them to missing values.  We want to use num_distinct_markers
; different colors, so we need num_distinct_markers sets of X and
; Y points.
  levels = (/0.,5.,25.,30.,40.,50.,60.,70.,80.,90.,100/) ; bin settings (bin0 = < 0., ; bin1 = 0.:4.999, etc.)
  nlevels = dimsizes(levels)
  colors  = span_color_rgba("wh-bl-gr-ye-re",nlevels+1)
  ;Given the number of desired color values, 
  ;return an array of RGB triplets or RGBA quadruplets 
  ;that nicely span the given color map.

  num_distinct_markers = nlevels+1        ; number of distinct markers
  lat_new = new((/num_distinct_markers,npts/),float,-999)
  lon_new = new((/num_distinct_markers,npts/),float,-999)

; Group the points according to which range they fall in. At the
; same time, create the label that we will use later in the labelbar
;       
  do i = 0, num_distinct_markers-1
    if (i.eq.0) then
      indexes = ind(R.lt.levels(0))
    end if
    if (i.eq.num_distinct_markers-1) then
      indexes = ind(R.ge.max(levels))
    end if
    if (i.gt.0.and.i.lt.num_distinct_markers-1) then       
      indexes = ind(R.ge.levels(i-1).and.R.lt.levels(i))
    end if
;; Now that we have the set of indexes whose values fall within 
; the given range, take the corresponding lat/lon values and store
; them, so later we can color this set of markers with the appropriate
; color.
;
    if (.not.any(ismissing(indexes))) then 
      npts_range = dimsizes(indexes)   ; # of points in this range.
      lat_new(i,0:npts_range-1) = lat(indexes)
      lon_new(i,0:npts_range-1) = lon(indexes)
    end if
    delete(indexes)            ; Necessary b/c "indexes" may be a different
    ; size next time.
  end do

;===========================================================================
; Begin plotting section.
;
  wks = gsn_open_wks("png",picture_path+"prep_total")              ; send graphics to PNG file

;---Set up some map resources.
  res                       = True
  ;res@gsnMaximize           = True             ; Maximize plot in frame.
  res@gsnDraw               = False            ; Will draw later
  res@gsnFrame              = False            ; Don't advance the frame
  res@gsnLeftString         = ""
  res@gsnRightString        = ""
  res@tmYROn                =False
  res@tmXTOn                =False
  res@gsnAddCyclic          =False
  res@tmBorderThicknessF  =4    ; border 6x as thick,图框粗细
  res@tmXBMajorThicknessF =4
  res@tmXBMinorThicknessF =4
  res@tmYLMajorThicknessF =4 
  res@tmYLMinorThicknessF =4 
  res@gsnAddCyclic          =False
  ; res@pmTickMarkDisplayMode = "Always"
  res@tmEqualizeXYSizes   =True 
  res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号

  ; res@vpWidthF=1.0
  ; res@vpHeightF=0.6

;-----------地图属性-----------------
  mpres                       =res
  mpres@mpMinLatF             = 24.
  mpres@mpMaxLatF             = 36.
  mpres@mpMinLonF             = 86.
  mpres@mpMaxLonF             = 104.
  mpres@gsnMajorLonSpacing    =3. 
  mpres@gsnMinorLonSpacing    =1.    
  mpres@gsnMajorLatSpacing    =3. 
  mpres@gsnMinorLatSpacing    =1
  mpres@mpDataBaseVersion       = "MediumRes"
  ; mpres@mpDataSetName           = "Earth..4"
  ; mpres@mpOutlineSpecifiers     = (/"China:states"/)
  mpres@mpFillOn              =False
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
  mpres1@cnLevels           =(/3000,4000,5000/)
  mpres1@cnLineColors       =(/"gray0","gray40","gray70"/)
  mpres1@cnLineThicknessF   = 8.
  mpres1@cnLineLabelsOn     = False            ; turn off line labels
  mpres1@cnInfoLabelOn      = False            ; turn off info label

;-----------------------------绘制地形图-------------------------

; First plot is original topo data drawn as line contours over a map.
  ;topo_line_plot=gsn_csm_contour_map(wks,elev,mpres)
  elev = smth9_Wrap(elev,0.5,0.25,False)
  elev = smth9_Wrap(elev,0.5,0.25,False)
  elev = smth9_Wrap(elev,0.5,0.25,False)

  map1 = gsn_csm_contour_map(wks,elev,mpres1)
  map2 = gsn_csm_map(wks,mpres)
  lon_range = 4. ;切除了过小的闭合等值线，使not busy
  lat_range = 4.
  map2 = add_isolines_by_range(wks,map1,map2,lon_range,lat_range)


;---------------------------------------绘制降水散点图--------------
;Create logical variables to hold the marker resources.
  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.
;
; Loop through each grouping of markers, and draw them one set at
; a time, assigning the proper color and size with gsn_marker.
;
  ; base_size = 0.01
  pmid = new(num_distinct_markers,graphic)
  do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
      gsres@gsMarkerColor      = colors(i,:)
      gsres@gsMarkerSizeF      = 0.015
      ; gsres@gsMarkerSizeF      = base_size * (i+1)/6.
      ; gsres@gsMarkerThicknessF = 0.5*(i+1)
      pmid(i) = gsn_add_polymarker(wks,map2,lon_new(i,:),lat_new(i,:),gsres)
    end if
  end do

;-------------------------绘制图例--------------------------
 ;---Draw labelbar and advance frame.
 add_labelbar(wks,map2,levels,colors) ;降水散点图图例

 ;高原三条等高线图例
  resLg = True
  resLg@lgItemType = "Lines" 
  resLg@lgMonoDashIndex = True
  resLg@lgPerimOn    =False
  resLg@lgDashIndex  = 0
  resLg@lgLineColors = (/"gray0","gray40","gray70"/)
  resLg@lgLineThicknessF = 8
  resLg@vpWidthF  = 0.15
  resLg@vpHeightF = 0.15
  resLg@lgPerimFill = 0
  resLg@lgPerimFillColor = "Background"
  resLg@lgLabelFontHeightF = 0.06
  resLg@lgTitleString = ""
  lbid = gsn_create_legend(wks, 3, (/"3 km",\
  "4 km","5 km"/), resLg)
 
; 将图例放置在图中
  amres = True
  amres@amParallelPosF = -0.38 ;左右，正右左负
  amres@amOrthogonalPosF = 0.34 ;上下，正下负上
  dumLg = gsn_add_annotation(map2, lbid, amres)

; shp="/mnt/d/data/TPshp/DBATP/DBATP_Line.shp"  
;  lnres                  = True     
;  lnres@gsLineColor      ="orange"
;  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
;  shp_plot               = gsn_add_shapefile_polylines(wks,map2,shp,lnres) 

 
  draw(map2)
  frame(wks)
end     
 

