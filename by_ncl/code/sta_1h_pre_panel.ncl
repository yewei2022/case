;读取高原站点1h降水量，组图
;文件读取方式 http://www.ncl.ucar.edu/Applications/list_io.shtml
;散点图参考 http://www.ncl.ucar.edu/Applications/polyg.shtml polyg_8_lbar.ncl
;  ncl /mnt/d/ncl_related/code/sta_1h_pre_panel.ncl


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



; 加彩色散点的图例,
; 在panel中加公共图例(gsn_labelbar_ndc)和
;在子图中(gsn_create_labelbar+gsn_add_annotation)
undef("add_labelbar")
  procedure add_labelbar(wks,levels,colors)
  local lbres, labels
  begin
  nlevels = dimsizes(levels)
  ; Draw a labelbar
  lbres                    = True
  lbres@vpWidthF           = 0.35           ; width
  lbres@vpHeightF          = 0.05            ; height
  lbres@lbPerimOn          = False            ; Turn off perimeter.
  lbres@lbOrientation      = "Horizontal"     ; Default is vertical.
  lbres@lbLabelAlignment   = "InteriorEdges"  ; Default is "BoxCenters",标签如何与标签栏对齐
  lbres@lbFillColors       = colors           ; Colors for boxes.
  lbres@lbMonoFillPattern  = True             ; Fill them all solid.
  lbres@lbLabelFontHeightF = 0.01            ; label font height
  lbres@lbBoxLinesOn       =True
  lbres@lbBoxEndCapStyle   = "TriangleBothEnds"
  lbres@lbTitleString      ="(mm)"
  lbres@lbTitleFontHeightF = 0.01
  lbres@lbTitleJust        ="BottomRight"
  lbres@lbTitlePosition    ="Bottom"

  labels = sprintf("%2.1f",levels)

  gsn_labelbar_ndc(wks,nlevels+1,labels,0.3,0.25,lbres)
  frame(wks)                                   ; now advance the frame

  end

; ; 一张图的图例 
; undef("add_labelbar")
;   procedure add_labelbar(wks,levels,colors)
;   local lbres, labels
;   begin
;   nlevels = dimsizes(levels)
;   ; Draw a labelbar
;   lbres                    = True
;   lbres@vpWidthF           = 0.6           ; width
;   lbres@vpHeightF          = 0.12            ; height
;   lbres@lbPerimOn          = False            ; Turn off perimeter.
;   lbres@lbOrientation      = "Horizontal"     ; Default is vertical.
;   lbres@lbLabelAlignment   = "InteriorEdges"  ; Default is "BoxCenters",标签如何与标签栏对齐
;   lbres@lbFillColors       = colors           ; Colors for boxes.
;   lbres@lbMonoFillPattern  = True             ; Fill them all solid.
;   lbres@lbLabelFontHeightF = 0.015            ; label font height
;   lbres@lbBoxLinesOn       =True
;   ; lbres@lbBoxEndCapStyle   = "TriangleBothEnds"
;   lbres@lbTitleString      ="(mm)"
;   lbres@lbTitleFontHeightF = 0.015
;   lbres@lbTitleJust        ="BottomRight"
;   lbres@lbTitlePosition    ="Bottom"
;   lbres@lbLabelStride      =1

;   labels = sprintf("%2.0f",levels)

;   gsn_labelbar_ndc(wks,nlevels+1,labels,0.2,0.17,lbres)
;   frame(wks)                                   ; now advance the frame

;   end

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
      indexes = ind(R.gt.max(levels))
    end if
    if (i.gt.0.and.i.lt.num_distinct_markers-1) then    ;中间的范围，因为格式不一样，要分别设置
      indexes = ind(R.gt.levels(i-1).and.R.le.levels(i))
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
  gsres@gsMarkerSizeF      = 0.005  ;daily 0.015 hourly0.005
  do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
      gsres@gsMarkerColor      = colors(i,:)
      pmid(i) = gsn_add_polymarker(wks,mapplot,lon_new(i,:),lat_new(i,:),gsres)      
    end if
  end do
  mapplot@$unique_string("pmid")$ = pmid
  end



;-----------------------------main code------------------------------------
begin

;读站点数据
  filepath= "/mnt/d/data/precipitation/1h/tp_1h_pre.txt" 
  pic_path="/mnt/d/ncl_related/pictures/rain/sta_1h/"
  f       = asciiread(filepath,-1,"string")
  sta     = stringtofloat(str_get_field(f(1:),1,",")) 
  lat     = stringtofloat(str_get_field(f(1:),2,",")) 
  lon     = stringtofloat(str_get_field(f(1:),3,",")) 
  npts  = numAsciiRow(filepath)-1
  R =new((/24,npts/),float)
  R@_FillValue = 32766.
  do i=0,23
    ii=i+4
    R(i,:)= stringtofloat(str_get_field(f(1:),ii,",")) 
  end do
  R=R*0.1
  R_total=dim_sum_n_Wrap(R, 0)

  ; print(sta)
  ; exit
;   ;注意：这里的txt文件最后一行数据下面需要空白行，才能读出文件的正确行数
;   ;大多数txt文件都满足要求，但还是需要核实一下
;   R=where(R.gt.0.1, R, R@_FillValue) ;剔除了小于等于0.1mm的降水
;   ; print(num(.not.ismissing(R)))
;   ; exit
;   ; index=ind(R.ge.25.) 
;   ; print(R(index))
;   ; exit 


  levels = (/0.5,1,1.5,2,2.5,3,3.5,4,4.5,5/) ; hour pre
  ; levels = (/0.,5.,25.,30,40,50,60,70,80,90,100./) ; total pre

  nlevels = dimsizes(levels)
  colors  = span_color_rgba("wh-bl-gr-ye-re",nlevels+1)
  num_distinct_markers = nlevels+1        ; number of distinct markers
  lat_new = new((/24,num_distinct_markers,npts/),float,-999)
  lon_new = new((/24,num_distinct_markers,npts/),float,-999)
  do i= 0,23
    ref   =mark_lat_lon(num_distinct_markers,R(i,:),levels,lat,lon,npts)
    lat_new(i,:,:)=ref[0]
    lon_new(i,:,:)=ref[1]
  end do

  ref   =mark_lat_lon(num_distinct_markers,R_total,levels,lat,lon,npts)
  lat_new_total=ref[0]
  lon_new_total=ref[1]

;===========================================================================
; Begin plotting section.
  wks = gsn_open_wks("ps",pic_path+"sta_1h_pre")              ; send graphics to PNG file

;---Set up some map resources.
  res                       = True
  ; res@gsnMaximize           = True             ; Maximize plot in frame.
  res@gsnDraw               = False            ; Will draw later
  res@gsnFrame              = False            ; Don't advance the frame
  res@gsnLeftString         = ""
  res@gsnRightString        = ""
  res@tmYROn                =False
  res@tmXTOn                =False

  ; res@vpWidthF=1.0
  ; res@vpHeightF=0.6

  ;-----------地图属性-----------------
  mpres                       =res
  mpres@mpMinLatF             = 25.
  mpres@mpMaxLatF             = 40.
  mpres@mpMinLonF             = 85
  mpres@mpMaxLonF             = 105
  mpres@gsnMajorLonSpacing    =5. 
  mpres@gsnMinorLonSpacing    =1.    
  mpres@gsnMajorLatSpacing    =5. 
  mpres@gsnMinorLatSpacing    =1

  mpres@mpFillOn=False
  mpres@mpDataBaseVersion       = "MediumRes"
  mpres@mpDataSetName           =  "/mnt/d/ChinaMap/Earth..4"
  mpres@mpOutlineSpecifiers     = (/"China:states"/)
  ;边界线
  mpres@mpOutlineBoundarySets   = "National"
  mpres@mpNationalLineColor     = "black"
  mpres@mpProvincialLineColor   = "black"
  mpres@mpGeophysicalLineColor  = "black"

  ;高原轮廓线
  shp_line="/mnt/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 4.0           ; 6x thickness 


; 绘制各时刻小时降水-------------------------
  map=new(24,graphic)
  plot_topoline=new(24,graphic)
  ;1 绘制前两张是过程累计降水，日降雨量
  do iplot=0,23
    map(iplot) = gsn_csm_map(wks,mpres)
    plot_rain_scatter(wks,colors,num_distinct_markers,\
                      lat_new(iplot,:,:),lon_new(iplot,:,:),map(iplot))
    plot_topoline(iplot) = gsn_add_shapefile_polylines(wks,map(iplot),shp_line,lnres) 
  end do 
;组图
  time_dir= "/mnt/d/data/total_1h_pre.nc"  
  data_time    = addfile(time_dir, "r")
  time    = data_time->time(13:36)
  YYYYMMDDHH=cd_calendar(time,-3)
  pres=True
  pres@gsnPanelCenter        =True ;对齐为默认值
  pres@gsnPanelFigureStrings =YYYYMMDDHH
  pres@gsnPanelFigureStringsFontHeightF=0.01
  pres@amJust                ="topLeft" ;图注放在左上方
  pres=True
  pres@gsnFrame        = False  ; don't advance frame yet，以便后面使用gsn_(legend)_ndc 
  ; pres@gsnPanelBottom = 0.1   ;底部增加空白
  gsn_panel(wks,map,(/4,6/), pres) 
; 绘制降水散点图图例
  add_labelbar(wks,levels,colors) 


; ; 绘制日累计降水
;   map = gsn_csm_map(wks,mpres)
;   plot_rain_scatter(wks,colors,num_distinct_markers,lat_new_total,lon_new_total,map)
;   plot_topoline= gsn_add_shapefile_polylines(wks,map,shp_line,lnres) 
; draw(map)
; add_labelbar(wks,levels,colors) 

end     
 

