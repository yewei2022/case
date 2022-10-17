;读取4个时次micaps地面相同站点降水资料，相加，画日降水散点图
;时间改三处，79,120，256
;文件读取方式 http://www.ncl.ucar.edu/Applications/list_io.shtml
;散点图参考 http://www.ncl.ucar.edu/Applications/polyg.shtml polyg_8_lbar.ncl
;调色标位置大小，画布大小vpwidth，坐标标签大小


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

undef("add_labelbar")
procedure add_labelbar(wks,plot,levels,colors)
local lbres, labels
begin
  nlevels = dimsizes(levels)
;----------------------------------------------------------------------
; Draw a labelbar
;----------------------------------------------------------------------
  lbres                    = True
  lbres@vpWidthF           = 0.05           ; width
  lbres@vpHeightF          = 0.6            ; height
  lbres@lbPerimOn          = False            ; Turn off perimeter.
  ;lbres@lbOrientation      = "Horizontal"     ; Default is vertical.
  lbres@lbLabelAlignment   = "InteriorEdges"  ; Default is "BoxCenters",标签如何与标签栏对齐
  lbres@lbFillColors       = colors           ; Colors for boxes.
  lbres@lbMonoFillPattern  = True             ; Fill them all solid.
  lbres@lbLabelFontHeightF = 0.02            ; label font height
  lbres@lbBoxLinesOn       =True
  lbres@lbBoxEndCapStyle   = "TriangleBothEnds"

  labels = sprintf("%2.0f",levels)

  lbid = gsn_create_labelbar(wks,nlevels+1,labels,lbres)
   
;
; Now, create some annotation resources indicating how we want to
; attach the labelbar to the plot. Here, we are using the top right
; corner of the labelbar as the point which we are going to position
; it, and then we use amParallelPosF and amOrthogonalPosF to indicate
; where we want to place it.
;
; amParallelPosF/amOrthogonalPosF                                       ;
;    0.0/ 0.0  -  annotation in dead center of plot
;    0.5/ 0.5  -  annotation at bottom right of plot
;    0.5/-0.5  -  annotation at top right of plot
;   -0.5/-0.5  -  annotation at top left of plot
;   -0.5/ 0.5  -  annotation at bottom left of plot
;                                                                 
  amres                  = True
  amres@amJust           = "TopCenter"
  amres@amParallelPosF   =  0.56   ; keep labelbar centered
  amres@amOrthogonalPosF =  -0.50  ; move down and outside of plot

;
; Give both annotation id and labelbar id unique names.
;
; Attaching them to plot with unique names ensures that
; labelbar "lives" outside this procedure.
;
  tmpid1 = "anno"+unique_string("id")
  tmpid2 = "lbar"+unique_string("id")

  plot@$tmpid1$ = gsn_add_annotation(plot,lbid,amres)
  plot@$tmpid2$ = lbid
end

;------------------------main code---------------------------------
begin

;-------Options--------

  levels = (/0,5,10,20,30,40,50,60,70,80,90/) ; bin settings (bin0 = < 0., ; bin1 = 0.:4.999, etc.)
  nlevels = dimsizes(levels)
  colors  = span_color_rgba("wh-bl-gr-ye-re",nlevels+1)

;1.读取数据
  filepath = "/cygdrive/D/ncl_related/data/precipitation/6h/"  ;参数文件路径
  time_num="2512_2612"
  picture_path="/cygdrive/D/ncl_related/pictures/rain/"
  f=asciiread(filepath+time_num+".txt",-1,"string")
  station = stringtofloat(str_get_field(f,2," ")) ;从数组f中获取站号
  lat = stringtofloat(str_get_field(f,3," ")) ;获取纬度值lon
  lon = stringtofloat(str_get_field(f,4," ")) ;获取经度值lat
  R1  = stringtofloat(str_get_field(f,5," ")) ;获取第1个时次站点数据
  R2  = stringtofloat(str_get_field(f,6," ")) ;获取第2个时次站点数据
  R3  = stringtofloat(str_get_field(f,7," ")) ;获取第3个时次站点数据
  R4  = stringtofloat(str_get_field(f,8," ")) ;获取第4个时次站点数据

  ;风暴数据,用于标记风暴位置或路径
 
  fiTY = "/cygdrive/D/ncl_related/data/TC/bio042008.txt"
 
  ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
  nrow = numAsciiRow(fiTY)
 
  YYYYMMDDHH = new(nrow, "string")
  latt = new(nrow, "float")
  lonn = new(nrow, "float")
 
  data = asciiread(fiTY, -1, "string")
  YYYYMMDDHH = str_get_field(data, 3, ",")
  latt = stringtofloat(str_get_field(data, 7, ",")) *0.1
  lonn = stringtofloat(str_get_field(data, 8, ",")) *0.1

  DateChar = stringtochar(YYYYMMDDHH)
  MM = chartostring(DateChar(:,5:6))
  DD = chartostring(DateChar(:,7:8))
  HH = chartostring(DateChar(:,9:10))

 
;相同站点降水数据相加
 ;设置缺测值
 R1@_FillValue = 9999.0
 R2@_FillValue = 9999.0
 R3@_FillValue = 9999.0
 R4@_FillValue = 9999.0


 npts=3275 ;2512-2612 3275;  2612-2712是3280；2712-2812 3162
 R=new(npts,float)
 R=R1+R2+R3+R4
 R@_FillValue =  9.96921e+36 

 ;print(R(607))
 ;exit

 ;print(lon(607))
 ;exit

;------------------------------ 
; Create X and Y arrays to hold the points for each range and initialize
; them to missing values.  We want to use num_distinct_markers
; different colors, so we need num_distinct_markers sets of X and
; Y points.
;
  num_distinct_markers = nlevels+1        ; number of distinct markers
  lat_new = new((/num_distinct_markers,npts/),float,-999)
  lon_new = new((/num_distinct_markers,npts/),float,-999)
;
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
  wks = gsn_open_wks("png",picture_path+"micaps_"+time_num)              ; send graphics to PNG file

;---Set up some map resources.
  mpres                       = True
  mpres@mpFillOn              = False
  ;mpres@gsnMaximize           = True             ; Maximize plot in frame.
  mpres@gsnDraw               = False            ; Will draw later
  mpres@gsnFrame              = False            ; Don't advance the frame
  mpres@gsnLeftString         = ""
  mpres@gsnRightString        = ""
  mpres@tmYROn                =False
  mpres@tmXTOn                =False
  ;mpres@pmTickMarkDisplayMode = "Always"

  mpres@mpMinLatF             = 15.
  mpres@mpMaxLatF             = 45.
  mpres@mpMinLonF             = 70
  mpres@mpMaxLonF             = 110
  mpres@gsnMajorLonSpacing    =5. 
  mpres@gsnMinorLonSpacing    =1.    
  mpres@gsnMajorLatSpacing    =5. 
  mpres@gsnMinorLatSpacing    =1


  ;设置坐标轴标签粗细
  mpres@tmXBLabelFont=22
  mpres@tmYLLabelFont=22
  mpres@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
  mpres@tmYLLabelFontHeightF =0.02

  mpres@vpWidthF=1.0
  mpres@vpHeightF=0.6

  mpres@mpDataBaseVersion       = "MediumRes"
  mpres@mpDataSetName           = "Earth..4"
  mpres@mpOutlineSpecifiers     = (/"China:states"/)
  ;边界线
  mpres@mpOutlineBoundarySets           = "National"
  mpres@mpNationalLineColor             = "black"
  mpres@mpProvincialLineColor           = "black"
  mpres@mpGeophysicalLineColor          = "black"
  mpres@mpNationalLineThicknessF        = 1
  mpres@mpProvincialLineThicknessF      = 1
  mpres@mpGeophysicalLineThicknessF     = 1

;画地图
  map = gsn_csm_map(wks,mpres)

;--Create logical variables to hold the marker resources.
  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.
;
; Loop through each grouping of markers, and draw them one set at
; a time, assigning the proper color and size with gsn_marker.
;
  base_size = 0.01
  pmid = new(num_distinct_markers,graphic)
  do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
      gsres@gsMarkerColor      = colors(i,:)
      gsres@gsMarkerSizeF      = base_size * (i+1)/3.
      gsres@gsMarkerThicknessF = 0.7*(i+1)
      pmid(i) = gsn_add_polymarker(wks,map,lon_new(i,:),lat_new(i,:),gsres)
    end if
  end do

;---Draw labelbar and advance frame.
  add_labelbar(wks,map,levels,colors)

;某日风暴路径绘制
  resLine =True
  resLine@gsLineThicknessF = 4
  resLine@gsLineColor      = "black"

  txres = True                     ; text mods desired
  txres@txFontHeightF   = 0.03                     ; font smaller. default big
  txres@txFontThicknessF=4
  txres@txFontColor     ="black"


  resTx = True
  resTx@txFontHeightF = 0.015
  resTx@txFontColor = "black"
  resTx@txJust = "CenterLeft"

  nnrow = 5
  dumLine =new(nnrow, graphic)
  lon1=new(nnrow,float)
  lat1=new(nnrow,float)
  DD1=new(nnrow,string)
  HH1=new(nnrow,string)
  ;存取26日18时，27日经纬度数组
  ii = 0
  do i=0,nrow-1
    if (((DD(i) .eq. "25").and.(HH(i) .eq. "12"))\
      .or.((DD(i) .eq. "25").and.(HH(i) .eq. "18"))\
      .or.((DD(i) .eq. "26").and.(HH(i) .eq. "00"))\
      .or.((DD(i) .eq. "26").and.(HH(i) .eq. "06"))\
      .or.((DD(i) .eq. "26").and.(HH(i) .eq. "12"))) then
      lon1(ii)=lonn(i)
      lat1(ii)=latt(i)
      DD1(ii)=DD(i)
      HH1(ii)=HH(i)
      ii = ii+1
    end if 
  end do

  ;画线
  do i = 0, nnrow-2
    xx = (/ lon1(i), lon1(i+1)/)
    yy = (/ lat1(i), lat1(i+1)/)
    dumLine(i) = gsn_add_polyline(wks, map, xx, yy, resLine)
  end do 

  ;画点; 标注日期

  dumDate = new(nnrow,graphic)
  do i = 0, nnrow-1
    dum = gsn_add_text(wks,map,"~F37~p",lon1(i), lat1(i),txres) 
    dumDate = gsn_add_text(wks,map, DD1(i)+HH1(i), lon1(i)+1, lat1(i), resTx)
  end do

shp="/cygdrive/D/ncl_related/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,map,shp,lnres) 
 
  draw(map)
  frame(wks)
end     
 

