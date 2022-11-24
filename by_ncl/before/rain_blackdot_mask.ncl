;读取4个时次micaps地面相同站点降水资料，相加，画日降水散点图
;时间改三处，23,43，216
;文件读取方式 http://www.ncl.ucar.edu/Applications/list_io.shtml
;散点图参考 http://www.ncl.ucar.edu/Applications/polyg.shtml polyg_8_lbar.ncl
;调色标位置大小，画布大小vpwidth，坐标标签大小


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


begin
;-----------------读地形文件----------------------------------------------------
   data1_dir  ="/cygdrive/D/ncl_related/data/other/"     ;地形数据路径
   setfileoption("nc","MissingToFillValue",False)        ;设置缺省值
   elev= read_elev_data(data1_dir+"ETOPO5.DAT")

;--------------------读站点数据--------
;1.读取数据
  filepath = "/cygdrive/D/ncl_related/data/micaps_rain/"  ;参数文件路径
  time_num="2612_2712"
  picture_path="/cygdrive/D/ncl_related/pictures/rain/"
  f=asciiread(filepath+time_num+".txt",-1,"string")
  station = stringtofloat(str_get_field(f,2," ")) ;从数组f中获取站号
  lat = stringtofloat(str_get_field(f,3," ")) ;获取纬度值lon
  lon = stringtofloat(str_get_field(f,4," ")) ;获取经度值lat
  R1  = stringtofloat(str_get_field(f,5," ")) ;获取第1个时次站点数据
  R2  = stringtofloat(str_get_field(f,6," ")) ;获取第2个时次站点数据
  R3  = stringtofloat(str_get_field(f,7," ")) ;获取第3个时次站点数据
  R4  = stringtofloat(str_get_field(f,8," ")) ;获取第4个时次站点数据

 
;相同站点降水数据相加
 ;设置缺测值
 R1@_FillValue = 9999.0
 R2@_FillValue = 9999.0
 R3@_FillValue = 9999.0
 R4@_FillValue = 9999.0
 ;注意：这里的txt文件最后一行数据后面回车
 ;也就是下面有一行空白，才能读出文件的正确行数
 npts = numAsciiRow(filepath+time_num+".txt") ;获取行数-站点数
 ;npts=3280 ;2512-2612 3275;  2612-2712是3280；2712-2812 3162
 R_old=new(npts,float)
 R_old=R1+R2+R3+R4
 R_old@_FillValue =  9.96921e+36 

  ;---Open shapefile and read Mississippi River Basin lat/lon values.
  f1       = addfile("/cygdrive/D/ncl_related/data/TPshp/DBATP/DBATP_Polygon.shp", "r")
  mrb_lon = f1->x
  mrb_lat = f1->y

;---Start with data all missing
    R = new(dimsizes(R_old),typeof(R_old),R_old@_FillValue)

 do i=0,npts-1
        if(gc_inout(lat(i),lon(i),mrb_lat,mrb_lon)) then
          R(i) = R_old(i)
        else
          R(i)=R_old@_FillValue
        end if
      end do


;------------------------------ 
; Create X and Y arrays to hold the points for each range and initialize
; them to missing values.  We want to use num_distinct_markers
; different colors, so we need num_distinct_markers sets of X and
; Y points.
  levels = (/0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90./) ; bin settings (bin0 = < 0., ; bin1 = 0.:4.999, etc.)
  nlevels = dimsizes(levels)

  num_distinct_markers = nlevels+1        ; number of distinct markers
  lat_new = new((/num_distinct_markers,npts/),float,-999)
  lon_new = new((/num_distinct_markers,npts/),float,-999)
;
; Group the points according to which range they fall in. At the
; same time, create the label that we will use later in the labelbar

  labels = new(num_distinct_markers,string)  ; Labels for legend.       
  do i = 0, num_distinct_markers-1
    if (i.eq.0) then
      indexes = ind(R.lt.levels(0))
      labels(i) = "x < " + levels(0)
    end if
    if (i.eq.num_distinct_markers-1) then
      indexes = ind(R.ge.max(levels))
      labels(i) = "x >= " + max(levels)
    end if
    if (i.gt.0.and.i.lt.num_distinct_markers-1) then       
      indexes = ind(R.ge.levels(i-1).and.R.lt.levels(i))
      labels(i) = levels(i-1) + " <= x < " + levels(i)
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
  wks = gsn_open_wks("png",picture_path+"micaps_blackdot_"+time_num)              ; send graphics to PNG file

;---Set up some map resources.
  res                       = True
  ;res@gsnMaximize           = True             ; Maximize plot in frame.
  res@gsnDraw               = False            ; Will draw later
  res@gsnFrame              = False            ; Don't advance the frame
  res@gsnLeftString         = ""
  res@gsnRightString        = ""
  res@tmYROn                =False
  res@tmXTOn                =False
  ;res@pmTickMarkDisplayMode = "Always"
  ;设置坐标轴标签粗细
  res@tmXBLabelFont=22
  res@tmYLLabelFont=22
  res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
  res@tmYLLabelFontHeightF =0.02

  ; res@vpWidthF=1.0
  ; res@vpHeightF=0.6

;-----------地图属性-----------------
  mpres                       =res
  mpres@mpMinLatF             = 15.
  mpres@mpMaxLatF             = 45.
  mpres@mpMinLonF             = 70
  mpres@mpMaxLonF             = 110
  mpres@gsnMajorLonSpacing    =5. 
  mpres@gsnMinorLonSpacing    =1.    
  mpres@gsnMajorLatSpacing    =5. 
  mpres@gsnMinorLatSpacing    =1
  mpres@mpDataBaseVersion       = "MediumRes"
  mpres@mpDataSetName           = "Earth..4"
  mpres@mpOutlineSpecifiers     = (/"China:states"/)
  mpres@mpOutlineOn = True
  mpres@mpOutlineBoundarySets           = "National"
  mpres@mpNationalLineColor             = "black"
  mpres@mpProvincialLineColor           = "black"
  mpres@mpGeophysicalLineColor          = "black"
  mpres@mpNationalLineThicknessF        = 1
  mpres@mpProvincialLineThicknessF      = 1
  mpres@mpGeophysicalLineThicknessF     = 1

;---高原地形填图
  topores                    = res
  topores@cnFillOn           = True             ; turn off contour fill
  topores@cnFillPalette      = "OceanLakeLandSnow"
  topores@cnLinesOn          = False
  topores@cnLevelSelectionMode="ManualLevels" 
  topores@cnLevelSpacingF    = 1000
  topores@cnMinLevelValF     =0            
  topores@cnFillMode         = "RasterFill"     ; much faster than AreaFill
  ;topores@cnLineColor        =(/"black","blue","red"/)
  ;topores@cnLineThicknessF   =10.
  ;topores@cnLevelSelectionMode = "ExplicitLevels" 
  ;topores@cnLevels = (/3000,4000,5000/)
  ;topores@cnLevelSpacingF    = 125              ; 125,NCL picks 2000
  topores@lbBoxLinesOn       = True ; Label bar
  topores@lbLabelAutoStride  = True            ; let NCL determine label spacing
  topores@cnLineLabelsOn     = False            ; turn off line labels
  topores@cnInfoLabelOn      = False            ; turn off info label
  topores@pmLabelBarHeightF  =0.08 
  topores@lbBoxEndCapStyle   ="TriangleBothEnds"
  topores@pmLabelBarOrthogonalPosF  = 0.15 ;色标离图的距离
  topores@lbTitleString                  ="elevation/m"
  topores@lbTitleJust                    ="BottomRight"
  topores@lbTitlePosition                ="Bottom"

  ;---高原3km,4km，5km等值线
topolineres                    = res
topolineres@cnFillOn           = False             ; turn off contour fill
topolineres@cnLinesOn          = True
topolineres@cnLevelSelectionMode="ManualLevels" 
topolineres@cnLevelSpacingF    = 5000
topolineres@cnMinLevelValF     = 3000  
topolineres@cnLevelSpacingF    = 1000
topolineres@cnLineThicknessF   = 8.
topolineres@cnLineLabelsOn     = False            ; turn off line labels
topolineres@cnInfoLabelOn      = False            ; turn off info label
;绘制地图-------------------------
map = gsn_csm_map(wks,mpres)
; topo_plot = gsn_csm_contour(wks,elev,topores)
topo_line_plot=gsn_csm_contour(wks,elev,topolineres)
; overlay(map, topo_plot)
overlay(map, topo_line_plot)

;--绘制散点图---------------------------
  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.
  gsres@gsMarkerColor = "black"
;
; Loop through each grouping of markers, and draw them one set at
; a time, assigning the proper color and size with gsn_marker.

  size=new(num_distinct_markers,float) ;将size存起来制作图例
  thick=new(num_distinct_markers,float)
  base_size = 0.01
  do i = 0, num_distinct_markers-1
    size(i)=base_size * (i+1)/5.
    thick(i)=0.5*(i+1)
  end do
  pmid = new(num_distinct_markers,graphic)
  do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
      gsres@gsMarkerSizeF      = base_size * (i+1)/5.
      gsres@gsMarkerThicknessF = 0.5*(i+1)
      pmid(i) = gsn_add_polymarker(wks,map,lon_new(i,:),lat_new(i,:),gsres)
    end if
  end do
  ;print(thick)
  ;exit

; 绘制降水散点图图例
  resLg = True 
  resLg@lgItemType = "Markers"
  resLg@lgMonoMarkerIndex = True
  ;resLg@lgMarkerColors = colours
  resLg@lgMonoMarkerColor=True  
  resLg@lgMarkerIndex = 16
  resLg@lgMonoLineThickness=False  ;这两个必须关掉，才可以设置后面的
  resLg@lgMonoMarkerSize=False
  resLg@lgMarkerSizes = size(1:)
  resLg@lgMarkerThicknesses=thick(1:)
  resLg@vpWidthF = 0.2
  resLg@vpHeightF = 0.4
  resLg@lgPerimFill = 0;底色白色0，透明为-1
  resLg@lgPerimFillColor = "Background" 
  resLg@lgLabelFontHeightF = 0.06
  resLg@lgTitleFontHeightF = 0.015
  resLg@lgTitleString = "Precipitation/mm"
  resLg@lgBoxMinorExtentF=0.3 ;设置文本离boxline左边的空白

 
  lbid = gsn_create_legend(wks,num_distinct_markers-1, labels(1:), resLg)
 
; 将图例放置在图中
  amres = True
  amres@amParallelPosF = -0.395 ;左右，正右左负
  amres@amOrthogonalPosF = 0.168 ;上下，正下负上
  dumLg = gsn_add_annotation(map, lbid, amres)


;---------------------------某日风暴路径绘制--------------------------------------

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


  resLine =True
  resLine@gsLineThicknessF = 4
  resLine@gsLineColor      = "black"

  txres = True                     ; text mods desired
  txres@txFontHeightF   = 0.03                     ; font smaller. default big
  txres@txFontThicknessF=6
  txres@txFontColor     ="red"


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
    if (((DD(i) .eq. "26").and.(HH(i) .eq. "12"))\
      .or.((DD(i) .eq. "26").and.(HH(i) .eq. "18"))\
      .or.((DD(i) .eq. "27").and.(HH(i) .eq. "00"))\
      .or.((DD(i) .eq. "27").and.(HH(i) .eq. "06"))\
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

;shp="/cygdrive/D/ncl_related/data/TPshp/DBATP/DBATP_Line.shp"  
;  lnres                  = True     
;  lnres@gsLineColor      ="orange"
;  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
;  shp_plot               = gsn_add_shapefile_polylines(wks,map,shp,lnres) 

 
  draw(map)
  frame(wks)
end     
 

