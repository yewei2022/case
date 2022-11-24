;---------------------------------------------------
; Plot FY2D TBB at 200810262130_0000 (UTC)
;地形绘制参考http://www.ncl.ucar.edu/Applications/topo.shtml   topo_1.ncl
;labelbar的标标题设置参考 http://www.ncl.ucar.edu/Document/Graphics/Resources/lb.shtml#lbTitlePosition
;-改三处时间50，90,98，--------------------------------------------------
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
   ;printVarSummary(elev)
   ;exit

;-----------------读TBB并处理----------------------------------
 data_dir       ="/cygdrive/D/ncl_related/data/TBB/"
 file_name      ="FY2D_TBB_IR1_OTG_20081026_1230"   
 pictures_dir   ="/cygdrive/D/ncl_related/pictures/TBB/" 
  
 nlat           =1201
 nlon           =1201
 setfileoption("bin", "ReadByteOrder", "LittleEndian")
 tbb            = fbindirread (data_dir+file_name+".grd", 0,(/nlat,nlon/), "integer")
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
 ;printVarSummary(tbb)
 ;exit

;---------------------------读降水并处理--------------------------------------------------------------------


;读取数据

  picture_path="/cygdrive/D/ncl_related/pictures/rain/"    ;图片存放
  filepath = "/cygdrive/D/ncl_related/data/micaps_rain/"  ;数据
  time_num="08102702"
  f=asciiread(filepath+time_num+".000",-1,"string")
  station = stringtofloat(str_get_field(f(3::),1," ")) ;从数组f中获取站号
  lon = stringtofloat(str_get_field(f(3::),2," ")) ;获取经度值
  lat = stringtofloat(str_get_field(f(3::),3," ")) ;获取纬度值
  R  = stringtofloat(str_get_field(f(3::),13," ")) ;获取6h降水


 npts=3519  ;时间LST2702=utc2618=2612之后6h,降水量2702(3519),2714(3831),2720(3617)
 R@_FillValue =  9999

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



;---- Begin to plot.
wks = gsn_open_wks("png",pictures_dir+"dot_"+file_name)




res = True
 res@gsnDraw                        = False
 res@gsnFrame                       = False
 ;res@gsnMaximize                    = True
 ;res@gsnPaperOrientation            = "portrait"
 ;res@gsnAddCyclic                   = False         ; regional data, don't add longitude cyclic point
 ;res@pmTickMarkDisplayMode          = "Always" ;关掉不需要的刻度,次刻度线和右上刻度线
 res@gsnLeftString=""
 res@gsnRightString=""
 res@tmYROn=False
 res@tmXTOn=False

;设置坐标轴标签粗细
  res@tmXBLabelFont=22
  res@tmYLLabelFont=22
  res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
  res@tmYLLabelFontHeightF =0.02

  res@vpWidthF=1.0  ;画布大小
  res@vpHeightF=0.6


;---- Map Set

 ;绘制省界
 mpres                                =res
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
 mpres@mpMinLonF                    = 70
 mpres@mpMaxLonF                    = 110
 mpres@mpMinLatF                    = 15
 mpres@mpMaxLatF                    = 45
 mpres@gsnMajorLonSpacing    =5. 
 mpres@gsnMinorLonSpacing    =1.    
 mpres@gsnMajorLatSpacing    =5. 
 mpres@gsnMinorLatSpacing    =1

;--散点属性------------Create logical variables to hold the marker resources.
  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.

;tbb填图
 tbbres                                = res
 tbbres@cnFillOn                       = True
 tbbres@cnFillMode                     = "RasterFill" ;矩形栅格填充
 tbbres@cnLinesOn                      = False
 tbbres@pmLabelBarOrthogonalPosF       = 0.02 ;色标条与图之间的距离
 tbbres@pmLabelBarParallelPosF         = 0.45 ;把lb上下移动，上正下负
 tbbres@pmLabelBarWidthF               = 0.08     ;色标变窄
 tbbres@pmLabelBarHeightF              = 0.6
 tbbres@lbBoxLinesOn                   = True ; 色标之间分界线
 tbbres@lbOrientation                  = "vertical"            ; 垂直坐标
 tbbres@lbBoxEndCapStyle               = "TriangleBothEnds"
 tbbres@lbTitleOn                      = True
 tbbres@lbTitleString                  = "TBB/~S~o~N~C"
 tbbres@lbTitleJust                    = "TopRight"
 tbbres@lbTitlePosition                = "Bottom"
 tbbres@lbTitleFontHeightF             = 0.02

;----颜色设置
  tbbres@cnFillPalette= "MPL_cool" 
  tbbres@cnLevelSelectionMode = "ExplicitLevels"
  tbbres@cnLevels = (/-32,-40,-52,-60,-70,-80/) 
  ;tbbcolors=ispan(0, 125, 18)
  tbbres@cnFillColors=(/0,18,36,54,90,108,-1/)

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

;绘图
mapplot = gsn_csm_map(wks, mpres)
tbb_plot = gsn_csm_contour(wks,tbb,tbbres)
topo_plot = gsn_csm_contour(wks,elev,topores)
overlay(mapplot, topo_plot)
overlay(mapplot, tbb_plot)

;--------------------画降水散点图---------------------------------
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
      pmid(i) = gsn_add_polymarker(wks,mapplot,lon_new(i,:),lat_new(i,:),gsres)
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
  dumLg = gsn_add_annotation(mapplot, lbid, amres)

;某时刻风暴位置标记

 ;风暴数据,用于标记风暴位置
 
 fiTY = "/cygdrive/D/ncl_related/data/TC/bio042008.txt"

 ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
 nrow = numAsciiRow(fiTY)
 
 YYYYMMDDHH = new(nrow, "string")
 lat1 = new(nrow, "float")
 lon1 = new(nrow, "float")
 
 cdata = asciiread(fiTY, -1, "string")
 YYYYMMDDHH = str_get_field(cdata, 3, ",")
 lat1 = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
 lon1 = stringtofloat(str_get_field(cdata, 8, ",")) *0.1

 DateChar = stringtochar(YYYYMMDDHH)
 MM = chartostring(DateChar(:,5:6))
 DD = chartostring(DateChar(:,7:8))
 HH = chartostring(DateChar(:,9:10))
 txres               = True                     ; text mods desired
 txres@txFontHeightF = 0.03                     ; font smaller. default big
 txres@txFontThicknessF=6
 txres@txFontColor="red"
 do i = 0, nrow-1
   if (HH(i) .eq. "12".and.DD(i) .eq. "26") then
      dum = gsn_add_text(wks,mapplot,"~F37~p",lon1(i), lat1(i),txres) 
    end if
 end do
;shp="/cygdrive/D/ncl_related/data/TPshp/DBATP/DBATP_Line.shp"  
;  lnres                  = True     
;  lnres@gsLineColor      ="orange"
;  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
;  shp_plot               = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres) 
        
draw(mapplot)
frame(wks)

end 
