;TBB填色，其后6h降水(剔除了<0.1mm以及高原以外站点降水)
; 挑选高原站点 即mask
;高原轮廓线
;labelbar的标标题设置参考 http://www.ncl.ucar.edu/Document/Graphics/Resources/lb.shtml#lbTitlePosition
;-改三处时间15,50,最后风暴位置标记
;代码路径  ncl /mnt/h/d/ncl_related/code/FY2D_awx_6hdot_mask.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin
;-----------------读TBB并处理----------------------------------
 data_dir       ="/mnt/h/d/data/TBB/"
 day_str="27"
 hh_str="06"
 rain_day="27"
 rain_hh="12"

 file_name      ="FY2D_TBB_IR1_OTG_200810"+day_str+"_"+hh_str+"30"   
 pictures_dir   ="/mnt/h/d/ncl_related/pictures/TBB/" 
  
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
  picture_path="/mnt/h/d/ncl_related/pictures/rain/"    ;图片存放
  filepath = "/mnt/h/d/data/precipitation/6h/"  ;数据
  time_num="200810"+rain_day+rain_hh
  prefix="SURF_WEA_PRE_6HOUR_"
  f=asciiread(filepath+prefix+time_num+".TXT",-1,"string")
  station = stringtofloat(str_get_field(f(1::),5," ")) ;从数组f中获取站号(第5列),第二行开始读取
  lon = stringtofloat(str_get_field(f(1::),7," ")) ;获取经度值，第7列
  lat = stringtofloat(str_get_field(f(1::),6," ")) ;获取纬度值
  R1  = stringtofloat(str_get_field(f(1::),9," ")) ;获取6h降水

 npts = numAsciiRow(filepath+prefix+time_num+".TXT")-1 ;获取站点数
 R1@_FillValue =  9999.

 ;---Open shapefile and read Mississippi River Basin lat/lon values.
  f1      = addfile("/mnt/h/d/data/TPshp/DBATP/DBATP_Polygon.shp", "r")
  mrb_lon = f1->x
  mrb_lat = f1->y

;---Start with data all missing
  R = new(dimsizes(R1),typeof(R1),R1@_FillValue)

 do i=0,npts-1
   R(i)=where(gc_inout(lat(i),lon(i),mrb_lat,mrb_lon), R1(i), R1@_FillValue)
 end do
 R=where(R.ge.0.1, R, R@_FillValue) ;剔除了0.1mm的降水

;------------------------------ 
; Create X and Y arrays to hold the points for each range and initialize
; them to missing values.  We want to use num_distinct_markers
; different colors, so we need num_distinct_markers sets of X and
; Y points.

  levels = (/0.,5.,10.,20.,30.,40.,50./) ; bin settings (bin0 = < 0., ; bin1 = 0.:4.999, etc.)
  nlevels = dimsizes(levels)
  num_distinct_markers = nlevels+1        ; number of distinct markers
  lat_new = new((/num_distinct_markers,npts/),float,-999);设置成二维，标记*站点数
  lon_new = new((/num_distinct_markers,npts/),float,-999)
;
; Group the points according to which range they fall in. At the
; same time, create the label that we will use later in the labelbar

  labels = new(num_distinct_markers,string)  ; Labels for legend. ;这里只是在设置图例而已  
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



;---- Begin to plot.
wks = gsn_open_wks("x11",pictures_dir+day_str+hh_str)

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

;设置坐标轴标签粗细，图例在左时候设置的
  ; res@vpWidthF=1.0  ;画布大小
  ; res@vpHeightF=0.6 


;---- Map Set

 ;绘制省界
 mpres                                =res
 mpres@mpDataSetName                  = "/mnt/h/d/ChinaMap/Earth..4"
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
 mpres@tmEqualizeXYSizes=True
 mpres@tmXBLabelFontHeightF = 0.02 ; x轴标签字号

 ;地图边界
 mpres@mpMinLonF                    = 80
 mpres@mpMaxLonF                    = 105
 mpres@mpMinLatF                    = 15
 mpres@mpMaxLatF                    = 40
 mpres@gsnMajorLonSpacing    =5. 
 ; mpres@gsnMinorLonSpacing    =1.    
 mpres@gsnMajorLatSpacing    =5. 
 ; mpres@gsnMinorLatSpacing    =1
 mpres@pmTickMarkDisplayMode ="Always" ;坐标轴标签上添加符号度


;tbb填图
 tbbres                               = res
 tbbres@gsnLeftString                =day_str+hh_str
 tbbres@gsnLeftStringFontHeightF     =0.025
 tbbres@gsnLeftStringOrthogonalPosF  =-0.12
 tbbres@gsnLeftStringParallelPosF    =0.01
 tbbres@txPerimOn = True ;标签背景框
 tbbres@txPerimThicknessF=2.5
 tbbres@txBackgroundFillColor ="white"

 tbbres@cnFillOn                       = True
 tbbres@cnFillMode                     = "RasterFill" ;矩形栅格填充
 tbbres@cnLinesOn                      = False
 ;----颜色设置
  tbbres@cnLevelSelectionMode = "ExplicitLevels"
  tbbres@cnLevels = (/-32,-40,-52,-60,-70,-80/) 
  ; tbbres@cnFillPalette= "MPL_cool" 
  ; tbbres@cnFillColors=(/0,18,36,54,90,108,-1/)
  tbbres@cnFillPalette= "perc2_9lev" 
  tbbres@cnFillColors=(/8,7,6,5,4,2,-1/)


 ;----------------------------色标在右的设置-------------------------------------
 ; tbbres@lbOrientation                  = "vertical"            ; 垂直坐标
 ; tbbres@pmLabelBarOrthogonalPosF       = 0.02 ;色标条与图之间的距离
 ; tbbres@pmLabelBarParallelPosF         = 0.45 ;把lb上下移动，上正下负
 ; tbbres@pmLabelBarWidthF               = 0.08     ;色标变窄
 ; tbbres@pmLabelBarHeightF              = 0.6
  ; tbbres@lbTitleFontHeightF             = 0.02
 ;-----------------------------------------------------------------

;---------------------------------------------------------------- 
 tbbres@pmLabelBarOrthogonalPosF       = 0.02 ;上下移动
 tbbres@pmLabelBarParallelPosF         = 0.5 ;把lb左右移动，上正下负
 tbbres@pmLabelBarWidthF               = 0.55     
 tbbres@pmLabelBarHeightF              = 0.08
 tbbres@lbTitleFontHeightF             = 0.024
 tbbres@lbLabelFontHeightF             = 0.022
;------------------------------------------------------------
 tbbres@lbBoxLinesOn                   = True ; 色标之间分界线
 tbbres@lbBoxEndCapStyle               = "TriangleBothEnds"
 tbbres@lbTitleOn                      = True
 tbbres@lbTitleString                  = "TBB/~S~o~N~C"
 tbbres@lbTitleJust                    = "TopRight"
 tbbres@lbTitlePosition                = "Bottom"



;绘图
mapplot = gsn_csm_map(wks, mpres)
tbb_plot = gsn_csm_contour(wks,tbb,tbbres)
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
  resLg@lgTitleString = "Precipitation(mm)"
  resLg@lgBoxMinorExtentF=0.2  ;0.3 ;设置文本离boxline左边的空白
 
  lbid = gsn_create_legend(wks,num_distinct_markers-1, labels(1:), resLg)
 
; 将图例放置在图中
  amres = True
  amres@amParallelPosF = 0.75  ;-0.384 原在左的设置;左右，正右左负
  amres@amOrthogonalPosF =0.    ;0.168 原在左的设置;上下，正下负上

  dumLg = gsn_add_annotation(mapplot, lbid, amres)

;某时刻风暴位置标记

 ;风暴数据,用于标记风暴位置
 
 fiTY = "/mnt/h/d/data/TC/bio042008.txt"

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
   if (HH(i) .eq. hh_str.and.DD(i) .eq. day_str) then
      dum = gsn_add_text(wks,mapplot,"~F37~p",lon1(i), lat1(i),txres) 
    end if
 end do

;高原轮廓线
shp="/mnt/h/d/data/TPshp/DBATP/DBATP_Line.shp"  
 lnres                  = True     
 lnres@gsLineColor      ="orange"
 lnres@gsLineThicknessF = 8.0           ; 6x thickness  
 shp_plot               = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres) 
        
draw(mapplot)
frame(wks)

end 
