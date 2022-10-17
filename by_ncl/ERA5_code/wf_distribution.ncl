; 某时刻 wf 分布图
; ncl /mnt/h/d/ncl_related/ERA5_code/wf_distribution.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

undef("mark_topoline")
  procedure mark_topoline(wks,map)
  begin
  shp="/mnt/h/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="black"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,map,shp,lnres) 
  map@$unique_string("shp_plot")$ = shp_plot
  end

begin

;1.读取数据
  ;1.1 变量数据
  pic_dir  ="/mnt/h/d/ncl_related/pictures/terrain/ERA5/wf/"
  data_dir = "/mnt/h/d/data/ERA5/"    
  data2    = addfile(data_dir+"ERA5_wind.nc", "r") 
  timeset_ymdh=data2->time
  path_format = "%Y%N%D%H"   ; YYYYMMMDDHH
  time  = cd_string(timeset_ymdh,path_format)
  need_time="200810"+"2612"
  time_ind=str_match_ind(time, need_time)

  u_in     =short2flt(data2 ->u100(time_ind,{10:35},{85:105}))
  v_in     =short2flt(data2 ->v100(time_ind,{10:35},{85:105}))

  ; ; examine time series;two methods-two var type
  ; time    =data1->time
  ; YYYYMMDDHH_f=cd_calendar(time,-3)


  ;1.2 地形数据
  topo        = addfile("/mnt/h/d/data/other/ETOPO2v2c_f4.nc", "r")
  elev1       = topo->z({15:50},{70:110})
  elev2        = tofloat(elev1)
  copy_VarCoords(elev1, elev2) 
  ; printVarSummary(elev2)
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
  elev2!0    = "lat"
  elev2!1    = "lon"
  elev2&lat  = lat_topo
  elev2&lon  = lon_topo


  ;将水平速度插值成和地形相同分辨率
  elev=elev2({10:35},{85:105})
  ; printVarSummary(elev_test) 600*600
  lon_in =u_in&longitude
  lat_in =u_in&latitude
  lon_out=elev&lon
  lat_out=elev&lat
  u     =linint2_Wrap(lon_in, lat_in, u_in, False, lon_out, lat_out, 0)
  delete(u_in)
  v     =linint2_Wrap(lon_in, lat_in, v_in, False, lon_out, lat_out, 0)
  delete(v_in) 

;2 数据处理 
;; 地形梯度
  nlat=dimsizes(elev&lat)
  mlon=dimsizes(elev&lon)
  lat    =elev&lat
  lon    =elev&lon

  dlon   = (lon(2)-lon(1))*0.0174533 ; convert to radians 0.017=3.14/180
  dHdX   = new ( dimsizes(elev),typeof(elev)); pre-allocate space
  copy_VarCoords(elev, dHdX)
  do nl=0,nlat-1                      ; loop over each latitude
  dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
  dHdX(nl:nl,:) = center_finite_diff_n (elev(nl:nl,:), dX , False,0,1) ; 1指的经度那一维
  end do

  dY   =(lat(2)-lat(1))*6378388.*0.0174533 ;标量
  dHdY = new ( dimsizes(elev),typeof(elev)); pre-allocate space
  copy_VarCoords(elev, dHdY)
  dHdY= center_finite_diff_n (elev, dY , False,0,0)

;2.3 计算 wf
  ; 取站点 现在文章里的图 但是不太合理 位置稍微偏离一点点 dhdy就有很大的变化
  ; 算速度,会自动输出离设定位置最近的点的值

  lat_cona=27.98 ;错那
  lon_cona=91.95
  lat_bomi=29.87;波密
  lon_bomi=95.77
  lat_zayu=28.65 ;察隅
  lon_zayu=97.47

  wf_all=v*dHdY
  copy_VarCoords(v, wf_all)

wks = gsn_open_wks("png", pic_dir+need_time)
 
res = True
res@gsnDraw     = False
res@gsnFrame    = False
res@gsnMaximize = False
res@gsnLeftString  = ""
res@gsnRightString = ""
res@tmYROn         =False
res@tmXTOn         =False
;mpres@pmTickMarkDisplayMode  = "Always"
res@tmBorderThicknessF  =2.5    ; border 6x as thick,图框粗细
res@tmXBMajorThicknessF =2.5
res@tmXBMinorThicknessF =2.5
res@tmYLMajorThicknessF =2.5 
res@tmYLMinorThicknessF =2.5   
; ; ;设置坐标轴标签粗细
; res@tmEqualizeXYSizes   =True 
; ; 现在一个坐标轴设置好后，
; ;(tmXBLabelFontHeightF,(tmXBMajorLengthF,(tmXBMinorLengthF
; ; 其他坐标轴就不用重复设置了
; res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号


;设置地图
mpres = res
mpres@mpFillOn                = False
mpres@mpMinLatF               = 27.5
mpres@mpMaxLatF               = 28.5
mpres@mpMinLonF               = 91.5
mpres@mpMaxLonF               = 92.5
; mpres@gsnMajorLatSpacing    =5. 
; mpres@gsnMinorLatSpacing    =1
; mpres@gsnMajorLonSpacing    =5. 
; mpres@gsnMinorLonSpacing    =1. 
mpres@mpDataSetName              = "/mnt/h/d/ChinaMap/Earth..4"
mpres@mpDataBaseVersion          = "MediumRes"
mpres@mpAreaMaskingOn            = True
mpres@mpMaskAreaSpecifiers       = (/"China"/) ;这个和上个两个命令好像是没啥用
mpres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)
mpres@pmTickMarkDisplayMode ="Always" ;坐标轴标签上添加符号度
mpres@gsnAddCyclic=False

;-------------填色--------------
cRes = mpres
cRes@cnFillOn = True ;cn means ContourPlot
cRes@cnLinesOn = False   
cRes@cnSpanFillPalette = True
cRes@cnFillPalette = "cmp_flux"
cRes@cnLineLabelsOn = False  
cRes@cnInfoLabelOn = False ; set false will not set an information label 
cRes@cnLevelSelectionMode = "ManualLevels" 
cRes@cnMinLevelValF       = -0.5             ; set min contour level
cRes@cnMaxLevelValF       = 0.5               ; set max contour level
cRes@cnLevelSpacingF      = 0.1              ; set contour spacing
;cRes@cnLineThicknessF = 4
;cRes@cnLevelSelectionMode = "ExplicitLevels"
;cRes@cnLevels = (/0, 5, 10, 15,20/)
;cRes@cnFillColors=(/120,150,180,210,240,250/)

cRes@lbBoxLinesOn = True ; Label bar
cRes@lbLabelAutoStride    = True            ; let NCL determine label spacing
; cRes@lbLabelFontHeightF   = 0.015
; cRes@lbOrientation        = "Vertical"         ; vertical label bar
; cRes@pmLabelBarWidthF = 0.03     
cRes@pmLabelBarHeightF = 0.1    
cRes@pmLabelBarOrthogonalPosF  = 0.15 ;色标离图的距离
cRes@lbBoxEndCapStyle = "TriangleBothEnds"


plot = gsn_csm_contour_map(wks,wf_all, cRes)
; mark_topoline(wks,plot)

;站点标记
  mkres               = True                     ; text mods desired
  mkres@txFontColor   ="black"
  mkres@txFontHeightF = 0.045
  mkres@txFontThicknessF=3
  dum_sta = gsn_add_text(wks,plot,"~F35~z",lon_cona,lat_cona,mkres) 


draw(plot)
frame(wks)

end