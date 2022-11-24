;使用ERA5数据
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;代码路径  ncl /mnt/h/d/ncl_related/ERA5_code/400hPa_hgt_tmp_barb.ncl

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
  txres@txFontHeightF = 0.03                     ; font smaller. default big
  txres@txFontThicknessF=6
  txres@txFontColor="red"
  dum_tc = gsn_add_text(wks,map,"~F37~p",stdlon,stdlat,txres) 
  map@$unique_string("dum_tc")$ = dum_tc
  end

undef("mark_topoline")
  procedure mark_topoline(wks,map)
  begin
  shp="/mnt/h/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,map,shp,lnres) 
  map@$unique_string("shp_plot")$ = shp_plot
  end

undef("mark_area")
  procedure mark_area(wks,map)
  begin
  ypts = (/33.0,  33.0,  28.0, 28.0, 33.0/)
  xpts = (/89.0,  98.0,  98.0, 89.0, 89.0/)
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 7.0                       ; thickness of lines
  dum_box = new(4,graphic)  
  ; draw each line separately. Each line must contain two points.
  do i = 0 , 3
    dum_box(i)=gsn_add_polyline(wks,map,xpts(i:i+1),ypts(i:i+1),resp)      
  end do 
  map@$unique_string("dum_box")$ = dum_box
  end


begin

; ；算风暴移动距离
; gcdist = gc_latlon(23.4,90.4,25.0,90.6,2,4)
; print(gcdist)
; exit
; 2008102600,182N,881E;2008102606,187N,884E;
; 2008102612,198N,888E;2008102618,210N,896E;
; 2008102700,234N,904E;2008102706,250N,906E

data_dir="/mnt/h/d/data/ERA5/"   
pic_dir ="/mnt/h/d/ncl_related/pictures/500/ERA5/"
data    = addfile(data_dir+"ERA5.nc", "r")
level   = data->level ;hPa
time    =data->time
YYYYMMDDHH_f=cd_calendar(time,-3)

day="26"
hh="12"
need_time="200810"+day+hh
time_locat = ind(YYYYMMDDHH_f.eq.need_time) ;找该时间所处位置
path_time=tostring(need_time)

pres=400
u = short2flt(data ->u(time_locat,{pres},:,:))  ;单位：m/s
v = short2flt(data->v(time_locat,{pres},:,:))
hgtData = short2flt(data->z(time_locat,{pres},:,:))  ; hgt 的单位：gpm
tmpData =short2flt(data->t(time_locat,{pres},:,: ))  ;单位：K

hgt = hgtData/9.8/10.
copy_VarCoords(hgtData, hgt)
tmp=tmpData-273.15
copy_VarCoords(tmpData, tmp)


;风暴数据,用于标记风暴位置或路径
fiTY = "/mnt/h/d/data/TC/bio042008.txt"
nrow = numAsciiRow(fiTY)
YYYYMMDDHH = new(nrow, "string")
lat = new(nrow, "float")
lon = new(nrow, "float")
cdata = asciiread(fiTY, -1, "string")
YYYYMMDDHH = str_get_field(cdata, 3, ",")
lat_tc = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
lon_tc = stringtofloat(str_get_field(cdata, 8, ",")) *0.1
DateChar = stringtochar(YYYYMMDDHH)
MM = chartostring(DateChar(:,5:6))
DD = chartostring(DateChar(:,7:8))
HH = chartostring(DateChar(:,9:10))


wks = gsn_open_wks("x11", pic_dir+tostring(pres)+path_time)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize=True
res@gsnLeftString = ""
res@gsnRightString = ""
res@tmYROn=False
res@tmXTOn=False
res@tmBorderThicknessF  =2.5    ; border 6x as thick,图框粗细
res@tmXBMajorThicknessF =2.5
res@tmXBMinorThicknessF =2.5
res@tmYLMajorThicknessF =2.5 
res@tmYLMinorThicknessF =2.5  

; ; ;坐标轴标签大小
; res@tmEqualizeXYSizes   =True 
; ; 现在一个坐标轴设置好后，
; ;(tmXBLabelFontHeightF,(tmXBMajorLengthF,(tmXBMinorLengthF
; ; 其他坐标轴就不用重复设置了
; res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号

;设置地图
mpres = res
mpres@mpFillOn                = False
mpres@mpMinLatF               = 15.
mpres@mpMaxLatF               = 35.
mpres@mpMinLonF               = 75.
mpres@mpMaxLonF               = 100.
mpres@gsnMajorLatSpacing    =5. 
mpres@gsnMinorLatSpacing    =1
mpres@gsnMajorLonSpacing    =5. 
mpres@gsnMinorLonSpacing    =1.    
mpres@mpDataSetName              = "/mnt/h/d/ChinaMap/Earth..4"
mpres@mpDataBaseVersion          = "MediumRes"
mpres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)
; mpres@pmTickMarkDisplayMode ="Always" ;坐标轴标签上添加符号度
; mpres@tmXBLabels=(/"75~S~o~N~E","80~S~o~N~E","85~S~o~N~E","90~S~o~N~E",\
;                   "95~S~o~N~E","100~S~o~N~E","105~S~o~N~E"/)
; mpres@tmYLLabels=(/"20~S~o~N~N","25~S~o~N~N","30~S~o~N~N","35~S~o~N~N"/)


mpres@tmXBLabels=(/"75~S~o~N~E","80~S~o~N~E","85~S~o~N~E","90~S~o~N~E",\
                  "95~S~o~N~E","100~S~o~N~E","105~S~o~N~E"/)
mpres@tmYLLabels=(/"15~S~o~N~N","20~S~o~N~N","25~S~o~N~N","30~S~o~N~N","35~S~o~N~N"/)


;--------------风barb图--------------
windRes = res
windRes@vcGlyphStyle    = "WindBarb"
windRes@vcWindBarbScaleFactorF=2 ;将m/s 转换为节；一长杆代表5m/s
windRes@vcWindBarbTickLengthF =0.4 ;指示风速的杆子长短，默认是0.3
windRes@vcWindBarbTickSpacingF=0.125 ;默认分杆之间的间隔为0.125
windRes@vcWindBarbLineThicknessF=3
windRes@vcMinDistanceF = 0.04
windRes@vcRefMagnitudeF = 5.;量级
windRes@vcRefLengthF    = 0.03 ;指示方向的杆子的长短
windRes@vcRefAnnoOn   = True
windRes@vcRefAnnoFontHeightF =0.015
windRes@vcRefAnnoZone = 1
windRes@vcRefAnnoSide = "Top"
windRes@vcRefAnnoJust = "TopRight"
windRes@vcRefAnnoString1 = "5 m/s"
windRes@vcRefAnnoString2On = False
windRes@vcRefAnnoOrthogonalPosF=0
windRes@vcRefAnnoPerimSpaceF=0.08
windRes@vcRefAnnoArrowSpaceF=1
windRes@gsnLeftStringFontHeightF=0.025
windRes@gsnLeftString=tostring(pres)+"hPa "+day+"/"+hh
windRes@gsnLeftStringOrthogonalPosF=-0.11
windRes@gsnLeftStringParallelPosF =0.01
windRes@txPerimOn=True
windRes@txPerimThicknessF=2.5
windRes@txBackgroundFillColor ="white"



;--------------叠加高度场--------------
hgtRes = res
hgtRes@cnFillOn  = False  ;cn means ContourPlot
hgtRes@cnLinesOn = True   
;hgtRes@cnMonoLineColor = False
hgtRes@cnLineColor="blue"
hgtRes@cnLineLabelsOn = True  
hgtRes@cnLineLabelFontHeightF=0.015
hgtRes@cnInfoLabelOn  = False ; set false will not set an information label
hgtRes@cnLevelSelectionMode = "ManualLevels" 
hgtRes@cnMinLevelValF = 552
hgtRes@cnMaxLevelValF = 592
hgtRes@cnLevelSpacingF  = 4
hgtRes@cnLineThicknessF = 6

;--------------叠加温度等值线--------------
tmRes = res
tmRes@cnFillOn = False  ;cn means ContourPlot
tmRes@cnLinesOn = True   
tmRes@cnLineColor="red"
tmRes@cnLineDashPattern="1"
;tmRes@cnMonoLineColor = False
tmRes@cnLineLabelsOn = True  
tmRes@cnLineLabelFontHeightF=0.015
tmRes@cnLineLabelFontColor="red"
tmRes@cnInfoLabelOn = False ; set false will not set an information label
tmRes@cnLineThicknessF = 6
tmRes@cnLevelSelectionMode = "ManualLevels" 
tmRes@cnMinLevelValF =-24
tmRes@cnMaxLevelValF = 8
tmRes@cnLevelSpacingF = 4 

tmRes@gsnLeftStringFontHeightF=0.025
tmRes@gsnLeftString="(a)"
tmRes@gsnLeftStringOrthogonalPosF=-0.11
tmRes@gsnLeftStringParallelPosF =0.01
tmRes@txPerimOn=True
tmRes@txPerimThicknessF=2.5
tmRes@txBackgroundFillColor ="white"

; ;只叠加0摄氏度的
; tmRes@cnLevelSelectionMode = "ExplicitLevels"
; tmRes@cnLevels           =(/0/)


;绘图
; hgt = smth9_Wrap(hgt,0.5,0.25,False)
; hgt = smth9_Wrap(hgt,0.5,0.25,False)
; tmp = smth9_Wrap(tmp,0.5,0.25,False)
; tmp = smth9_Wrap(tmp,0.5,0.25,False)
mapplot = gsn_csm_map(wks, mpres)
; contour_height = gsn_csm_contour(wks,hgt, hgtRes)
; contour_tmp = gsn_csm_contour(wks,tmp, tmRes)
windPlot = gsn_csm_vector(wks, u , v , windRes)
overlay(mapplot, windPlot)
; overlay(mapplot, contour_height)
; overlay(mapplot, contour_tmp)
mark_tc(wks,YYYYMMDDHH,path_time,mapplot,nrow,lon_tc,lat_tc)
mark_topoline(wks,mapplot)
; mark_area(wks,mapplot) 

; ;----------------------绘制站点位置-------------------------------
;   lat_sta=new(3,float)
;   lon_sta=new(3,float)
;   lat_sta!0="lat"
;   lat_sta@units="degrees_north"
;   lon_sta!0="lon"
;   lon_sta@units="degrees_east"
;   lat_sta(0)=27.98 ;错那
;   lon_sta(0)=91.95
;   lat_sta(1)=29.87;波密
;   lon_sta(1)=95.77
;   lat_sta(2)=28.65 ;察隅
;   lon_sta(2)=97.47
;   ;站点标记
;   mkres               = True                     ; text mods desired
;   mkres@txFontColor   ="red"
;   mkres@txFontHeightF = 0.045
;   mkres@txFontThicknessF=3
;   dum_sta = gsn_add_text(wks,mapplot,"~F35~z",lon_sta,lat_sta,mkres) 

draw(mapplot)
frame(wks)
end