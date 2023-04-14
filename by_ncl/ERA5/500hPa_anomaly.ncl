;使用ERA5数据，500hPa高度场求距平，
;参考函数 计算距平 http://www.ncl.ucar.edu/Document/Functions/Built-in/dim_rmvmean_n.shtml
;计算平均值 http://www.ncl.ucar.edu/Document/Functions/Built-in/dim_avg_n.shtml
;时间改1处 43
;代码路径  ncl /mnt/h/d/ncl_related/ERA5_code/500hPa_anomaly.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

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

;1. 数据读入
  pic_dir  = "/mnt/h/d/ncl_related/pictures/500/anomaly/" 
  data_dir = "/mnt/h/d/data/ERA5/geopotential/"
  nc_files = systemfunc("ls "+data_dir+"era5.geopotential.*.nc"); 
  f        = addfiles(nc_files,"r") ; f[0-29],30个文件
  ListSetType(f,"cat"); 指定f的连接类型位cat,因为join不行，因为第一维time365/366长度不等
  time0    = f[:]->time
  DDHH     = cd_string(time0, "%D%H")
  day_str = "2700"

  time_loc= str_match_ind(DDHH, day_str)
  ; print(DDHH(time_loc))
  ; exit
  hgt1    = short2flt(f[:]->z(time_loc,{500},:,:)) ;为了复制hgt坐标变量

;2. 计算距平
  hgt1         = hgt1/98. ;转换为位势什米
  hgt_RmvTime  = dim_rmvmean_n_Wrap(hgt1,0)
  time1        = hgt_RmvTime&time
  iymd         = cd_calendar(time1, -2)
  ; print(iymd(29)) ;29-第30年的数据
  ; exit
  hgt_ano          = hgt_RmvTime(29,:,:)


;风暴数据,用于标记风暴位置或路径
fiTY = "/mnt/h/d/data/TC/bio042008.txt"
nrow = numAsciiRow(fiTY)
YYYYMMDDHH = new(nrow, "string")
lat = new(nrow, "float")
lon = new(nrow, "float")
cdata = asciiread(fiTY, -1, "string")
lat_tc = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
lon_tc = stringtofloat(str_get_field(cdata, 8, ",")) *0.1
YYYYMMDDHH = str_get_field(cdata, 3, ",")
DateChar = stringtochar(YYYYMMDDHH)
MM = chartostring(DateChar(:,5:6))
DD = chartostring(DateChar(:,7:8))
HH = chartostring(DateChar(:,9:10))
DDHH_TC=DD+HH


wks = gsn_open_wks("ps", pic_dir+"anormaly_ERA5_"+day_str)


res = True
res@gsnDraw     = False
res@gsnFrame    = False
res@gsnMaximize =True
res@gsnLeftString  = ""
res@gsnRightString = ""
res@tmYROn         =False
res@tmXTOn         =False
;mpres@pmTickMarkDisplayMode  = "Always"
;mpres@tmXBMajorLengthF = 0.005

; ; ;设置坐标轴标签粗细
; res@tmEqualizeXYSizes   =True 
; res@tmXBLabelFontHeightF = 0.015 ; x轴标签字号


;设置地图
mpres = res
mpres@mpFillOn                = False
mpres@mpMinLatF               = 20.
mpres@mpMaxLatF               = 35.
mpres@mpMinLonF               = 75.
mpres@mpMaxLonF               = 100.
mpres@gsnMajorLonSpacing      = 5. 
mpres@gsnMinorLonSpacing      = 1.    
mpres@gsnMajorLatSpacing      = 5. 
mpres@gsnMinorLatSpacing      = 1.
mpres@mpDataBaseVersion       = "MediumRes"
mpres@mpDataSetName           = "/mnt/h/d/ChinaMap/Earth..4"
mpres@mpOutlineSpecifiers     = (/"China:states"/)

mpres@tmXBLabels=(/"75~S~o~N~E","80~S~o~N~E","85~S~o~N~E","90~S~o~N~E",\
                  "95~S~o~N~E","100~S~o~N~E","105~S~o~N~E"/)
mpres@tmYLLabels=(/"20~S~o~N~N","25~S~o~N~N","30~S~o~N~N","35~S~o~N~N"/)

;--------------叠加高度场--------------

heightRes = res
heightRes@cnFillOn        = False  ;cn means ContourPlot
heightRes@cnLinesOn       = True   
;heightRes@cnMonoLineColor= False
heightRes@cnLineColor     ="blue"
;heightRes@cnSpanLinePalette = True
;heightRes@cnLinePalette     = "rainbow"
heightRes@cnLineLabelsOn     = True  
heightRes@cnInfoLabelOn      = False ; set false will not set an information label
heightRes@cnLevelSelectionMode = "ManualLevels" 
;heightRes@cnMinLevelValF      = 552
;heightRes@cnMaxLevelValF      = 588
heightRes@cnLevelSpacingF      = 3 
heightRes@cnLineThicknessF     = 6
;heightRes@gsnSpreadColorStart = 1
;heightRes@gsnSpreadColorStart = 17
;heightRes@gsnSpreadColors     = True ; 6.1.0


;--------------叠加距平场--------------

anores                    = res
anores@cnFillOn           = False  ;cn means ContourPlot
anores@cnLinesOn          = True   
;anores@cnMonoLineColor   = False
anores@cnLineColor        = "blue"
;anores@cnSpanLinePalette = True
;anores@cnLinePalette     = "rainbow"
anores@cnLineLabelsOn     = True 
anores@cnLineLabelFontHeightF=0.015 
anores@gsnContourNegLineDashPattern = 15
anores@cnInfoLabelOn        = False ; set false will not set an information label
anores@cnLevelSelectionMode = "ManualLevels" 
anores@cnMinLevelValF      = -7
;anores@cnMaxLevelValF      = 588
anores@cnLevelSpacingF      = 1 
anores@cnLineThicknessF     = 6


anores@gsnLeftStringFontHeightF=0.025
anores@gsnLeftString="(c)"
anores@gsnLeftStringOrthogonalPosF=-0.11
anores@gsnLeftStringParallelPosF =0.01
anores@txPerimOn=True
anores@txPerimThicknessF=2.5
anores@txBackgroundFillColor ="white"

;绘图
hgt_ano = smth9_Wrap(hgt_ano,0.5,0.25,False)
;hgt_avg = smth9_Wrap(hgt_avg,0.5,0.25,False)


mapplot = gsn_csm_map(wks, mpres)
;contour_avg = gsn_csm_contour(wks,hgt_avg, heightRes)
contour_ano = gsn_csm_contour(wks,hgt_ano, anores)
;overlay(mapplot, contour_avg)
overlay(mapplot, contour_ano)
mark_area(wks,mapplot)
mark_tc(wks,DDHH_TC,day_str,mapplot,nrow,lon_tc,lat_tc)

;高原轮廓线
shp="/mnt/h/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres) 

draw(mapplot)
frame(wks)
end