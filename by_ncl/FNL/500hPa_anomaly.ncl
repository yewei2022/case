;使用NCAR数据，500hPa高度场求距平，
;参考函数 计算距平 http://www.ncl.ucar.edu/Document/Functions/Built-in/dim_rmvmean_n.shtml
;计算平均值 http://www.ncl.ucar.edu/Document/Functions/Built-in/dim_avg_n.shtml
;时间改1处 43
;代码路径  ncl /mnt/e/d/ncl_related/code/500hPa_anomaly.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 


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
  pic_dir  = "/mnt/e/d/ncl_related/pictures/500/anomaly/" 
  data_dir = "/mnt/e/d/data/hgt/"
  nc_files = systemfunc("ls "+data_dir+"hgt.*.nc"); 
  f        = addfiles(nc_files,"r") ; f[0-29],30个文件
  ListSetType(f,"cat"); 指定f的连接类型位cat,因为join不行，因为第一维time365/366长度不等
  time0    = f[:]->time
  utc_date = cd_calendar(time0, 0)
  month   = toint(utc_date(:,1))
  month1  = tostring(month)    
  day     = toint(utc_date(:,2))
  day1    = tostring(day)  
  MMDD    = month1+day1
  day_str = "1027"
  time_loc= str_match_ind(MMDD, day_str)
; print(time_loc)
; exit
  hgt1    = f[:]->hgt(time_loc,{500},{0:50},{65:110}) ;为了复制hgt坐标变量
  delete([/time0,utc_date,month,day,time_loc/])

;2. 计算距平
  hgt1         = hgt1/10. ;转换为位势什米
  hgt_RmvTime  = dim_rmvmean_n_Wrap(hgt1,0)
  time1        = hgt_RmvTime&time
  iymd         = cd_calendar(time1, -2)
  ; print(iymd(29)) ;29-第30年的数据
  ; exit
  hgt_ano          = hgt_RmvTime(29,:,:)

;3 计算均值
;  hgt_AvgTime = dim_avg_n_Wrap(hgt1, 0) 
;  hgt_avg     = hgt_AvgTime({500},:,:)



wks = gsn_open_wks("ps", pic_dir+"anormaly_"+day_str)


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

; ;设置坐标轴标签粗细
; res@tmXBLabelFont=22
; res@tmYLLabelFont=22
; res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
; res@tmYLLabelFontHeightF = 0.02


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
mpres@mpDataSetName           = "/mnt/e/d/ChinaMap/Earth..4"
mpres@mpOutlineSpecifiers     = (/"China:states"/)



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
anores@gsnContourNegLineDashPattern = 16
anores@cnInfoLabelOn        = False ; set false will not set an information label
anores@cnLevelSelectionMode = "ManualLevels" 
;anores@cnMinLevelValF      = 552
;anores@cnMaxLevelValF      = 588
anores@cnLevelSpacingF      = 0.5 
anores@cnLineThicknessF     = 6
;anores@gsnSpreadColorStart = 1
;anores@gsnSpreadColorStart = 17
;anores@gsnSpreadColors     = True ; 6.1.0




;绘图
hgt_ano = smth9_Wrap(hgt_ano,0.5,0.25,False)
;hgt_avg = smth9_Wrap(hgt_avg,0.5,0.25,False)


mapplot = gsn_csm_map(wks, mpres)
;contour_avg = gsn_csm_contour(wks,hgt_avg, heightRes)
contour_ano = gsn_csm_contour(wks,hgt_ano, anores)
;overlay(mapplot, contour_avg)
overlay(mapplot, contour_ano)
mark_area(wks,mapplot)

;高原轮廓线
shp="/mnt/e/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres) 

draw(mapplot)
frame(wks)
end