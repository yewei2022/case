;500hPa高度场求距平，某时刻高层绝对涡度分布以及风暴位置标记
;参考函数 计算距平 http://www.ncl.ucar.edu/Document/Functions/Built-in/dim_rmvmean_n.shtml
;计算平均值 http://www.ncl.ucar.edu/Document/Functions/Built-in/dim_avg_n.shtml
;时间该两处，29、72
;代码路径  ncl /mnt/d/ncl_related/code/500h_anrm_200_avo.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin

;1. 30年同一天的数据读入
  hgtfile_dir = "/mnt/d/data/hgt/"
  pic_dir  = "/mnt/d/ncl_related/pictures/500/anomaly/" 

  nc_files = systemfunc("ls "+hgtfile_dir+"hgt.*.nc"); 
  f        = addfiles(nc_files,"r") ; f[0-29],30个文件
  ListSetType(f,"cat"); 指定f的连接类型位cat,因为join不行，因为第一维time365/366长度不等

  ;hgt0    = f[0]->hgt ;为了复制hgt坐标变量
  hgt1    = new((/30,17,73,144/), float)
  ;copy_VarCoords_not_n(hgt0,hgt1,0) ; 除了time 365 以外的坐标信息
  do i= 0,29
    time0    = f[i]->time
    utc_date = cd_calendar(time0, 0)
    ; Store return information into more meaningful variables.
    ;year   = tointeger(utc_date(:,0))    ; Convert to integer for
    month  = tointeger(utc_date(:,1))    
    day    = tointeger(utc_date(:,2)) 
    time_loc= ind((month.eq.10.and.day.eq.26))
    hgt0         = f[i]->hgt
    hgt1(i,:,:,:)= hgt0(time_loc,:,:,:)
    delete([/time0,utc_date,month,day,time_loc,hgt0/])
  end do
;2. 计算距平
  hgt1         = hgt1/10. ;转换为位势什米
  hgt_RmvTime  = dim_rmvmean_n_Wrap(hgt1,0)
  ;time1        = hgt_RmvTime&time
  ;iymd         = cd_calendar(time1, -2)
  ;print(iymd(29)) ;20081026
  ;exit
  hgt_ano          = hgt_RmvTime(29,{500},:,:)
;3 计算均值
;  hgt_AvgTime = dim_avg_n_Wrap(hgt1, 0) 
;  hgt_avg     = hgt_AvgTime({500},:,:)
;4 27日00时的绝对涡度
  year   = 2008
  month  = 10
  day    = 27
  hour   = 00
  minute = 00
  second = 0
  ;设置时间格式读取数据####################################
  units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
  opt = 0
  timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
  file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
  file_time  = cd_string(timeset_ymdh,file_format)
  path_format = "%Y%N%D%H"   ; 台风路径文件中的时间形式
  path_time  = cd_string(timeset_ymdh,path_format)
  ;数据读取和处理
  fnl_dir  ="/mnt/d/data/fnl/"
  data      = addfile(fnl_dir+"fnl_"+file_time+".nc", "r")
  lon_in    = data->lon_0
  lat_in    = data->lat_0
  avo_in    = data->ABSV_P0_L100_GLL0 ({20000},:,:);26层
  ;---------------------------网格点插值--------------------
  lon_out=fspan(55,125,701)
  lat_out=fspan(0, 60, 601)
  avo_out=linint2_Wrap(lon_in, lat_in(::-1), avo_in(::-1,:), True, lon_out, lat_out, 0)
  delete(avo_in)
  avo =avo_out*1.e+5
  copy_VarCoords(avo_out, avo) 

wks = gsn_open_wks("png", pic_dir+"1026_anrm_2700_avo")


res = True
res@gsnDraw     = False
res@gsnFrame    = False
res@gsnMaximize =True
res@gsnLeftString  = ""
res@gsnRightString = ""
res@tmYROn         =False
res@tmXTOn         =False
;mpres@pmTickMarkDisplayMode  = "Always"
res@tmBorderThicknessF  =4    ; border 6x as thick,图框粗细
res@tmXBMajorThicknessF =4
res@tmXBMinorThicknessF =4
res@tmYLMajorThicknessF =4 
res@tmYLMinorThicknessF =4   

; ;设置坐标轴标签粗细
  res@tmEqualizeXYSizes   =True 
  ; 现在一个坐标轴设置好后，
  ;(tmXBLabelFontHeightF,(tmXBMajorLengthF,(tmXBMinorLengthF
  ; 其他坐标轴就不用重复设置了
  res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号


;设置地图
mpres = res
mpres@mpFillOn                = False
mpres@mpMinLatF               = 20.
mpres@mpMaxLatF               = 40.
mpres@mpMinLonF               = 70.
mpres@mpMaxLonF               = 105.
mpres@gsnMajorLatSpacing    =5. 
mpres@gsnMinorLatSpacing    =1
mpres@gsnMajorLonSpacing    =5. 
mpres@gsnMinorLonSpacing    =1. 
mpres@mpDataSetName              = "/mnt/d/ChinaMap/Earth..4"
mpres@mpDataBaseVersion          = "MediumRes"
mpres@mpAreaMaskingOn            = True
mpres@mpMaskAreaSpecifiers       = (/"China"/) ;这个和上个两个命令好像是没啥用
mpres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)



;--------------叠加气候场--------------

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
anores@cnLineLabelFontHeightF=0.02  
anores@gsnContourNegLineDashPattern = 16
anores@cnInfoLabelOn        = False ; set false will not set an information label
anores@cnLevelSelectionMode = "ManualLevels" 
;anores@cnMinLevelValF      = 552
;anores@cnMaxLevelValF      = 588
anores@cnLevelSpacingF      = 1 
anores@cnLineThicknessF     = 6


;--------------绝对涡度--------------

cRes = res
cRes@cnFillOn = True ;cn means ContourPlot
cRes@cnLinesOn = False   
cRes@cnSpanFillPalette = True
cRes@cnFillPalette ="BlueYellowRed"
cRes@cnLineLabelsOn = False  
cRes@cnInfoLabelOn = False ; set false will not set an information label 
; cRes@cnLevelSelectionMode = "ManualLevels" 
; cRes@cnMinLevelValF       = -6             ; set min contour level
; cRes@cnMaxLevelValF       = 6               ; set max contour level
; cRes@cnLevelSpacingF      = 2              ; set contour spacing
cRes@cnLevelSelectionMode   = "ExplicitLevels" 
cRes@cnLevels      =(/-6,-5,-4,-3,-2,-1,0/)
cRes@cnFillColors  =(/5,15,35,45,65,85,115,-1/)

cRes@lbBoxLinesOn = True ; Label bar
cRes@lbLabelAutoStride    = True            ; let NCL determine label spacing
; cRes@lbLabelFontHeightF   = 0.015
; cRes@lbOrientation        = "Vertical"         ; vertical label bar
; cRes@pmLabelBarWidthF = 0.03     
cRes@pmLabelBarHeightF = 0.1    
cRes@pmLabelBarOrthogonalPosF  = 0.15 ;色标离图的距离
cRes@lbBoxEndCapStyle = "TriangleBothEnds"




;绘图
; hgt_ano = smth9_Wrap(hgt_ano,0.5,0.25,False)

mapplot = gsn_csm_map(wks, mpres)
;contour_avg = gsn_csm_contour(wks,hgt_avg, heightRes)
contour_ano = gsn_csm_contour(wks,hgt_ano, anores)
contour_avo = gsn_csm_contour(wks,avo,cRes)
;overlay(mapplot, contour_avg)
overlay(mapplot, contour_avo)
overlay(mapplot, contour_ano)

;高原轮廓线
shp="/mnt/d/data/TPshp/DBATP/DBATP_Line.shp"  
  lnres                  = True     
  lnres@gsLineColor      ="orange"
  lnres@gsLineThicknessF = 8.0           ; 6x thickness  
  shp_plot               = gsn_add_shapefile_polylines(wks,mapplot,shp,lnres) 

;读取风暴中心位置并绘图   
  fiTY = "/mnt/d/data/TC/bio042008.txt"
  ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
  nrow = numAsciiRow(fiTY)
  YYYYMMDDHH = new(nrow, "string")
  lat = new(nrow, "float")
  lon = new(nrow, "float")
  cdata = asciiread(fiTY, -1, "string")
  YYYYMMDDHH = str_get_field(cdata, 3, ",")
  lat = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
  lon = stringtofloat(str_get_field(cdata, 8, ",")) *0.1
  DateChar = stringtochar(YYYYMMDDHH)
  MM = chartostring(DateChar(:,5:6))
  DD = chartostring(DateChar(:,7:8))
  HH = chartostring(DateChar(:,9:10))

  tctime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
  stdlon=0.
  stdlat=0.
  do i = 0, nrow-1
    if (i.eq.tctime_ind) then ;最方便
      stdlon=lon(i)
      stdlat=lat(i)
    end if
  end do
  ; 某时刻风暴位置标记---------------------------------
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.03                     ; font smaller. default big
  txres@txFontThicknessF=6
  txres@txFontColor="red"
  dum = gsn_add_text(wks,mapplot,"~F37~p",stdlon,stdlat,txres)
  ;4 贴上强降水区域的方框
  ypts = (/32.0,  32.0,  28.0, 28.0, 32.0/)
  xpts = (/89.0,  98.0,  98.0, 89.0, 89.0/)

  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 7.0                       ; thickness of lines
  dum_box = new(4,graphic)  
  ; draw each line separately. Each line must contain two points.
  do i = 0 , 3
    dum_box(i)=gsn_add_polyline(wks,mapplot,xpts(i:i+1),ypts(i:i+1),resp)      
  end do  

draw(mapplot)
frame(wks)
end