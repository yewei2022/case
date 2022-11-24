;  4.7 绘制台风路径 + 地形图
; load "./offical_sample_script/download/pvs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
; load "topo_2.ncl"
undef("set_res")

procedure set_all_res(res,mpres,windRes,heightRes,tmRes)
  local res,mpres,winRes,heightRes,tmRes
  begin
  res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@gsnLeftString = ""
    res@gsnRightString = ""
    ;res@vpHeightF = 1
    ;res@vpwidthF = 1

    ;设置地图
  mpres = res
    mpres@mpFillOn                = True ; fillon 后面修改填充颜色
    mpres@mpMinLatF               = 15.
    mpres@mpMaxLatF               = 26.
    mpres@mpMinLonF               = 81.
    mpres@mpMaxLonF               = 95.
    mpres@mpDataBaseVersion       = "MediumRes"
    mpres@mpDataSetName           = "Earth..4"
    ; mpres@mpOutlineSpecifiers     = (/"Oceans"/)
    mpres@pmTickMarkDisplayMode  = "Always"
    mpres@tmXBMajorLengthF = 0.005
    mpres@mpOceanFillColor       = "cyan"            ; array index in color map
    mpres@mpLandFillColor        = 180         ; for land, ocean, and inland h20
    mpres@mpInlandWaterFillColor = "cyan"
    ;mpres@tiMainString = "500hPa_2900"
    ;mpres@gsnDraw = False
    ; mapplot = gsn_csm_map(wks, mpres)

  ;--------------风矢量图--------------
  windRes = res
    ;windRes@vcGlyphStyle            = "FillArrow"   ; turn on curly vectors
    windRes@vcGlyphStyle    = "LineArrow"
    windRes@vcLineArrowThicknessF = 3
    windRes@vcMinDistanceF = 0.05
    windRes@vcRefMagnitudeF = 20.
    windRes@vcRefLengthF    = 0.02
    windRes@vcRefAnnoOn = True
    windRes@vcRefAnnoFontHeightF =0.015
    windRes@vcRefAnnoZone = 1
    windRes@vcRefAnnoSide = "Top"
    windRes@vcRefAnnoJust = "TopRight"
    windRes@vcRefAnnoString1 = "20 m/s"
    windRes@vcRefAnnoString2On = False
    windRes@vcLineArrowHeadMaxSizeF = 0.01
    windRes@vcLineArrowHeadMinSizeF = 0.01
  ;--------------叠加高度场--------------
  heightRes = res

    heightRes@cnFillOn = False  ;cn means ContourPlot
    heightRes@cnLinesOn = True   
    ;heightRes@cnMonoLineColor = False
    heightRes@cnLineColor="blue"
    ;heightRes@cnSpanLinePalette = True
    ;heightRes@cnLinePalette = "rainbow"

    heightRes@cnLineLabelsOn = True  
    heightRes@cnInfoLabelOn = False ; set false will not set an information label
    heightRes@cnLevelSelectionMode = "ManualLevels" 
    ;heightRes@cnMinLevelValF = 552
    ;heightRes@cnMaxLevelValF = 588
    heightRes@cnLevelSpacingF = 4 
    heightRes@cnLineThicknessF = 4
    ;heightRes@gsnSpreadColorStart = 1
    ;heightRes@gsnSpreadColorStart =17
    ;heightRes@gsnSpreadColors = True ; 6.1.0
    heightRes@lbBoxLinesOn = False ; Label bar
    heightRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
    heightRes@lbOrientation        = "Vertical"         ; vertical label bar

  ;--------------叠加温度场--------------
  tmRes = res
    tmRes@cnFillOn = False  ;cn means ContourPlot
    tmRes@cnLinesOn = True   
    tmRes@cnLineColor="red"
    tmRes@cnLineDashPattern="1"
    ;tmRes@cnMonoLineColor = False
    ;tmRes@cnSpanLinePalette = True
    ;tmRes@cnLinePalette = "rainbow"

    tmRes@cnLineLabelsOn = False  
    tmRes@cnInfoLabelOn = False ; set false will not set an information label
    tmRes@cnLevelSelectionMode = "ManualLevels" 
    tmRes@cnLevelSpacingF = 4 
    tmRes@cnLineThicknessF = 2
    tmRes@lbBoxLinesOn = False ; Label bar
    tmRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
    tmRes@lbOrientation        = "Vertical"         ; vertical label bar
end

; ----main code----
begin
;;;>>>手动设置一个时间 用各个数据生成一个带unit的 timeset_ymdh 并且转换为一个字符串 方便文件中读入
  year   = 2008
  month  = 10
  day    = 27
  hour   = 00
  minute = 00
  second = 0
  ;####################################
  units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
  opt = 0
  timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
  file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
  file_time  = cd_string(timeset_ymdh,file_format)
  path_format = "%Y%N%D%H"   ; 路径文件中的形式
  path_time  = cd_string(timeset_ymdh,path_format)

;;;路径数据
  fiTY = "~/practice/lx_ex_data/ye/bio042008.txt"

  ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
  nrow = numAsciiRow(fiTY)
  YYYYMMDDHH = new(nrow, "string")
  MMDDHH     = new(nrow, "string")
  lat        = new(nrow, "float")
  lon        = new(nrow, "float")
  cdata = asciiread(fiTY, -1, "string")
  YYYYMMDDHH = str_get_field(cdata, 3, ",")
  lat = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
  lon = stringtofloat(str_get_field(cdata, 8, ",")) *0.1
  MMDDHH = str_get_cols(YYYYMMDDHH, 5, 10) ; 在图中标识的MMDDHH

;;;>>>plot
  wks = gsn_open_wks("x11", "track"+file_time)  ;用x11方便检查
  ;   gsn_define_colormap(wks,"rainbow")
  res       = True
  mpres     = True
  windRes    = True
  heightRes = True
  tmRes     = True
  set_all_res(res,mpres,windRes,heightRes,tmRes) ; 用一个过程设置所有res,在过程中修改

  mapplot = gsn_csm_map(wks, mpres)
;;;>>>查找字符串，绘制路径
  righttime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
;>>>风暴位置标记
;>>>lines将台风各位置连接
    res_lines                      = True   
    res_lines@gsLineThicknessF     = 3.           ; 3x as thick
    res_lines@gsLineDashPattern    = 2
    dot  = new(nrow,graphic)    ; Make sure each gsn_add_polyxxx call
    line = new(nrow,graphic)    ; is assigned to a unique variable.
    line2= new(nrow,graphic)

    ; Loop through each date and add polylines to the plot.
    do i = 0,nrow-2
        ;  res_lines@gsLineColor           = cols(0)
        xx=(/lon(i),lon(i+1)/)
        yy=(/lat(i),lat(i+1)/)
        line(i) = gsn_add_polyline(wks,mapplot,xx,yy,res_lines)
    end do
;>>>line2将台风位置与text链接
    res_line2 = res_lines
    res_line2@gsLineDashPattern    = 2 ; 修改线形
    res_line2@gsLineThicknessF     = 1.; 修改粗细
    markerres               = True                     ; text mods desired
    markerres@txFontHeightF = 0.015                     ; 台风标记大小
    markerres@txFontColor="black"
    txres               = True                     ; text mods desired
    txres@txFontHeightF = 0.007                     ; 时间text大小
    txres@txFontColor="black"
    ; 设置一个量，标识标签和台风位置的间隔
    interm = 2.2
    latinterm = -0.8
    do i = 0, nrow-1
        ; if (i.ne.nrow-1)
        xx=(/lon(i),lon(i)+interm/)
        yy=(/lat(i),lat(i)+latinterm+0.15/)
        line2(i) = gsn_add_polyline(wks,mapplot,xx,yy,res_line2)
        ; end if
        ;   if (i.eq.righttime_ind) then
        ;     continue
        ;   end if
        dum  = gsn_add_text(wks,mapplot,"~F37~p",lon(i), lat(i),markerres) 
        text = gsn_add_text(wks,mapplot,MMDDHH(i),lon(i)+interm, lat(i)+latinterm,txres)
        interm    = interm *(-1)
        latinterm = latinterm *(-1)
    end do
    ; i = righttime_ind
    ;   txres@txFontColor="red"
    ;   dum  = gsn_add_text(wks,mapplot,"~F37~p",lon(i), lat(i),txres) 
    ;       text = gsn_add_text(wks,mapplot,path_time,lon(i)+2.5, lat(i),txres)
    ;;;叠加地形图


draw(mapplot)
frame(wks)
end
;;;绘图要素另写函数