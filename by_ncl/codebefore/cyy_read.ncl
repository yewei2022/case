;  4.5 处理台风数据 测试读入数据
;    中心 + 等高线 + 风场资料 合成图
; load "./offical_sample_script/download/pvs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
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
    mpres@mpFillOn                = False
    mpres@mpMinLatF               = 0.
    mpres@mpMaxLatF               = 60.
    mpres@mpMinLonF               = 50.
    mpres@mpMaxLonF               = 140.
    mpres@mpDataBaseVersion       = "MediumRes"
    mpres@mpDataSetName           = "Earth..4"
    mpres@mpOutlineSpecifiers     = (/"China:states"/)
    mpres@pmTickMarkDisplayMode  = "Always"
    mpres@tmXBMajorLengthF = 0.005
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
  ;###################################
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
;;;>>>数据读入
  ;场数据
  path_rawdata = "./practice/lx_ex_data/ye/" ; 修改路径
  f_1 = addfile(path_rawdata+"fnl_"+file_time+".nc", "r") ; 将字符串写进文件名中，方便修改 fnl_20081027_00_00.nc
  printFileVarSummary (f_1,"HGT_P0_L100_GLL0"); 注意原始数据中有些变量如 LAND_P0_L1_GLL0 有时间属性 
  test_levels = f_1->$"HGT_P0_L100_GLL0"$
  print(test_levels&lv_ISBL0)
  uData   = f_1->UGRD_P0_L100_GLL0(13,:,:)  ;单位：m/s
  vData   = f_1->VGRD_P0_L100_GLL0(13,:,:)
  hgtData = f_1->HGT_P0_L100_GLL0(13,:,:)  ;13={50000}  hgt 的单位：gpm
  tmpData = f_1->TMP_P0_L100_GLL0 (13,:,: )  ;单位：K

  hgt = hgtData/10 ; 修改为位势十米
  copy_VarCoords(hgtData, hgt)
  tmp=tmpData-273.15
  copy_VarCoords(tmpData, tmp)
  ; 平滑处理
  hgt = smth9_Wrap(hgt,0.5,0.25,False)
  hgt = smth9_Wrap(hgt,0.5,0.25,False)
  tmp = smth9_Wrap(tmp,0.5,0.25,False)
  tmp = smth9_Wrap(tmp,0.5,0.25,False)
  ;路径数据
  fiTY = "./practice/lx_ex_data/ye/bio042008.txt"

  ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
  nrow = numAsciiRow(fiTY)
  YYYYMMDDHH = new(nrow, "string")
  lat        = new(nrow, "float")
  lon        = new(nrow, "float")
  cdata = asciiread(fiTY, -1, "string")
  YYYYMMDDHH = str_get_field(cdata, 3, ",")
  lat = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
  lon = stringtofloat(str_get_field(cdata, 8, ",")) *0.1

  DateChar = stringtochar(YYYYMMDDHH) ; 路径数据字符串
;;;>>>plot
  wks = gsn_open_wks("x11", path_rawdata+"500hPa"+file_time)  ;用x11方便检查
  gsn_define_colormap(wks,"rainbow")
  res       = True
  mpres     = True
  windRes    = True
  heightRes = True
  tmRes     = True
  set_all_res(res,mpres,windRes,heightRes,tmRes) ; 用一个过程设置所有res,在过程中修改

  mapplot = gsn_csm_map(wks, mpres)

  contour_height = gsn_csm_contour(wks,hgt, heightRes)
  contour_tmp = gsn_csm_contour(wks,tmp, tmRes)
  windPlot = gsn_csm_vector(wks, uData , vData , windRes)
  ;叠加风矢量和位势高度
  overlay(mapplot, contour_height)
  overlay(mapplot, contour_tmp)
  overlay(mapplot, windPlot)

;;;>>>end plot
;;;>>>查找字符串，绘制路径
  righttime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
;风暴位置标记
  pres                  = True               ; polyline resources
  pres@gsLineColor      = "black"
  pres@txFontHeightF = 0.15        
  line = gsn_add_polyline(wks,mapplot,lat,lon,pres)      ; draw the traj
    
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.009                     ; font smaller. default big
  txres@txFontColor="black"
do i = 0, nrow-1
      if (i.eq.righttime_ind) then
        continue
      end if
      dum  = gsn_add_text(wks,mapplot,"~F37~p",lon(i), lat(i),txres) 
      ; text = gsn_add_text(wks,mapplot,"00",lon(i), lat(i),txres)
end do
i = righttime_ind
  txres@txFontColor="red"
  dum  = gsn_add_text(wks,mapplot,"~F37~p",lon(i), lat(i),txres) 
      text = gsn_add_text(wks,mapplot,path_time,lon(i), lat(i)+0.15,txres)
draw(mapplot)
frame(wks)
end
;;;绘图要素另写函数
