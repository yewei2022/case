; Figure6: (a) 2013Äê6ÔÂ3ÈÕ-5ÈÕ12-24Ê±LST¹Û²âÀÛ»ý½µË® (É¢µã)
;                  (b) Ä£ÄâµÄºó7ÌìÆ½¾ùµÄ12-24h(LST)ÀÛ»ý½µË®(ÌîÉ«)£»
;                  (c) 2013Äê6ÔÂ3ÈÕ-5ÈÕ12-24Ê±LST¹Û²â³öÏÖ¼«¶ËÐ¡Ê±½µË®µÄÕ¾µã(ÐÇ±ê)ºÍCTRLÄ£Äâºó7Ìì³öÏÖ>60mm¡¤h-1½µË®µÄÆµ´Î·Ö²¼(ÌîÉ«)£»
;                  (d)2013Äê6ÔÂ3ÈÕ-5ÈÕ00-24Ê±LST¹Ø¼üÇøÇøÓòÆ½¾ùµÄÐ¡Ê±½µË®(À¶ÏßÊÇ¹Û²â£¬ºìÏßÊÇCTRLÄ£Äâ)



     ;---To plot multiple lines, you must put them into a mulidimensional array.

     ;  data_xy     = new((/2,20/),float)

     ;     filepath1 = "./air_quality_obs_20170506_09-0507_12.txt"  ;´ò¿ªÎÄ¼þ
         
     ;     air_data = readAsciiTable(filepath1, 301, "float", 3) ;301ÁÐ£¬ÊÇ20¸öÐ¡Ê±µÄ²»Í¬ÖÖÀàµÄÅ¨¶È
     ;                                                           ;ÕâÀïÁÐÊý²»¶ÔµÄ»°¶ÁÊý¾Ý»á¶Á´í
     ;     printVarSummary(air_data)
     
     ; data_xy(1,:) = air_data(6,1::15) ;obs 6µÄÕ¾ºÅÊÇ13510£¬»¨¶¼Ê¦·¶; 8µÄÕ¾ºÅÊÇ13530£¬¾ÅÁúÕòÕòÁú
     ; data_xy@_FillValue = -999

     ; 用来查找局部最大最小值！！！
     
      filename = (/"2017-05-06_09:00:00","2017-05-06_10:00:00","2017-05-06_11:00:00",\
                        "2017-05-06_12:00:00","2017-05-06_13:00:00","2017-05-06_14:00:00","2017-05-06_15:00:00", \
                        "2017-05-06_16:00:00","2017-05-06_17:00:00","2017-05-06_18:00:00","2017-05-06_19:00:00", \
                        "2017-05-06_20:00:00","2017-05-06_21:00:00","2017-05-06_22:00:00","2017-05-06_23:00:00",\
                        "2017-05-07_00:00:00","2017-05-07_01:00:00","2017-05-07_02:00:00","2017-05-07_03:00:00",\
                        "2017-05-07_04:00:00","2017-05-07_05:00:00","2017-05-07_06:00:00"/)
     
      ; obs station max hourly rainfall 两个站点的降水序列， 0610UTC-0706UTC
     G1059_huadu_10_06UTC = (/0,0,0,0,0,0,0,1,100.1,116.3,70.5,0.4,5.8,10.4,7.3,9.8,1.8,0.1,0.1,0,0/)
     G3322_zengcheng_10_06UTC = (/0,0,0,0,0,0,0,0,0,0,2.7,16.3,184.4,150.3,47.9,9.7,5.5,7.3,9.4,1.3,0.6/)

        ijres = True
        ijres@returnInt = True  ;whether return integer
        ijres@useTime   = 0     ;which time in the file should be used when extracting XLAT/XLONG arrays

; ;;;;;;;;;;;;;;;;;;;;;;;;;;;storm 1;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
      lonStart = 113.4
      lonEnd = 113.7
      latStart = 23.1
      latEnd = 23.4

       ijres = True
       ijres@returnInt = True  ;whether return integer
       ijres@useTime   = 0     ;which time in the file should be used when extracting XLAT/XLONG arrays

     ; first calculate max point location, max time is 19-20UTC
     a              = addfiles( "../wrfout_d01_" + "2017-05-06_15:00:00", "r")
     times          = wrf_user_getvar(a, "times", -1) 
     lat         = wrf_user_getvar(a, "XLAT",0)
     lon         = wrf_user_getvar(a, "XLONG",0)
     rainnc = wrf_user_getvar(a,"RAINNC",0)
     ntimes         = dimsizes(times) 
     printVarSummary(rainnc)

     a_next              = addfiles( "../wrfout_d01_" + "2017-05-06_16:00:00", "r") 
     rainnc_next = wrf_user_getvar(a_next,"RAINNC",0)

     ; calculate hourly rainfall and its location
     hourly_rain = rainnc_next-rainnc
     hourly_rain@lat2d = lat
     hourly_rain@lon2d = lon

     loc  = wrf_user_ll_to_xy(a,(/lonStart,lonEnd/),(/latStart,latEnd/),ijres) ; (0,0)-start_x, (0,1)-end_x, (1,0)-start_y, (1,1)-end_y
                             ; a region where max rainfall appears 关注的核心区域，最大值在这个里面出现
     hourly_rain_regional := hourly_rain(loc(1,0):loc(1,1),loc(0,0):loc(0,1))

     ;---Convert to 1D
     hourly_rain_1D      = ndtooned(hourly_rain(loc(1,0):loc(1,1),loc(0,0):loc(0,1)))
     dsizes_hourly_rain = dimsizes(hourly_rain(loc(1,0):loc(1,1),loc(0,0):loc(0,1)))
     ;---Resolve the 1D indices back to their original 3D array. get the max rain location in the region;最大值在这个区域中的坐标位置
     max_hourly_rain_indices  = ind_resolve(maxind(hourly_rain_1D),dsizes_hourly_rain)

     print(hourly_rain(loc(1,0)+max_hourly_rain_indices(0,0),loc(0,0)+max_hourly_rain_indices(0,1)))

     ; storage hourly rain series of a point 存放这个点的小时降水序列
     storm1_max_hourly_rain = new((/dimsizes(filename)-1/),float)
     do it = 0, dimsizes(filename)-1-1
   
 
     a              = addfiles( "../wrfout_d01_" + filename(it), "r")
     times          = wrf_user_getvar(a, "times", -1) 
     lat         = wrf_user_getvar(a, "XLAT",0)
     lon         = wrf_user_getvar(a, "XLONG",0)
     rainnc = wrf_user_getvar(a,"RAINNC",0)
     ntimes         = dimsizes(times) 
     printVarSummary(rainnc)

     a_next              = addfiles( "../wrfout_d01_" + filename(it+1), "r") 
     rainnc_next = wrf_user_getvar(a_next,"RAINNC",0)

     ; calculate hourly rainfall
     hourly_rain = rainnc_next-rainnc
     hourly_rain@lat2d = lat
     hourly_rain@lon2d = lon 

     printVarSummary(hourly_rain)

                                           ;regional y start + indice y            ;regional x start + indice x
 
     storm1_max_hourly_rain(it) = hourly_rain(loc(1,0)+max_hourly_rain_indices(0,0),loc(0,0)+max_hourly_rain_indices(0,1))

     end do
     

;;;;;;;;;;;;;;;;;;;;;;;;;;;storm 2;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   

     lonStart = 112.8
     lonEnd = 113.0
     latStart = 23.4
     latEnd = 23.6

     ; first calculate max point location, max time is 19-20UTC
     a              = addfiles( "../wrfout_d01_" + "2017-05-06_21:00:00", "r")
     times          = wrf_user_getvar(a, "times", -1) 
     lat         = wrf_user_getvar(a, "XLAT",0)
     lon         = wrf_user_getvar(a, "XLONG",0)
     rainnc = wrf_user_getvar(a,"RAINNC",0)
     ntimes         = dimsizes(times) 
     printVarSummary(rainnc)

     a_next              = addfiles( "../wrfout_d01_" + "2017-05-06_22:00:00", "r") 
     rainnc_next = wrf_user_getvar(a_next,"RAINNC",0)

     ; calculate hourly rainfall and its location
     hourly_rain := rainnc_next-rainnc
     hourly_rain@lat2d = lat
     hourly_rain@lon2d = lon


     loc  := wrf_user_ll_to_xy(a,(/lonStart,lonEnd/),(/latStart,latEnd/),ijres) ; (0,0)-start_x, (0,1)-end_x, (1,0)-start_y, (1,1)-end_y
            print(max(hourly_rain(loc(1,0):loc(1,1),loc(0,0):loc(0,1))))


                    ; a region where max rainfall appears 关注的核心区域，最大值在这个里面出现
     hourly_rain_regional := hourly_rain(loc(1,0):loc(1,1),loc(0,0):loc(0,1))

     ;---Convert to 1D
     hourly_rain_1D      := ndtooned(hourly_rain(loc(1,0):loc(1,1),loc(0,0):loc(0,1)))
     dsizes_hourly_rain := dimsizes(hourly_rain(loc(1,0):loc(1,1),loc(0,0):loc(0,1)))
     ;---Resolve the 1D indices back to their original 3D array. get the max rain location in the region;最大值在这个区域中的坐标位置
     max_hourly_rain_indices  := ind_resolve(maxind(hourly_rain_1D),dsizes_hourly_rain)
     print(hourly_rain(loc(1,0)+max_hourly_rain_indices(0,0),loc(0,0)+max_hourly_rain_indices(0,1)))

     ; storage hourly rain series of a point 存放这个点的小时降水序列
     storm2_max_hourly_rain = new((/dimsizes(filename)-1/),float)
     do it = 0, dimsizes(filename)-1-1
   
 
     a              = addfiles( "../wrfout_d01_" + filename(it), "r")
     times          = wrf_user_getvar(a, "times", -1) 
     lat         = wrf_user_getvar(a, "XLAT",0)
     lon         = wrf_user_getvar(a, "XLONG",0)
     rainnc = wrf_user_getvar(a,"RAINNC",0)
     ntimes         = dimsizes(times) 
     printVarSummary(rainnc)

     a_next              = addfiles( "../wrfout_d01_" + filename(it+1), "r") 
     rainnc_next = wrf_user_getvar(a_next,"RAINNC",0)

     ; calculate hourly rainfall
     hourly_rain := rainnc_next-rainnc
     hourly_rain@lat2d = lat
     hourly_rain@lon2d = lon 

    

                                           ;regional y start + indice y            ;regional x start + indice x
 
     storm2_max_hourly_rain(it) = hourly_rain(loc(1,0)+max_hourly_rain_indices(0,0),loc(0,0)+max_hourly_rain_indices(0,1))

     end do
    
     
       type = "ps"     ;¿ªÊ¼Ö¸¶¨ÎÄ¼þÊä³öµÄÀàÐÍ
       wks = gsn_open_wks(type ,"./max_hourly_rainfall_xy") 
      
       res4                       = True             ; plot mods desired
     ;  res3@gsnDraw = False
       res4@gsnFrame = False
       res4@vpXF                 = 0.53
       res4@vpYF                 = 0.517
       res4@vpWidthF             = 0.41
       res4@vpHeightF            = 0.35
       res4@gsnLeftString   = "  hourly rainfall"
       res4@gsnLeftStringFontHeightF=14
       res4@gsnLeftStringOrthogonalPosF=0.002
       ;res4@tmXBLabelsOn=True
       ;res4@tmYLLabelsOn=True
       res4@tmXBMajorLengthF=0.003
       res4@tmYLMajorLengthF=0.006
        res4@tmXBMode          = "Explicit"         ; explicit labels
        res4@tmXBValues        = ispan(1,24,1)
        res4@tmXBLabels        = (/"17"," ","19"," ","21"," ","23"," ","01"," ","03"," ","05"," ","07"," ","09"," ","11"," ","13"," ","15"," "/) 
        ;res@tmXBLabelAngleF = 90
        ;res@tmXBLabelStride = 100
        res4@tmXTOn = False 
        res4@tmXTMode          = "Automatic"         ; explicit labels
        res4@tiXAxisString= "Time (LST)"
        res4@tmXBLabelFontHeightF = 0.01
     ;   res3@gsnYRefLine = 0
     ;   res3@gsnXYBarChart = False
      ;  res3@gsnAboveYRefLineColor = "blue"
        res4@tiYAxisString = "mm/h"
        res4@tiYAxisFontHeightF = 0.012
        res4@tmYLLabelDeltaF=-0.6
      ;  res4@tiYAxisConstantSpacingF
      ;  res4@tiYAxisSide="Right"
       ; res3@gsnXYBarChartBarWidth = 0.6
        res4@trYMinF = 0
        res4@trYMaxF = 180
        
       res4@xyDashPattern  = 0                                         ; Make curves all solid
       res4@xyMarkLineMode = "MarkLines"               ; Markers *and* lines
       res4@xyLineColors = (/"blue","blue","red","red"/) 
       res4@xyMarkers      = (/16,5,16,5/)                                ; 3 different markers
       res4@xyMarkerColors = (/"blue","blue","red","red"/)           ; 3 different colors 
       res4@xyMarkerSizeF=0.004
       res4@lgLegendOn = True
       
       res4@pmLegendDisplayMode    = "Always"            ; turn on legend
       res4@pmLegendSide           = "Bottom"               ; Change location of 
     ;  res3@pmLegendParallelPosF   = .5                  ; move units right
     ;  res3@pmLegendOrthogonalPosF = 0.05             ; move units down
       res4@pmLegendParallelPosF   = .8                 ; move units right
       res4@pmLegendOrthogonalPosF = -1.1           ; move units down
       res4@pmLegendWidthF         = 0.06                ; Change width and
       res4@pmLegendHeightF        = 0.032                ; height of legend.
      ; res3@pmLegend
       res4@lgPerimOn              = True               ; turn off box around
       res4@lgPerimColor         = "gray"
       res4@lgLabelFontHeightF     = .007                 ; label font height
       res4@xyExplicitLegendLabels = (/" G1059 "," storm1 "," G3322 "," storm2 "/)         ; create explicit labels
       
       series = new((/4,21/),float)

       series(0,:) = G1059_huadu_10_06UTC
       series(1,:) = storm1_max_hourly_rain
       series(2,:) = G3322_zengcheng_10_06UTC
       series(3,:) = storm2_max_hourly_rain
       
       plot_xy= gsn_csm_xy(wks, fspan(0.5,21.5,21), series(:,0:20),res4)  ;Creates and draws an XY plot.
       
       
       draw(plot_xy)
       
     ; txres3               = True                     ; text mods desired
     ; txres3@txFontHeightF = 0.015                    ; font smaller. default big
     ; dum4 = gsn_add_text(wks,plot_xy,"(c)",114.7,21.7,txres2)    
     
     ;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     frame(wks)


