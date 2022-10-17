; read ERA5，plot the classic skew-T with Temperature, Dewpoint, and wind barbs
; writte by yewei on Thur Nov 18,2021
;  ncl /mnt/d/ncl_related/ERA5_code/skewt_from_ERA5.ncl

; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/skewt_func.ncl"
load "/mnt/d/ncl_related/ERA5_code/skewt_func_dafis.ncl"

;--------------------main code-------------
begin

;1 read-----------------------------------------------------
  pic_dir="/mnt/d/ncl_related/pictures/skewt_tlogp/ERA5/"    
  data_dir= "/mnt/d/data/ERA5/"   
  data    = addfile(data_dir+"ERA5.nc", "r")  
  time    =data->time
  YYYYMMDDHH_f=cd_calendar(time,-3)

  need_time=2008102700
  time_locat= ind(YYYYMMDDHH_f.eq.need_time) ;找该时间所处位置
  path_time=tostring(need_time)

  u  = short2flt(data ->u(time_locat,{700:150},{29.75},{95.75})) ;bomi
  v  = short2flt(data ->v(time_locat,{700:150},{29.75},{95.75}))
  t= short2flt(data ->t(time_locat,{700:150},{29.75},{95.75}))
  rh =short2flt(data->r(time_locat,{700:150},{29.75},{95.75}))  

  ; u  = short2flt(data ->u(time_locat,{700:150},{28.0},{92.0})) ;cuona
  ; v  = short2flt(data ->v(time_locat,{700:150},{28.0},{92.0}))
  ; t= short2flt(data ->t(time_locat,{700:150},{28.0},{92.0}))
  ; rh =short2flt(data->r(time_locat,{700:150},{28.0},{92.0}))  

  ; u  = short2flt(data ->u(time_locat,{700:150},{28.75},{97.5})) ;chayu
  ; v  = short2flt(data ->v(time_locat,{700:150},{28.75},{97.5}))
  ; t= short2flt(data ->t(time_locat,{700:150},{28.75},{97.5}))
  ; rh =short2flt(data->r(time_locat,{700:150},{28.75},{97.5}))  

  ; u  = short2flt(data ->u(time_locat,{700:150},{29.5},{94.5})) ;linzhi
  ; v  = short2flt(data ->v(time_locat,{700:150},{29.5},{94.5}))
  ; t= short2flt(data ->t(time_locat,{700:150},{29.5},{94.5}))
  ; rh =short2flt(data->r(time_locat,{700:150},{29.5},{94.5}))  


  p  = data->level({700:150})  ;hPa

; 2 calculate meteo var
  td = dewtemp_trh(t,rh)
  copy_VarCoords(t, td)
  t=t-273.15
  td=td-273.15
  t@units="C"
  td@units="C"
  wspd = wind_speed(u,v)*2.5
  wdir = wind_direction(u,v,0)
  zall = stdatmus_p2tdz(p)
  z=zall(2,:)
  copy_VarCoords(p, z)


  wks  = gsn_open_wks ("ps", pic_dir+"ERA5_skewt_"+path_time)  ; send graphics to PNG file

; ; --- Create background skew-T and plot sounding----------------
  skewtOpts                 = True
  skewtOpts@tiMainString    = "Linzhi"+"(closest to Bomi)"+"  200810"+path_time
  skewtOpts@DrawMixRatio    = True 
  skewtOpts@vpWidthF  = 0.70
  skewtOpts@vpHeightF = 0.85
  ;skewtOpts@vpXF     = 0.01 
  skewtOpts@tiMainFontHeightF    = 0.018
  skewtOpts@tiMainPosition = "Right"
  skewtOpts@DrawFahrenheit  = False   ; default is True  
  skewtOpts@lineThicknessMixRatio     = 5

  dataOpts                  = True   ; options describing data and ploting
  dataOpts@Wthin            = 1      ; plot every n-th wind barb
  dataOpts@DrawWindBarbThk  = 5.0  
  dataOpts@PrintZ           = False    ; do not print Z
  ; dataOpts@ThermoInfo = False
  ; dataOpts@colWindP        = "blue"
  dataOpts@colTemperature  = "red"
  dataOpts@colDewPt        = "chartreuse3"
  dataOpts@colCape         = "chocolate3"
  dataOpts@colWindP        = "blue" 
  dataOpts@lineThicknessDewPt = 4.5
  dataOpts@lineThicknessTemperature = 4.5
  dataOpts@linePatternCape    = 17
  dataOpts@xpWind             = 45
  dataOpts@colTemperature  = "red"
  dataOpts@colDewPt        = "chartreuse3"
  dataOpts@colCape         = "chocolate3"
  dataOpts@colWindP        = "blue" 
  dataOpts@lineThicknessDewPt = 4.5
  dataOpts@lineThicknessTemperature = 4.5
  dataOpts@linePatternCape    = 17
  dataOpts@DrawWindBarbThk  = 5.0  ;无效了

    
  

  skewt_bkgd = skewT_BackGround (wks, skewtOpts)
  skewt_data = skewT_PlotData   (wks, skewt_bkgd,p,t,td,z \
                                    , wspd,wdir, dataOpts)
  draw (skewt_bkgd)
  draw (skewt_data)
  frame(wks)

end