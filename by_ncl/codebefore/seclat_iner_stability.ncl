 mask_start_time = get_cpu_time() 

;27日00时过台风中心纬向剖面：垂直速度矢量，惯性稳定度等值线，相当位温等值线
;换时次时，需改动等值线间隔66，时间180，风暴中心经纬度168,以及最后标记的位置
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;---------------------------------------绘图属性设置------------------------------------------
undef("set_all_res")
procedure set_all_res(res,etores,cnres,eqtres)
 local res,etores,cnres,eqtres
 begin
 
  ;--------------------------------------------------
   res = True
   res@gsnDraw       = False
   res@gsnFrame      = False
   res@gsnMaximize   = True
   res@tmXTOn        = False ;关刻度
   res@tmYROn        = False
   res@gsnLeftString = "" ;左上角字符
   res@gsnRightString=""
   res@tiYAxisString ="P/hPa"
   res@tiYAxisFont   =22
   res@tmYRMode      = "Automatic" ;关掉右边的height（km）字样
   ;res@gsnMajorLonSpacing    =5. 
   ;res@gsnMinorLonSpacing    =1.    
   ;res@gsnMajorLatSpacing    =5. 
   ;res@gsnMinorLatSpacing    =1



  ;设置坐标轴
    res@tmXBMode       ="Explicit"
    res@tmXBValues     =(/85,90,95/)
    res@tmXBLabels     =(/"85E","90E","95E"/)
    res@tmXBMinorValues=(/86,87,88,89,91,92,93,94/)
    res@tmXBLabelFont=22
    res@tmYLLabelFont=22
    res@tmXBLabelFontHeightF =0.025 ; x轴标签字号
    res@tmYLLabelFontHeightF =0.025



  ;--------------------惯性稳定度等值线---------------------------
   cnres=res
   cnres@cnLinesOn         =True
   cnres@cnFillOn          =False                  ; turn on color
   cnres@cnLineColor       ="blue"
   cnres@cnLineLabelsOn    =True    
   cnres@cnLineLabelFontColor  ="blue"
   cnres@cnLineLabelPerimOn=True
   cnres@cnLineLabelFormat="0@;*.2f"
   cnres@cnInfoLabelOn         = False ; set false will not set an information label
   ;cnres@cnFillPalette        = "NEO_div_vegetation_c"            ; set color map
   cnres@cnLevelSelectionMode  = "ManualLevels" 
   ;cnres@cnLabelMasking = True
   ;cnres@cnLineLabelPlacementMode = "Constant" ;使标签是线的一部分
   ;cnres@cnLineDashSegLenF    = 0.3
   ;cnres@cnLineLabelBackgroundColor = -1 ;设置标签背景透明
   cnres@cnLineThicknessF = 6
   cnres@gsnContourNegLineDashPattern=15
   ;cnres@cnMinLevelValF       = -8             ; set min contour level
   ;cnres@cnMaxLevelValF       = 8               ; set max contour level
   cnres@cnLevelSpacingF       = 0.2              ; set contour spacing,2600-0.1,2700-0.2
   ;cnres@lbBoxLinesOn         = True ; Label bar
   ;cnres@lbLabelAutoStride    = True            ; let NCL determine label spacing
   ;cnres@lbOrientation        = "vertical"            ; 垂直坐标
   ;cnres@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
   ;cnres@lbBoxEndCapStyle     = "TriangleBothEnds"
   ;cnres@pmLabelBarWidthF     = 0.05     ;色标变窄
   ;cnres@pmLabelBarHeightF    =0.05
   ;cnres@cnFillDrawOrder       = "PreDraw"

  ;---------vector-----------------
    ;垂直速度矢量叠加在contour上，属性设置合并
    cnres@vcMapDirection        = False ;绘制垂直剖面图要设置为False
    cnres@vcGlyphStyle          = "LineArrow"
    cnres@vcLineArrowThicknessF = 4
    cnres@vcMinDistanceF        = 0.025
    cnres@vcRefMagnitudeF       = 15 ;设置单位量级
    cnres@vcRefLengthF          = 0.02
    cnres@vcLineArrowHeadMaxSizeF = 0.01
    cnres@vcLineArrowHeadMinSizeF = 0.01
    cnres@vcRefAnnoOn           = True ;anno代表矢量注释
    cnres@vcRefAnnoFontHeightF  = 0.015
    cnres@vcRefAnnoZone         = 1
    cnres@vcRefAnnoSide         = "Top"
    cnres@vcRefAnnoJust         = "TopRight"
    ;cnres@vcRefAnnoString1 = "5 g cm~S~-1~N~hPa~S~-1~Ns~S~-1~N)"
    cnres@vcRefAnnoString1      = "15"
    cnres@vcRefAnnoString2On    = False


 ;相当位温 
 eqtres=res
    eqtres@cnFillOn    = False  ;cn means ContourPlot
    eqtres@cnLinesOn   = True 
    ;res@cnMonoLineColor = False
    eqtres@cnLineColor ="black"
    eqtres@cnLabelMasking = True
    eqtres@cnLineLabelPlacementMode = "Constant" ;使标签是线的一部分,标签沿线等距分布；默认Randomized
    eqtres@cnLineDashSegLenF    = 0.3  ;标签之间的距离
    eqtres@cnLineLabelBackgroundColor = -1 ;设置标签背景透明
    eqtres@cnInfoLabelOn    = False ; set false will not set an information label
    eqtres@cnLineThicknessF = 4
    eqtres@cnLevelSelectionMode = "ManualLevels" 
    eqtres@cnLevelSpacingF  = 2   


  ;-----------地形---------------------------
   etores=True
    etores@gsnMaximize      = True
    etores@gsnDraw          = False           ; don't draw the plots yet
    etores@gsnFrame         = False           ; don't advance the frame
    etores@trYMaxF          =1000
    etores@gsnYRefLine      =1000.          ; create a reference line 
    etores@gsnBelowYRefLineColor = "black"  ; above ref line fill red
end


undef("read_elev_data")
function read_elev_data(topo_file)
local nlat, nlon, topo_file, lat, lon
begin
  ;---Read data as a straight binary file
  nlat = 2160
  nlon = 4320
  setfileoption("bin","ReadByteOrder","BigEndian")
  elev= cbinread(topo_file,(/nlat,nlon/),"short")

  ;---Create 1D coordinate arrays
  lat       = fspan(90,-90,nlat)
  lon       = fspan(0,360,nlon)
  lat!0     = "lat"
  lon!0     = "lon"
  lat@units = "degrees_north"
  lon@units = "degrees_east"
  lat&lat   = lat
  lon&lon   = lon

  ;---Attach the coordinate arrays
  elev!0    = "lat"
  elev!1    = "lon"
  elev&lat  = lat
  elev&lon  = lon
 
  return(elev)
end


begin

;----------------------------------------------读地形文件----------------------------
   ;lon1       fspan(0,360,4320)  ;gsn_xy    
   lat1     =fspan(90,-90,2160)  ;用于gsn_xy 地形图
   ;datadir="/cygdrive/D/ncl_related/data/"
   datadir  ="/cygdrive/D/ncl_related/data/other/"   ;地形数据路径
   setfileoption("nc","MissingToFillValue",False)    ;设置缺省值
   elev= read_elev_data(datadir+"ETOPO5.DAT")
   ;printVarSummary(elev)
   elev0 = 1013.25*(1-elev*0.0065/288.15)^5.25145        ;将地形数据的单位：m转化为hPa
   copy_VarCoords(elev,elev0) 
   ;printVarSummary(elev0)

;----------------------------------------------选择经纬度-------------------------------
  ; stdlon=88.1  ;2600
  ; stdlat=18.2
  stdlon=90.4    ;2700 
  stdlat=23.4   
  ;2008102600,182N,881E;2008102606,187N,884E;
  ;2008102612,198N,888E;2008102618,210N,896E;
  ;2008102700,234N,904E;2008102706,250N,906E,
    elev1=elev0(:,{stdlon})       ;固定某一经度的地形数据    
    elev1@_FillValue= -9.96921e+36
  strlon=85 ;绘图横坐标范围
  endlon=95

;---------------------------------------------读变量-----------------------------------------------------
  ;手动设置一个时间 用各个数据生成一个带unit的 timeset_ymdh 并且转换为一个字符串 方便文件中读入
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


   path_rawdata= "/cygdrive/D/ncl_related/data/fnl/"  ; 变量数据路径
   data       = addfile(path_rawdata+"fnl_"+file_time+".nc", "r") 
   rh_in      = data->RH_P0_L100_GLL0   ;lv_ISBL4/lat_0/lon_0,21层
   tmp_in = data->TMP_P0_L100_GLL0(5:25,:,: ) ;单位：K,lv_ISBL0/lat_0/lon_0,26层，选择lv_ISBL4对应的21层
   u_in       = data->UGRD_P0_L100_GLL0(5:25,:,:) 
   v_in       = data->VGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
   vvel_in    = data->VVEL_P0_L100_GLL0(:,:,: ) ;Pa/s,lv_ISBL4, lat_0, lon_0 
   pres1      =data->lv_ISBL4  ;Pa
   lon_in     =data->lon_0
   lat_in     =data->lat_0
   pres       =pres1/100
   pres@units ="hPa"


;---------------------------网格点插值--------------------
lon_out=fspan(70,110,401)
lat_out=fspan(0, 50, 501)
   u      =linint2_Wrap(lon_in, lat_in(::-1), u_in(:,::-1,:), True, lon_out, lat_out, 0)
   delete(u_in)
   v      =linint2_Wrap(lon_in, lat_in(::-1), v_in(:,::-1,:), True, lon_out, lat_out, 0)
   delete(v_in)
   vvel   =linint2_Wrap(lon_in, lat_in(::-1), vvel_in(:,::-1,:), True, lon_out, lat_out, 0)
   delete(vvel_in)
   rh     =linint2_Wrap(lon_in, lat_in(::-1), rh_in(:,::-1,:), True, lon_out, lat_out, 0)
   delete(rh_in)
   tmp=linint2_Wrap(lon_in, lat_in(::-1), tmp_in(:,::-1,:), True, lon_out, lat_out, 0)
   delete(tmp_in)

;-------------------------计算各层散度----------------------------------------------------------
  div=uv2dvF_Wrap(u(:,:,:),v(:,:,:)) ;用这个球面间谐波函数计算散度更精确
   ;di = uv2dv_cfd (u,v,lat,lon, 3) 
   div=div*1.e5
   h_ldiv=div(:,{stdlat},:)

;----------------------------------设定垂直速度矢量-------------------------------------------------
  ;设定纬向风速，将Pa/s转换成m/s   
   w1 = omega_to_w(vvel, conform(vvel,pres1, 0), tmp)  ; w[*][*][*] (m/s)，要求温度K，气压Pa
   w =w1*100.
   copy_VarCoords(w1, w)
   h_lw = w(:,{stdlat},:)
   h_lw!0="level" ;有这步没这步画出来的图都不一样
   h_lw&level=pres
  ;设定经向风速v经度
   h_lu = v(:,{stdlat},:)
   h_lu!0="level"
   h_lu&level=pres

 ;----------------------------计算相当位温---------------------------------------------
 eqt         = pot_temp_equiv(pres1,tmp, rh, 0, "rh")  ;可直接使用这个函数计算,用Pa来算
 copy_VarCoords(tmp, eqt)
 h_leqt      =eqt(:,{stdlat},:)
 ;h_leqt     =dim_avg_n_Wrap(eqt(:,:,{92:98}), 2)
 h_leqt!0    ="level"
 h_leqt&level=pres

;----------------------------计算惯性稳定度---------------------------------------------
 ;1.设置三维的距离变量,地转参数f
  r= new ( dimsizes(u),typeof(u)); pre-allocate space
  copy_VarCoords(u, r)
  f= new(dimsizes(u), typeof(u))
  copy_VarCoords(u, f)
  x= new ( dimsizes(u),typeof(u))
  copy_VarCoords(u, x) 
  y= new(dimsizes(u), typeof(u))
  copy_VarCoords(u, y)
 ;2.给f、x，y，r 变量赋值
  nlat=dimsizes(lat_out)
  mlon=dimsizes(lon_out)
  dY   =6378388.*0.0174533 ;一个纬度的单位弧长
  dX= new(dimsizes(u), float)
  copy_VarCoords(u, dX)
  ;dlon = (lon_out(2)-lon_out(1))*0.0174533 ; convert to radians 0.017=3.14/180
  ;不对，上面这个式子代表0.1个经度的弧长
  ;应该改成1个经度的
  dlon=0.0174533 ;R*1°/180°
  do nl=0,nlat-1                      ; loop over each latitude
    f(:,nl,:) = 2*7.292*sin(lat_out(nl)*3.14159/180.0)*0.00001  ;f=2wsin0
    y(:,nl,:) = dY*lat_out(nl)
    dX(:,nl,:)= 6378388.*cos(0.0174533*lat_out(nl))*dlon  ; constant at this latitude，R*COS()*单位经度弧长
  end do
  do nl=0,mlon-1                      ; loop over each latitude
    x(:,:,nl) =dX(:,:,nl)*lon_out(nl)
  end do
 ;3.计算每个格点相对风暴中心的x,y,r
  X=x
  Y=y
  X0= x(:,{stdlat},{stdlon})
  Y0= y(:,{stdlat},{stdlon})

  X=x-conform(X, X0, 0) 
  Y=y-conform(Y, Y0, 0)

  r=sqrt(X^2+Y^2)
  r@units="m"

  
 ;4计算垂直涡度
  vr = uv2vr_cfd (u,v,lat_out,lon_out, 2) 
  copy_VarCoords(u, vr)

 ;5.计算切向风
  rinv = 1. / where(r.ne.0,r,r@_FillValue)
  Vt=v*X*rinv-u*Y*rinv
  copy_VarCoords(u, Vt)

 ;6. 计算惯性稳定度
  iner_stabl= (vr+f)*(2.*Vt*rinv+f)
  copy_VarCoords(u, iner_stabl)
  iner_stabl=iner_stabl*1.e7
  h_l_iner_stabl  =iner_stabl(:,{stdlat},:)
  h_l_iner_stabl!0="level"
  h_l_iner_stabl&level=pres



pic_dir ="/cygdrive/D/ncl_related/pictures/storm/sec_iner_stabl/lat/"
wks = gsn_open_wks("png",pic_dir+stdlat+"N"+"_iner_stabl_eqt_"+file_time) 


;------------绘图---------------------

res     = True
 etores = True
 cnres  = True
 eqtres = True
 set_all_res(res,etores,cnres,eqtres)

;print(min(h_lvvel))
;print(max(h_lrh))
;exit


plot_vector= gsn_csm_pres_hgt_vector(wks,h_l_iner_stabl({150:1000},{strlon:endlon}),\
                                     h_lu({150:1000},{strlon:endlon}),h_lw({150:1000},{strlon:endlon}),cnres)
plot_eqt  = gsn_csm_pres_hgt(wks, h_leqt({150:1000},{strlon:endlon}), eqtres ) 
plot_hgt = gsn_csm_xy(wks,lat1,elev1,etores) 
overlay(plot_vector,plot_eqt)
overlay(plot_vector,plot_hgt)
draw(plot_vector)

;某时刻风暴位置标记
 txres               = True                     ; text mods desired
 txres@txFontHeightF = 0.03                     ; font smaller. default big
 txres@txFontThicknessF=6
 txres@txFontColor   ="red"
 ;dum = gsn_add_text(wks,plot_vector,"~F37~p",1000,stdlat,txres) 
;page coordinates are normalized. go from 0->1 not 0->8" etc. 第一个为左右，左小右大，第二个为上下，下小上大
gsn_text_ndc(wks,"~F37~p",0.61,0.18,txres)  
;2600-0.435;2700-0.61
frame(wks)
end

mask_end_time = get_cpu_time()
print(mask_end_time)