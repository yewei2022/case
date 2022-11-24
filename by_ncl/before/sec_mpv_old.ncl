;计算相对位温，垂直涡度，绝对动量，MPV,画区域平均沿92E的剖面图
;old与new的区别在于：
;old计算相当位温是自己写公式计算的，new则用函数
;old里面mpv*e5，是错的
; These files are loaded by default in NCL V6.2.0 and newer
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"


;;-------------------绘图属性----------------------------
undef("set_all_res")
procedure set_all_res(res,res2,MPVres,eqtres,vores,Mres,txres)
 local res,res2,MPVres,eqtres,vores,Mres,txres
 begin
 res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@tmXTOn = False ;关刻度
    res@tmYROn = False
    ;res@gsnLeftString = "hPa" ;左上角字符
    res@tiYAxisString="P/hPa"
    res@tiYAxisFont=22
    ;res@tiYAxisFontHeightF=0.5

    ;设置坐标轴标签粗细
    res@gsnMaximize=True
    res@tmXBLabelFont=22
    res@tmYLLabelFont=22
    res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
    res@tmYLLabelFontHeightF =0.02


 ;-------------------湿位涡剖面属性----------------------------
 MPVres=res
    MPVres@cnFillOn = True  ;cn means ContourPlot
    MPVres@cnLinesOn = True   
    ;res@cnMonoLineColor = False
    MPVres@cnLineColor="black"
    ;res@cnSpanLinePalette = True
    MPVres@cnFillPalette = "WhiteBlue"
    ;res@cnLineLabelsOn = False 
    ;res@cnLineLabelFormat="0@;*.1f" 
    ;res@cnLineLabelInterval=0.05
    MPVres@cnInfoLabelOn = False ; set false will not set an information label
    MPVres@cnLevelSelectionMode = "ExplicitLevels" 
    MPVres@cnLevels=(/-0.3,-0.2,-0.1,-0.05,0/)
    MPVres@cnFillColors=(/240,200,160,120,80,0/)
    ;res@cnLevelSelectionMode = "ManualLevels" 
    ;res@cnMinLevelValF = 552
    ;res@cnMaxLevelValF = 588
    ;res@cnLevelSpacingF = 4 
    ;res@cnLineThicknessF = 4
    ;res@gsnSpreadColorStart = 1
    ;res@gsnSpreadColorStart =17
    ;res@gsnSpreadColors = True ; 6.1.0
    MPVres@lbBoxLinesOn = True ; Label bar
    MPVres@lbLabelAutoStride     = True            ; let NCL determine label spacing
    ;res@lbOrientation        = "Vertical"         ; vertical label bar
    ;res@pmLabelBarWidthF = 0.05     ;色标变窄
    MPVres@pmLabelBarHeightF =0.05
    MPVres@lbBoxEndCapStyle = "TriangleBothEnds"
    MPVres@gsnContourZeroLineThicknessF = 2.            ; doubles thickness of zero contour
    MPVres@gsnContourNegLineDashPattern = 1         ; sets negative contours to dash pattern 1
    MPVres@tmYRMode             = "Automatic"          ; turn off special labels on right axis  
    MPVres@cnFillDrawOrder      = "Predraw"  ;调整画图的顺序 

 
 eqtres=res
    eqtres@cnFillOn = False  ;cn means ContourPlot
    eqtres@cnLinesOn = True 
    ;res@cnMonoLineColor = False
    eqtres@cnLineColor="black"
    eqtres@cnLabelMasking = True
    eqtres@cnLineLabelPlacementMode = "Constant"
    eqtres@cnLineDashSegLenF    = 0.3
    eqtres@cnLineLabelBackgroundColor = -1
    eqtres@cnInfoLabelOn = False ; set false will not set an information label
    eqtres@cnLineThicknessF = 6
    eqtres@cnLevelSelectionMode = "ManualLevels" 
    eqtres@cnLevelSpacingF = 2   
    eqtres@tmYRMode             = "Automatic"          ; turn off special labels on right axis   

 ;地形填充属性
 res2=True
    ;res2@gsnBelowYRefLineColor = "black" 
    ;res2@gsnXYBelowFillColors="black"
    res2@gsnMaximize      = True
    res2@trYMaxF=1000
    res2@gsnYRefLine      = 1000.          ; create a reference line 
    res2@gsnBelowYRefLineColor = "black"              ; above ref line fill red
    res2@gsnDraw        = False             ; don't draw the plots yet
    res2@gsnFrame       = False        ; don't advance the frame

 vores=res
    vores@cnFillOn = True  ;cn means ContourPlot
    vores@cnLinesOn = False  
    vores@cnInfoLabelOn = False ; set false will not set an information label 
    vores@cnFillPalette = "BlueWhiteOrangeRed"
    ;vores@cnLevelSelectionMode = "ExplicitLevels" 
    ;vores@cnLevels=(/-0.3,-0.2,-0.1,-0.05,0/)
    ;vores@cnFillColors=(/240,200,160,120,80,0/)
    vores@cnLevelSelectionMode = "ManualLevels" 
    ;res@cnLevelSpacingF = 4 
    ;res@cnLineThicknessF = 4
    vores@lbBoxLinesOn = True ; Label bar
    vores@lbLabelAutoStride     = True            ; let NCL determine label spacing
    ;res@lbOrientation        = "Vertical"         ; vertical label bar
    ;res@pmLabelBarWidthF = 0.05     ;色标变窄
    vores@pmLabelBarHeightF =0.05
    vores@lbBoxEndCapStyle = "TriangleBothEnds"
    vores@tmYRMode             = "Automatic"          ; turn off special labels on right axis 
    vores@cnFillDrawOrder      = "Predraw"  ;调整画图的顺序

 Mres=res
    Mres@cnFillOn = False  ;cn means ContourPlot
    Mres@cnLinesOn = True   
    ;res@cnMonoLineColor = False
    Mres@cnLineColor="red"
    Mres@cnLineDashPattern= 14
    Mres@cnLineLabelsOn=True
    Mres@cnLineLabelFontColor="red"
    Mres@cnLabelMasking = True
    Mres@cnLineLabelPlacementMode = "Constant"
    Mres@cnLineDashSegLenF    = 0.3
    Mres@cnLineLabelBackgroundColor = -1
    Mres@cnInfoLabelOn = False ; set false will not set an information label
    Mres@cnLineThicknessF = 3
    Mres@cnLevelSpacingF=10
    Mres@cnLevelSelectionMode = "ManualLevels"   
    Mres@tmYRMode             = "Automatic"          ; turn off special labels on right axis
    Mres@cnLineDrawOrder      = "Predraw"  ;调整画图的顺序

 txres               = True                     ; text mods desired
    txres@txFontHeightF = 0.03                     ; font smaller. default big
    txres@txFontColor   ="red"



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

;-------------main code----------------------------------

begin
;----------------------------读数据----------------------------------------------------
 pic_dir ="/cygdrive/D/ncl_related/pictures/MPV/"

;读地形文件----------------------------
   ;lon1       fspan(0,360,4320)  ;gsn_xy    
   lat1     =fspan(90,-90,2160)  ;用于gsn_xy 地形图
   data1_dir  ="/cygdrive/D/ncl_related/data/other/"       ;数据、图片所在路径
   setfileoption("nc","MissingToFillValue",False)    ;设置缺省值
   elev= read_elev_data(data1_dir+"ETOPO5.DAT")
   ;printVarSummary(elev)
   elev0 = 1013.25*(1-elev*0.0065/288.15)^5.25145        ;将地形数据的单位：m转化为hPa
   copy_VarCoords(elev,elev0) 
   ;printVarSummary(elev0)
   stdlon=92;选择固定某一角度
    elev1=elev0(:,{stdlon})       ;固定某一经度的地形数据
    elev1@_FillValue= -9.96921e+36

;读变量-----------------------------------------------------
  ;;;>>>手动设置一个时间 用各个数据生成一个带unit的 timeset_ymdh 并且转换为一个字符串 方便文件中读入
   year   = 2008
   month  = 10
   day    = 26
   hour   = 12
   minute = 00
   second = 0
  ;设置时间格式读取数据####################################
   units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
   opt = 0
   timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
   file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)
   data2_dir= "/cygdrive/D/ncl_related/data/fnl/"  ; 物理量数据的路径  
   data   = addfile(data2_dir+"fnl_"+file_time+".nc", "r")                        

    uData  = data ->UGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s lv_ISBL0=26层
    vData  = data->VGRD_P0_L100_GLL0(5:25,:,:) ;26
    tmpdata=data->TMP_P0_L100_GLL0 (5:25,:,: ) ;单位：K 26层
    avo    =data->ABSV_P0_L100_GLL0 (5:25,:,:);26层
    rh     =data->RH_P0_L100_GLL0 ( :, :, :) ;lv_ISBL4=21层 %
    pres1   =data->lv_ISBL4  ;Pa
    lon    =data->lon_0
    lat    =data->lat_0
    pres   =pres1/100
    pres@units="hPa"

;printVarSummary(uData)
;print(pres)
;exit 
;-----------------------------------------------------------



wks = gsn_open_wks("png",pic_dir+"sec_eqt_mpv_mpv1"+file_time) 
;gsn_define_colormap(wks,"BlueWhiteOrangeRed") 



;----------------------------计算相当位温-----------------------------------
 ;eqt= pot_temp_equiv(pres1,tmpdata, rh, 0, "rh")  ;可直接使用这个函数计算

  tmp    =tmpdata-273.15
  copy_VarCoords(rh, tmp)
  es=new(dimsizes(tmp),typeof(tmp))
  qs=new(dimsizes(tmp),typeof(tmp))
  q=new(dimsizes(tmp),typeof(tmp))
  e=new(dimsizes(tmp),typeof(tmp))
  tlcl=new(dimsizes(tmp),typeof(tmp))
  theta=new(dimsizes(tmp),typeof(tmp))
  eqt=new(dimsizes(tmp),typeof(tmp))

  es    =(6.112*exp((17.67*tmp)/(tmp+243.5))) ;饱和水汽压
  qs    =0.62197*es/((conform(es,pres,0))-0.378*es) ;饱和比湿 
  ;由于pres一维，qs三维，为进行正常计算，用conform
  ;https://www.ncl.ucar.edu/Document/Functions/Built-in/conform.shtml
  q     =rh*qs/100 ;比湿 
  e     =(conform(es,pres,0))*q/(0.62197+q)+1e-10 ;水汽压
  tlcl  =55.0+2840.0/(3.5*log(tmp+273.16)-log(e)-4.805) ;抬升凝结温度
  theta =(tmp+273.16)*((1000/(conform(es,pres,0)))^(0.2854*(1.0-0.28*q)));位温
  eqt   =theta*exp(((3376./tlcl)-2.54)*q*(1.0+0.81*q)) ;相当位温

  copy_VarCoords(tmp, eqt)
  ;h_leqt  =eqt(:,:,{stdlon})
  h_leqt = dim_avg_n_Wrap(eqt(:,:,{92:98}), 2)

  h_leqt!0="level"
  h_leqt&level=pres


  
;-----------------------------计算湿位涡----------------------------    
 g      =9.8
 dTdp   = center_finite_diff_n(eqt,pres1,False,0,0)
 dudp   = center_finite_diff_n(uData,pres1,False,0,0)
 dvdp   = center_finite_diff_n(vData,pres1,False,0,0)

 nlat=dimsizes(lat)
 mlon=dimsizes(lon)

 dlon   = (lon(2)-lon(1))*0.0174533 ; convert to radians 0.017=3.14/180
 dTdX   = new ( dimsizes(eqt),typeof(eqt)); pre-allocate space
 do nl=0,nlat-1                      ; loop over each latitude
  dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
  dTdX(:,nl:nl,:) = center_finite_diff_n (eqt(:,nl:nl,:), dX , True,0,2)
 end do

 dY   =6378388.*0.0174533
 dTdY = new ( dimsizes(eqt),typeof(eqt)); pre-allocate space
 do nl=0,mlon-1                      ; loop over each latitude
  dTdY(:,:,nl:nl) = center_finite_diff_n (eqt(:,:,nl:nl), dY , False,0,1)
 end do

 mpv1=-g*avo*dTdp 
 mpv2=g*(dvdp*dTdX-dudp*dTdY)
 mpv =mpv2
 copy_VarCoords(tmp, mpv)

 h_lmpv=dim_avg_n_Wrap(mpv(:,:,{92:97}), 2)*1.e+5

 ;h_lmpv      =mpv(:,:,{stdlon})*1.e+5
 copy_VarCoords(tmp, h_lmpv)
 h_lmpv!0    ="level"
 h_lmpv&level=pres


;----------------------------计算垂直涡度------------------
 vr = uv2vr_cfd (uData,vData,lat,lon, 2) 
  copy_VarCoords(rh, vr)
  ;h_lvr  =vr(:,:,{stdlon})*1.e+5
  h_lvr = 1.e+5*dim_avg_n_Wrap(vr(:,:,{92:97}), 2)
  copy_VarCoords(rh, h_lvr)
  h_lvr!0="level"
  h_lvr&level=pres
  h_lvr!1="lat"
  h_lvr&lat=lat
  ;printVarSummary(h_lvr)
  ;exit

;----------------------------计算绝对动量-----------------
 M= new ( dimsizes(uData),typeof(uData)); pre-allocate space
  copy_VarCoords(uData, M)
  f= new(dimsizes(lat), float)
  copy_VarCoords(lat, f)
  y= new(dimsizes(lat), float)
  copy_VarCoords(lat, y)
  do nl=0,nlat-1                      ; loop over each latitude
    f(nl)=2*7.292*sin(lat(nl)*3.14159/180.0)*0.00001  ;f=2wsin0
    y(nl)=dY*lat(nl) 
  end do
  M=uData-(conform(uData,f,1))*(conform(uData,y,1))
  ;h_lM  =M(:,:,{stdlon})
  h_lM = dim_avg_n_Wrap(M(:,:,{92:97}), 2)

  copy_VarCoords(uData, h_lM)
  h_lM!0="level"
  h_lM&level=pres

;------------绘图---------------------

res    = True
 res2   = True
 MPVres = True
 eqtres = True
 vores  = True
 Mres   = True
 txres  = True
 set_all_res(res,res2,MPVres,eqtres,vores,Mres,txres)

plot_eqt   = gsn_csm_pres_hgt(wks, h_leqt({100:700},{25:35}), eqtres ) 
plot_vr= gsn_csm_pres_hgt(wks, h_lvr({100:700},{25:35}), vores ) 
;plot_MPV  = gsn_csm_pres_hgt(wks, h_lmpv({100:700},{25:35}), MPVres ) 
plot_M = gsn_csm_pres_hgt(wks, h_lM({100:700},{25:35}), Mres ) 
;plot_hgt = gsn_csm_xy(wks,lat1,elev1,res2) 
overlay(plot_vr,plot_eqt)
overlay(plot_vr,plot_M)
;overlay(plot_vr,plot_hgt)


;dum = gsn_add_text(wks,plot_vr,"~F37~p",lat({28}),pres1({100000}),txres) 


draw(plot_vr)
frame(wks)
end