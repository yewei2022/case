;和mpv_level相比，改了MPV1绘图属性
; These files are loaded by default in NCL V6.2.0 and newer
;修改气压层在L28、200/201,读数据和出图

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;--------------------------read data--------------------------------------------

begin

;读变量-----------------------------------------------------
  ;;;>>>手动设置一个时间 用各个数据生成一个带unit的 timeset_ymdh 并且转换为一个字符串 方便文件中读入
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
   file_format  = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time    = cd_string(timeset_ymdh,file_format)
   data_dir="/cygdrive/D/ncl_related/data/fnl/"   

   pic_dir ="/cygdrive/D/ncl_related/pictures/MPV/600_mpv1/" ;注意改层次

   data   = addfile(data_dir+"fnl_"+file_time+".nc", "r")  
    uData  = data ->UGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s lv_ISBL0=26层
    vData  = data->VGRD_P0_L100_GLL0(5:25,:,:) ;26
    tmpdata=data->TMP_P0_L100_GLL0 (5:25,:,: ) ;单位：K 26层
    avo    =data->ABSV_P0_L100_GLL0 (5:25,:,:)
    rh     =data->RH_P0_L100_GLL0 ( :, :, :) ;lv_ISBL4=21层 %
    pres1   =data->lv_ISBL4  ;Pa
    lon    =data->lon_0
    lat    =data->lat_0
    pres   =pres1/100
    pres@units="hPa"
    tmp    =tmpdata-273.15
    copy_VarCoords(rh, tmp)

;读rashmi路径
;风暴数据,用于标记风暴位置或路径
 
 fiTY = "/cygdrive/D/ncl_related/data/TC/bio042008.txt"
 ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
 nrow = numAsciiRow(fiTY) 
 YYYYMMDDHH = new(nrow, "string")
 latr = new(nrow, "float")
 lonr = new(nrow, "float")
 
 cdata = asciiread(fiTY, -1, "string")
 YYYYMMDDHH = str_get_field(cdata, 3, ",")
 latr = stringtofloat(str_get_field(cdata, 7, ",")) *0.1
 lonr = stringtofloat(str_get_field(cdata, 8, ",")) *0.1

 DateChar = stringtochar(YYYYMMDDHH)
 MM = chartostring(DateChar(:,5:6))
 DD = chartostring(DateChar(:,7:8))
 HH = chartostring(DateChar(:,9:10))



;printVarSummary(tmpdata)
;print(pres)
;exit 

;----------------------------计算相当位温-----------------------------------
 eqt= pot_temp_equiv(pres1,tmpdata, rh, 0, "rh")  ;可直接使用这个函数计算,用Pa来算
 copy_VarCoords(tmp, eqt)


;-------------------------------------------------------------------


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
 mpv =mpv1+mpv2
 copy_VarCoords(tmp, mpv)
 mpv =mpv*1.e+6
 copy_VarCoords(tmp, mpv1)
 mpv1 =mpv1*1.e+6


;----------------------------------------------------------------




wks = gsn_open_wks("png", pic_dir+file_time+"MPV1")

;--------------------属性设置-----------
res = True
 res@gsnDraw = False
 res@gsnFrame = False
 res@gsnLeftString = ""
 res@gsnRightString = ""
 res@tmYROn=False
 res@tmXTOn=False
 ;res@pmTickMarkDisplayMode   = "Always" ;不显示次刻度线
 ;res@tmXBMajorLengthF = 0.005  ;主刻度线长度,0.005太短了
 ;res@gsnMaximize=True
 ;res@vpWidthF=0.6
 ;res@vpHeightF=0.4

;设置坐标轴标签粗细
;res@gsnMaximize=True
res@tmXBLabelFont=22
res@tmYLLabelFont=22
res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
res@tmYLLabelFontHeightF =0.02



mpres=res
 mpres@mpFillOn                = False
 mpres@mpMinLatF               = 10.
 mpres@mpMaxLatF               = 40.
 mpres@mpMinLonF               = 70.
 mpres@mpMaxLonF               = 110.
 mpres@mpDataBaseVersion       = "Mediumres"
 mpres@mpDataSetName           = "Earth..4"
 mpres@mpOutlineSpecifiers     = (/"China:states"/)


mpvres=res
 mpvres@cnFillOn = True  ;cn means ContourPlot
 mpvres@cnLinesOn = True   
 mpvres@cnLineColor="black"
 mpvres@cnFillPalette = "BlueWhiteOrangeRed"
 mpvres@cnLineLabelsOn = True 
 mpvres@cnLineLabelFontColor="red"
 ;mpvres@cnLineLabelFormat="0@;*.1f" 
 mpvres@cnLineLabelInterval=0.5
 mpvres@cnInfoLabelOn = False ; set false will not set an information label
 ;mpvres@cnLevelSelectionMode = "ExplicitLevels" 
 ;mpvres@cnLevels=(/-3,-2.5,-2,-1.5,0,/)
 ;mpvres@cnFillColors=(/240,220,190,160,120,0/)
 mpvres@cnLevelSelectionMode = "ManualLevels" 
 mpvres@cnMinLevelValF = -3
 mpvres@cnMaxLevelValF = 3
 mpvres@cnLevelSpacingF = 0.5 
 ;mpvres@cnLineThicknessF = 4
 mpvres@gsnSpreadColorStart = 45
 mpvres@gsnSpreadColorEnd=220
 ;mpvres@gsnSpreadColors = True ; 6.1.0
 mpvres@lbLabelFontHeightF=0.015
 mpvres@lbBoxLinesOn = True ; Label bar
 ;mpvres@lbLabelAutoStride     = True            ; let NCL determine label spacing
 mpvres@lbOrientation        = "Vertical"         ; vertical label bar
 mpvres@pmLabelBarWidthF = 0.05     ;色标变窄
 mpvres@pmLabelBarHeightF=0.5
 mpvres@pmLabelBarOrthogonalPosF  = 0.025 ;色标离图的距离
 mpvres@lbBoxEndCapStyle = "TriangleBothEnds"

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


;----------------绘图----------------------
mpv= smth9_Wrap(mpv,0.5,0.25,False)
do i=0,3
 eqt= smth9_Wrap(eqt,0.5,0.25,False)
end do 
mapplot= gsn_csm_map(wks, mpres)
plot1  = gsn_csm_contour(wks,mpv1({60000},:,:),mpvres)
plot2  = gsn_csm_contour(wks,eqt({60000},:,:),eqtres)
overlay(mapplot, plot1)
overlay(mapplot, plot2)

;某时刻风暴位置标记
 txres               = True                     ; text mods desired
 txres@txFontHeightF = 0.03                     ; font smaller. default big
 txres@txFontColor   ="red"
 do i = 0, nrow-1
   if (HH(i) .eq. "00".and.DD(i) .eq."27") then
      dum = gsn_add_text(wks,mapplot,"~F37~p",lonr(i), latr(i),txres) 
    end if
 end do
draw(mapplot)
frame(wks)
end