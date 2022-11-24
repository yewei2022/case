;0305 绘制96E垂直速度等值线，散度填色垂直剖面
;计算等压面上的比湿 mixhum_ptrh，
;网址 http://www.ncl.ucar.edu/Document/Functions/Built-in/mixhum_ptrh.shtml
;将Pa/s转换成m/s，用于与v合成
;网址 https://www.ncl.ucar.edu/Document/Functions/Contributed/omega_to_w.shtml
;根据气压计算高度http://www.ncl.ucar.edu/Document/Functions/Built-in/stdatmus_p2tdz.shtml
;气压纬度矢量图gsn_csm_pres_hgt_vector 
;http://www.ncl.ucar.edu/Document/Graphics/Interfaces/gsn_csm_pres_hgt_vector.shtml

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;---------------------------------------绘图属性设置------------------------------------------
undef("set_all_res")
procedure set_all_res(res,cnres,vvelres,etores)
 local res,cnres,vvelres,etores
 begin
 
  ;--------------------------------------------------
   res = True
   res@gsnDraw = False
   res@gsnFrame = False
   res@tmXTOn = False ;关刻度
   res@tmYROn = False
   res@tmYRMode = "Automatic"  ; turn off special labels on right axis   
   res@gsnLeftString = "" ;左上角字符
   res@gsnRightString=""
   res@tiYAxisString="P/hPa"
   res@tiYAxisFont=22

  ;设置坐标轴标签粗细
    ;res@gsnMaximize=True
    res@tmXBLabelFont=22
    res@tmYLLabelFont=22
    res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
    res@tmYLLabelFontHeightF =0.02

  ;--------------------散度---------------------------
   cnres=res
   cnres@cnFillOn          = True                  ; turn on color
   cnres@cnInfoLabelOn = False ; set false will not set an information label
   cnres@cnFillPalette     = "NEO_div_vegetation_c"            ; set color map
   cnres@cnLinesOn         = False                 ; no contour lines
   cnres@cnLevelSelectionMode = "ManualLevels" 
   cnres@cnMinLevelValF       = -8             ; set min contour level
   cnres@cnMaxLevelValF       = 8               ; set max contour level
   cnres@cnLevelSpacingF      = 2              ; set contour spacing
   cnres@lbBoxLinesOn      = True ; Label bar
   cnres@lbLabelAutoStride     = True            ; let NCL determine label spacing
   ;cnres@lbOrientation        = "vertical"            ; 垂直坐标
   cnres@pmLabelBarOrthogonalPosF  = 0.08
   cnres@lbBoxEndCapStyle  = "TriangleBothEnds"
   ;cnres@pmLabelBarWidthF = 0.05     ;色标变窄
   cnres@pmLabelBarHeightF =0.05
   cnres@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
   cnres@cnFillDrawOrder   = "PreDraw"


  ;---------垂直速度等值线-----------------
   vvelres=res
   vvelres@cnFillOn = False  ;cn means ContourPlot
   vvelres@cnLinesOn = True   
   vvelres@cnLineColor="black"
   vvelres@cnLineThicknessF =7  
   vvelres@cnLineLabelsOn = True 
   vvelres@cnLineLabelFormat="0@;*.2f" 
   vvelres@cnLineLabelFontHeightF=0.015
   vvelres@cnLineLabelPlacementMode = "Constant"
   vvelres@cnLineDashSegLenF    = 0.3
   vvelres@cnLineLabelBackgroundColor = -1  
   vvelres@cnLabelMasking = True
   vvelres@cnLevelSpacingF=0.2
   vvelres@cnInfoLabelOn = False ; set false will not set an information label
   vvelres@gsnContourZeroLineThicknessF = 8            ; doubles thickness of zero contour
   vvelres@gsnContourNegLineDashPattern = 1         ; sets negative contours to dash pattern 1
   vvelres@tmYRMode             = "Automatic"          ; turn off special labels on right axis 


  ;-----------地形---------------------------
   etores=True
    etores@gsnMaximize      = True
    etores@trYMaxF=1000
    etores@gsnYRefLine      = 1000.          ; create a reference line 
    etores@gsnBelowYRefLineColor = "black"              ; above ref line fill red
    etores@gsnDraw        = False             ; don't draw the plots yet
    etores@gsnFrame       = False        ; don't advance the frame

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
 
 data1_dir  ="/cygdrive/D/ncl_related/data/fnl/" ;fnl数据路径
 data2_dir  ="/cygdrive/D/ncl_related/data/other/"    ;地形数据路径
 pic_dir  ="/cygdrive/D/ncl_related/pictures/vapour/sec/" ;图片存储路径

;读地形文件----------------------------
   ;lon1       fspan(0,360,4320)  ;gsn_xy    
   lat1     =fspan(90,-90,2160)  ;用于gsn_xy 地形图
   setfileoption("nc","MissingToFillValue",False)    ;设置缺省值
   elev= read_elev_data(data2_dir+"ETOPO5.DAT")
   ;printVarSummary(elev)
   elev0 = 1013.25*(1-elev*0.0065/288.15)^5.25145        ;将地形数据的单位：m转化为hPa
   copy_VarCoords(elev,elev0) 
   ;printVarSummary(elev0)

  stdlon=96;选择固定某一角度
    elev1=elev0(:,{stdlon})       ;固定某一经度的地形数据
    elev1@_FillValue= -9.96921e+36

; 数据处理-----------------------------------------------------
  ;;;>>>手动设置一个时间 用各个数据生成一个带unit的 timeset_ymdh 并且转换为一个字符串 方便文件中读入
   year   = 2008
   month  = 10
   day    = 26
   hour   = 18
   minute = 00
   second = 0
  ;设置时间格式读取数据####################################
   units  = "hours since 1-1-1 00:00:0.0"  ; 时间设置，注意单位
   opt = 0
   timeset_ymdh = cd_inv_calendar(year,month,day,hour,minute,second,units,opt)
   file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)

  ;计算等压面上的比湿
   data   = addfile(data1_dir+"fnl_"+file_time+".nc", "r") 
   rh     = data->RH_P0_L100_GLL0   ;lv_ISBL4/lat_0/lon_0，21层
   tmpdata= data->TMP_P0_L100_GLL0(5:25,:,: ) ;单位：K,lv_ISBL0/lat_0/lon_0,26层，选择对应的21层
   u      = data->UGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
   v      = data->VGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
   vvel   = data->VVEL_P0_L100_GLL0(:,:,: ) ;Pa/s,lv_ISBL4, lat_0, lon_0 
   pres1  =data->lv_ISBL4  ;Pa,21层
   lon    =data->lon_0
   lat    =data->lat_0
   pres   =pres1*0.01
   pres@units="hPa"

   h_lvvel  =vvel(:,:,{stdlon})



  ;计算各层散度
  div=uv2dvF_Wrap(u(:,::-1,:),v(:,::-1,:)) ;用这个球面间谐波函数计算散度更精确
  ;di = uv2dv_cfd (u,v,lat1,lon1, 3) 
  div=div*1.e5
  h_ldiv=div(:,:,{stdlon})






;dir="/cygdrive/d/ncl_related/data/"
wks = gsn_open_wks("png",pic_dir+"w_div"+file_time)


;------------绘图---------------------

res    = True
 cnres  = True
 vvelres  =True
 etores=True
 set_all_res(res,cnres,vvelres,etores)


;print(min(h_lvvel))
;print(max(h_lrh))
;exit
plot_div=gsn_csm_pres_hgt(wks,h_ldiv(:,{25:35}), cnres ) 
plot_vvel  = gsn_csm_pres_hgt(wks, h_lvvel(:,{25:35}), vvelres ) 
plot_hgt = gsn_csm_xy(wks,lat1,elev1,etores)
overlay(plot_div,plot_vvel) 
overlay(plot_div,plot_hgt)


draw(plot_div)
frame(wks)
end