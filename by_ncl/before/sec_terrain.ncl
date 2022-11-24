;0305 绘制96E地形剖面

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;---------------------------------------绘图属性设置------------------------------------------
undef("set_all_res")
procedure set_all_res(res,tmpres,etores)
 local res,tmpres,etores
 begin
 
  ;--------------------------------------------------
   res = True
   res@gsnDraw = False
   res@gsnFrame = False
   res@gsnMaximize   =True
   res@tmXTOn = False ;关刻度
   res@tmYROn = False
   tmpres@tmYRMode          = "Automatic"          ; turn off special labels on right axis   
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



  ;---------温度-------------
   tmpres=res
   tmpres@cnFillOn          = True                  ; turn on color
   tmpres@cnFillPalette     = "BlueRedGray"              ; set color map
   tmpres@cnLinesOn         = False                 ; no contour lines
   tmpres@cnLevelSelectionMode="ManualLevels" 
   tmpres@lbBoxLinesOn      = True ; Label bar
   ;tmpres@lbOrientation    = "vertical"            ; 垂直坐标
   tmpres@pmLabelBarOrthogonalPosF  = 0.08
   tmpres@lbBoxEndCapStyle  = "TriangleBothEnds"
   ;tmpres@pmLabelBarWidthF = 0.05     ;色标变窄
   tmpres@pmLabelBarHeightF =0.05
   tmpres@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
   tmpres@cnFillDrawOrder   = "PreDraw"


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
 pic_dir  ="/cygdrive/D/ncl_related/pictures/terrain/" ;图片存储路径

;读地形文件----------------------------
   ;lon1       fspan(0,360,4320)  ;gsn_xy    
   lat1     =fspan(90,-90,2160)  ;用于gsn_xy 地形图
   setfileoption("nc","MissingToFillValue",False)    ;设置缺省值
   elev= read_elev_data(data2_dir+"ETOPO5.DAT")
   ;printVarSummary(elev)
   elev0 = 1013.25*(1-elev*0.0065/288.15)^5.25145        ;将地形数据的单位：m转化为hPa
   copy_VarCoords(elev,elev0) 
   ;printVarSummary(elev0)

  stdlon=95;选择固定某一角度
    elev1=elev0(:,{stdlon})       ;固定某一经度的地形数据
    elev1@_FillValue= -9.96921e+36

; 数据处理-----------------------------------------------------
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
   file_format = "%Y%N%D_%H_%M" ; 文件夹中的形式YYYYMMDD_HH_MM
   file_time  = cd_string(timeset_ymdh,file_format)

  ;计算等压面上的比湿
   data   = addfile(data1_dir+"fnl_"+file_time+".nc", "r") 
   rh     = data->RH_P0_L100_GLL0   ;lv_ISBL4/lat_0/lon_0，21层
   tmpdata= data->TMP_P0_L100_GLL0(5:25,:,: ) ;单位：K,lv_ISBL0/lat_0/lon_0,26层，选择对应的21层
   pres1  =data->lv_ISBL4  ;Pa,21层
   pres   =pres1*0.01
   pres@units="hPa"

   h_ltmp  =tmpdata(:,:,{stdlon})
   h_ltmp!0="level"
   h_ltmp&level=pres







;dir="/cygdrive/d/ncl_related/data/"
wks = gsn_open_wks("png",pic_dir+file_time+"_terrain_tmp_"+stdlon+"E")


;------------绘图---------------------

res    = True
tmpres  = True
 etores=True
 set_all_res(res,tmpres,etores)


;print(min(h_lvvel))
;print(max(h_ltmp))
;exit


plot_tmp   = gsn_csm_pres_hgt(wks,h_ltmp(:,{25:35}), tmpres) 
plot_hgt = gsn_csm_xy(wks,lat1,elev1,etores) 
overlay(plot_tmp,plot_hgt)

draw(plot_tmp)
frame(wks)
end