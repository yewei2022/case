;用fnl数据算整层水汽通量,~散度,关于通量和~散度的填色可切换：加粗地图省界国界边界线
;比湿q:kg/kg 速度u:m/s 气压:Pa，g:m/s2; vibeta(q*u/g)→量级100kg/(m*s)
;添加累计降水量大于25mm的站点
;参考函数：垂直积分 http://www.ncl.ucar.edu/Document/Functions/Built-in/vibeta.shtml
;改时间只改一处19
;代码路径  ncl /mnt/h/d/ncl_related/code/vertical_integral_vf.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"

begin

;手动设置一个时间,需要改两处
;用各个数据生成一个带unit的 timeset_ymdh 
;转换为一个字符串 方便文件中读入
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

;数据读取，再分析数据
data_dir="/mnt/h/d/data/fnl/"
  pic_dir="/mnt/h/d/ncl_related/pictures/vapour/plane/"
  data   = addfile(data_dir+"fnl_"+file_time+".nc", "r") 
  rh     = data->$"RH_P0_L100_GLL0"$   ;lv_ISBL4/lat_0/lon_0，21层
  tmpdata= data->TMP_P0_L100_GLL0(5:25,:,: ) ;单位：K,lv_ISBL0/lat_0/lon_0,26层，选择对应的21层
  u      = data->UGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
  v      = data->VGRD_P0_L100_GLL0(5:25,:,:)  ;单位：m/s,lv_ISBL0, lat_0, lon_0
  pres1  =data->lv_ISBL4  ;Pa,21层
  psfc   =data->PRES_P0_L1_GLL0
  lon1    =data->lon_0
  lat1    =data->lat_0
;fnl数据处理-----------------------
  pres   =pres1*0.01 ;把Pa变成hPa,用于计算比湿
  pres@units="hPa"
  q      = mixhum_ptrh(conform(tmpdata,pres,0),tmpdata,rh,-2) ;每层都算了
   ; specific humidity (g/kg)，要求温度K,气压hPa
  q=q*0.001
  copy_VarCoords(rh,q)
  q@units="kg/kg"
  qu = new(dimsizes(rh),float)
  qv = new(dimsizes(rh),float)
  g=9.8
  qu= q*u/g
  qv= q*v/g
  copy_VarCoords(rh,qu)
  copy_VarCoords(rh,qv)
  linlog=1
  pbot=110000
  ptop=10000     
  vint_qu=vibeta(pres1(::-1),qu(lat_0|:,lon_0|:,lv_ISBL4|::-1),linlog,psfc,pbot,ptop)  ;kg/kg，pa
  copy_VarCoords(psfc,vint_qu)
  vint_qv=vibeta(pres1(::-1),qv(lat_0|:,lon_0|:,lv_ISBL4|::-1),linlog,psfc,pbot,ptop) 
  copy_VarCoords(psfc,vint_qv)

  ;------------------------整层水汽通量散度大小--------------------------
  mm1=uv2dv_cfd(qu,qv,lat1,lon1,2)
  copy_VarCoords(rh,mm1) 
  mm=vibeta(pres1(::-1),mm1(lat_0|:,lon_0|:,lv_ISBL4|::-1),linlog,psfc,pbot,ptop)  
  mm=mm*1.e4 ;单位、量级1e-4kg/(m**2*s)
  copy_VarCoords_not_n(rh,mm,0) 
  ;printVarSummary(mm)
  ;exit

  ;-------------------整层水汽通量数值---------------------------------------------
  ;mm=sqrt(vint_qu^2+vint_qv^2)  ;单位、量级100kg/(m*s)
  ;copy_VarCoords(psfc, mm)
  ;--------------------------------------------------

;数据读取，风暴数据,用于标记风暴位置或路径 
  file_tc = "/mnt/h/d/data/TC/bio042008.txt"
  ; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
  nrow = numAsciiRow(file_tc)
  YYYYMMDDHH = new(nrow, "string")
  lat2 = new(nrow, "float")
  lon2 = new(nrow, "float")
  ;以字符串形式读入
  data3 = asciiread(file_tc, -1, "string")
  YYYYMMDDHH = str_get_field(data3, 3, ",")
  lat2 = stringtofloat(str_get_field(data3, 7, ",")) *0.1
  lon2 = stringtofloat(str_get_field(data3, 8, ",")) *0.1
  DateChar = stringtochar(YYYYMMDDHH)
  MM = chartostring(DateChar(:,5:6))
  DD = chartostring(DateChar(:,7:8))
  HH = chartostring(DateChar(:,9:10))
  tctime_ind = str_match_ind(YYYYMMDDHH, path_time) ; 返回YYYYMMDDHH中和path_time（文件时间）匹配的位置ind
  stdlon=0.
  stdlat=0.
  do i = 0, nrow-1
    if (i.eq.tctime_ind) then ;最方便
    stdlon=lon2(i)
    stdlat=lat2(i)
    end if
  end do

;--------------------读站点数据--------------------
  filepath = "/mnt/h/d/data/precipitation/6h/"  ;参数文件路径
  name="rain12_py"
  f       = asciiread(filepath+name+".txt",-1,"string")
  sta     = stringtointeger(str_get_field(f,1," "))  ;读取站号
  lat3    = stringtofloat(str_get_field(f,2," ")) ;获取纬度值lon
  lon3    = stringtofloat(str_get_field(f,3," ")) ;获取经度值lat
  R1      = stringtofloat(str_get_field(f,6," ")) ;获取第2618时次站点数据
  R2      = stringtofloat(str_get_field(f,7," ")) ;获取第2700个时次站点数据
  R3      = stringtofloat(str_get_field(f,8," ")) ;获取第2706个时次站点数据
  R4      = stringtofloat(str_get_field(f,9," ")) ;获取第2712个时次站点数据
  npts = numAsciiRow(filepath+name+".txt") ;获取行数-站点数
  ;注意：这里的txt文件最后一行数据下面需要空白行，才能读出文件的正确行数
  ;大多数txt文件都满足要求，但还是需要核实一下
  R_old=new(npts,float)
  R_old=R1+R2+R3+R4
  ;设置缺测值
  R_old@_FillValue = 9999.0

  ;---Open shapefile and read Mississippi River Basin lat/lon values.
  f1      = addfile("/mnt/h/d/data/TPshp/DBATP/DBATP_Polygon.shp", "r")
  mrb_lon = f1->x
  mrb_lat = f1->y
  ;---Start with data all missing,mask掉高原以外区域
  R = new(dimsizes(R_old),typeof(R_old),R_old@_FillValue)
  do i=0,npts-1
    R(i)=where(gc_inout(lat3(i),lon3(i),mrb_lat,mrb_lon), R_old(i), R_old@_FillValue)
  end do
  lat_new = new(npts,float,-999)
  lon_new = new(npts,float,-999)
  indexes=ind(R.ge.25.) ;剔除了25mm以下的降水
  if (.not.any(ismissing(indexes))) then 
    npts_range = dimsizes(indexes)   ; # of points in this range.
    lat_new(0:npts_range-1) = lat3(indexes)
    lon_new(0:npts_range-1) = lon3(indexes)
  end if
  delete(indexes)     


wks = gsn_open_wks("x11", pic_dir+"vapdiv_"+file_time)

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize=True
res@gsnLeftString = ""
res@gsnRightString = ""
res@tmYROn=False
res@tmXTOn=False
;res@pmTickMarkDisplayMode   = "Always" ;不显示次刻度线
res@tmBorderThicknessF  =4    ; border 6x as thick,图框粗细
res@tmXBMajorThicknessF =4
res@tmXBMinorThicknessF =4
res@tmYLMajorThicknessF =4 
res@tmYLMinorThicknessF =4 


;设置坐标轴标签粗细，图例在左时候设置的
res@tmEqualizeXYSizes   =True 
; 上个命令，(tmXBLabelFontHeightF,(tmXBMajorLengthF,(tmXBMinorLengthF设置好后
; 其他坐标轴就不用重复设置了
res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号

;设置地图
cnres                         = res                    ; plot mods desired
cnres@mpFillOn                = False
cnres@mpMinLatF               = 10.
cnres@mpMaxLatF               = 60.
cnres@mpMinLonF               = 70.
cnres@mpMaxLonF               = 135.
cnres@gsnMajorLonSpacing    =10. 
cnres@gsnMinorLonSpacing    =5.    
cnres@gsnMajorLatSpacing    =10. 
cnres@gsnMinorLatSpacing    =5
cnres@mpDataSetName              = "/mnt/h/d/ChinaMap/Earth..4"
cnres@mpDataBaseVersion          = "MediumRes" ; or "Ncarg4_1"
cnres@mpAreaMaskingOn            = True
cnres@mpMaskAreaSpecifiers       = (/"China"/) ;这和上两个命令好像是没啥用
cnres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)

cnres@mpOutlineBoundarySets      = "National"
cnres@mpNationalLineThicknessF   = 2
cnres@mpProvincialLineThicknessF = 2
cnres@mpGeophysicalLineThicknessF= 2

;-------------------水汽通量散度填色-----------------------------------------
cnres@cnFillOn        = False                    ; turn on color
cnres@cnLinesOn       = False                   ; turn off contour lines
cnres@cnLineLabelsOn  = False                   ; tuen off line labels
cnres@cnInfoLabelOn   = False 
; cnres@cnFillPalette       = "BlueWhiteOrangeRed"  ;水汽通量散度填色
; cnres@cnLevelSelectionMode= "ManualLevels" 
; cnres@cnMinLevelValF      =-10
; cnres@cnMaxLevelValF      = 10
; cnres@cnLevelSpacingF     = 1

; ;---------------------水汽通量填色------------------------------------------
; ;cnres@cnFillPalette         = "MPL_PuBuGn" ;青色系,水汽通量填色
; ;cnres@cnLevelSelectionMode = "ManualLevels" ;先用这个试探数值范围，再统一
; ;cnres@cnLevelSelectionMode  ="ExplicitLevels" 
; ;cnres@cnLevels              =(/100,200,300,400,500,600,700,800,900,1000/)
; ;cnres@cnFillColors          =(/20,30,40,50,60,70,80,90,100,110,120/)
; ;-----------------------------------------------------------------------

; cnres@lbLabelBarOn          = True
; cnres@lbLabelStride         = 2         ; 两个level只标一个标签
; ; cnres@lbLabelAutoStride   = True            ; let NCL determine label spacing
; cnres@lbLabelFontHeightF    = 0.02
; cnres@lbBoxLinesOn          = True; Label bar
; ;cnres@lbLeftMarginF        =-0.5   ;离左边图的间隔
; cnres@lbOrientation         = "Vertical"         ; vertical label bar
; cnres@pmLabelBarWidthF      = 0.05     ;色标变窄
; ;cnres@pmLabelBarHeightF    = 0.5    ;色标变短
; cnres@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
; cnres@lbBoxEndCapStyle      = "TriangleBothEnds"
; cnres@cnFillDrawOrder       ="PreDraw"


; ;--------------水汽通量矢量--------------
; vRes = res
; vRes@vcGlyphStyle            = "CurlyVector"   ; turn on curly vectors
; ;vRes@vcGlyphStyle           = "LineArrow"
; vRes@vcLineArrowThicknessF   = 3
; vRes@vcMinDistanceF          = 0.015
; vRes@vcLineArrowHeadMaxSizeF = 0.008
; vRes@vcLineArrowHeadMinSizeF = 0.008
; vRes@vcRefMagnitudeF      = 400.
; vRes@vcRefLengthF         = 0.02
; vRes@vcRefAnnoOn          = True
; vRes@vcRefAnnoFontHeightF = 0.015
; vRes@vcRefAnnoZone = 1.5
; vRes@vcRefAnnoSide = "Top"
; vRes@vcRefAnnoJust = "TopRight" ;默认值
; ;vRes@vcRefAnnoString1 = "10 gcm~S~-1~N~hPa~S~-1~N~s~S~-1~N~"
; vRes@vcRefAnnoString1 = "400"
; vRes@vcRefAnnoString2On = False
; ;vRes@vcRefAnnoParallelPosF=0.33
; ;vRes@vcRefAnnoOrthogonalPosF=-0.1 ;调整上下位置
; ;vRes@vcRefAnnoBackgroundColor=-1 ;将矢量注释的填充色设为透明
; ;vRes@vcRefAnnoPerimOn=False ;将矢量注释的边框去掉


;------------------------绘图---------------------------------

; mm=smth9(mm,0.5,0.25,False)
; vplot= gsn_csm_vector(wks, vint_qu , vint_qv , vRes)
plot=gsn_csm_contour_map(wks, mm, cnres)
; overlay(plot,vplot)

; ;绘制≥25mm的站点位置--------------
; gsres               = True                     ; text mods desired
; gsres@txFontColor   ="red"
; gsres@txFontHeightF = 0.03
; pmid = gsn_add_text(wks,plot,"~F35~z",lon_new,lat_new,gsres)

; ; 贴上强降水区域的方框
; ;************************************************
;  ypts = (/32.0,  32.0,  28.0, 28.0, 32.0/)
;  xpts = (/89.0,  98.0,  98.0, 89.0, 89.0/)

;  resp                  = True                      ; polyline mods desired
;  resp@gsLineColor      = "black"                     ; color of lines
;  resp@gsLineThicknessF = 6.0                       ; thickness of lines
;  dum_box = new(4,graphic)  
; ; draw each line separately. Each line must contain two points.
;   do i = 0 , 3
;     dum_box(i)=gsn_add_polyline(wks,plot,xpts(i:i+1),ypts(i:i+1),resp)      
;   end do 

; ;某时刻风暴位置标记
; txres               = True                     ; text mods desired
; txres@txFontHeightF = 0.03                     ; font smaller. default big
; txres@txFontThicknessF=6
; txres@txFontColor="red"
; dum_storm = gsn_add_text(wks,plot,"~F37~p",stdlon,stdlat,txres) 


; ;高原轮廓线
; shp="/mnt/h/d/data/TPshp/DBATP/DBATP_Line.shp"  
;   lnres                  = True     
;   lnres@gsLineColor      ="orange"
;   lnres@gsLineThicknessF = 8.0           ; 6x thickness  
;   shp_plot               = gsn_add_shapefile_polylines(wks,plot,shp,lnres) 

draw(plot)
frame(wks)
end