;1种数据，用fnl数据算整层水汽通量组图，还没改完：
;比湿q:kg/kg 速度u:m/s 气压:hPa，g:m/s2; vibeta(q*u/g)→量级100kg/(m*s)
;参考函数：垂直积分 http://www.ncl.ucar.edu/Document/Functions/Built-in/vibeta.shtml
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"

begin

;数据读取
  data_dir="/cygdrive/D/ncl_related/data/fnl/2600_2812_12/"
  pic_dir="/cygdrive/D/ncl_related/pictures/vapour/"
  nc_files = systemfunc("ls "+data_dir+"fnl_*_*_*.nc")
  f = addfiles(nc_files,"r") ; f[0-12]
  ListSetType(f,"join"); 打开多个文件并且返回文件句柄(是列表)

  rh     = f[:]->$"RH_P0_L100_GLL0"$   ;lv_ISBL4/lat_0/lon_0，21层
  tmpdata1= f[:]->$"TMP_P0_L100_GLL0"$ ;单位：K,lv_ISBL0/lat_0/lon_0,26层，选择对应的21层
  tmpdata= tmpdata1(:,5:25,:,: ) ;单位：K,lv_ISBL0/lat_0/lon_0,26层，选择对应的21层

  u1      = f[:]->$"UGRD_P0_L100_GLL0"$  ;单位：m/s,lv_ISBL0, lat_0, lon_0
  u      = u1(:,5:25,:,: )  ;单位：m/s,lv_ISBL0, lat_0, lon_0

  v1      = f[:]->$"VGRD_P0_L100_GLL0"$  ;单位：m/s,lv_ISBL0, lat_0, lon_0
  v      = v1(:,5:25,:,: )  ;单位：m/s,lv_ISBL0, lat_0, lon_0

  pres1  =f[:]->$"lv_ISBL4"$  ;Pa,21层
  psfc   =f[:]->$"PRES_P0_L1_GLL0"$

;数据处理
  pres   =pres1*0.01 ;把Pa变成hPa
  pres@units="hPa"
  q      = mixhum_ptrh(conform(tmpdata,pres,1),tmpdata,rh,-2) ;每层都算了
   ; specific humidity (g/kg)，要求温度K,气压hPa
  q=q*0.001
  copy_VarCoords(rh,q)
  q@units="kg/kg"
  qu = new(dimsizes(rh),float)
  qv = new(dimsizes(rh),float)
  g=9.8
  qu= q*u
  qv= q*v
  copy_VarCoords(rh,qu)
  copy_VarCoords(rh,qv)
linlog=1
pbot=110000
ptop=10000     
vint_qu=vibeta(pres1(::-1),qu(:,lat_0|:,lon_0|:,lv_ISBL4|::-1),linlog,psfc,pbot,ptop)/g  
copy_VarCoords(psfc,vint_qu)
vint_qv=vibeta(pres1(::-1),qv(:,lat_0|:,lon_0|:,lv_ISBL4|::-1),linlog,psfc,pbot,ptop)/g  
copy_VarCoords(psfc,vint_qv)
mm=sqrt(vint_qu^2+vint_qv^2)
copy_VarMeta(vint_qu,mm)

printVarSummary(mm)
exit

;max1=min(vapdiv({25:30},{90:100}))
;print(max1)
;exit


;风暴数据,用于标记风暴位置或路径 
fiTY = "/cygdrive/D/ncl_related/data/TC/bio042008.txt"
; 获取文本文件的行数，相应的还有numAsciiCol函数用于获取列数
nrow = numAsciiRow(fiTY)
YYYYMMDDHH = new(nrow, "string")
lat = new(nrow, "float")
lon = new(nrow, "float")
;以字符串形式读入
data3 = asciiread(fiTY, -1, "string")
YYYYMMDDHH = str_get_field(data3, 3, ",")
lat = stringtofloat(str_get_field(data3, 7, ",")) *0.1
lon = stringtofloat(str_get_field(data3, 8, ",")) *0.1
DateChar = stringtochar(YYYYMMDDHH)
MM = chartostring(DateChar(:,5:6))
DD = chartostring(DateChar(:,7:8))
HH = chartostring(DateChar(:,9:10))


wks = gsn_open_wks("png", pic_dir+day+hour)
gsn_define_colormap(wks,"MPL_PuBuGn") 

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize=True
res@gsnLeftString = ""
res@gsnRightString = ""
res@tmYROn=False
res@tmXTOn=False
;res@pmTickMarkDisplayMode   = "Always" ;不显示次刻度线
;res@tmXBMajorLengthF = 0.005  ;主刻度线长度
;res@vpWidthF=0.6
;res@vpHeightF=0.4



;设置坐标轴标签粗细
res@tmXBLabelFont=22
res@tmYLLabelFont=22
res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
res@tmYLLabelFontHeightF =0.02




;设置地图
cnres                         = res                    ; plot mods desired
cnres@mpFillOn                = False
cnres@mpMinLatF               = 16.
cnres@mpMaxLatF               = 36.
cnres@mpMinLonF               = 80.
cnres@mpMaxLonF               = 105.
cnres@gsnMajorLonSpacing    =5. 
cnres@gsnMinorLonSpacing    =1.    
cnres@gsnMajorLatSpacing    =4. 
cnres@gsnMinorLatSpacing    =1
cnres@mpDataBaseVersion       = "MediumRes"
cnres@mpDataSetName           = "Earth..4"
cnres@mpOutlineSpecifiers     = (/"China:states"/)
  ;边界线
   cnres@mpOutlineBoundarySets           = "National"
   cnres@mpNationalLineColor                        = "black"
   cnres@mpProvincialLineColor                = "black"
   cnres@mpGeophysicalLineColor          = "black"
   cnres@mpNationalLineThicknessF         = 3
   cnres@mpProvincialLineThicknessF  = 3
   cnres@mpGeophysicalLineThicknessF = 3

;水汽通量填色
cnres@cnFillOn        = True                    ; turn on color
cnres@cnLinesOn       = True                   ; turn off contour lines
cnres@cnLineLabelsOn  = False                   ; tuen off line labels
cnres@cnInfoLabelOn   =False 
cnres@cnFillPalette         = "MPL_PuBuGn"
;cnres@cnLevelSelectionMode = "ManualLevels" ;先用这个试探数值范围，再统一
cnres@cnLevelSelectionMode  ="ExplicitLevels" 
cnres@cnLevels              =(/100,200,300,400,500,600,700,800,900,1000/)
cnres@cnFillColors          =(/20,30,40,50,60,70,80,90,100,110,120/)

cnres@lbLabelBarOn          = True
;cnres@lbLabelStride        = 2         ; 两个level只标一个标签
;cRes@lbLabelAutoStride     = True            ; let NCL determine label spacing
cnres@lbBoxLinesOn          = True; Label bar
;cnres@lbLeftMarginF        =-0.5   ;离左边图的间隔
cnres@lbOrientation         = "Vertical"         ; vertical label bar
cnres@pmLabelBarWidthF      = 0.05     ;色标变窄
;cnres@pmLabelBarHeightF    = 0.5    ;色标变短
cnres@pmLabelBarOrthogonalPosF  = 0.02 ;色标离图的距离
cnres@lbBoxEndCapStyle      = "TriangleBothEnds"
cnres@cnFillDrawOrder       ="PreDraw"


;--------------水汽通量矢量--------------


vRes = res
vRes@vcGlyphStyle            = "CurlyVector"   ; turn on curly vectors
;vRes@vcGlyphStyle    = "LineArrow"
vRes@vcLineArrowThicknessF = 3
vRes@vcMinDistanceF = 0.015
vRes@vcLineArrowHeadMaxSizeF = 0.01
vRes@vcLineArrowHeadMinSizeF = 0.01
vRes@vcRefMagnitudeF = 100.
vRes@vcRefLengthF    = 0.02
vRes@vcRefAnnoOn = True
vRes@vcRefAnnoFontHeightF =0.015
vRes@vcRefAnnoZone = 1.5
vRes@vcRefAnnoSide = "Top"
vRes@vcRefAnnoJust = "TopRight" ;默认值
;vRes@vcRefAnnoString1 = "10 gcm~S~-1~N~hPa~S~-1~N~s~S~-1~N~"
vRes@vcRefAnnoString1 = "100"
vRes@vcRefAnnoString2On = False
;vRes@vcRefAnnoParallelPosF=0.33
;vRes@vcRefAnnoOrthogonalPosF=-0.1 ;调整上下位置
;vRes@vcRefAnnoBackgroundColor=-1 ;将矢量注释的填充色设为透明
;vRes@vcRefAnnoPerimOn=False ;将矢量注释的边框去掉



;contour on contour
;plot = gsn_csm_contour_map_overlay(wks,vf,hgt,cRes, heightRes)
vplot= gsn_csm_vector(wks, vint_qu , vint_qv , vRes)
mm=smth9(mm,0.5,0.25,False)
plot=gsn_csm_contour_map(wks, mm, cnres)
;vectorplot= gsn_csm_vector(wks, qu , qv , vRes)
;heightplot=gsn_csm_contour(wks,hgt, heightRes)
;mapplot=gsn_csm_map(wks, mpres)


;叠加矢量和等值线
overlay(plot,vplot)
;overlay(plot,heightplot)
;overlay(mapplot,vectorplot)

;某时刻风暴位置标记
txres               = True                     ; text mods desired
txres@txFontHeightF = 0.03                     ; font smaller. default big
txres@txFontColor="red"
do i = 0, nrow-1
   if ((HH(i) .eq."00").and.(DD(i) .eq."28")) then
      dum = gsn_add_text(wks,plot,"~F37~p",lon(i), lat(i),txres) 
    end if
end do
draw(plot)
frame(wks)
end