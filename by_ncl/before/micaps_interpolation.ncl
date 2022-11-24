;读取4个时次micaps地面相同站点降水资料，将缺测值赋值为0，相加再插值，画日降水分布图
;文件读取方式 http://www.ncl.ucar.edu/Applications/list_io.shtml
;插值函数http://www.ncl.ucar.edu/Document/Functions/interp.shtml
;问题，插值后数据误差，data({28},{92})应该是97，但是变成83，还能接受

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;1.读取数据


  filepath = "/cygdrive/D/ncl_related/data/micaps_rain/rain_common3.txt"  ;参数文件路径
  picture_path="/cygdrive/D/ncl_related/pictures/rain/"
  f=asciiread(filepath,-1,"string")
  station = stringtofloat(str_get_field(f,2," ")) ;从数组lines中获取站号
  lat = stringtofloat(str_get_field(f,3," ")) ;从数组lines中获取纬度值lon
  lon = stringtofloat(str_get_field(f,4," ")) ;从数组lines中获取经度值lat
  R1  = stringtofloat(str_get_field(f,5," ")) ;从数组lines中获取第1个时次站点数据
  R2  = stringtofloat(str_get_field(f,6," ")) ;从数组lines中获取第1个时次站点数据
  R3  = stringtofloat(str_get_field(f,7," ")) ;从数组lines中获取第1个时次站点数据
  R4  = stringtofloat(str_get_field(f,8," ")) ;从数组lines中获取第1个时次站点数据



;相同站点数据相加
 ;设置缺测值
 R1@_FillValue = 9999.0
 R2@_FillValue = 9999.0
 R3@_FillValue = 9999.0
 R4@_FillValue = 9999.0

 valid1 = ind(ismissing(R1)) ; 找到数组a中所有缺测值的位置
 R1(valid1) = 0 ; 将缺测值所在位置元素赋值为0
 valid2 = ind(ismissing(R2)) ; 找到数组a中所有缺测值的位置
 R2(valid2) = 0 ; 将缺测值所在位置元素赋值为0
 valid3= ind(ismissing(R3)) ; 找到数组a中所有缺测值的位置
 R3(valid3) = 0 ; 将缺测值所在位置元素赋值为0
 valid4 = ind(ismissing(R4)) ; 找到数组a中所有缺测值的位置
 R4(valid4) = 0 ; 将缺测值所在位置元素赋值为0


 R=new(3322,float)
 R=R1+R2+R3+R4

 ;print(R(607))
 ;exit

 ;print(lon(607))
 ;exit

;创建存放插值后生成数据的数组 在北纬0-50，东经70-130这个矩形框内插值

  olon = new(241,"float");
  olat = new(201,"float");
  data = new((/201,241/),"float")

  olon=fspan(70, 130, 241)
  olat=fspan(0, 50, 201)


;2.设置数组属性，为了符合netcdf规定的数据格式，使函数能够识别经纬度
     olon!0          = "lon"
     olon@long_name  = "lon"
     olon@units      = "degrees_east"
     olon&lon        = olon

     olat!0          = "lat"
     olat@long_name  = "lat"
     olat@units      = "degrees_north"
     olat&lat        = olat

;3.调用插值函数

  R@_FillValue =  9.96921e+36  
  rscan = (/10,2,0.5/)   ;连续的有效半径大小，最大为10，依次递减
  data  = obj_anal_ic_Wrap(lon,lat,R,olon,olat,rscan, False)  ;Cressman插值

;4.---Attach the coordinate arrays
  data!0    = "lat"
  data!1    = "lon"
  data&lat  = olat
  data&lon  = olon


  ;printVarSummary(data1)

  ;printMinMax(data, 0)
  ;print(data({28},{92}))
  ;exit
  ; fill_value in data
  ;oned_data = ndtooned(R)
  ;data_check = dimsizes(ind(ismissing(oned_data).eq.True)) ; 
  ;print("fill_values' number in data1:")
  ;print(data_check) ;将缺测赋值为0后是1，未处理缺测值是3255



wks = gsn_open_wks("png", picture_path+"rain_27_impr")

;-------------------属性设置----------------------------------------------------
  res          = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@gsnLeftString  = "" ;左方字符串
  res@gsnRightString = ""
  res@gsnAddCyclic  = False
  ;设置坐标轴标签粗细
  res@gsnMaximize=True
  res@tmXBLabelFont=22
  res@tmYLLabelFont=22
  res@tmXBLabelFontHeightF = 0.02 ; x轴标签字号
  res@tmYLLabelFontHeightF =0.02
;地图
  res@mpFillOn       = False
  res@mpMinLatF               = 15.
  res@mpMaxLatF               = 35.
  res@mpMinLonF               = 80.
  res@mpMaxLonF               = 105.
  res@mpDataBaseVersion       = "MediumRes"
  res@mpDataSetName           = "Earth..4"
  res@mpOutlineSpecifiers     = (/"China:states"/)
  res@pmTickMarkDisplayMode  = "Always"
  res@tmXBMajorLengthF = 0.005 ;底部主刻度线长度

  ;抹掉中国行政区以外位置，为啥无效？
  res@mpAreaMaskingOn = True   ;使能填充覆盖
  res@mpMaskAreaSpecifiers = (/"China:states"/)   ;China:states

  ;等值线
  ;画等值线图
  res = res
  res@cnFillOn = True  ;cn means ContourPlot
  res@cnLinesOn = True   
  res@cnMonoLineColor = False ;线条不使用同一个颜色
  res@cnSpanLinePalette = True
  res@cnFillPalette = "prcp_1"
  ;res@cnLinePalette = "prcp_1" ;线条的颜色
  res@cnLineLabelsOn = False   
  res@cnInfoLabelOn = False ; set false will not set an information label
  res@cnLevelSelectionMode = "ManualLevels" 
  res@cnLineThicknessF = 2
  res@cnLevelSpacingF = 10 

  res@lbBoxLinesOn = True ; Label bar
  res@lbLabelAutoStride     = True            ; let NCL determine label spacing
  res@lbOrientation        = "Vertical"         ; vertical label bar 
  res@pmLabelBarWidthF = 0.05     ;色标变窄
  ;res@pmLabelBarHeightF = 0.02    ;色标变细
  res@lbBoxEndCapStyle = "TriangleBothEnds"

;绘图
  data= smth9_Wrap(data,0.5,0.25,False)
  plot = gsn_csm_contour_map(wks,data,res) 


;绘制三个站点位置
lat2=new(3,double)
lon2=new(3,double)
lat2!0="lat"
lat2@units="degrees_north"
lon2!0="lon"
lon2@units="degrees_east"




lat2(0)=olat({28.00}) ;错那
lon2(0)=olon({92.00})
lat2(1)=olat({29.50}) ;波密
lon2(1)=olon({95.50})
lat2(2)=olat({28.50}) ;察隅
lon2(2)=olon({97.25})


;站点标记
mkres               = True                     ; text mods desired
mkres@gsMarkerColor="black"
mkres@gsMarkerIndex=16
mkres@gsMarkerSizeF=0.02
mkres@gsMarkerThicknessF=0.4

;do i =0,2
dum = gsn_add_polymarker(wks,plot,lon2,lat2,mkres) 
;end do 


draw(plot)
frame(wks)
end
