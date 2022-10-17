load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
; 利用中国地面降水日值0.5°×0.5°格点的nc数据绘图
; 参考函数 http://www.ncl.ucar.edu/Document/Functions/Built-in/system.shtml
; 参考函数 http://www.ncl.ucar.edu/Document/Functions/Built-in/filedimdef.shtml
; 参考函数 http://www.ncl.ucar.edu/Document/Functions/Contributed/readAsciiTable.shtml

begin
;-------------------读取数据----------------------------------------------------

 
  dir    = "/cygdrive/D/ncl_related/data/cma_0.5_0.5_rain/nc/26_28/"  
  fRoot  = "SURF_CLI_CHN_PRE_DAY_GRID_0.5_" 
  pic_dir="/cygdrive/D/ncl_related/pictures/rain/"
  nc_files = systemfunc("ls "+dir+fRoot+"*.nc")
  f        = addfiles(nc_files,"r") ; f[0-3]
  ListSetType(f,"join"); 打开多个文件并且返回文件句柄(是列表)
  p        = f[:]->$"op"$ 

  ;printVarSummary(p)
  ;exit


wks = gsn_open_wks("png", pic_dir+"cma")

;---Set up some map resources.
  mpres                       = True
  mpres@mpFillOn              = False ;必须设置，不然中国以外的地方为缺测值会被填色
  mpres@gsnMaximize           = True             ; Maximize plot in frame.
  mpres@gsnDraw               = False            ; Will draw later
  mpres@gsnFrame              = False            ; Don't advance the frame
  mpres@gsnLeftString         = ""
  mpres@gsnRightString        = ""
  mpres@tmYROn                =False
  mpres@tmXTOn                =False
  mpres@gsnAddCyclic          =False
  ;mpres@pmTickMarkDisplayMode= "Always"

  mpres@mpMinLatF             = 16.
  mpres@mpMaxLatF             = 36.
  mpres@mpMinLonF             = 80
  mpres@mpMaxLonF             = 105
  mpres@gsnMajorLonSpacing    =5. 
  mpres@gsnMinorLonSpacing    =1.    
  mpres@gsnMajorLatSpacing    =4. 
  mpres@gsnMinorLatSpacing    =1


  ;设置坐标轴标签粗细
  mpres@tmXBLabelFont         =22
  mpres@tmYLLabelFont         =22
  mpres@tmXBLabelFontHeightF  =0.02 ; x轴标签字号
  mpres@tmYLLabelFontHeightF  =0.02


  mpres@mpDataBaseVersion       = "MediumRes"
  mpres@mpDataSetName           = "Earth..4"
  mpres@mpOutlineSpecifiers     = (/"China:states"/)
  ;边界线
  mpres@mpOutlineBoundarySets   = "National"
  mpres@mpNationalLineColor     = "black"
  mpres@mpProvincialLineColor   = "black"
  mpres@mpGeophysicalLineColor  = "black"
  mpres@mpNationalLineThicknessF    = 3
  mpres@mpProvincialLineThicknessF  = 3
  mpres@mpGeophysicalLineThicknessF = 3


;画等值线图
preRes                  = mpres
preRes@cnFillOn         = True  ;cn means ContourPlot
preRes@cnLinesOn        = True   
preRes@cnMonoLineColor  = False ;线条不使用同一个颜色
preRes@cnFillPalette    = "wh-bl-gr-ye-re"
preRes@cnLineLabelsOn   = False   
preRes@cnInfoLabelOn    = False ; set false will not set an information label
;统一子图的等值线，以共用色标
preRes@cnLevelSelectionMode = "ManualLevels" 
preRes@cnMinLevelValF   = 10
preRes@cnMaxLevelValF   = 120
preRes@cnLevelSpacingF  = 10
preRes@cnLineThicknessF = 2
preRes@lbLabelBarOn     = False ;关闭每幅图的色标

;组图属性
pres=True
pres@lbBoxLinesOn          = True ; Label bar
pres@lbLabelAutoStride     = True            ; let NCL determine label spacing
;preRes@lbOrientation      = "Vertical"         ; vertical label bar 色标垂直
pres@pmLabelBarWidthF      = 0.8     ;色标变窄
;pres@pmLabelBarHeightF    = 0.02    ;色标变细
pres@lbLabelFontHeightF    =0.01

pres@lbBoxEndCapStyle      = "TriangleBothEnds"
pres@gsnPanelCenter        =True ;对齐为默认值
pres@gsnPanelLabelBar      =True
pres@gsnPanelFigureStrings =(/"20081026","20081027","20081028"/)
pres@gsnPanelFigureStringsFontHeightF=0.01
pres@amJust                ="topLeft" ;图注放在左上方



;绘图
plot=new(3, graphic)
p   = smth9_Wrap(p,0.5,0.25,False)
p   = smth9_Wrap(p,0.5,0.25,False)
do iplot     =0,2
  plot(iplot)= gsn_csm_contour_map(wks,p(iplot,:,:),preRes)  
end do
gsn_panel(wks, plot, (/1,3/), pres) 
;draw(plot) ;不管怎样都会出两张图，但是关掉这个的话出的是两张组图
;frame(wks) ;这个关掉就OK     

 end                                
                                 
                                 

