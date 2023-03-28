
;  ncl /mnt/d/case/code/by_ncl/FNL/prep_6h_bar.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; bar_10.ncl
;
; Concepts illustrated:
;   - Drawing horizontal filled bars using gsnXYBarChart
;   - Changing the aspect ratio of a bar plot
;   - Setting the minimum/maximum value of the X and Y axis in a bar plot
;   - Explicitly setting tickmarks and labels on the left Y axis
;   - Turning off tickmarks, but not the labels
;   - Overlaying XY plots on each other
;   - Drawing grid lines on an XY plot
;*************************************************
; This script requires NCL V6.4.0 or later to 
; run. See bar_old_10.ncl for an older way of
; creating this plot.
;*************************************************
; This file is loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

pic_dir ="/mnt/d/case/pic/rain/"
filepath= "/mnt/d/case/data/precipitation/6h/prep_6h_3sta.txt"  ;参数文件路径
f1      = asciiread(filepath,-1,"string")
cona     = stringtointeger(str_get_field(f1,1," "))  ;读取站号
bome     = stringtofloat(str_get_field(f1,2," ")) 
zayu    = stringtofloat(str_get_field(f1,3," ")) 
stdarr = new((/8,3/),"float")
stdarr(:,0)=cona
stdarr(:,1)=bome
stdarr(:,2)=zayu

wks = gsn_open_wks("eps",pic_dir+"prep_6h_bar")  ; send graphics to PNG file  
gsn_define_colormap(wks,"wgne15")
colors=(/14,8,4/) ;红 绿 蓝


; ;========若是柱子在时刻之间，8个数据=========================
; sres@trXMinF = 1
; sres@trXMaxF = 9
; sres@tmXBValues        = ispan(1, 9, 1)
; sres@tmXBLabels = (/"26/00","26/06","26/12","26/18","27/00","27/06","27/12","27/18","28/00"/)
; xpos1=fspan(1.35,8.36,8) ;位置调整，先定头，再定尾


;===============若是柱子在时刻上===============================
sres = True
sres@vpWidthF = 0.5
sres@vpHeightF = 0.4
sres@tmXTOn=False
sres@tmYROn=False
sres@vpXF = .15
sres@trXMinF = 0.4
sres@trXMaxF = 8.6
sres@trYMinF = 0
sres@trYMaxF = 60
sres@gsnDraw = True
sres@gsnFrame = False
sres@gsnXYBarChart = True
sres@gsnXYBarChartBarWidth = 0.25           ; change bar widths
sres@tmXBMode          = "Explicit"         ; explicit labels
sres@tmXBValues        = ispan(1, 8, 1)
sres@tmXBLabels = (/"26/06","26/12","26/18","27/00","27/06","27/12","27/18","28/00"/)
sres@gsnRightString = ""
sres@tmEqualizeXYSizes=True
sres@tmXBLabelFontHeightF=0.018
sres@tmXBLabelAngleF=45
sres@gsnLeftStringFontHeightF=0.018
; sres@gsnLeftString="(d)"
sres@gsnLeftStringOrthogonalPosF=-0.15
sres@gsnLeftStringParallelPosF =0.018
sres@txPerimOn=True

sres@tiYAxisString = "Precipitation (mm)"

xpos1=fspan(0.751,7.758,8) ;位置调整，先定头，再定尾
intep=sres@gsnXYBarChartBarWidth
  
sres@gsnXYBarChartColors = colors(0)
plot1 = gsn_csm_xy(wks,xpos1,stdarr(:,0),sres)   ; draw each time series
sres@gsnXYBarChartColors = colors(1)         ; seperately, not
plot2 = gsn_csm_xy(wks,xpos1+intep,stdarr(:,1),sres)   ; advancing the frame        ; but tweaking where
sres@gsnXYBarChartColors = colors(2)         ; seperately, not
plot3 = gsn_csm_xy(wks,xpos1+intep*2.,stdarr(:,2),sres)   ; advancing the frame        ; but tweaking where

;==============填充图案=============================
resp =True
resp@gsFillIndex  = 14
resp@gsFillScaleF = 1
resp@gsFillColor  = "black"

cona_x=new((/8,4/),float)
cona_y=new((/8,4/),float)
xpos11=xpos1-intep/2.
do i=0,7  
  cona_x(i,:)=(/xpos11(i),xpos11(i)+intep,xpos11(i)+intep,xpos11(i)/)
  cona_y(i,:)=(/0,0,stdarr(i,0),stdarr(i,0)/)
end do
dum_cona=new(8, graphic)
do i=0,7
  dum_cona(i) = gsn_add_polygon(wks,plot1,cona_x(i,:),cona_y(i,:),resp) 
end do  

bome_x=new((/2,4/),float)
bome_y=new((/2,4/),float)
xpos11=xpos1+intep/2.
;波密只有ind=1,2时填充
do i=2,3  
  bome_x(i-2,:)=(/xpos11(i),xpos11(i)+intep,xpos11(i)+intep,xpos11(i)/)
  bome_y(i-2,:)=(/0,0,stdarr(i,1),stdarr(i,1)/)
end do

dum_bome=new(2, graphic)
do i=0,1
  dum_bome(i) = gsn_add_polygon(wks,plot2,bome_x(i,:),bome_y(i,:),resp) 
end do  


;Í¼Àý
lbres                    = True          ; labelbar only resources
lbres@vpWidthF           = 0.15           ; labelbar width
lbres@vpHeightF          = 0.05           ; labelbar height
lbres@lbBoxMajorExtentF  = 0.32          ; puts space between color boxes
lbres@lbMonoFillPattern  = True          ; Solid fill pattern
lbres@lbLabelFontHeightF = 0.018        ; font height. default is small
lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
lbres@lbPerimOn          = False
lbres@lgPerimColor      = "white"
lbres@lbFillColors       = colors(0)
gsn_labelbar_ndc(wks,1,"Cona",0.18,0.34,lbres) ; draw left labelbar column
lbres@lbFillColors       = colors(1)
gsn_labelbar_ndc(wks,1,"Bomi",0.33,0.34,lbres) ; draw left labelbar column
lbres@lbFillColors       = colors(2)
gsn_labelbar_ndc(wks,1,"Zayu",0.48,0.34,lbres) ; draw left labelbar column

draw(plot1)
draw(plot2)
frame(wks)      
end