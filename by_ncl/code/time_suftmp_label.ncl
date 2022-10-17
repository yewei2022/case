;读取fortran提取的micaps地面温度时间演变
;为了标注出降雪时间段，将数据分段处理
;出现问题：由于分段导致，折线断掉，之后处理：手动连线
; ncl /mnt/h/d/ncl_related/code/time_suftmp_label.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"


begin

;1.读取数据
  filepath = "/mnt/h/d/data/tmp/"  ;参数文件路径
  picture_path="/mnt/h/d/ncl_related/pictures/tmp/"
  data_cuo = asciiread(filepath+"surtmp55690.txt",-1,"string") ;以字符串形式读取参数文件入数组
  tmp_cuo  = stringtointeger(str_get_field(data_cuo,2," ")) ;从数组中获取温度值
  data_bo = asciiread(filepath+"surtmp56227.txt",-1,"string") ;以字符串形式读取参数文件入数组
  tmp_bo  = stringtointeger(str_get_field(data_bo,2," ")) ;从数组中获取温度值
  data_cha = asciiread(filepath+"surtmp56434.txt",-1,"string") ;以字符串形式读取参数文件入数组
  tmp_cha  = stringtointeger(str_get_field(data_cha,2," ")) ;从数组中获取温度值


;2.将变量整合到一个数组,并设置时间属性
  ;画线
  data      = new((/3,9/),integer)
  data(0,:) = tmp_cuo
  data(1,:) = tmp_cha
  data(2,:) = tmp_bo

  ;将数据分段,画点
  data_split  = new((/7,9/),integer,-999.)

  data_split(0,0:1)  = tmp_cuo(0:1) ;圈
  data_split(0,8)    = tmp_cuo(8) 
  data_split(1,2:7)  = tmp_cuo(2:7) ;花
  data_split(2,0)    = tmp_cha(0) ;圈
  data_split(3,1:8)  = tmp_cha(1:8) ;点
  data_split(4,0:1)  = tmp_bo(0:1) ;圈
  data_split(5,2)    = tmp_bo(2) ;点
  data_split(5,5:8)  = tmp_bo(5:8) 
  data_split(6,3:4)  = tmp_bo(3:4) ;花




  ;设置时间属性，从nc文件中获取
  nc_files = systemfunc("ls "+filepath+"fnl_*_*_*.nc")
  f = addfiles(nc_files,"r") ; f[0-8],9个文件
  ListSetType(f,"join"); 打开多个文件并且返回文件句柄(是列表),指定如何返回文件变量
  f_time =  f[:]->TMP_P0_L100_GLL0@initial_time; 注意原始数据中有些变量如 LAND_P0_L1_GLL0 有时间属性
  ;print(f_time) ; 是每个文件的具体时间,具有单位, initial_time : 10/27/2008 (00:00)
  ;exit
  tmpdata= f[:]->$"TMP_P0_L100_GLL0"$ 
  tmpdata!0 = "time"

  ;;;>>>时间格式的修改
  ; tmpdata&time = f_time
  time_set = grib_stime2itime(f_time) ; 将grib格式时间由字符串转化为整型
  time_set_yy = time_set/1000000
  time_set_mm = time_set/10000%100
  time_set_dd = time_set/100%100
  time_set_hh = time_set%100
  time_set_min = time_set_hh
  time_set_sec = time_set_hh
  time_set_min = 0
  time_set_sec = 0
  ; print(time_set_hh)
  ; print(time_set)
  ; tmp的时间格式重新写为 ju格式
  tmpdata&time = cd_inv_calendar(time_set_yy,time_set_mm,\
  time_set_dd,time_set_hh,time_set_min,time_set_sec, "hour since 2008-10-26 00:00:00",0)


 data!1="time"
 data&time=tmpdata&time



wks = gsn_open_wks("ps", picture_path+"micap_suf_tmp_label")
gsn_define_colormap(wks,"wgne15")
colors=(/14,8,4/) ;红 绿 蓝

 ;------------------------------绘图属性设置------------------------------------------
res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize=True
  res@tmXTOn = False ;关刻度
  res@tmYROn = False
  res@tiYAxisString="T (~S~o~N~C)"
  res@tmBorderThicknessF  =3   ; border 6x as thick,图框粗细
  res@tmXBMajorThicknessF =3
  res@tmXBMinorThicknessF =3
  res@tmYLMajorThicknessF =3
  res@tmYLMinorThicknessF =3  
  res@tmXBLabelAngleF = 48.   ;角度
  res@tmEqualizeXYSizes =True
  res@tmXBLabelFontHeightF = 0.015 ; 坐标轴标签字号

  res@gsnLeftStringFontHeightF=0.015
  res@gsnLeftString="(b)"
  res@gsnLeftStringOrthogonalPosF=-0.1
  res@gsnLeftStringParallelPosF =0.01
  res@txPerimOn=True
  res@txPerimThicknessF=2.5

  res@trYMaxF=14
  res@trYMinF=-18


  ;设置画图区域大小
  res@vpWidthF=0.5
  res@vpHeightF=0.4

  Tres=res
  Tres@xyMarkLineMode = "Markers"               ; Markers *and* lines
  Tres@xyMarkers =(/4,3,4,16,4,16,3/)
  Tres@xyMarkerThicknessF=5
  Tres@xyMarkerSizeF=10

  ;时间轴设置
  trestick = True
  trestick@ttmFormat = "%D/%H" ; 时间轴的格式
  trestick@ttmAxis   = "XB"
  trestick@ttmNumTicks = 9 ; 显示所有9个点



;绘点
time_axis_labels(data&time,Tres,trestick) ; call the formatting procedure 
plot  = gsn_csm_xy (wks,data&time,data_split,Tres) ; create plot

;绘线
  nline=3
  nrow=9
  dumLine = new((/nline,nrow/), graphic) 
  x= new(nrow, typeof(data&time))
  y= new(nrow, typeof(data&time))

  ;属性设置
  resLine = True
  resLine@gsLineThicknessF = 8
  resLine@tfPolyDrawOrder="PreDraw"

  ;开始画线
  x=data&time
  do j=0,nline-1
    y=data(j,:)
    if (j.eq.0) then
      resLine@gsLineColor = colors(0)
    end if
    if (j.eq.1) then
      resLine@gsLineColor = colors(2)
    end if
    if (j.eq.2) then
      resLine@gsLineColor = colors(1)
    end if 
  do i = 0, nrow-2
    xx = (/ x(i), x(i+1)/)
    yy = (/ y(i), y(i+1)/)
    dumLine(j,i) = gsn_add_polyline(wks, plot, xx, yy, resLine)
  end do
  delete(y)
  end do
 
; ; 绘制图例
;   lgcolours = (/"white","white","white","blue","red","purple"/)

;   resLg = True
;   resLg@lgMonoItemType=False
;   resLg@lgItemTypes = (/"Markers","Markers","Markers","Lines","Lines","Lines"/)
 
;   resLg@lgMarkerColors=(/"black","black","black","black","black","black"/)
;   resLg@lgMarkerIndexes=(/4,3,16,0,0,0/)
;   resLg@lgMarkerSizeF = 0.02  
;   resLg@lgMarkerThicknessF=5
;   resLg@lgMonoDashIndex=True ;一定要把这个打开！不然线型不能一致
;   ;resLg@lgLineDashIndexes=(/0,0,0,0,0,0/) ;这个好像打开也没用
;   resLg@lgLineColors = lgcolours
;   resLg@lgLineThicknessF = 8 

;   resLg@vpWidthF = 0.25
;   resLg@vpHeightF = 0.2
 
;   resLg@lgLabelFontHeightF = 0.15
 
;   lbid = gsn_create_legend(wks, 6, (/"No precipitation",\
;   "Snowfall","Rainfall","Cuona","Chayu","Bomi"/), resLg)
 
; ; 将图例放置在图中
;   amres = True
;   amres@amParallelPosF = 0.34 ;左右，正右左负
;   amres@amOrthogonalPosF = -0.35 ;上下，正下负上
;   dumLg = gsn_add_annotation(plot, lbid, amres)

draw(plot)
frame(wks)
end
