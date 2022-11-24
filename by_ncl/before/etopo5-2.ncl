load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;����shpae�ļ�ֱ�ӻ�����ظ�ԭ
begin
  datadir ="/home/acer/TN/eps/"      ;�ļ�����·��
  datadir1 ="/home/acer/data/shp/"      ;�ļ�����·��
  datadir2 ="/home/acer/data/daily/"
   setfileoption("nc","MissingToFillValue",False)    ;����ȱʡֵ
  wks = gsn_open_wks("eps",datadir+"topo2") 
;  elev = read_elev_data(datadir1+"ETOPO5.DAT")
;  printVarSummary(elev)
  in = addfile(datadir2+"hgt.2011.nc","r")    ;��ȡnc�ļ�
  hgt1=in->hgt(0,0,:,:)
   printVarSummary(hgt1)
;  elev1 = where(elev.lt.3000.,elev@_FillValue,elev)
;  copy_VarCoords(elev,elev1)     
;  printVarSummary(elev1)

  res                    = True
  
  res@mpMinLatF            = 0      ; ���û�ͼ�ľ�γ��Χ
  res@mpMaxLatF            = 60.
  res@mpMinLonF            = 60.
  res@mpMaxLonF            = 140.
  
  res@gsnDraw        = False             ; don't draw the plots yet
  res@gsnFrame       = False        ; don't advance the frame
  res@gsnMaximize        = True             ; maximize plot in frame
  res@cnFillOn           = True             ; turn on contour fill
  res@gsnAddCyclic = True
; res@cnLinesOn          = False            ; turn off contour lines
  res@cnLineLabelsOn     = False            ; turn off line labels
  res@cnInfoLabelOn      = False            ; turn off info label
  res@lbBoxLinesOn       = False            ; turn off labelbar box lines
  res@gsnAddCyclic       = True           ; don't add longitude cyclic point
  res@mpFillOn           = False            ; turn off map fill
  res@tiMainString       = "ETOPO5.DAT"     ; main title
  res@pmLabelBarWidthF   = 0.8              ; default is too short
   
  
;  res2=True
;  res2@gsnDraw        = False             ; don't draw the plots yet
;  res2@gsnFrame       = False        ; don't advance the frame
;  res2@cnLinesOn = False

;  res2@cnFillOn           = True             ; turn on contour fill
;  res2@cnFillColors   =(/"black","black","black","black","black","black","black","black","black","black"/)
;    res2@cnLinesOn            = False             ; turn off contour lines
;  res2@gsnAddCyclic       = False            ; don't add longitude cyclic point
;  res2@lbLabelBarOn    = False                  ; will draw a panel label bar instead
;  res2@cnLineLabelsOn       =  False      ; Turn off contour lines
  plot = gsn_csm_contour_map(wks,hgt1,res)  ;map��ͼֻ����Ϊ����
 
  lnres=True
  lnres@gsLineColor      = "white"
  lnres@gsLineThicknessF = 6.0
  map = gsn_add_shapefile_polylines(wks,plot,datadir1+"Qinghai-Tibet_Plateau.shp",lnres)

;  plot_hgt=gsn_csm_contour(wks,map,lnres)
;  overlay(plot,map)
  draw(plot)
  frame(wks)
end