load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"      ;可以转化lon为-180-180
;----------------------------------------------------------------------
; This function reads a binary file containing elevation data and
; generates the necessary lat/lon coordinate arrays for plotting later.
; The information on the binary file is provided at the beginning of
; this script.
;
; The binary file was downloaded from:
;     http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO5/TOPO/ETOPO5/
;----------------------------------------------------------------------
undef("read_elev_data")
function read_elev_data(topo_file)
local nlat, nlon, topo_file, lat, lon
begin
;---Read data as a straight binary file
  nlat = 2160
  nlon = 4320
  setfileoption("bin","ReadByteOrder","BigEndian")
  elev = cbinread(topo_file,(/nlat,nlon/),"short")

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

;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
begin
 lon       = fspan(0,360,4320)  ;gsn_xy
  datadir ="/home/acer/TN/eps/"      ;文件所在路径
  datadir1 ="/home/acer/data/"      ;文件所在路径
  datadir2 ="/home/acer/data/daily/"
   setfileoption("nc","MissingToFillValue",False)    ;设置缺省值
  wks = gsn_open_wks("eps",datadir+"topo4") 
  stdlat=50;选择固定某一角度
  
  elev= read_elev_data(datadir1+"ETOPO5.DAT")
  printVarSummary(elev)
  in = addfile(datadir2+"uwnd.2011.nc","r")    ;读取nc文件
;  uwnd1=lonFlip(in->uwnd(:,:,:,:)) ;将经度范围转化为-180到180
  uwnd1=in->uwnd(:,:,:,:)
  uave=dim_avg_n_Wrap(uwnd1(206:211,:,{stdlat},:),0)
  printVarSummary(uave)
    
  elev0 = 1013.25*(1-elev*0.0065/288.15)^5.25145        ;将地形数据的单位：m转化为hPa
  copy_VarCoords(elev,elev0) 
  printVarSummary(elev0)

;  elev1=elev0({stdlat},:)       ;固定某一纬度的地形数据
;  printVarSummary(elev1)
;  elev1@_FillValue= -9.96921e+36
;  elev2=where(elev1.gt.1000.and.elev1.lt.2000,elev1@_FillValue,elev1)
;  copy_VarCoords(elev1,elev2)
;  printVarSummary(elev2)

 ;lat  = fspan(90,-90,73)
 ;lon  = fspan(0,357.5,144)
 ;elev1 = area_hi2lores_Wrap(elev0&lon,elev0&lat,elev0,True,1,lon,lat,False) 
  elev2=elev0({stdlat},:)       ;固定某一纬度的地形数据
    printVarSummary(elev2)
  elev2@_FillValue= -9.96921e+36
  elev3=where(elev2.lt.1000,elev2@_FillValue,elev2)
  copy_VarCoords(elev2,elev3)


;  printVarSummary(elev3)
;  uave1=uave(1,:) 
;   elev4=conform_dims(dimsizes(uave),elev3,1)  
; copy_VarCoords(uave,elev4)
 printVarSummary(elev3)    
  
;开始作图

  res                    = True 

 res@gsnDraw        = False             ; don't draw the plots yet
  res@gsnFrame       = False        ; don't advance the frame
  res@gsnMaximize        = True             ; maximize plot in frame
  res@cnFillOn           = True             ; turn on contour fill
  res@cnFillDrawOrder      = "Predraw"  ;调整画图的顺序
  res2=True
   ; res2@gsnBelowYRefLineColor = "black" 
   ;res2@gsnXYBelowFillColors="black"
    res2@gsnMaximize      = True
   res2@trYMaxF=1000
  res2@gsnYRefLine      = 1000.          ; create a reference line 
  res2@gsnBelowYRefLineColor = "black"              ; above ref line fill red
    res2@gsnDraw        = False             ; don't draw the plots yet
  res2@gsnFrame       = False        ; don't advance the frame  
;  res2@gsnBelowYRefLineColor = "black"             ; below ref line fill blue
 ; res2@cnMissingValFillColor="black" 
 
  ;res2@cnLinesOn = False

;  res2@cnFillOn           = True             ; turn on contour fill
;  res2@cnFillColors   =(/"black","black","black","black","black","black","black","black","black","black"/)
;  res2@cnLinesOn            = False             ; turn off contour lines
;  res2@cnLineLabelsOn       =  False      ; Turn off contour lines





  plot = gsn_csm_pres_hgt(wks,uave,res)  ;作为基地  
  plot_hgt = gsn_csm_xy(wks,lon,elev2,res2)  
  ;plot_hgt=gsn_contour(wks,aa,res2) 
   ;plot = gsn_csm_contour(wks,aa,res)
  overlay(plot,plot_hgt)
  draw(plot)
  frame(wks)
end