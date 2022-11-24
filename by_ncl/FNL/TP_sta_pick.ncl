;mask掉高原以外站点，输出高原站号，写入文件
;文件读取方式参考 http://www.ncl.ucar.edu/Applications/list_io.shtml
;写入文件参考 http://www.ncl.ucar.edu/Applications/write_ascii.shtml 
;代码路径  ncl /mnt/h/d/ncl_related/code/TP_sta_pick.ncl


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;--------------------读站点数据--------
  filepath = "/mnt/h/d/data/precipitation/6h/"  ;参数文件路径
  name    ="rain12_py"
  f1       = asciiread(filepath+name+".txt",-1,"string")
  sta_all = stringtointeger(str_get_field(f1(1:),1," "))  ;读取站号
  lat     = stringtofloat(str_get_field(f1(1:),2," ")) ;获取纬度值lon
  lon     = stringtofloat(str_get_field(f1(1:),3," ")) ;获取经度值lat
; print(sta_all)

 ;设置缺测值，mask掉高原以外的站点
  sta_all@_FillValue = 9999
;   ;注意：这里的txt文件最后一行数据下面需要空白行，才能读出文件的正确行数
;   ;大多数txt文件都满足要求，但还是需要核实一下
  npts = numAsciiRow(filepath+name+".txt")-1 ;获取行数-站点数

  ;---Open shapefile and read Mississippi River Basin lat/lon values.
  ;---注意这里不是画轮廓的DBATP_Line.shp,但是这两个画出来又好像是一样的
  f2      = addfile("/mnt/h/d/data/TPshp/DBATP/DBATP_Polygon.shp", "r")
  mrb_lon = f2->x
  mrb_lat = f2->y
  ;---Start with data all missing
  sta = new(dimsizes(sta_all),typeof(sta_all),sta_all@_FillValue)
  do i=0,npts-1
    sta(i)=where(gc_inout(lat(i),lon(i),mrb_lat,mrb_lon), sta_all(i), sta_all@_FillValue)
  end do

;输出高原站数和站号
  nmsg = num(.not.ismissing(sta))   ; count number of missing
  ; TP_sta=new(nmsg,integer)
  ista = ind(.not.ismissing(sta))
  tp_sta = sta(ista)
  print(tp_sta)

  tp_sta_path="/mnt/h/d/data/precipitation/"
  asciiwrite(tp_sta_path+"TP_sta.txt",tp_sta)   

end     
 

