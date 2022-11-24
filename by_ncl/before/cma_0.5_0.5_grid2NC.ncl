load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
; 把中国地面降水日值0.5°×0.5°格点数据txt转换为nc
; 参考函数 http://www.ncl.ucar.edu/Document/Functions/Built-in/system.shtml
; 参考函数 http://www.ncl.ucar.edu/Document/Functions/Built-in/filedimdef.shtml
; 参考函数 http://www.ncl.ucar.edu/Document/Functions/Contributed/readAsciiTable.shtml

begin
;-------------------转换数据----------------------------------------------------
  tStart = 20081025
  tLast  = 20081028  
  tJump  = 1
 
  diri   = "/cygdrive/D/ncl_related/data/cma_0.5_0.5_rain/txt/"
  diro   = "/cygdrive/D/ncl_related/data/cma_0.5_0.5_rain/nc/"  
  fRoot  = "SURF_CLI_CHN_PRE_DAY_GRID_0.5_" 
  ncols      = 128
  nrows      = 72
  xllcorner  = 72.0  ;lon
  yllcorner  = 18.0  ;lat
  lat        = fspan (yllcorner,yllcorner+(nrows-1)*0.5,nrows)
  lon        = fspan (xllcorner,xllcorner+(ncols-1)*0.5,ncols)
  lat@units      = "degrees_north"
  lon@units      = "degrees-east"
  ;dim_names = (/"lat","lon"/)
  ;dim_sizes = (/nrows,ncols/)
  ;dimUnlim = (/False,False/)
  do time = tStart, tLast, tJump
     fili = fRoot + time + ".txt"    ; name is of form SURF_CLI_CHN_PRE_DAY_GRID_0.5_20081025.txt and similar
     filo = fRoot + time + ".nc"
     system("/bin/rm -f diro+filo") ; remove the nc files
     ncdf = addfile(diro+filo,"c") ;creat nc file
     filedimdef(ncdf,"time",-1,True)
     ;filedimdef( ncdf, dim_names, dim_sizes, dimUnlim )
     ip       = readAsciiTable(diri+fili,ncols,"float",6)
     op       = ip(::-1,:)                       
     op@_FillValue	= -9999.000
     op!0    	= "lat"                      
     op!1   	= "lon"
     op&lat  	= lat                    
     op&lon  	= lon
     op@long_name	= "daily precipitation"
     op@units    	= "mm"
     ncdf->op  	= op
     delete(ncdf)
     delete(ip) ;一定要记得删掉之前的变量
     delete(op)
  end do 
 end                                
                                 
                                 

