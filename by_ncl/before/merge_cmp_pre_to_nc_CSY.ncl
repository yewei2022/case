;　将ＣＭＯＲＰＨ融合降水产品二进制文件转换为ＮＣ文件
; 关于CMORPH_PRE文件的说明
; 最西经度：７０，最东经度：１４０
; 最南纬度：１５，最北纬度：５９
; 空间分辨率：0.1　度
; 投影方式：等经纬度
; 特征值：－９９９表示缺测

;------ｇｒａｄｓ参考ｃｔｌ文件------
; DSET ^SEVP_CLI_CHN_MERGE_FY2_PRE_HOUR_GRID_0.10-%y4%m2%d2%h2.grd
; *
; UNDEF -999.0
; *
; OPTIONS   little_endian  template
; *
; TITLE  China Hourly Merged Precipitation Analysis
; *
; xdef  700 linear  70.05  0.10
; *
; ydef  440 linear  15.05  0.10 
; *
; ZDEF     1 LEVELS 1  
; *
; TDEF 9999 LINEAR 00Z01Aug2010 1hr 
; *
; VARS 2                           
; crain      1 00  CH01   combined analysis (mm/Hour)
; gsamp      1 00  CH02   gauge numbers
; ENDVARS

;------读取降水文件并添加属性------
undef("sy_define_read_cmorph")
function sy_define_read_cmorph(file_path)
local file_path
begin
    pre = new((/440,700/), "float", 0)
    pre = fbindirread(file_path, 0, (/440,700/), "float")

    ;设置维度信息
    pre!0               =   "lat"
    pre&lat             =   fspan(15.05, 58.95, 440)
    pre&lat@long_name   =   "latitude"
    pre&lat@units       =   "degrees_north"

    pre!1               =   "lon"
    pre&lon             =   fspan(70.05, 139.95, 700)
    pre&lon@long_name   =   "longitude"
    pre&lon@units       =   "degrees_east"

    pre@long_name       =   "China Hourly Merged Precipitation Analysis"
    pre@units           =   "mm/Hour"
    pre@_FillValue      =   -999

    return pre
end 



begin   
    ; ==============================================================
    ;   用户设置  
    ; ==============================================================
    diri        =   "../../DATA/merge_cmp_pre/"
    file_name   =   "*.grd"   
    N   =   440
    M   =   700

    tinvl   =   1   ; 时间间隔

    ; ==============================================================
    ;   读入数据
    ; ==============================================================
    files   =   systemfunc("ls " + diri + file_name)
    ntime   =   dimsizes(files)
    
    pre     =   new((/ntime, N, M/), "float", -999)
    do i = 0, ntime - 1
        pre(i,:,:) = sy_define_read_cmorph(files(i))
    end do

    ; ***********************************************
    ;   添加时间维度属性
    ; ***********************************************
    pre!0       =   "time"
    pre&time    =   ispan(0, (ntime - 1) * tinvl, tinvl)
    pre&time@units  =   "hours since 2017-08-09 01:00:0.0"

    ; ==============================================================
    ;   写入数据
    ; ==============================================================
    foName  =   "MERGE_CMP_PRE_HOUR.nc"
    system("rm -rf " + diri + foName)
    fo  =   addfile(diri + foName, "c")
    fo->pre = pre
end 