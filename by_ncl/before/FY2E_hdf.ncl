;---------------------------------------------------
; Plot FY2E TBB at 20130822_0000 (UTC)
;---------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

;---- Read FY2 satellite data TBB.
        f         = addfile("FY2E_TBB_IR1_NOM_20130822_0000.hdf","r")
        tbb = f->FY2E_TBB_Hourly_Product

        LVR = tbb@LowerValidRange ; The valid range of data is [LVR,UVR].
        UVR = tbb@UpperValidRange

;---- Read FY2 LatLon lookup table.
        f0                         = data_dir+"NOM_ITG_2288_2288(0E0N)_LE.dat"
        dims                       = (/2288,2288/)
        type                       = "float"
        ; Band 1 is lon: data1
        data1                      = fbindirread(f0,0,dims,type)
        ; Band 2 is lat: data2
        data2                      = fbindirread(f0,1,dims,type)
        
        lonCenter                 = 104.5                ; footpoint of FY-2D
        latCenter                 = 0
        lon2d                     = data1+lonCenter
        lat2d                     = data2+latCenter
        lat2d@units               = "degrees_north"
        lon2d@units               = "degrees_east"

;---- Add meta data for TBB.
        tbb@lat2d                 = lat2d
        tbb@lon2d                 = lon2d
        tbb@units                 = "degree Kelvin"
        tbb@long_name         = "Temperature of Bright Blackbody"
        tbb@coordinates = "lat2d lon2d"
        tbb@_FillValue  = -999.

;---- Replace data outside the valid range with missing value.
        tbb_1d                   = ndtooned(tbb) ;二维变一维
        ind_notValid             = ind(tbb_1d.lt.LVR.or.tbb_1d.gt.UVR)
        tbb_1d(ind_notValid)     = tbb@_FillValue
        tbb                      = onedtond(tbb_1d,dimsizes(tbb))

        tbb = tbb-273.15        ;(convert degree Kelvin to degree Celsius)
        tbb@units = "degree Celsius"
        printMinMax(tbb,0)
        printVarSummary(tbb)



;---- Begin to plot.
        wks = gsn_open_wks("png","TBB_plot")
        res = True

        res@gsnDraw                                 = False
        res@gsnFrame                                = False
        res@gsnMaximize                                = True
        res@gsnPaperOrientation                = "portrait"
        res@gsnAddCyclic                         = False         ; regional data
        ;tfDoNDCOverlay should be set False (default) when lat2d and lon2d are used.
        ;res@tfDoNDCOverlay                        = True                        
        res@pmTickMarkDisplayMode          = "Always"

        ;---- Map Set
        res@mpDataSetName                                = "Earth..4"
        res@mpDataBaseVersion                        = "MediumRes"
        res@mpOutlineSpecifiers                        = (/"China","China:Provinces"/)
        res@mpOutlineBoundarySets           = "NoBoundaries"
        res@mpNationalLineColor                        = "black"
        res@mpProvincialLineColor                = "black"
        res@mpGeophysicalLineColor          = "black"
        res@mpNationalLineThicknessF         = 3
        res@mpProvincialLineThicknessF  = 3
        res@mpGeophysicalLineThicknessF = 3

        res@mpMinLonF                                = 105
        res@mpMaxLonF                                = 128
        res@mpMinLatF                                = 18
        res@mpMaxLatF                                = 35

        res@trGridType                                = "TriangularMesh"
        res@cnFillOn                                 = True
        res@cnFillMode                                 = "RasterFill"
        res@cnLinesOn                                 = False
        res@lbOrientation        = "vertical"            ; 垂直坐标
        res@pmLabelBarOrthogonalPosF  = 0.08
        res@lbBoxEndCapStyle  = "TriangleBothEnds"
        res@pmLabelBarWidthF = 0.05     ;色标变窄


        res@gsnLeftString                        = "TBB"
        res@gsnRightString                        = "~S~o~N~C"

        ;---- Color set
        cmap  = read_colormap_file("precip3_16lev")
        cmap1 = cmap(4::,:)
        cmap1(0,:) = cmap(0,:)
        cmap1 = cmap1(::-1,:)
        res@cnFillPalette = cmap1
        
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels = (/-30,-40,-50,-60,-70,-80/)
        res@lbOrientation = "vertical"

        plot = gsn_csm_contour_map_ce(wks,tbb,res)
        
        draw(plot)
        frame(wks)

end 
