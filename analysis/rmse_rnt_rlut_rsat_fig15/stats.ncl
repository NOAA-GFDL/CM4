

NUMSTAT = 5
DTR = acos(0.0)/90.

function compute_statistics ( x[*][*][*]:numeric, o[*][*][*]:numeric, w[*][*]:numeric )
local dims, stats, avgo
begin
  dims = dimsizes(x)
  stats = new((/dims(0),NUMSTAT/),typeof(x))
  stats(:,0) = wgt_areaave2(x,w,0) ; average
  avgo = wgt_areaave2(o,w,0)
  stats(:,1) = stats(:,0) - avgo      ; bias
  stats(:,2) = wgt_arearmse2(x,o,w,0) ; rms
  stats(:,3) = pattern_cor(x,o,w,0)   ; correlation (centered)
  stats(:,4) = pattern_cor(x,o,w,1)   ; correlation (uncentered)
  return stats
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
function quartile_values ( x[*][*]:numeric )
local dimx,frac,np,stat,xs,is,fstat,istat,dstat
begin
  dimx = dimsizes(x)
  frac = (/ 0.25, 0.50, 0.75 /)
  np = dimx(0)
  dimx(0) = 5
  stat = new(dimx,typeof(x))

  xs = x
  is = dim_pqsort_n(xs,2,0)
  stat(0,:) = xs(0,:)       ; minimum
  stat(4,:) = xs(np-1,:) ; maximum
  do i = 1, dimsizes(frac)
    fstat = frac(i-1)*tofloat(np)
    istat = toint(fstat)
    dstat = fstat-istat
    stat(i,:) = xs(istat-1,:)*(1.0-dstat) + xs(istat,:)*dstat
  end do

  return stat
end

function get_model_name (fpath[1]:string)
local parse,model
begin
 ;parse = str_split(str_get_cols(fpath,str_index_of_substr(fpath,"/",-1)+1,-1),"_")
  parse = str_split(fpath,"_")
  model = parse(2)
  return model
end

function isAM4model(fpath[1]:string)
local parse
begin
  parse = str_split(str_get_cols(fpath,str_index_of_substr(fpath,"/",-1)+1,-1),"_")
  ; return True if model=AM4 or exper=amip
  if (parse(2) .eq. "GFDL-AM4" .or. parse(3) .eq. "amip") then
    return True
  else
    return False
  end if
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

procedure append_array_value(x[1],A[1]:logical,att[1]:string,append[1]:logical)
local a1
begin
  if (append) then
    a1 = array_append_record(A@$att$,x,0)
  else
    a1 = array_append_record(x,A@$att$,0)
  end if
  delete(A@$att$)
  A@$att$ = a1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; input: list of filenames for runs from the same model
; output: first run member (usually r1i1p1)

function find_first_run (files[*]:string)
local rip,run,n,parse,r,i,p,indi,indp
begin
  run = new(dimsizes(files),integer)   ; r*
  do n = 0, dimsizes(files)-1
    parse = str_split(files(n),"_")
    rip = parse(4)
    if (str_get_cols(rip,0,0) .ne. "r") then
      print("ERROR: invalid rip: "+rip+", file = "+files(n))
      status_exit(1)
    end if
    if (rip .eq. "r1i1p1") then
      return files(n)
    end if
    indi = str_index_of_substr(rip,"i",1)
    indp = str_index_of_substr(rip,"p",1)
    r = toint(str_get_cols(rip,1,indi-1))
    i = toint(str_get_cols(rip,indi+1,indp-1))
    p = toint(str_get_cols(rip,indp+1,-1))
    run(n) = r*10000 + i*100 + p
  end do
  return files(minind(run))
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

procedure mask_data_array( data:numeric, lmask[*][*]:logical, fill[1]:numeric )
local rank
begin
  rank = dimsizes(dimsizes(data))
  if (rank .eq. 2) then
    data = where(lmask, data, fill)
  else
    data = where(conform(data,lmask,(/rank-2,rank-1/)), data, fill)
  end if
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

load "/home/Huan.Guo/awg/xanadu/utils/CM4_paper_plots_BW/graphics/define_regions.ncl"
load "/home/Huan.Guo/awg/xanadu/utils/CM4_paper_plots_BW/graphics/box_whisker.ncl"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

begin

  if (.not.isvar("regionFile")) then
    regionFile = "$BW_PACKAGE_ROOT/regions.txt"
  end if

  if (.not.isvar("var")) then
    var = "pr"
  end if
  if (.not.isvar("plotSSN")) then
    plotSSN = True
  end if
  if (.not.isvar("averageOn")) then
    averageOn = False
  end if
  if (.not.isvar("maskReg")) then
    maskReg = "global"
  end if
  if (.not.isvar("removePslAvgOn")) then
    removePslAvgOn = False
  end if

  cpath = var+"/cmip5/obsgrid";
  mpath = var+"/cm4/obsgrid";
  opath = var+"/obs";

  if (plotSSN) then
    cfiles = systemfunc("/bin/ls "+cpath+"/*-ssnclim.nc")
    mfiles = systemfunc("/bin/ls "+mpath+"/*-ssnclim.nc")
    ofile = systemfunc("/bin/ls "+opath+"/*-ssnclim.nc")
    numTime = 5
  else
    cfiles = systemfunc("/bin/ls "+cpath+"/*-clim.nc")
    mfiles = systemfunc("/bin/ls "+mpath+"/*-clim.nc")
    ofile = systemfunc("/bin/ls "+opath+"/*-clim.nc")
    numTime = 14 ; 12 months + 2 cyclic points
  end if

  if (maskReg .eq. "land" .or. maskReg .eq. "ocean") then
    maskfile = systemfunc("/bin/ls "+mpath+"/sftlf*.nc")
    if (dimsizes(maskfile) .gt. 1) then
      print("ERROR: more than one landsea mask file found")
      status_exit(1)
    end if
    fm = addfile(maskfile,"r")
    sftlf = fm->sftlf
    delete(fm)
  else
    if (maskReg .ne. "global") then
      print("WARNING: invalid masking region, defaulting to global")
    end if
  end if

  if (var .eq. "psl") then
    maskfile = systemfunc("/bin/ls "+opath+"/orog*.nc")
    if (dimsizes(maskfile) .gt. 1) then
      print("ERROR: more than one orography file found")
      status_exit(1)
    end if
    fm = addfile(maskfile,"r")
    orog = fm->orog
    maskPSL = True
    print("Masking PSL")
  else
    maskPSL = False
  end if

  ; set up geographic regions
  REGIONS = define_regions(regionFile)
  if (REGIONS@numReg .eq. 0) then
    print("ERROR: number of regions is zero")
    status_exit(1)
  end if
  numReg = REGIONS@numReg

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; read the observed data
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  fo = addfile(ofile,"r")
  odata = fo->$var$
  area = conform(odata(0,:,:),cos(odata&$odata!1$*DTR),0)
  copy_VarCoords(odata(0,:,:),area)

  ; remove global mean from PSL
  if (var .eq. "psl" .and. removePslAvgOn) then
    odata = odata - conform(odata,wgt_areaave2(odata,area,0),0)
  end if

  ; average of observed data
  obsavg = new(numTime,typeof(odata))
  if (plotSSN) then
    obsavg = wgt_areaave2(odata,area,0)
  else
    obsavg(1:numTime-2) = wgt_areaave2(odata,area,0)
    obsavg(0) = obsavg(numTime-2) ; cyclic point
    obsavg(numTime-1) = obsavg(1) ; cyclic point
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; find all model runs
  ; group them by model name
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  model = new(dimsizes(cfiles),string)
  model_runs = True
  do m = 0, dimsizes(cfiles)-1
    fc = addfile(cfiles(m),"r")
    cfilename = str_get_cols(cfiles(m),str_index_of_substr(cfiles(m),"/",-1)+1,-1)
    model(m) = get_model_name(cfilename)
    ; count number of runs for each model
    if (isatt(model_runs,model(m))) then
      fnames = array_append_record(model_runs@$model(m)$,cfilename,0)
      delete(model_runs@$model(m)$)
      model_runs@$model(m)$ = fnames   ;model_runs@$model(m)$ + 1
      delete(fnames)
    else
      model_runs@$model(m)$ = cfilename
      ;model_runs@$model(m)$ = 1
    end if
    delete(fc)
  end do

  oneRunPerModel = True
  numRuns = dimsizes(cfiles)
  if (oneRunPerModel) then
    model_names = getvaratts(model_runs)
    numRuns = dimsizes(model_names)
  end if
  dims = (/ numRuns, numReg, numTime, NUMSTAT /)
  stats = new(dims,float)
  delete(dims)

  do m = 0, numRuns-1  ;dimsizes(cfiles)-1
    if (oneRunPerModel) then
      cfile = cpath + "/" + find_first_run(model_runs@$model_names(m)$)
      print("model: "+model_names(m)+", file: "+cfile)
    else
      cfile = cfiles(m)
    end if
    fc = addfile(cfile,"r")
    cdata = fc->$var$

    ; remove global mean from PSL
    if (var .eq. "psl" .and. removePslAvgOn) then
      cdata = cdata - conform(cdata,wgt_areaave2(cdata,area,0),0)
    end if

   ; scale data if necessary
    if (var .eq. "pr") then
      cdata = cdata * 86400.
    end if
    delete(fc)

    ; compute statistics for each region
    do reg = 0, numReg-1
      ; define the geographic region (lon/lat)
      regName  = REGIONS@id(reg)
      regTitle = REGIONS@title(reg)
      minLon   = REGIONS@lonbeg(reg)
      maxLon   = REGIONS@lonend(reg)
      minLat   = REGIONS@latbeg(reg)
      maxLat   = REGIONS@latend(reg)

      ; shift data longitude for this region
      cdata1 = my_lonPivot(cdata,minLon)
      odata1 = my_lonPivot(odata,minLon)
      area1  = my_lonPivot(area,minLon)
      if (maskReg .eq. "land") then
        mask1  = where(my_lonPivot(sftlf,minLon) .gt. 0.99, True, False)
        mask_data_array( cdata1, mask1, cdata1@_FillValue )
        mask_data_array( odata1, mask1, odata1@_FillValue )
        mask_data_array(  area1, mask1, 0.0 )
      else if (maskReg .eq. "ocean") then
        mask1  = where(my_lonPivot(sftlf,minLon) .lt. 0.01, True, False)
        mask_data_array( cdata1, mask1, cdata1@_FillValue )
        mask_data_array( odata1, mask1, odata1@_FillValue )
        mask_data_array(  area1, mask1, 0.0 )
        mask1  = my_lonPivot(sftlf,minLon)
      else if (maskPSL) then
        mask1  = where(my_lonPivot(orog,minLon) .lt. 500., True, False)
        mask_data_array( cdata1, mask1, cdata1@_FillValue )
        mask_data_array( odata1, mask1, odata1@_FillValue )
        mask_data_array(  area1, mask1, 0.0 )
      end if
      end if
      end if

      if (plotSSN) then
        stats(m,reg,:,:) = compute_statistics( cdata1(:,{minLat:maxLat},{minLon:maxLon}), \
                                               odata1(:,{minLat:maxLat},{minLon:maxLon}), \
                                               area1({minLat:maxLat},{minLon:maxLon}) )
      else
        stats(m,reg,1:numTime-2,:) = compute_statistics( cdata1(:,{minLat:maxLat},{minLon:maxLon}), \
                                                         odata1(:,{minLat:maxLat},{minLon:maxLon}), \
                                                         area1({minLat:maxLat},{minLon:maxLon}) )
      end if
      delete([/cdata1,odata1,area1/])
      if (isvar("mask1")) then
        delete(mask1)
      end if
    end do
    delete(cdata)
  end do

  if (.not.plotSSN) then
    stats(:,:,0,:) = stats(:,:,numTime-2,:) ; cyclic point
    stats(:,:,numTime-1,:) = stats(:,:,1,:) ; cyclic point
  end if

  ; setup blue markers for "CM3" model
  if (oneRunPerModel) then
    blue = ind(model_names .eq. "GFDL-CM3")
    purple = ind(model_names .eq. "GFDL-CM2p1")
   ;green = array_append_record(ind(model_names .eq. "GFDL-ESM2G"), \
   ;                            ind(model_names .eq. "GFDL-ESM2M"),0)
  else
    blue = ind(model .eq. "GFDL-CM3")
    purple = ind(model .eq. "GFDL-CM2p1")
   ;green = array_append_record(ind(model .eq. "GFDL-ESM2G"), \
   ;                            ind(model .eq. "GFDL-ESM2M"),0)
  end if
  
  ;;;;;;;;;;;;;;;;;
  ;;  CM4 model  ;;
  ;;;;;;;;;;;;;;;;;
  dims = (/ dimsizes(mfiles), numReg, numTime, NUMSTAT /)
  mstats = new(dims,float)

  do m = 0, dimsizes(mfiles)-1
    fm = addfile(mfiles(m),"r")
    mdata = fm->$var$

    ; remove global mean from PSL
    if (var .eq. "psl" .and. removePslAvgOn) then
      mdata = mdata - conform(mdata,wgt_areaave2(mdata,area,0),0)
    end if

    ; scale precip
    if (var .eq. "pr") then
      mdata = mdata * 86400.
    end if

    if (m .eq. 0 .and. isfilevaratt(fm,var,"long_name")) then
      plotTitle = fm->$var$@long_name
      print( " plotTitle " + plotTitle ) 
    end if
    delete(fm)

    ; compute statistics for each region
    do reg = 0, numReg-1
      ; define the geographic region (lon/lat)
      regName  = REGIONS@id(reg)
      regTitle = REGIONS@title(reg)
      minLon   = REGIONS@lonbeg(reg)
      maxLon   = REGIONS@lonend(reg)
      minLat   = REGIONS@latbeg(reg)
      maxLat   = REGIONS@latend(reg)

      ; shift data longitude for this region
      mdata1 = my_lonPivot(mdata,minLon)
      odata1 = my_lonPivot(odata,minLon)
      area1  = my_lonPivot(area,minLon)

      if (maskReg .eq. "land") then
        mask1  = where(my_lonPivot(sftlf,minLon) .gt. 0.99, True, False)
        mask_data_array( mdata1, mask1, mdata1@_FillValue )
        mask_data_array( odata1, mask1, odata1@_FillValue )
        mask_data_array(  area1, mask1, 0.0 )
      else if (maskReg .eq. "ocean") then
        mask1  = where(my_lonPivot(sftlf,minLon) .lt. 0.01, True, False)
        mask_data_array( mdata1, mask1, mdata1@_FillValue )
        mask_data_array( odata1, mask1, odata1@_FillValue )
        mask_data_array(  area1, mask1, 0.0 )
      else if (maskPSL) then
        mask1  = where(my_lonPivot(orog,minLon) .lt. 500., True, False)
        mask_data_array( mdata1, mask1, mdata1@_FillValue )
        mask_data_array( odata1, mask1, odata1@_FillValue )
        mask_data_array(  area1, mask1, 0.0 )
      end if
      end if
      end if

      if (plotSSN) then
        mstats(m,reg,:,:) = compute_statistics( mdata1(:,{minLat:maxLat},{minLon:maxLon}), \
                                                odata1(:,{minLat:maxLat},{minLon:maxLon}), \
                                                area1({minLat:maxLat},{minLon:maxLon}) )
      else
        mstats(m,reg,1:numTime-2,:) = compute_statistics( mdata1(:,{minLat:maxLat},{minLon:maxLon}), \
                                                          odata1(:,{minLat:maxLat},{minLon:maxLon}), \
                                                          area1({minLat:maxLat},{minLon:maxLon}) )
      end if
      delete([/mdata1,odata1,area1/])
      if (isvar("mask1")) then
        delete(mask1)
      end if
    end do
  end do

  if (.not.plotSSN) then
    mstats(:,:,0,:) = mstats(:,:,numTime-2,:) ; cyclic point
    mstats(:,:,numTime-1,:) = mstats(:,:,1,:) ; cyclic point
  end if


  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;  plotting section  ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; create informational label
  if (oneRunPerModel) then
    infolabel = "Number of CMIP5 models = "+numRuns
  else
    atts = getvaratts(model_runs)
    numModels = dimsizes(atts)
    infolabel = "Number of CMIP5 models = "+numModels+", Total number of runs = "+numRuns
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; box-whisker seasonal plots ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (plotSSN) then
    resbw = True
    if (isvar("plotTitle")) then
      resbw@plotTitle = plotTitle
    end if
    if (.not.ismissing(blue(0))) then
      resbw@blue = blue
    end if
    if (.not.ismissing(purple(0))) then
      resbw@purple = purple
    end if
    resbw@var = var
    resbw@info = infolabel
    resbw@regID = REGIONS@id
    resbw@regTitle = REGIONS@title
    resbw@averageOn = averageOn
    resbw@maskedRegionName = maskReg

    ; set the index for the AM4 model
    am4flag = new(dimsizes(mfiles),logical)
    am4flag = False
    do i = 0, dimsizes(mfiles)-1
      if (isAM4model(mfiles(i))) then
         am4flag(i) = True
      end if
    end do
    if (num(am4flag) .eq. 1) then
      resbw@am4 = ind(am4flag)
    end if

    print( " before box_whisker " +  REGIONS@title )
    box_whisker_ssn( stats, mstats, resbw)
    exit
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;; monthly climatology plots ;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  xstats = conform(stats(:,0,:,0),fspan(0,13,14),1)
  yaxis_labels = (/ "Average", "Bias", "RMSE", "Corr (centered)", "Corr (uncentered)" /)
  shade_colors = (/ "grey85", "grey65" /)

  istat0 = 1
  if (averageOn) then
    istat0 = 0
  end if
  
  res = True
  res@gsnFrame = False
  res@gsnDraw = False
  res@xyMonoLineColor = True
  res@xyLineThicknessF = 1
  res@xyLineColor = shade_colors(0)
  ; define x-axis for monthly means
  res@tmXBMode = "Explicit"
  res@tmXBValues = (/1,2,3,4,5,6,7,8,9,10,11,12/)
  res@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
  res@trXMinF  =   0.5 
  res@trXMaxF  =  12.5
  if (isvar("plotTitle")) then
   ;res@tiMainString = plotTitle
    res@gsnLeftString = plotTitle
  end if

  respl = True
  if (dimsizes(mfiles) .eq. 1) then
    respl@gsLineThicknessF = 3.0
  else if (dimsizes(mfiles) .eq. 2) then
    respl@gsLineThicknessF = 2.0
  else
    respl@gsLineThicknessF = 1.0
  end if
  end if


  resplo = True
  resplo@gsLineColor = "black"
  resplo@gsLineThicknessF = 2.0

  respl2 = True
  respl2@gsLineThicknessF = 1.0

  respg = True
  respg@tfPolyDrawOrder = "Predraw"

  xpoly = new(2*numTime,float)
  ypoly = new(2*numTime,float)
  xpoly(0:numTime-1) = xstats(0,:)
  xpoly(numTime:2*numTime-1:-1) = xstats(0,:)

  do reg = 0, numReg-1
    regName  = REGIONS@id(reg)
    regTitle = REGIONS@title(reg)

    if (maskReg .eq. "global") then
      wks   = gsn_open_wks ("ps",var+".plot."+regName)
      res@gsnRightString = regTitle
    else
      wks   = gsn_open_wks ("ps",var+".plot."+regName+"."+maskReg)
      res@gsnRightString = regTitle + " ("+maskReg+")"
    end if

    do istat = istat0, 4
      res@tiYAxisString = yaxis_labels(istat)
      quartiles = quartile_values(stats(:,reg,:,istat)) ; quartiles(nstat,month) nstat=(min,lqrt,med,uqrt,max)

      ; set nice limits from min/max
      ymin = min(quartiles(0:4:4,:))
      ymax = max(quartiles(0:4:4,:))
      ymin = min((/ymin,min(mstats(:,reg,:,istat))/))
      ymax = max((/ymax,max(mstats(:,reg,:,istat))/))
      minmax = nice_mnmxintvl(ymin,ymax,20,True)
      res@trYMinF  =  minmax(0)
      res@trYMaxF  =  minmax(1)

      ; initial legend labels & colors
      LEG = True
      LEG@labels = "CM4"
      LEG@colors = "red"

      ; plot the extremes
      plot = gsn_csm_xy( wks, xstats(0:1,:), quartiles(0:4:4,:), res ) 

      ; shading for cmip5 models
      do n = 0, 1
        ypoly(0:numTime-1) = quartiles(n,:)
        ypoly(numTime:2*numTime-1:-1) = quartiles(4-n,:)
        respg@gsFillColor = shade_colors(n)
        str = unique_string("poly")
        plot@$str$ = gsn_add_polygon(wks, plot, xpoly, ypoly, respg)
      end do

      ; observation (average only)
      if (istat .eq. 0) then
        str = unique_string("obs")
        plot@$str$ = gsn_add_polyline(wks, plot, xstats(0,:), obsavg, resplo)
      end if

      ; GFDL/CMIP5 model
      if (.not.ismissing(purple(0))) then
        respl2@gsLineColor = "purple"
        do n = 0, dimsizes(purple)-1
          str = unique_string("cm2p1_")
          plot@$str$ = gsn_add_polyline(wks, plot, xstats(purple(n),:), stats(purple(n),reg,:,istat), respl2)
        end do
        append_array_value("CM2p1",LEG,"labels",False)
        append_array_value("purple",LEG,"colors",False)
      end if
      if (.not.ismissing(blue(0))) then
        respl2@gsLineColor = "blue"
        do n = 0, dimsizes(blue)-1
          str = unique_string("cm3_")
          plot@$str$ = gsn_add_polyline(wks, plot, xstats(blue(n),:), stats(blue(n),reg,:,istat), respl2)
        end do
        append_array_value("CM3",LEG,"labels",False)
        append_array_value("blue",LEG,"colors",False)
      end if

      ; CM4 model
      do n = 0, dimsizes(mfiles)-1
        respl@gsLineColor = "red"
        if (isAM4model(mfiles(n))) then
          respl@gsLineColor = "darkgreen"
          append_array_value("AM4",LEG,"labels",True)
          append_array_value("darkgreen",LEG,"colors",True)
        end if
        str = unique_string("cm4_")
        plot@$str$ = gsn_add_polyline(wks, plot, xstats(n,:), mstats(n,reg,:,istat), respl)
      end do

      draw(plot)

      ; draw informational label
      txres = True
      txres@txFontHeightF = 0.013
      txres@txFontColor = "grey85"
      gsn_text_ndc(wks, infolabel, 0.50, 0.05, txres)

      ; create a legend
      lgres = True
      lgres@lgAutoManage = False
      lgres@lgLineColors = LEG@colors
      lgres@lgLineThicknessF = 1.0
      lgres@vpWidthF = 0.14
      lgres@vpHeightF = 0.07
      lgres@lgPerimOn = False
      lgres@lgLabelFontHeightF = 0.05
      lgres@lgMonoDashIndex = True
      lgres@lgLabelOffsetF = .05

      xpos = 0.62
      ypos = 0.77
      gsn_legend_ndc(wks, dimsizes(LEG@labels), LEG@labels, xpos, ypos, lgres)
      delete(LEG)
        

      frame(wks)
    end do
  end do

end
