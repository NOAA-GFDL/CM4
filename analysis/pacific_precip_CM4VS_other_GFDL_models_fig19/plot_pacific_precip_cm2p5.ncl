
load "$BW_PACKAGE_ROOT/simple_legend.ncl"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

NUMSTAT = 5

function correlation(x,y,type)
local mx,my,mxx,myy,mxy,corr
begin
  mxx = dim_avg_Wrap(x*x)
  myy = dim_avg_Wrap(y*y)
  mxy = dim_avg_Wrap(x*y)
  if (type .eq. 0) then
    mx = dim_avg_Wrap(x)
    my = dim_avg_Wrap(y)
    corr = (mxy-mx*my)/(sqrt(mxx-mx*mx) * sqrt(myy-my*my)) ; centered (Pearson)
  else
    corr = mxy/sqrt(mxx * myy) ; uncentered
  end if
  if (dimsizes(dimsizes(x)) .eq. 3) then
    copy_VarCoords(x(:,:,0),corr)
  end if
  return corr
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; dimensions of x & o: [ssn,lat,lon]
function compute_zonal_statistics ( x[*][*][*]:numeric, o[*][*][*]:numeric )
local dims, stats, avgo
begin
  dims = dimsizes(x)
  stats = new((/NUMSTAT,dims(0),dims(1)/),typeof(x))
  stats(0,:,:) = dim_avg_Wrap(x) ; average
  avgo = dim_avg_Wrap(o)
  stats(1,:,:) = stats(0,:,:) - avgo  ; bias
  stats(2,:,:) = dim_rmsd_Wrap(x,o)   ; rms
  stats(3,:,:) = correlation(x,o,0)   ; correlation (centered)
  stats(4,:,:) = correlation(x,o,1)   ; correlation (uncentered)
  return stats
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mask_data_array( data:numeric, lmask[*][*]:logical )
local rank, data1, fill
begin
  rank = dimsizes(dimsizes(data))
  data1 = data
  fill = data1@_FillValue
  if (rank .eq. 2) then
    data1 = where(lmask, data, fill)
  else
    data1 = where(conform(data,lmask,(/rank-2,rank-1/)), data, fill)
  end if
  return data1
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function get_model_name (fpath[1]:string)
local parse,model
begin
 ;parse = str_split(str_get_cols(fpath,str_index_of_substr(fpath,"/",-1)+1,-1),"_")
  parse = str_split(fpath,"_")
  model = parse(2)
  return model
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

begin

  var = "pr"
  cpath = var+"/cmip5/obsgrid";
  mpath = var+"/cm4/obsgrid";
  opath = var+"/obs";

; data files
  cfiles = systemfunc("/bin/ls "+cpath+"/*-ssnclim.nc")
  mfiles = systemfunc("/bin/ls "+mpath+"/*-ssnclim.nc")
  ofile = systemfunc("/bin/ls "+opath+"/*-ssnclim.nc")
  numTime = 5

  ; compute statistics for east pacific region
  regTitle = "150W-90W"
  minLon   = 210
  maxLon   = 270
  minLat   = -14
  maxLat   = 14

; land-sea mask
  maskfile = systemfunc("/bin/ls "+mpath+"/sftlf*.nc")
  if (dimsizes(maskfile) .gt. 1) then
    print("ERROR: more than one landsea mask file found")
    status_exit(1)
  end if
  fm = addfile(maskfile,"r")
  sftlf = fm->sftlf
  delete(fm)

  ; number latitude points in region
  latSlice = sftlf({minLat:maxLat},0)
  latitudes = latSlice&$latSlice!0$
  numLat = dimsizes(latitudes)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; read the observed data
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  fo = addfile(ofile,"r")
  odata = fo->$var$
 ;area = conform(odata(0,:,:),cos(odata&$odata!1$*DTR),0)
 ;copy_VarCoords(odata(0,:,:),area)
  obsavg = dim_avg_Wrap(odata(:,{minLat:maxLat},{minLon:maxLon}))

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
  dims = (/ numRuns, NUMSTAT, numTime, numLat /)
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

   ; scale data if necessary
    if (var .eq. "pr") then
      cdata = cdata * 86400.
    end if
    delete(fc)

    mask1  = where(sftlf .lt. 0.01, True, False)
    cdata1 = mask_data_array(cdata, mask1)
    odata1 = mask_data_array(odata, mask1)
  
    stats(m,:,:,:) = compute_zonal_statistics( cdata1(:,{minLat:maxLat},{minLon:maxLon}), \
                                               odata1(:,{minLat:maxLat},{minLon:maxLon}) )
    delete([/cdata1,odata1/])
  end do

  ; setup blue markers for "CM3" model
  if (oneRunPerModel) then
    blue = ind(model_names .eq. "GFDL-CM3")
    purple = ind(model_names .eq. "GFDL-CM2p1")
   ;purple = array_append_record(ind(model_names .eq. "GFDL-ESM2G"), \
   ;                            ind(model_names .eq. "GFDL-ESM2M"),0)
  else
    blue = ind(model .eq. "GFDL-CM3")
    purple = ind(model .eq. "GFDL-CM2p1")
   ;purple = array_append_record(ind(model .eq. "GFDL-ESM2G"), \
   ;                            ind(model .eq. "GFDL-ESM2M"),0)
  end if

  ;;;;;;;;;;;;;;;;;
  ;;  CM4 model  ;;
  ;;;;;;;;;;;;;;;;;
  dims = (/ dimsizes(mfiles), NUMSTAT, numTime, numLat /)
  dims_avg = (/ numLat /)

  mstats = new(dims,float)
  mstats_avg = new(dims_avg, float)

  do m = 0, dimsizes(mfiles)-1
    fm = addfile(mfiles(m),"r")
    mdata = fm->$var$
    if (var .eq. "pr") then
      mdata = mdata * 86400.
    end if
    if (m .eq. 0 .and. isfilevaratt(fm,var,"long_name")) then
      plotTitle = fm->$var$@long_name
    end if
    delete(fm)

    mask1  = where(sftlf .lt. 0.01, True, False)
    mdata1 = mask_data_array(mdata, mask1)
    odata1 = mask_data_array(odata, mask1)
  
    mstats(m,:,:,:) = compute_zonal_statistics( mdata1(:,{minLat:maxLat},{minLon:maxLon}), \
                                                odata1(:,{minLat:maxLat},{minLon:maxLon}) )

        alist = [/latitudes,mstats(m,0,0,:)/]
	print("==============================")
	print_table([/"LAT","PRE"/],"%15s %15s")
	print_table(alist,"%15.5f %15.5f")

    delete([/mdata1,odata1/])
  end do

; calculate the ensemble average
  mstats_avg = dim_avg_n(mstats(:,0,0,:),0) 
	alist = [/latitudes,mstats_avg(:)/]
	print("***********Ensemble Average ***********")
	print_table([/"LAT","PRE"/],"%15s %15s")
	print_table(alist,"%15.5f %15.5f")

  

  ;;;;;;;;;;;;;;;;;;;
  ;;  CM2.5 model  ;;
  ;;;;;;;;;;;;;;;;;;;
  fname = "/home/h1g/awg/xanadu/utils/CM4_paper_plots_BW/delworth_et_al_2012_fig8_cm25.dat"
  filin = asciiread(fname,-1,"string")
  data0 = str_split_csv(filin,"        ",3)
  dims0 = dimsizes(data0)
  data1 = tofloat(data0(1:dims0(0)-1,:))
  dat = data1(:,1)
  lat = data1(:,0)
  lat!0 = "lat"
  lat&lat = lat
  dat!0 = "lat"
  dat&lat = lat

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

  yLat = conform(stats(:,0,0,:),latitudes,1)
  yLat@units = "degrees_N"
  yLat@long_name = "Latitude"
  axis_labels = (/ "Avg", "Bias", "RMSE", "Corr", "uCorr" /)
  season_labels = (/ "Annual", "DJF", "MAM", "JJA", "SON" /)

  ; line colors
  legend_colors = (/ "black",     "blue", "firebrick4", "purple", "red",  "darkgreen"/)
  legend_labels = (/ "GPCP-v2.3", "CM3",  "CM2p5",      "CM2p1",  "CM4",  "AM4"/)
  shade_colors = (/ "grey85", "grey65" /)

  res = True
  res@gsnFrame = False
  res@gsnDraw = False
  res@xyMonoLineColor = True
  res@xyMonoDashPattern = True
  res@xyLineThicknessF = 1
  res@xyLineColor = shade_colors(0)
  res@tiYAxisString = "Latitude"
  res@gsnLeftString = "Precipitation ("+regTitle+")"

  respl = True
  if (dimsizes(mfiles) .eq. 1) then
    respl@gsLineThicknessF = 3.0
  else if (dimsizes(mfiles) .eq. 2) then
    respl@gsLineThicknessF = 2.0
  else
    respl@gsLineThicknessF = 1.0
  end if
  end if

  resb = True
  resb@gsLineThicknessF = 1.0

  reso = True
  reso@gsLineColor = legend_colors(ind(legend_labels .eq. "GPCP-v2.3"))
  reso@gsLineThicknessF = 2.0

  respg = True
  respg@tfPolyDrawOrder = "Predraw"

  ; shaded polygons for cmip6 models
  xpoly = new(2*numLat,float)
  ypoly = new(2*numLat,float)
  ypoly(0:numLat-1) = latitudes
  ypoly(numLat:2*numLat-1:-1) = latitudes


  do istat = 0, 0
    wks   = gsn_open_wks ("ps",var+".plot."+str_lower(axis_labels(istat)))
    res@tiXAxisString = axis_labels(istat) + " (mm/day)"

    do isea = 0, 0
      res@gsnRightString = season_labels(isea)

      quartiles = quartile_values(stats(:,istat,isea,:)) ; quartiles(qstat,lat)

      ; set nice limits from min/max
      ymin = min(quartiles(0:4:4,:))
      ymax = max(quartiles(0:4:4,:))
      ymin = min((/ymin,min(mstats(:,istat,isea,:))/))
      ymax = max((/ymax,max(mstats(:,istat,isea,:))/))
      minmax = nice_mnmxintvl(ymin,ymax,20,True)
      res@trXMinF  =  minmax(0)
      res@trXMaxF  =  minmax(1)

      ; plot the extremes
      plot = gsn_csm_xy( wks, quartiles(0:4:4,:), yLat(0:1,:), res )

      ; get the plot limits
      getvalues plot
        "trXMinF" : trXMinF
        "trXMaxF" : trXMaxF
        "trYMinF" : trYMinF
        "trYMaxF" : trYMaxF
      end getvalues

      ; shading for cmip5 models
      do n = 0, 1
        xpoly(0:numLat-1) = quartiles(n,:)
        xpoly(numLat:2*numLat-1:-1) = quartiles(4-n,:)
        respg@gsFillColor = shade_colors(n)
        str = unique_string("poly")
        plot@$str$ = gsn_add_polygon(wks, plot, xpoly, ypoly, respg)
      end do

      ; observation (average only)
      if (istat .eq. 0) then
        str = unique_string("obs")
        plot@$str$ = gsn_add_polyline(wks, plot, obsavg(isea,:), yLat(0,:), reso)
      end if

      ; CM2p1 model
      if (.not.ismissing(purple(0))) then
        resb@gsLineColor = legend_colors(ind(legend_labels .eq. "CM2p1"))
        do n = 0, dimsizes(purple)-1
          str = unique_string("cm2p1_")
          plot@$str$ = gsn_add_polyline(wks, plot, stats(purple(n),istat,isea,:), yLat(purple(n),:), resb)
        end do
      end if

      ; CM2.5 model
      if (istat .eq. 0  .and. isea .eq. 0 ) then 
        resb@gsLineColor = legend_colors(ind(legend_labels .eq. "CM2p5"))
        str = unique_string("cm2p5_")
        plot@$str$ = gsn_add_polyline(wks, plot, dat({minLat:maxLat}), lat({minLat:maxLat}), resb)
      end if

      ; CM3 model
      if (.not.ismissing(blue(0))) then
        resb@gsLineColor = legend_colors(ind(legend_labels .eq. "CM3"))
        do n = 0, dimsizes(blue)-1
          str = unique_string("cm3_")
          plot@$str$ = gsn_add_polyline(wks, plot, stats(blue(n),istat,isea,:), yLat(blue(n),:), resb)
        end do
      end if

      ; CM4 model
      do n = 0, dimsizes(mfiles)-1
        respl@gsLineColor = legend_colors(ind(legend_labels .eq. "CM4"))
        if (isAM4model(mfiles(n))) then
          respl@gsLineColor = legend_colors(ind(legend_labels .eq. "AM4"))
        end if
        str = unique_string("cm4_")
        plot@$str$ = gsn_add_polyline(wks, plot, mstats(n,istat,isea,:), yLat(n,:), respl)
      end do

      ;;;;;;;;;;;;;;;;;;;;;;;;
      ;  add/attach a legend
      ;;;;;;;;;;;;;;;;;;;;;;;;
      ls = 1
      if (istat .eq. 0) then
        ls = 0
      end if
      xleg = trXMinF + 0.67*(trXMaxF-trXMinF)
      yleg = trYMinF + 0.55*(trYMaxF-trYMinF)
      lgres = True
      lgres@lineColors = legend_colors(ls::-1)
      lgres@lineDashSegLen = .050 ; not same as plot seg length? (.12)
      lgres@fontHeight = .012
      legend = simple_legend (wks,plot,legend_labels(ls::-1),xleg,yleg,lgres)
      delete(lgres@lineColors)


      draw(plot)

      ; draw informational label
      txres = True
      txres@txFontHeightF = 0.013
      txres@txFontColor = "grey85"
      gsn_text_ndc(wks, infolabel, 0.50, 0.05, txres)

      frame(wks)
    end do
  end do
end

