
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"



;;; NCL script used to calculate relative vorticity and other necessary dynamical terms
;;; for Rencurrel and Rose (2019). Output files available as "rel_vort_$model.cam.h0.clim.nc"
;;; at http://thredds.atmos.albany.edu:8080/thredds/catalog/cesm/aquaplanet_som/catalog.html


begin
	model=(/"CO2","ctrl","N1_1PW","N1_2PW","N1_3PW","N1_4PW","N2_1PW","N2_2PW","N2_3PW","N2_4PW","N3_1PW","N3_2PW","N3_3PW","N3_4PW","N4_1PW","N4_2PW","N4_3PW","N4_4PW","N6_1PW","N6_2PW","N6_3PW","N6_4PW","N8_1PW","N8_2PW","N8_3PW","N8_4PW"/)

	do t= 0,24
	   sfile=addfile("/network/daes-scratch/data/rose_scr/mr311693/OHT_repository/QAqu_"+model(t)+".cam.h0.clim.nc","r")
	   nfile=addfile("/network/rit/home/mr311693/ncl_scripts/HC_data/output/rel_vort_"+model(t)+".cam.h0.clim.nc","c")
	   lat=sfile->lat
	   lev=sfile->lev
	   lon=sfile->lon
	   u=sfile->U
	   v=sfile->V

   	   rel_vort=uv2vr_cfd(u,v,lat,lon,2)

   	   rv=rm_single_dims(rel_vort)


	   omega=sfile->OMEGA

	   temp=sfile->T

	   p= lev*100.

	  ; w_convert= omega_to_w(omega, p, temp)

       	   rgas = 287.058            ; J/(kg-K) => m2/(s2 K)
           g    = 9.80665            ; m/s2
	   rho  = p/(rgas*temp)         ; density => kg/m3
	   w    = -omega/(rho*g)     ; array operation







   	   pi=3.14159265359
    	   f = 2.*2.*pi/ ( 60.*60.*24. ) *sin ( pi/180.*lat(:) )
   	   r=dim_avg_n_Wrap(rv, 2)

   	   f1=doubletofloat(f)
	   f_test=onedtond(f,(/26,96/))
   	   R_0=-r/f_test
   
	   nfile->lat=lat
	   nfile->lon=lon
	   nfile->lev=lev
	   nfile->U=u
	   nfile->V=v
           nfile->rel_vort=rv
	   nfile->coriolis=f_test
	   nfile->R_0=R_0
	   nfile->W=w
	end do 
	
end
