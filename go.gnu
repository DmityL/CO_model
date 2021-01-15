#clear
reset
#set term postscript eps enhanced color solid 'Arial' 16
#set output "output.eps"
#set term windows
#set out

# Множитель e^(-tau), где tau - оптичеcкая толщина
tgau(v,t,v0,d) = exp(-abs(t)*exp(-(v-v0)**2/d**2)) 

load 'const_param.gnu'

do for [clump = 6:10] {
	print(sprintf("Clump%03d",clump))

	# Загружаем начальное приближение и константы
	load sprintf("_spect_comb/Clump%03d.dat.dat",clump)



	# Яркостные температуры на выходе из каждого слоя  для перехода 2->1
	temp(v,r,d_v,T0) =                 T0*(1.0/(exp(T0/tx1) - 1) - 1.0/(exp(T0/tbg) - 1))*(1-tgau(v+d_v,t1*r,v1,d1))
	#temp2(v,r,d_v,T0) =  temp(v,r,d_v,T0)*tgau(v+d_v,t2*r,v2,d2) + T0*(1.0/(exp(T0/tx2) - 1) - 1.0/(exp(T0/tbg) - 1))*(1-tgau(v+d_v,t2*r,v2,d2))
	#temp3(v,r,d_v,T0) = temp2(v,r,d_v,T0)*tgau(v+d_v,t3*r,v3,d3) + T0*(1.0/(exp(T0/tx3) - 1) - 1.0/(exp(T0/tbg) - 1))*(1-tgau(v+d_v,t3*r,v3,d3))
	#temp4(v,r,d_v,T0) = temp3(v,r,d_v,T0)*tgau(v+d_v,t4*r,v4,d4) + T0*(1.0/(exp(T0/tx4) - 1) - 1.0/(exp(T0/tbg) - 1))*(1-tgau(v+d_v,t4*r,v4,d4))

	# Профили перехода 1>0
	profile_12CO_10(v) = temp(v,1.0,0,5.5321);	# Профиль линии 12CO 1->0
	profile_13CO_10(v) = temp(v,1.0/rat,-35,5.288);	# Профиль линии 13CO 1->0

	# Профили перехода 2>1
	profile_12CO_21(v) = temp(v,1.0,-70,11.0639);	# Профиль линии 12CO 2->1
	profile_13CO_21(v) = temp(v,1.0/rat,-105,10.57738);	# Профиль линии 13CO 2->1


	# Общий профиль (сумма 4 линий)
	pro(x) =  profile_12CO_21(x)  + profile_13CO_21(x) + profile_12CO_10(x) +  profile_13CO_10(x);  # Общий профиль

	################## Fitting #######################
	FIT_LIMIT = 3.e-9
	FIT_START_LAMBDA = .1
	FIT_MAXITER = 50

	set fit errorvariables

	spectra = sprintf("_spect_comb/Clump%03d.dat.spec",clump)

	fit pro(x) spectra using 1:2:3  via tx1,t1,v1,d1

	set print sprintf("_spect_comb/Clump%03d.fit.dat",clump)
	print "V1	Tx1	t1	d1"
	print sprintf("%5.2f	%5.1f	%5.1f	%5.3f",v1,tx1,t1,d1)
	print sprintf("±%5.3f	±%3.1f	±%3.1f	±%5.3f",v1_err,tx1_err,t1_err,d1_err)
	set print

	#tx1=ftx(mtx1)
	#tx2=ftx(mtx2)
	#tx3=ftx(mtx3)
	#nu12=fnu(mnu12)
	#nu13=fnu(mnu13)


}