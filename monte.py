import numpy as np
from scipy.optimize import minimize
import emcee
import matplotlib.pyplot as plt
import sys  
import re
import emcee
import corner
from scipy.signal import lfilter, lfilter_zi, filtfilt, butter

plot = 1


tx_init=29.0
v_init=-19.00
t_init=74.0
d_init=0.899
tx2,t2,v2,d2 = [0,0,0,0]
tx3,t3,v3,d3 = [0,0,0,0]
N_fix = 0
rat=70.0
tbg=2.7

rx_dict = {
    'tbg': re.compile(r'tbg = (?P<tbg>.+)\r'),
	'rat': re.compile(r'rat = (?P<rat>.+)\r'),
	'tx1': re.compile(r'tx1 = (?P<tx1>.+)\r'),
	't1': re.compile(r't1 = (?P<t1>.+)\r'),
	'v1': re.compile(r'v1 = (?P<v1>.+)\r'),
	'd1': re.compile(r'd1 = (?P<d1>.+)\r'),
	'tx2': re.compile(r'tx2 = (?P<tx2>.+)\r'),
	't2': re.compile(r't2 = (?P<t2>.+)\r'),
	'v2': re.compile(r'v2 = (?P<v2>.+)\r'),
	'd2': re.compile(r'd2 = (?P<d2>.+)\r'),
	'tx3': re.compile(r'tx3 = (?P<tx3>.+)\r'),
	't3': re.compile(r't3 = (?P<t3>.+)\r'),
	'v3': re.compile(r'v3 = (?P<v3>.+)\r'),
	'd3': re.compile(r'd3 = (?P<d3>.+)\r')
}

def _parse_line(line):
    """
    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex
    """
    for key, rx in rx_dict.items():
        match = rx.search(line)
        if match:
            return key, match
    # if there are no matches
    return None, None



def tgau(v_,t_,v0_,d_):
	return np.exp(-np.fabs(t_)*np.exp(-(v_-v0_)**2.0/d_**2.0))

def temp(theta,v,r,d_v,T0):
	tx,t,v0,d = theta
	return T0*(1.0/(np.exp(T0/tx) - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t*r,v0,d))

def temp_n_radtran(theta,fix_params,n_fix,order,v,r,d_v,T0):
	tx1,t1,v1,d1 = theta
	if n_fix == 0:
		return T0*(1.0/(np.exp(T0/tx1)  - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t1* r,v1, d1))
	if n_fix == 1:
		tx2,t2,v2,d2 = fix_params
		if order == '21':
			tx1,t1,v1,d1 = fix_params
			tx2,t2,v2,d2 = theta
		T_first =									T0*(1.0/(np.exp(T0/tx1)  - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t1*r,v1,d1))
		return T_first*tgau(v+d_v,t2*r,v2,d2) +	T0*(1.0/(np.exp(T0/tx2)  - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t2*r,v2,d2))
	if n_fix == 2:
		tx2,t2,v2,d2,tx3,t3,v3,d3 = fix_params
		if order == '213':
			tx2,t2,v2,d2 = theta
			tx1,t1,v1,d1,tx3,t3,v3,d3 = fix_params
		if order == '312':
			tx3,t3,v3,d3 = theta
			tx1,t1,v1,d1,tx2,t2,v2,d2 = fix_params
		
		T_first =									 T0*(1.0/(np.exp(T0/tx1) - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t1*r,v1,d1))
		T_second = T_first*tgau(v+d_v,t2*r,v2,d2) + T0*(1.0/(np.exp(T0/tx2) - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t2*r,v2,d2))
		return    T_second*tgau(v+d_v,t3*r,v3,d3) + T0*(1.0/(np.exp(T0/tx3) - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t3*r,v3,d3))


def temp_n_simple(theta,fix_params,n_fix,v,r,d_v,T0):
	tx,t,v0,d = theta
	if n_fix == 0:
		return T0*(1.0/(np.exp(T0/tx)  - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t* r,v0,  d))
	if n_fix == 1:
		tx1,t1,v01,d1 = fix_params
		T_first = T0*(1.0/(np.exp(T0/tx)  - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t* r,v0, d))
		T_second = T0*(1.0/(np.exp(T0/tx1)  - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t1* r,v01, d1))
		return T_first+T_second
	if n_fix == 2:
		tx1,t1,v01,d1,tx2,t2,v02,d2 = fix_params
		T_first =  T0*(1.0/(np.exp(T0/tx)  - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t* r,v0, d))
		T_second = T0*(1.0/(np.exp(T0/tx1)  - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t1* r,v01, d1))
		T_third =  T0*(1.0/(np.exp(T0/tx2)  - 1.0) - 1.0/(np.exp(T0/tbg) - 1.0))*(1.0-tgau(v+d_v,t2* r,v02, d2))
		return T_first+T_second+T_third

	

def sum_lines(theta,v):
	return temp(theta,v,1.0,0.0,5.5321) +  temp(theta,v,1.0/rat,-35,5.288) + temp(theta,v,1.0,-70.0,11.0639) + temp(theta,v,1.0/rat,-105,10.57738);

def sum_lines_n(theta,fix_params,n_fix,order,v):
	return temp_n_radtran(theta,fix_params,n_fix,order,v,1.0,0.0,5.5321) +  temp_n_radtran(theta,fix_params,n_fix,order,v,1.0/rat,-35,5.288) + temp_n_radtran(theta,fix_params,n_fix,order,v,1.0,-70.0,11.0639) + temp_n_radtran(theta,fix_params,n_fix,order,v,1.0/rat,-105,10.57738);

def log_likelihood_temp(theta, x, y, yerr):
	model = sum_lines(theta,x)
	sigma2 = yerr ** 2 
	return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

def log_likelihood_temp_n(theta, fix_params, n_fix, order, x, y, yerr):
	model = sum_lines_n(theta,fix_params, n_fix,order,x)
	sigma2 = yerr ** 2 
	return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def get_velo_indices(v,vi,size_kms):
	logics1 = np.logical_and(v>=vi-size_kms/2.0, v<=vi+size_kms/2.0)
	logics2 = np.logical_and(v>=vi-size_kms/2.0+35, v<=vi+size_kms/2.0+35)
	logics3 = np.logical_and(v>=vi-size_kms/2.0+70, v<=vi+size_kms/2.0+70)
	logics4 = np.logical_and(v>=vi-size_kms/2.0+105, v<=vi+size_kms/2.0+105)
	ind = np.where(logics1 | logics2 | logics3 | logics4)
	return ind


for clump in range(18,19):
	N_fix = 0

	c = "%03d"%(clump)
	filelen = file_len("./_spect_comb/Clump"+c+".dat.dat")
	
	with open("./_spect_comb/Clump"+c+".dat.dat", 'r') as file_object:
		for i in range(filelen):
			line = file_object.readline()
			# at each line check for a match with a regex
			key, match = _parse_line(line)
			# extract school name
			if key == 'tx1':
				tx_init = float(match.group('tx1'))
			if key == 't1':
				t_init = float(match.group('t1'))
			if key == 'v1':
				v_init = float(match.group('v1'))
			if key == 'd1':
				d_init = float(match.group('d1'))
			if key == 'tx2':
				tx2 = float(match.group('tx2'))
				N_fix = 1
			if key == 't2':
				t2 = float(match.group('t2'))
			if key == 'v2':
				v2 = float(match.group('v2'))
			if key == 'd2':
				d2 = float(match.group('d2'))
			if key == 'tx3':
				tx3 = float(match.group('tx3'))
				N_fix = 2
			if key == 't3':
				t3 = float(match.group('t3'))
			if key == 'v3':
				v3 = float(match.group('v3'))
			if key == 'd3':
				d3 = float(match.group('d3'))

	print("Initial estimation:")
	print(tx_init,t_init,v_init,d_init)
	if N_fix > 0:
		print("Fixed parameters:")
		print(tx2,t2,v2,d2)	
		if N_fix > 1:
			print(tx3,t3,v3,d3)	



	spectra = np.genfromtxt("./_spect_comb/Clump"+c+".dat.spec",  delimiter='\t') #dtype=None,
	xi_full = spectra[:,0].copy()
	yi_full = spectra[:,1].copy()
	yi_full_unfilter = spectra[:,1].copy()
	yerri_full = spectra[:,2].copy()*2.0
	
	# Apply low-pass filter
	b, a = butter(3, 0.3)
	zi = lfilter_zi(b, a)
	z, _ = lfilter(b, a, yi_full, zi=zi*yi_full[0])
	z2, _ = lfilter(b, a, z, zi=zi*z[0])
	yi_full = filtfilt(b, a, yi_full)


	size_fit_kms = 6.0
	ind = get_velo_indices(xi_full,v_init,size_fit_kms)
	xi = np.array(xi_full)[ind] 
	yi = np.array(yi_full)[ind] 
	yerri = np.array(yerri_full)[ind] 

	ind2,ind3 = [0,0]
	xi2,yi2,yerri2,xi3,yi3,yerri3 = [0,0,0,0,0,0]
	if N_fix == 1:
		ind2 = get_velo_indices(xi_full,v2,size_fit_kms)
		xi2 = np.array(xi_full)[ind2] 
		yi2 = np.array(yi_full)[ind2] 
		yerri2 = np.array(yerri_full)[ind2] 
	if N_fix == 2:
		ind2 = get_velo_indices(xi_full,v2,size_fit_kms)
		xi2 = np.array(xi_full)[ind2] 
		yi2 = np.array(yi_full)[ind2] 
		yerri2 = np.array(yerri_full)[ind2] 
		ind3 = get_velo_indices(xi_full,v3,size_fit_kms)
		xi3 = np.array(xi_full)[ind3] 
		yi3 = np.array(yi_full)[ind3] 
		yerri3 = np.array(yerri_full)[ind3] 




	nll_temp = lambda *args: -log_likelihood_temp_n(*args)
	fix_params = [0]
	order = '12'
	if N_fix == 1: 
		fix_params = [tx2, t2, v2, d2]
	if N_fix == 2: 
		fix_params = [tx2, t2, v2, d2, tx3, t3, v3, d3]
	initial_temp = np.array([tx_init, t_init, v_init, d_init]) 
	soln_temp = minimize(nll_temp, initial_temp, args=(fix_params, N_fix, order, xi, yi, yerri))
	tx_est, t_est, v_est, d_est = soln_temp.x
	tx_est2, t_est2, v_est2, d_est2 = [0,0,0,0]
	tx_est3, t_est3, v_est3, d_est3 = [0,0,0,0]

	print("Maximum likelihood estimates (Comp 1):")
	print("Tx = {0:.3f}".format(tx_est))
	print("t = {0:.3f}".format(t_est))
	print("v = {0:.3f}".format(v_est))
	print("d = {0:.3f}".format(d_est))
	
	if N_fix == 1:
		fix_params = [tx_est, t_est, v_est, d_est]
		initial_temp = np.array([tx2, t2, v2, d2]) 
		order = '21'
		soln_temp2 = minimize(nll_temp, initial_temp, args=(fix_params, N_fix, order, xi2, yi2, yerri2))
		tx_est2, t_est2, v_est2, d_est2 = soln_temp2.x
		tx2,t2,v2,d2 = [tx_est2, t_est2, v_est2, d_est2] 
		print("Maximum likelihood estimates (Comp 2):")
		print("Tx = {0:.3f}".format(tx_est2))
		print("t = {0:.3f}".format(t_est2))
		print("v = {0:.3f}".format(v_est2))
		print("d = {0:.3f}".format(d_est2))

	if N_fix == 2:
		fix_params = [tx_est, t_est, v_est, d_est, tx3, t3, v3, d3]
		initial_temp = np.array([tx2, t2, v2, d2]) 
		order = '213'
		soln_temp2 = minimize(nll_temp, initial_temp, args=(fix_params, N_fix, order, xi2, yi2, yerri2))
		tx_est2, t_est2, v_est2, d_est2 = soln_temp2.x
		tx2,t2,v2,d2 = [tx_est2, t_est2, v_est2, d_est2] 
		print("Maximum likelihood estimates (Comp 2):")
		print("Tx = {0:.3f}".format(tx_est2))
		print("t = {0:.3f}".format(t_est2))
		print("v = {0:.3f}".format(v_est2))
		print("d = {0:.3f}".format(d_est2))

		fix_params = [tx_est, t_est, v_est, d_est, tx2, t2, v2, d2]
		initial_temp = np.array([tx3, t3, v3, d3]) 
		order = '312'
		soln_temp3 = minimize(nll_temp, initial_temp, args=(fix_params, N_fix, order, xi3, yi3, yerri3))
		tx_est3, t_est3, v_est3, d_est3 = soln_temp3.x
		tx3,t3,v3,d3 = [tx_est3, t_est3, v_est3, d_est3] 
		print("Maximum likelihood estimates (Comp 3):")
		print("Tx = {0:.3f}".format(tx_est3))
		print("t = {0:.3f}".format(t_est3))
		print("v = {0:.3f}".format(v_est3))
		print("d = {0:.3f}".format(d_est3))




	if plot == 1:
		if N_fix == 1: 
			fix_params = [tx2, t2, v2, d2]
		if N_fix == 2: 
			fix_params = [tx2, t2, v2, d2, tx3, t3, v3, d3]

		print("Saving spectra...")
		plt.clf()
		x0 = np.linspace(-35, 105, 700)
		plt.figure(figsize=(10,4))
		
		#plt.plot(x0, sum_lines_n(initial_temp,fix_params,N_fix,x0), label="Initial guess", lw=1)
		#xi_all = np.concatenate([xi, xi2])
		x_con = np.concatenate((xi, xi2, xi3), axis=None)
		y_con = np.concatenate((yi, yi2, yi3), axis=None)
		plt.plot(x_con, y_con, 'o', label="Data used for fit", markersize=3)	
		#plt.errorbar(xi, yi, yerr=yerri, fmt=".k", capsize=0, label="Data used for fit")
		plt.plot(xi_full, yi_full_unfilter, label="Obs. data", lw=1)
		plt.plot(xi_full, yi_full, label="Obs. data (filter)", lw=1)
		order = '12'
		plt.plot(x0, sum_lines_n([tx_est,t_est,v_est,d_est],fix_params,N_fix,order,x0), label="ML estimation", lw=2)
		plt.title('Spectra of Clump'+c)

		#plt.plot(x0, m_true * np.cos(x0) + b_true, "k", alpha=0.3, lw=3, label="truth")
		#plt.plot(x0, np.dot(np.vander(x0, 2), [m_ml, b_ml]), ":k", label="ML")
		plt.legend(fontsize=12)
		plt.xlim(-35,105)
		plt.xlabel("Vlsr (km/s)")
		plt.ylabel("Tmb (K)");
		plt.savefig('_pdf_plots/Clump'+c+'_spectra.pdf')
		plt.savefig('_pdf_plots/Clump'+c+'_spectra.eps')

	
	def log_prior_temp(theta):
		tx,t,v0,d = theta
		if (1.0 < tx < 100.0) and (0.0 < t < 150.0) and ( -25.0 < v0 < -5.0) and (0.1 < d < 3.0):
			return 0.0
		return -np.inf

	def log_probability_temp(theta, x, y, yerr):
		lp = log_prior_temp(theta)
		if not np.isfinite(lp):
			#print("Inf")
			return -np.inf
		#print(lp)
		return lp+log_likelihood_temp(theta, x, y, yerr)

	def log_probability_temp_n(theta, fix_params, n_fix, order, x, y, yerr):
		lp = log_prior_temp(theta)
		if not np.isfinite(lp):
			#print("Inf")
			return -np.inf
		#print(lp)
		return lp+log_likelihood_temp_n(theta, fix_params, n_fix, order, x, y, yerr)

	


	nwalk = 32
	steps = 3000
	
	fix_params = [0]
	if N_fix == 1:
		fix_params = [tx2, t2, v2, d2]

	pos = soln_temp.x + 1e-2 * np.random.randn(nwalk, 4)
	nwalkers, ndim = pos.shape

	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability_temp_n, args=(fix_params, N_fix, '', xi, yi, yerri))
	sampler.run_mcmc(pos, steps, progress=True);

	if plot == 1:

		print("Making steps plot")
		plt.clf()
		plt.title('Step plot of Clump'+c)
		fig, axes = plt.subplots(4, figsize=(10, 7), sharex=True)
		samples = sampler.get_chain()
		#labels = ["m", "b", "log(f)"]
		labels = ["Tex", "t", "v0", "d"]
		for i in range(ndim):
			ax = axes[i]
			ax.plot(samples[:, :, i], "k", alpha=0.3)
			ax.set_xlim(0, len(samples))
			ax.set_ylabel(labels[i])
			#ax.yaxis.set_label_coords(-0.1, 0.5)

		axes[-1].set_xlabel("step number");

		fig.savefig('_pdf_plots/Clump'+c+'_step.pdf')

	 
	flat_samples = sampler.get_chain(discard=50, thin=15, flat=True)
	tx1 = np.percentile(flat_samples[:, 0], 50)
	t1 = np.percentile(flat_samples[:, 1], 50)
	v1 = np.percentile(flat_samples[:, 2], 50)
	d1 = np.percentile(flat_samples[:, 3], 50)
	print(tx1,t1,v1,d1)

	if plot == 1:
		print("Making corner plot")
		fig2 = corner.corner(
			flat_samples, labels=labels, truths=[tx1, t1, v1, d1]
		);
		fig2.savefig('_pdf_plots/Clump'+c+'_corner.pdf')

	#plt.clf()
	#print("Making spectra plot")
	#x0 = np.linspace(-35, 105, 500)
	#inds = np.random.randint(len(flat_samples), size=100)
	#for ind in inds:
	#    sample = flat_samples[ind]
	#    plt.plot(x0, np.dot(np.vander(x0, 2), sample[:2]), "C1", alpha=0.1)
	#plt.errorbar(xi, yi, yerr=yerri, fmt=".k", capsize=0,label="Observed")
	#plt.plot(x0, sum_lines([tx_true,t_true,v_true,d_true],x0), label="True")
	#plt.legend(fontsize=14)
	#plt.xlim(-35, 100)
	#plt.xlabel("Vlsr")
	#plt.ylabel("Intensity");
	#
	#plt.savefig('Spectra.pdf')

	from IPython.display import display, Math
	from datetime import datetime
	now = datetime.now()
	dt_string = now.strftime("%d.%m.%Y %H:%M:%S")
	f = open("ML_estimate.txt", "a")
	f.write("\n============================================\n Clump "+c+": starting estimation at "+dt_string+" with steps = "+str(steps)+", Nwalkers = "+str(nwalk)+"\n")

	f.write("Maximum likelihood estimates (Comp 1):\n")
	f.write("Tx = {0:.3f}\n".format(tx_est))
	f.write("t = {0:.3f}\n".format(t_est))
	f.write("v = {0:.3f}\n".format(v_est))
	f.write("d = {0:.3f}\n".format(d_est))
	
	if tx_est2 > 0:
		f.write("Maximum likelihood estimates (Comp 2):\n")
		f.write("Tx = {0:.3f}\n".format(tx_est2))
		f.write("t = {0:.3f}\n".format(t_est2))
		f.write("v = {0:.3f}\n".format(v_est2))
		f.write("d = {0:.3f}\n".format(d_est2))
	if tx_est3 > 0:
		f.write("Maximum likelihood estimates (Comp 3):\n")
		f.write("Tx = {0:.3f}\n".format(tx_est3))
		f.write("t = {0:.3f}\n".format(t_est3))
		f.write("v = {0:.3f}\n".format(v_est3))
		f.write("d = {0:.3f}\n".format(d_est3))

	text = "Tex\tt\tv0\td\n"
	unc = ""
	for i in range(ndim):
		mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
		q = np.diff(mcmc)
		#txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
		#txt = txt.format(mcmc[1], q[0], q[1], labels[i])
		txt = "{0:.3f}"
		txt = txt.format(mcmc[1])
		txt_unc = "[-{0:.3f}, +{1:.3f}]"
		txt_unc = txt_unc.format(q[0],q[1])
		if i > 0:
			text += "\t"
			unc += "\t"
		text += txt
		unc += txt_unc
		#display(Math(txt))

	f.write(text)
	f.write("\n")
	f.write(unc+'\n')

	sampler.reset()



	if N_fix == 1:

		fix_params = [tx1,t1,v1,d1]

		pos = soln_temp2.x + 1e-2 * np.random.randn(nwalk, 4)
		nwalkers, ndim = pos.shape

		sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability_temp_n, args=(fix_params, N_fix, '21', xi, yi, yerri))
		sampler.run_mcmc(pos, steps, progress=True);

		flat_samples = sampler.get_chain(discard=50, thin=15, flat=True)
		tx2 = np.percentile(flat_samples[:, 0], 50)
		t2 = np.percentile(flat_samples[:, 1], 50)
		v2 = np.percentile(flat_samples[:, 2], 50)
		d2 = np.percentile(flat_samples[:, 3], 50)
		print(tx2,t2,v2,d2)

		if plot == 1:
			print("Making corner plot for comp 2")
			fig3 = corner.corner(
				flat_samples, labels=labels, truths=[tx2, t2, v2, d2]
			);
			fig3.savefig('_pdf_plots/Clump'+c+'_corner_comp2.pdf')

		text = "Comp 2\n"
		unc = ""
		for i in range(ndim):
			mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
			q = np.diff(mcmc)
			#txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
			#txt = txt.format(mcmc[1], q[0], q[1], labels[i])
			txt = "{0:.3f}"
			txt = txt.format(mcmc[1])
			txt_unc = "[-{0:.3f}, +{1:.3f}]"
			txt_unc = txt_unc.format(q[0],q[1])
			if i > 0:
				text += "\t"
				unc += "\t"
			text += txt
			unc += txt_unc

		sampler.reset()



	f.write(text)
	f.write("\n")
	f.write(unc+'\n')
	f.close()



	print("Saving Monte-Carlo spectra...")
	plt.clf()
	x0 = np.linspace(-35, 105, 700)
	plt.figure(figsize=(10,4))
	
	#plt.plot(x0, sum_lines_n(initial_temp,fix_params,N_fix,x0), label="Initial guess", lw=1)
	#xi_all = np.concatenate([xi, xi2])
	x_con = np.concatenate((xi, xi2, xi3), axis=None)
	y_con = np.concatenate((yi, yi2, yi3), axis=None)
	plt.plot(x_con, y_con, 'o', label="Data used for fit", markersize=3)	
	#plt.errorbar(xi, yi, yerr=yerri, fmt=".k", capsize=0, label="Data used for fit")
	plt.plot(xi_full, yi_full_unfilter, label="Obs. data", lw=1)
	plt.plot(xi_full, yi_full, label="Obs. data (filter)", lw=1)
	order = ''
	fix_params = [0]
	if N_fix == 1:
		fix_params = [tx2,t2,v2,d2]
	if N_fix == 2:
		fix_params = [tx2,t2,v2,d2,tx3,t3,v3,d3]
	plt.plot(x0, sum_lines_n([tx1,t1,v1,d1],fix_params,N_fix,order,x0), label="MC estimation", lw=2)
	plt.title('Spectra of Clump'+c)

	#plt.plot(x0, m_true * np.cos(x0) + b_true, "k", alpha=0.3, lw=3, label="truth")
	#plt.plot(x0, np.dot(np.vander(x0, 2), [m_ml, b_ml]), ":k", label="ML")
	plt.legend(fontsize=12)
	plt.xlim(-35,105)
	plt.xlabel("Vlsr (km/s)")
	plt.ylabel("Tmb (K)");
	plt.savefig('_pdf_plots/Clump'+c+'_spectra_MC.pdf')
	plt.savefig('_pdf_plots/Clump'+c+'_spectra_MC.eps')
