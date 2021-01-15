import numpy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from matplotlib.ticker import MaxNLocator
from matplotlib import rcParams
rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
rcParams['font.size']=14.0
rcParams['font.sans-serif']='Verdana'



v0 = -17.5
dv = 35


fig= plt.figure(figsize=(12,4))
ax= fig.add_axes([0.1,0.1,0.8,0.8])

#fig, ax = plt.subplots()

#fig, axs = plt.subplots(2, sharex=True, sharey=True, gridspec_kw={'hspace': 0})
#fig.suptitle('Sharing both axes')
#axs[0].plot(*numpy.loadtxt("05343+3605.spec",unpack=True,usecols=(0,1)), color='black', linewidth=1.0, linestyle='-')
#axs[0].plot(*numpy.loadtxt("spectra_model.spec",unpack=True,usecols=(0,1)), color='red', linewidth=1.0, linestyle='-')

# Hide x labels and tick labels for all but bottom plot.
#for ax in axs:
#    ax.label_outer()


plt.xlabel('Velocity (km s$^{-1}$)')
plt.ylabel('T$_{mb} (K)$')
#plt.ylim((25,250))




for i in range(1,146):
	i_str = "%03d" % i
	print("Working on Clump"+i_str)
	plt.cla()
	plt.xlim((-35,110))
	plt.plot(*numpy.loadtxt("./_spect_comb/Clump"+i_str+".dat.spec",unpack=True,usecols=(0,1)), color='black', linewidth=1.0, linestyle='-')
	plt.plot(*numpy.loadtxt("./_spect_comb/Clump"+i_str+".mod.spec",unpack=True,usecols=(0,1)), color='red', linewidth=1.0, linestyle='-')

	ax.axvline(x=v0+dv/2,linewidth=1, color='black')
	ax.axvline(x=v0+dv/2+dv,linewidth=1, color='black')
	ax.axvline(x=v0+dv/2+dv*2,linewidth=1, color='black')

	ax.text(v0+10, 5, r'$^{12}$CO(1-0)', fontsize=15,horizontalalignment='center')
	ax.text(v0+dv+10, 5, r'$^{13}$CO(1-0)', fontsize=15,horizontalalignment='center')
	ax.text(v0+dv*2+10, 5, r'$^{12}$CO(2-1)', fontsize=15,horizontalalignment='center')
	ax.text(v0+dv*3+10, 5, r'$^{13}$CO(2-1)', fontsize=15,horizontalalignment='center')

	#plt.xlim(left=-30,right=100)
	plt.savefig("init_model_plot/Clump"+i_str+'_profile.jpg',dpi=300,format='jpg',bbox_inches='tight')








