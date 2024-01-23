
using PyCall
using LinearAlgebra
using Printf
using StatsBase
using Statistics

flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Start")
plt=pyimport("matplotlib.pyplot")
pygui()
coil=[292.52 273.38 274.40 285.63 281.25 285.37]
typeA=[160 147 150 140 143 148]
typeD=[85 85 85 90 71]
typeC=[23 128 88 58 88]
paba=[7.6 10.5 10.8 10.6]
x=[1.5,3,4.5,6]
y=[mean(coil),mean(typeA),mean(typeD),mean(typeC)]
yerr=[std(coil),std(typeA),std(typeD),std(typeC)]

plt.errorbar(x,y,yerr,fmt="o", linewidth=1,markersize=4, capsize=4,color="b",label="Mean and STD for aligned chains")
plt.plot([9,9],[0,300],color="k",linestyle="--", linewidth=1)
plt.fill_between([0,18],[9,9],[14,14],alpha=.2, linewidth=0,color="r")
plt.errorbar(10,14,10,fmt="o",linewidth=2,capsize=5,color="g",markersize=0,capthick=2,label="Range for self assembled chains")
plt.errorbar(12,11.5,6.5,fmt="o",linewidth=2,capsize=5,color="g",markersize=0,capthick=2)
plt.errorbar(14,7.5,3,fmt="o",linewidth=2,capsize=5,color="g",markersize=0,capthick=2)
plt.xlim(left=1,right=15)
plt.xticks([1.5,3,4.5,6,10,12,14],labels=["CC","A","B","C","no","oaba","paba"])
plt.tick_params(labelsize=16) 
plt.ylim(top=300,bottom=0)
plt.ylabel("Persistence Length in nm",fontsize=16)
plt.text(1.2,10,"Experimental Value: 9-14 nm",fontsize=14)
plt.legend(loc="upper right",fontsize=14, borderpad=1,labelspacing=2,framealpha=1)
plt.tight_layout()
plt.savefig("PersitenceLength.pdf")
plt.show()
#plt.close()

flush(stdout) ##Only necessary if you run the program in Atom IDE
println("Calculations Finished")
