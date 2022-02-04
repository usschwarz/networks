#script to extract and compare the lenght of the actin filaments from different simulations
#used to generate data for Figure 3E and F
#author: Julia Jaeger


import matplotlib.pyplot as plt
import readdy
import numpy as np
import re
import glob
import seaborn as sns


# In[ ]:


def analyse_actin_lengths(file,finished):
    
    files_normal=sorted(glob.glob(file))
    lengths=[len(i) for i in files_normal]
    first_files=[]
    second_files=[]
    for i in files_normal:
        if len(i)==min(lengths):
            first_files.append(i)
        else:
            second_files.append(i)
    files_normal=first_files+second_files
    
    print(files_normal)
    n_files=int(len(files_normal))

    traj = readdy.Trajectory(files_normal[0])
    time, topology_records = traj.read_observable_topologies()
    time_b, counts = traj.read_observable_number_of_particles()

    if finished:
        l=1
    else:
        l=2
    for i in range(len(files_normal)-l): #because the last file is not finished
        print(files_normal[i+1])
        traj1 = readdy.Trajectory(files_normal[i+1])
        time1, topology_records_1 = traj1.read_observable_topologies()
        time1_b, counts_1 = traj1.read_observable_number_of_particles()


        time=np.concatenate((time, [time[-1]+i for i in time1]), axis=0)
        counts=np.concatenate((counts, counts_1), axis=0)
        topology_records=np.concatenate((topology_records,topology_records_1), axis=0)   
        
    length = np.zeros((203,len(topology_records)))
    sum_actin=np.zeros((len(topology_records)))
    plot_time=[time]
    av_length=[i/46. for i in counts.T[4]]
    

    for i in range(len(topology_records)-1):
        plot_time=np.append(plot_time,[time],axis=0)
        # for each time step
    t=0
    for topologies in topology_records:
        # gather average polymer length
        #print(len(topologies))
        i=0
        for top in topologies:
            #print(len(top.edges))
            length[i,t] = len(top.edges)+1
            sum_actin[t]+= len(top.edges)+1
            i+=1
        t+=1
        
    for i in range(len(sum_actin)):
        #subtract spectrin
        sum_actin[i]=(sum_actin[i]-157*39)/46.
        
    #plt.plot(time*0.00001,sum_actin)
    #plt.axhline(46*8,xmin=0,xmax=100*(n_files-1),color="k",linewidth=3)
    #plt.title("available actin: " + str(4926+46*6))
    #plt.show()

    return np.array(time),length,n_files,av_length


# In[ ]:


def plot_lengths(time,length,length2):
    
    new_length=np.array([length[0]])
    for i in range(45):
        new_length=np.append(new_length,[length[i+1]],axis=0)
        
    for i in range(46):
        new_length=np.append(new_length,[length2[i]],axis=0)
    
    final=np.array(new_length.T[-1])
   
    #plt.plot(time*0.00001,new_length.T)
    #plt.axhline(8,xmin=0,xmax=100*(n_files-1),color="k",linewidth=3)
    #plt.ylim(1,30)
    #plt.xlabel(r"time in $\mu s$")
    #plt.show()  
    
    sns.set_style('darkgrid')
    sns.histplot(final,bins=np.linspace(min(final)-0.5,max(final)-0.5,max(final)-min(final)+1))
    plt.show()
    
    return final
 


# In[ ]:


time_normal,length_normal,n_files_normal,av_lengths=analyse_actin_lengths("caps/616/caps*.h5",True)
time_normal_2,length_normal_2,n_files_normal_2,av_lengths_no=analyse_actin_lengths("caps/616/no_caps*.h5",True)
time_normal_3,length_normal_3,n_files_normal_3,av_lengths_3=analyse_actin_lengths("caps/616/less_adducin*.h5",True)
time_normal_4,length_normal_4,n_files_normal_4,av_lengths_4=analyse_actin_lengths("caps/616/tropomyosin*.h5",True)

np.save("caps/616/time.npy",time_normal)
np.save("caps/616/caps.npy",length_normal)
np.save("caps/616/av_lengths.npy",av_lengths)

np.save("caps/616/time_no.npy",time_normal_2)
np.save("caps/616/no_caps.npy",length_normal_2)
np.save("caps/616/av_lengths_no.npy",av_lengths_no)

np.save("caps/616/time_less_adducin.npy",time_normal_3)
np.save("caps/616/less_adducin.npy",length_normal_3)
np.save("caps/616/av_lengths_less_adducin.npy",av_lengths_3)

np.save("caps/616/time_tropomyosin.npy",time_normal_4)
np.save("caps/616/tropomyosin.npy",length_normal_4)
np.save("caps/616/av_lengths_tropomyosin.npy",av_lengths_4)


time_normal,length_normal,n_files_normal,av_lengths=analyse_actin_lengths("caps/1232/caps*.h5",True)
time_normal_2,length_normal_2,n_files_normal_2,av_lengths_no=analyse_actin_lengths("caps/1232/no_caps*.h5",True)
time_normal_3,length_normal_3,n_files_normal_3,av_lengths_3=analyse_actin_lengths("caps/1232/less_adducin*.h5",True)
time_normal_4,length_normal_4,n_files_normal_4,av_lengths_4=analyse_actin_lengths("caps/1232/tropomyosin*.h5",True)

np.save("caps/1232/time.npy",time_normal)
np.save("caps/1232/caps.npy",length_normal)
np.save("caps/1232/time_no.npy",time_normal_2)
np.save("caps/1232/no_caps.npy",length_normal_2)

np.save("caps/1232/av_lengths.npy",av_lengths)
np.save("caps/1232/av_lengths_no.npy",av_lengths_no)

np.save("caps/1232/time_less_adducin.npy",time_normal_3)
np.save("caps/1232/less_adducin.npy",length_normal_3)
np.save("caps/1232/av_lengths_less_adducin.npy",av_lengths_3)

np.save("caps/1232/time_tropomyosin.npy",time_normal_4)
np.save("caps/1232/tropomyosin.npy",length_normal_4)
np.save("caps/1232/av_lengths_tropomyosin.npy",av_lengths_4)


time_normal,length_normal,n_files_normal,av_lengths=analyse_actin_lengths("caps/1848/caps*.h5",True)
time_normal_2,length_normal_2,n_files_normal_2,av_lengths_no=analyse_actin_lengths("caps/1848/no_caps*.h5",True)
time_normal_3,length_normal_3,n_files_normal_3,av_lengths_3=analyse_actin_lengths("caps/1848/less_adducin*.h5",True)
time_normal_4,length_normal_4,n_files_normal_4,av_lengths_4=analyse_actin_lengths("caps/1848/tropomyosin*.h5",True)

np.save("caps/1848/time.npy",time_normal)
np.save("caps/1848/caps.npy",length_normal)
np.save("caps/1848/time_no.npy",time_normal_2)
np.save("caps/1848/no_caps.npy",length_normal_2)

np.save("caps/1848/av_lengths.npy",av_lengths)
np.save("caps/1848/av_lengths_no.npy",av_lengths_no)

np.save("caps/1848/time_less_adducin.npy",time_normal_3)
np.save("caps/1848/less_adducin.npy",length_normal_3)
np.save("caps/1848/av_lengths_less_adducin.npy",av_lengths_3)

np.save("caps/1848/time_tropomyosin.npy",time_normal_4)
np.save("caps/1848/tropomyosin.npy",length_normal_4)
np.save("caps/1848/av_lengths_tropomyosin.npy",av_lengths_4)


time_normal,length_normal,n_files_normal,av_lengths=analyse_actin_lengths("caps/2463/caps*.h5",True)
time_normal_2,length_normal_2,n_files_normal_2,av_lengths_no=analyse_actin_lengths("caps/2463/no_caps*.h5",True)
time_normal_3,length_normal_3,n_files_normal_3,av_lengths_3=analyse_actin_lengths("caps/2463/less_adducin*.h5",True)
time_normal_4,length_normal_4,n_files_normal_4,av_lengths_4=analyse_actin_lengths("caps/2463/tropomyosin*.h5",True)

np.save("caps/2463/time.npy",time_normal)
np.save("caps/2463/caps.npy",length_normal)
np.save("caps/2463/time_no.npy",time_normal_2)
np.save("caps/2463/no_caps.npy",length_normal_2)

np.save("caps/2463/av_lengths.npy",av_lengths)
np.save("caps/2463/av_lengths_no.npy",av_lengths_no)

np.save("caps/2463/time_less_adducin.npy",time_normal_3)
np.save("caps/2463/less_adducin.npy",length_normal_3)
np.save("caps/2463/av_lengths_less_adducin.npy",av_lengths_3)

np.save("caps/2463/time_tropomyosin.npy",time_normal_4)
np.save("caps/2463/tropomyosin.npy",length_normal_4)
np.save("caps/2463/av_lengths_tropomyosin.npy",av_lengths_4)


time_normal,length_normal,n_files_normal,av_lengths=analyse_actin_lengths("caps/3079/caps*.h5",True)
time_normal_2,length_normal_2,n_files_normal_2,av_lengths_no=analyse_actin_lengths("caps/3079/no_caps*.h5",True)
time_normal_3,length_normal_3,n_files_normal_3,av_lengths_3=analyse_actin_lengths("caps/3079/less_adducin*.h5",True)
time_normal_4,length_normal_4,n_files_normal_4,av_lengths_4=analyse_actin_lengths("caps/3079/tropomyosin*.h5",True)

np.save("caps/3079/time.npy",time_normal)
np.save("caps/3079/caps.npy",length_normal)
np.save("caps/3079/time_no.npy",time_normal_2)
np.save("caps/3079/no_caps.npy",length_normal_2)

np.save("caps/3079/av_lengths.npy",av_lengths)
np.save("caps/3079/av_lengths_no.npy",av_lengths_no)

np.save("caps/3079/time_less_adducin.npy",time_normal_3)
np.save("caps/3079/less_adducin.npy",length_normal_3)
np.save("caps/3079/av_lengths_less_adducin.npy",av_lengths_3)

np.save("caps/3079/time_tropomyosin.npy",time_normal_4)
np.save("caps/3079/tropomyosin.npy",length_normal_4)
np.save("caps/3079/av_lengths_tropomyosin.npy",av_lengths_4)



