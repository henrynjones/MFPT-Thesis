import numpy as np
import matplotlib.pyplot as plt

fine = 0.001*np.arange(1,16)
#course = 0.02*np.linspace(1,2,5)
#xes = np.concatenate((fine,course))
#xes = 0.2*((0.5)**np.arange(0,7))
ctimes = 2000*((0.2/0.001)**2)
#rtimes = 1000*((0.2/xes)**2)
#s_i,nc = np.max(xes)
#with min(xes) = 0.003125, and ctimes = ctimes*0.85, the max of ctime was 
#we were running 8 million steps
#now we are running 80 million steps
steps = int(ctimes)
del(ctimes)
#del(rtimes)
numpts = 1
numfunctcalls = 3000
histoarr = np.zeros((np.size(fine),numpts*numfunctcalls+2))
histoarr[:,0]=fine


"""
Want this function to take in epsilon

inputs: epsilon
computer: grid spacing and appropriate number of steps
then return the MFPT for this epsilon after scaling for the change in grid size


make our default grid spacing be 0.02, so calculate steps in term of the steps for this result.


"""
def walk(epsilonarr = np.ndarray, loopnum = int):
    numpts=1
    startarr = np.zeros(shape = (numpts,2)) 
    startarr += 0.5
    
    #this is np.max(ctimes)
    s_inc = 0.001
    basisvects = np.array([(0,s_inc),(0,-s_inc),(s_inc,0),(-s_inc,0)])

    #steps in out walk for now

    #boundaries for folding in
    boundaries = np.array([(0, 1), (0, 1)])
    csize = np.diff(boundaries, axis=1).ravel()

    #steps= insteps
    numpts= int((startarr.size)/2) #size counts both columns of our x,y matrix
    twodboundaries = np.tile(boundaries[:,0],(numpts,1)).T
    twodsize = np.tile(csize,(numpts,1)).T  
    #The old stuff
    #summarr=np.swapaxes(np.cumsum(basisvects[np.random.randint(0,int(basisvects.size/2),size=(numpts,steps))],axis=1),1,0) #this might have to be axis=1 now
    #print('now fold into circle')
    #trajectory_fold = (np.abs((summarr + startarr - twodboundaries.T + twodsize.T) % (2 * twodsize.T) - twodsize.T) + twodboundaries.T) #tranpose of the tranposes
    #summarr.shape

    #try at hyper-flattened stuff
    trajectory_fold = np.abs((np.swapaxes(np.cumsum(basisvects[np.random.randint(0,int(basisvects.size/2),size=(numpts,steps))],axis=1),1,0) + startarr - twodboundaries.T + twodsize.T) % (2 * twodsize.T) - twodsize.T) + twodboundaries.T
    #del(summarr)
    trajectory_fold = np.swapaxes(trajectory_fold,1,0)

    #may need to replace swapped_traj with trajectory_fold
    avg_per_eps = np.array([])
    print('now do searching')
    #changed from 'for epsilon in epsilonarr'
    for j in range(np.size(epsilonarr)):
        epsilon = epsilonarr[j]
        trickier = np.where((trajectory_fold[:,:,0].round(decimals=6) ==1.0) & (abs(trajectory_fold[:,:,1].round(decimals=6)-0.5)<=epsilon))
        #times = trickier[1]
        #locs = trickier[0]
        #del(trickier)
        #print(j)
        #j +=1
        #print(times
        #print(locs)
        eps_specific_misses =0
        inds_locs = np.unique(trickier[0], return_index=True,axis=0)[1]
        if int(np.size(inds_locs)) != int(np.size(startarr)/2):
            #make the last column the 'DID not absorb' counter
            print('Run did not absorb')
            histoarr[j,-1] +=1
            avg_per_eps = np.append(avg_per_eps,0)
            continue
        print('run was absorbed')
        times = trickier[1][inds_locs]
        del(trickier,inds_locs)
        #now we now times is just a single value
        histoarr[j,numpts*loopnum +1 : numpts*loopnum + numpts+1] = times
        #our within function averging which now doesn't matter
        #do we need atimes and avg_per_eps now?
        atimes = np.average(times.reshape(-1, numpts), axis=1)
        #del(times)
        avg_per_eps = np.append(avg_per_eps, atimes)
    del(trajectory_fold)
    
    # array([ 3847,  4171,  4174, ..., 90418, 90421, 90422
    #avgtime = np.average(avgtime.reshape(-1,epsilon_reps),axis=1)
    return avg_per_eps

#epsilon_reps = 5
#eps_repeated = np.repeat(epsilon_arr_narrow,epsilon_reps)
#create epsilon array
#range of epsilons to check
#epsilon_arr_narrow = np.linspace(0.0,0.2,11)
"""
Repeeat epsilon to get an build up our averages 

"""

#could be made to be s_inc
#fine = 0.001*np.arange(1,15)
#course = 0.02*np.linspace(1,2,5)
#xes = np.concatenate((fine,course))


#xes = 0.2*((0.5)**np.arange(0,7))
#ctimes = 2000*((0.2/xes)**2)
#epsilon_reps = 5
#eps_repeated = np.repeat(xes,epsilon_reps)
#initialize the array for averages outside function
avg_arr = np.zeros(np.size(fine))

#do ten loops to get a 1000 sample average
for i in range(numfunctcalls):
    avg_arr = avg_arr + walk(fine,i)
    print('outerloop')
    print(i)
avg_arr = avg_arr/numfunctcalls
"""
reshape to average over our epsilons and get final results
"""

#avg_arr = np.average(avg_arr.reshape(-1, epsilon_reps), axis=1)
newavg = np.array([],dtype = float)
for i in range(int(np.size(fine))):
    newavg = np.concatenate((newavg, np.array([np.average(histoarr[i,np.where((histoarr[i,1:-1]!= 0))])])))
"""
Had to change this now that we are switching things up
"""
#avg_arr = avg_arr.ravel()
#total = np.column_stack((fine,avg_arr))
newavg.tofile('C:\\Users\\h_jones\\Desktop\\w_out\\vm3000avg_finer5000_spring_breakf4.csv',sep = ',')

"""
Now lets make the histograms for the smallest four epsilson
"""
bins = steps/10
labels = bins* np.arange(11)
histarrfinal = np.empty((4,11))
for j in range(4):
    firstbin    = np.size(np.where(histoarr[j,1:-1] <=bins))
    secondbin   = np.size(np.where((histoarr[j,1:-1] >bins)     & (histoarr[j,1:-1]<= 2*bins)))
    thirdbin    = np.size(np.where((histoarr[j,1:-1] > 2*bins)  & (histoarr[j,1:-1]<= 3*bins)))
    fourthbin   = np.size(np.where((histoarr[j,1:-1] > 3*bins)  & (histoarr[j,1:-1]<= 4*bins)))
    fifthbin    = np.size(np.where((histoarr[j,1:-1] > 4*bins)  & (histoarr[j,1:-1]<= 5*bins)))
    sixthbin    = np.size(np.where((histoarr[j,1:-1] > 5*bins)  & (histoarr[j,1:-1]<= 6*bins)))
    seventhbin  = np.size(np.where((histoarr[j,1:-1] > 6*bins)  & (histoarr[j,1:-1]<= 7*bins)))
    eigthbin    = np.size(np.where((histoarr[j,1:-1] > 7*bins)  & (histoarr[j,1:-1]<= 8*bins)))
    ninthbin    = np.size(np.where((histoarr[j,1:-1] > 8*bins)  & (histoarr[j,1:-1]<= 9*bins)))
    tenthbin    = np.size(np.where((histoarr[j,1:-1] > 9*bins)  & (histoarr[j,1:-1]<= 10*bins)))
    incompletebin = histoarr[j,-1]
    histarrfinal[j] = np.array([firstbin,secondbin,thirdbin,fourthbin,fifthbin,sixthbin,seventhbin,eigthbin,ninthbin,tenthbin,incompletebin])

first = np.column_stack((labels,histarrfinal[0]))
first.tofile('C:\\Users\\h_jones\\Desktop\\w_out\\hist_0_001f4.csv',sep = ',')

second = np.column_stack((labels,histarrfinal[1]))
second.tofile('C:\\Users\\h_jones\\Desktop\\w_out\\hist_0_002f4.csv',sep = ',')

third = np.column_stack((labels,histarrfinal[2]))
third.tofile('C:\\Users\\h_jones\\Desktop\\w_out\\hist_0_003f4.csv',sep = ',')

fourth = np.column_stack((labels,histarrfinal[3]))
fourth.tofile('C:\\Users\\h_jones\\Desktop\\w_out\\hist_0_004f4.csv',sep = ',')


"""
Now save figure
"""
figbl = plt.figure(figsize = (9,9),dpi=1000)
axbl = plt.axes()
axbl.grid()

#axbl.scatter(finalxlist,finalylist, np.divide(time_array,60), c = 'r', s = 10.) #this one works
#axbl.plot_trisurf(newx,newy, avgtime)
#axbl.scatter(totalX,totalY,totaloverdue,c='k',s=10)
axbl.set_title('Epsilon and MFPT', pad=20)

# Set axes label
axbl.set_xlabel('Epsilon', labelpad=20)
axbl.set_ylabel('MFPT (steps)', labelpad=20)
#axbl.set_zlabel('time', labelpad=20)

axbl.scatter(fine, newavg)
plt.savefig('C:\\Users\\h_jones\\Desktop\\w_out\\3000calls_sping_break_finef4.png', bbox_inches = 'tight', pad_inches=0.5)
plt.close()