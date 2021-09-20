import sys

f1=sys.argv[1]  # input, Smarca_VHL.dat
f2=sys.argv[2]  # input, Smarca_Protac.dat
f3=sys.argv[3]  # output, CV1.dat

print(f1)
print(f2)
print(f3)

contacts=[0,0]  # want these to be the Smarca-VHL native contacts
mindist=[0,0]   # want these to be the minimum of either Smarca-VHL or Smarca-Protac min. distances
cv=[0,0]        # want these to be the combo cv


with open(f1) as f1:                            # deal with Smarca-VHL interactions
    s1=f1.readlines()
    if len(s1)>2:                               # i.e., if iter>1, thus invoked by '../westpa_scripts/runseg.sh'
        x1=s1[-2].strip().split()               # parent's values, not current one
        contacts[0]+=int(x1[1])                 # add to native contacts 
        mindist[0]=float(x1[3])                 # the min. distance
    x1=s1[-1].strip().split()                   # now the current values
    contacts[1]+=int(x1[1])                     # add to native contacts 
    mindist[1]=float(x1[3])                     # the min. distance

with open(f2) as f1:                            # deal with Smarca-Protac interactions
    s1=f1.readlines()                   
    if len(s1)>2:                               # i.e., if iter>1, thus invoked by '../westpa_scripts/runseg.sh'
        x1=s1[-2].strip().split()               # parent's values, not current one
        contacts[0]+=int(x1[1])                 # add to native contacts 
        mindist[0]=min(float(x1[3]),mindist[0]) # minimum of Smarca-VHL or Smarca-Protac min. distances
    x1=s1[-1].strip().split()                   # now the current values
    contacts[1]+=int(x1[1])                     # add to native contacts        
    mindist[1]=min(float(x1[3]),mindist[1])     # minimum of Smarca-VHL or Smarca-Protac min. distances

if len(s1)>2:                                   # i.e., if iter>1, thus invoked by '../westpa_scripts/runseg.sh'
    print("contacts",contacts[0],contacts[1])   # both current and parent values filled, and both required
    print("mindist",mindist[0],mindist[1])      # ditto
else:                                           # i.e., if iter==1, thus invoked by '../westpa_scripts/get_pcoord.sh
    print("contacts",contacts[1])               # only current value filled and required
    print("mindist",mindist[1])                 # ditto
    

if len(s1)>2:                                   # i.e., if iter>1, thus invoked by '../westpa_scripts/runseg.sh'
    if contacts[0]==0:
        cv[0]=mindist[0]                        # if parent had no contacts, the cv is the min. distance value
    else:                                       # otherwise, it is a negative number of contacts (so they can occur on same scale)
        cv[0]=-contacts[0]

if contacts[1]==0:
    cv[1]=mindist[1]                            # if current has no contacts, the cv is the min. distance value
else:
    cv[1]=-contacts[1]                          # otherwise, it is a negative number of contacts (so they can occur on same scale)

if len(s1)>2:                                   # i.e., if iter>1, thus invoked by '../westpa_scripts/runseg.sh'
    with open(f3,"w") as f:
        f.write(str(cv[0])+"\n")                # write parent and current value on different lines
        f.write(str(cv[1])+"\n")        
        print("CV1:[",cv[0],cv[1],"]")          # write to seg_logs
else:                                           # i.e., if iter==1, thus invoked by '../westpa_scripts/get_pcoord.sh'
    with open(f3,"w") as f:
        f.write(str(cv[1])+"\n")                # write only current value
        print("CV1:[",cv[1],"]")                # write to seg_logs
    


