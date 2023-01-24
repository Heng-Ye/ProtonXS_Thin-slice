import math

nthinslices=40
Emin=0
Emax=800
thinslicewidth=(Emax-Emin)/nthinslices

list_Emin=[]
list_Emax=[]
list_id=[]

#cnt=0
for i in range (nthinslices): 
    KE_min=-thinslicewidth/2+Emax-(i+0.5)*thinslicewidth
    KE_max=thinslicewidth/2+Emax-(i+0.5)*thinslicewidth
    ##int(ceil(Emax -KEff)/ΔE)
    list_Emin.append(KE_min)
    list_Emax.append(KE_max)
    list_id.append(i)

#    pl=i%2
#    if pl==0:
#         print('\\rowcolor{LightCyan}',i,' & ', int(KE_min), '-', int(KE_max), '\\\\')
#    if pl==1:
#         print(i,' & ', int(KE_min), '-', int(KE_max), '\\\\')


#print(list_Emin)


cnt=0
string_id=""
for i in range(len(list_Emin)): 
    pl=i%2
    if pl==0:
       string_id+=str(list_id[i])+' & '
       #print(i, ' & ')
    if pl==1:     
       #print(i, ' & ')
       string_id+=str(list_id[i])+' & '

print(string_id)
'''
cnt=0
for i in len(list_Emin): 
    pl=i%2
    if pl==0:
       print(list_Emin[i], ' & ')
    if pl==1:     
       print(list_Emin[i], ' & ')
'''
