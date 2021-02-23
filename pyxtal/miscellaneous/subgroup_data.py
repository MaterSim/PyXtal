from pyxtal.symmetry import Group
    
length1=0
length2=0
length3=0
length_greater=0
    
for g in range(1,231): 
    wyc = Group(g).get_max_t_subgroup()
    for re, h in zip(wyc['relations'], wyc['subgroup']):
        for pos in re:
            length=len(pos)
            if length==1:
                length1+=1
            
            elif length==2:
                length2+=1
            
            elif length==3:
                length3+=1
                #print(g, h, pos)   
            else:
                length_greater+=1
                print(g, h, pos)   

