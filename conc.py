

import csv
import numpy
import statistics
import itertools

def get_groups(filename):
    groups = {}
    with open(filename, 'r') as fp:
        reader = csv.reader(fp)
        next(reader)
        
        for line in reader:
            if (line[2] in groups) == False:
                groups.update({line[2]:0})
    fp.close()
    return(groups)


def group_count(filename, groups):
    with open(filename, 'r') as fp:
        reader = csv.reader(fp)
        next(reader)
        
        for line in reader:
            for sample in groups:
                if sample == line[2]:
                    groups[sample] = groups[sample]+1
    fp.close()
    return(groups)

            
def is_normalise(filename):
    print("\n***** Normalising metabolites to Internal Standard... *****\n")
    mat = []
    with open(filename, 'r') as fp:
        reader = csv.reader(fp)
        fline = next(reader)
        
        for line in reader:
            mat.append(line[3:])
            
    arr = numpy.array(mat)
    intarr = arr.astype(numpy.int)

    isnorm = []
    isnorm = intarr[:,0]/intarr[:,0].min()
    print(isnorm)
    finalarr = numpy.array(isnorm)

    for i in range(1, len(intarr[0])):
        val = intarr[:,i]/isnorm
        finalarr = numpy.vstack((finalarr, val))
    tr_finalarr = finalarr.transpose()
    tr_finalarr_list = tr_finalarr.tolist()
    fp.close()

    with open(filename, 'r') as fp:
        reader = csv.reader(fp)
        fline = next(reader)
        fin_list = []
        fin_list.append(fline)
        i = 0
        for line in reader:
            fin_list.append(line[:3]+tr_finalarr_list[i])
            i = i+1
    fp.close()
    return(fin_list, tr_finalarr)


def subtract_reg(metlist):
    print("***** Subtracting reagent blanks... *****\n")

    size = len(metlist) # calculates the number of samples!
    num_met = len(metlist[0])-3
    met_count = 3
    Sub_mat = [[None]*(size-1)]
    iscol = []
    Sub_mat = numpy.array(Sub_mat)
   
    for j in range(0, num_met): # cycling through metabolites
        j = j + met_count

        Reg_li = []
        Reg_avg = 0
        for i in range(1,size): # cycling through the number of samples
            if metlist[i][2] == 'R':
                Reg_li.append(metlist[i][j])
        Reg_avg = statistics.mean(Reg_li)

        Sub_li = []
        for i in range(1,size): # cycling through the number of samples
            sub_val = 0
            if(metlist[i][j] != 0):
                sub_val = metlist[i][j] - Reg_avg
            Sub_li.append(sub_val)

        Sub_li_arr = numpy.array([Sub_li])     
        Sub_mat = numpy.vstack((Sub_mat, Sub_li_arr))
        tr_Sub_mat = Sub_mat[1:].transpose()
        tr_Sub_mat_li = tr_Sub_mat.tolist()

    fin_sub_li = []
    fin_sub_li.append(metlist[0])

    for l in range(1, len(metlist)):
        temp = []
        init_li = []
        init_li = [metlist[l][0],metlist[l][1], metlist[l][2]]
        temp.append(init_li + tr_Sub_mat_li[l-1])
        fin_sub_li.extend(temp)

    return(fin_sub_li)


def stats(metlist, groups):
    for group in groups:
        print("In group ",group)
        size = len(metlist)
        num_met = len(metlist[0])-3
        sum = 0
        met_count = 3
        for j in range(0, num_met): # cycling through metabolites
            j = j + met_count
            sum = 0
            li = []
            for i in range(1,size): # cycling through the number of samples in each group
                if group == metlist[i][2]:
                    sum = sum + metlist[i][j]
                    li.append(metlist[i][j])
            try:
                print("CV of", metlist[0][j], ": ", statistics.stdev(li)/statistics.mean(li)*100)
            except ZeroDivisionError:
                print("\n Warning: The sum of the following list is 0 causing a Div by 0. Therefore skipping calculating stats for this group.")
                print(li)
                print("\n")
            


def linreg(metlist):
    print("\n\n***** Calculating linear regression model for each metabolite... *****\n")
#    print(metlist)
    size = len(metlist)
    num_met = len(metlist[0])-4
    sum = 0
    met_count = 4
    conc_li = []

    for j in range(0, num_met): # cycling through metabolites
        j = j + met_count
        sum = 0
        spike_li = []
        auc_li = []        
        li = []
        for i in range(1,size): # cycling through the number of samples in each group
            if metlist[i][2] == 'S':
                if metlist[i][j] != 0:
                    spike_li.append(metlist[i][1])
                    auc_li.append(metlist[i][j])                    
        spike_li = list(map(float, spike_li))
        fit = numpy.polyfit(spike_li, auc_li, 1)
        print("For ", metlist[0][j], " m = ",fit[0]," c = ",fit[1])
        conc_li.append((metlist[0][j], fit[0], fit[1]))
    return(conc_li)


def conc_cal(metlist, conc_li):
    print("\n\n***** Calculating concentrations... *****\n")

    size = len(metlist) # calculates the number of samples!
    num_met = len(metlist[0])-4
    met_count = 4
    Sub_mat = [[None]*(size-1)]
    iscol = []
    Sub_mat = numpy.array(Sub_mat)
   
    for j in range(0, num_met): # cycling through metabolites
        m = 0
        c = 0
        x = 0
        y = 0
        j = j + met_count

        Sub_li = []
        for i in range(1,size): # cycling through the number of samples
            for k in conc_li:
#                print(k, metlist[0][j])
                if k[0] == metlist[0][j]:                    
                    m = k[1]
                    c = k[2]
#                    print(m, c)
                    break
                
            y = metlist[i][j]

            if(metlist[i][j] != 0):
                x = (y-c)/m
            Sub_li.append(x)

        Sub_li_arr = numpy.array([Sub_li])     
        Sub_mat = numpy.vstack((Sub_mat, Sub_li_arr))
        tr_Sub_mat = Sub_mat[1:].transpose()
        tr_Sub_mat_li = tr_Sub_mat.tolist()

    fin_sub_li = []
    fin_sub_li.append(metlist[0])

    for l in range(1, len(metlist)):
        temp = []
        init_li = []
        init_li = [metlist[l][0],metlist[l][1], metlist[l][2], '']
        temp.append(init_li + tr_Sub_mat_li[l-1])
        fin_sub_li.extend(temp)

    return(fin_sub_li)


    
#############################################################################################

# main()            

filename = "input.csv"

groups = get_groups(filename)
groups = group_count(filename, groups)

# Calling the internal standard normalising method.
isnorm = is_normalise(filename)
for i in isnorm[0]:
    print(i)
print("\n")

print("\n CVs after internal standard normalisation")
stats(isnorm[0], groups)


# Calling the reagent blank subtraction method.
reg_sub_li = subtract_reg(isnorm[0])
for i in reg_sub_li:
    print(i)
print("\n")

# Calculating CVs of each metabolite in each group.
print("***** Reporting CVs of Peak Areas... *****\n")
print("Data has been IS normalised, Reagent Blank subtracted\n")
stats(reg_sub_li, groups)
'''
# Calculating the linear regression model for each metabolite.
conc_li = linreg(reg_sub_li)

# Calculating concentrations for metabolites
fin_conc_val = conc_cal(reg_sub_li, conc_li)
for i in fin_conc_val:
    print(i)
print("\n")

# Calculating CVs of each metabolite in each group.
print("***** Reporting CVs of Concentrations... *****\n")
print("Data has been IS normalised, Reagent Blank subtracted\n")
stats(fin_conc_val, groups)

'''
