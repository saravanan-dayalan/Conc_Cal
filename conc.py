

import csv
import numpy
import statistics
import itertools
import xlsxwriter

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


def stats(metlist, groups, workbook, text):
    worksheet = workbook.add_worksheet()
    worksheet.write(0, 0, text)
    worksheet.write(1, 0, "")
    out_count = 1
    for group in groups:
        out_count = out_count + 1
        op_text = "In group " + group
        worksheet.write(out_count, 0, op_text)
        out_count = out_count + 1
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
#                op_text = "CV of" + metlist[0][j] + ": " + str(statistics.stdev(li)/statistics.mean(li)*100)
                worksheet.write(out_count, 0, metlist[0][j])  
                worksheet.write(out_count, 1, statistics.stdev(li)/statistics.mean(li)*100)               
                #print("CV of", metlist[0][j], ": ", statistics.stdev(li)/statistics.mean(li)*100)
            except ZeroDivisionError:
                worksheet.write(out_count, 0, metlist[0][j])  
                worksheet.write(out_count, 1, "NA") 
                #print("\n Warning: sum = 0 causing a ZeroDivisionError. Therefore skipping calculating stats for this metabolite.")
                #print(li)
                #rint("\n")
            out_count = out_count + 1
            


def linreg(metlist):
    size = len(metlist)
    num_met = len(metlist[0])-3
    sum = 0
    met_count = 3
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
        #print("For ", metlist[0][j], " m = ",fit[0]," c = ",fit[1])
        conc_li.append((metlist[0][j], fit[0], fit[1]))
    return(conc_li)


def conc_cal(metlist, conc_li):
    size = len(metlist) # calculates the number of samples!
    num_met = len(metlist[0])-3
    met_count = 3
    Sub_mat = [[None]*(size-1)]
    iscol = []
    Sub_mat = numpy.array(Sub_mat)
   
    for j in range(0, num_met): # cycling through metabolites

        j = j + met_count

        Sub_li = []
        for i in range(1,size): # cycling through the number of samples
            for k in conc_li:
                m = 0
                c = 0
                x = 0
                y = 0              
                if k[0] == metlist[0][j]:                    
                    m = k[1]
                    c = k[2]
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
        init_li = [metlist[l][0],metlist[l][1], metlist[l][2]]
        temp.append(init_li + tr_Sub_mat_li[l-1])
        fin_sub_li.extend(temp)

    return(fin_sub_li)


def write_rawdata(workbook, filename):
    worksheet = workbook.add_worksheet()
    worksheet.write(0, 0, 'Raw Data')
    with open(filename, 'r') as fp:
        reader = csv.reader(fp)
        x = 1
        y = 0        
        for line in reader:     
            for j in range(0, len(line)):   
                if (x == 1):
                    worksheet.write(x, y, line[j])
                else:
                    if (y == 0 or y == 2):
                        worksheet.write(x, y, line[j])
                    else:
                        worksheet.write_number(x, y, float(line[j]))
                y = y + 1
            x = x + 1
            y = 0

                

def write_data(workbook, reader, text):
    worksheet = workbook.add_worksheet()
    worksheet.write(0, 0, text)
    x = 1
    y = 0        
    for line in reader:
        #print(line)
        for j in range(0, len(line)):   
            if (x == 1):
                worksheet.write(x, y, line[j])
            else:
                if (y == 0 or y == 2):
                    worksheet.write(x, y, line[j])
                else:
                    worksheet.write_number(x, y, float(line[j]))
            y = y + 1
        x = x + 1
        y = 0




#############################################################################################

# main()            

path = "/Users/sdayalan/Code/Proj1/Test_2/"

filename = path + "input.csv"
opfile = path + "out.xlsx"

workbook = xlsxwriter.Workbook(opfile)
write_rawdata(workbook, filename)

print("\nReading input file...")

groups = get_groups(filename)
groups = group_count(filename, groups)

# Calling the internal standard normalising method.
isnorm = is_normalise(filename)
text = "IS Normalised Data"
write_data(workbook, isnorm[0], text)

print("\n\nPerforming Internal Standard Normalisation...")

text = "CVs after Internal Standard Normalisation"
stats(isnorm[0], groups, workbook, text)

print("\n\nCalculating CVs...")

# Calling the reagent blank subtraction method.
reg_sub_li = subtract_reg(isnorm[0])
text = "IS Normalised, Reagent Blank Subtracted Data"
write_data(workbook, reg_sub_li, text)

print("\n\nPerforming Reagent Blank Subtraction...")

text = "CVs after Internal standard Normalisation and Reagent Blank Subtraction"
stats(reg_sub_li, groups, workbook, text)

print("\n\nCalculating CVs...")

# Calculating the linear regression model for each metabolite.
conc_li = linreg(reg_sub_li)

print("\n\nCalculating Linear Regression Models...")

# Calculating concentrations for metabolites
fin_conc_val = conc_cal(reg_sub_li, conc_li)

print("\n\nCalculating Concentrations...")

text = "IS Normalised, Reagent Blank Subtracted CONCENTRATION Data"
write_data(workbook, fin_conc_val, text)


text = "CVs of Concentrations after Internal standard Normalisation and Reagent Blank Subtraction"
stats(fin_conc_val, groups, workbook, text)

print("\n\nCalculating CVs...")

print("\n\nWriting all output to out.xlsx\n\n\n")

workbook.close()



