#Extracts IDR sequences from IUPRED results into individual and combined .txt files (FASTA format) in addition to a summary table.txt
#This version works on a group of IUPRED_results that are the sole files inside a single folder

"""
Created by Leo Dettori on 2021.02.12.
"""

#Function is defined
def Extract_IDRs(IUPRED_results_folder, tolerance, min_IDR_length, threshold, output = "input"):
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from subprocess import Popen, PIPE
    import os
    import numpy as np
    import re
    from math import log10, floor
    
    #Starting message:
    print("Extract_IDRs has started:\n\n")
    
    #In case output was not defined
    if output == "input":
        output = IUPRED_results_folder
        print("No output defined! Output was set to same destination as input!" + "\n")
        
    #Prepares to create a master folder to store all the results
    if output.endswith('/'):
        pass
    else:
        output = output + "/"
        
    #Updates output location
    output = output+"Extract_IDRs_Results"+'__len_'+str(min_IDR_length)+'_thr_'+str(threshold)+'_tol_'+str(tolerance)+"/"
    #Creates the master folder to store all the results
    os.mkdir(output)
        
    #Creates a list for all the files in the specified folder
    os.chdir(IUPRED_results_folder)
    files = [f for f in os.listdir('.') if f.endswith(".txt") | f.endswith(".result")]
    #print(files)

    #Creates Summary.txt to summarize run
    if output.endswith('/'):
        summary_output = output
        f4 = open(str(summary_output)+"Summary.txt", "a")
    else:
        summary_output = output + "/"
        f4 = open(str(summary_output)+"Summary.txt", "a")
    print("Summary.txt was successfully created!\n")
    print("# --------------------------------------------- #")
    f4.close()
   
    #Starts the loop through each IUPRED.result file
    for f in files:
        #Extracts the name of the current group
        this_name = f.split(".")[0]
        #Message for the user:
        print("\n"+str(this_name)+':'+'  min. IDR length = '+str(min_IDR_length)+'  threshold = '+str(threshold)+'  tolerance = '+str(tolerance)+"\n")
        print("Running...")
        #Creates a folder within the master folder for the current group's results to be stored in an orderly manner
        if output.endswith('/'):
            this_output = output+str(this_name)+'__len_'+str(min_IDR_length)+'_thr_'+str(threshold)+'_tol_'+str(tolerance)
            os.mkdir(this_output)

        else:
            this_output = output+ '/' +str(this_name)+'__len_'+str(min_IDR_length)+'_thr_'+str(threshold)+'_tol_'+str(tolerance)
            os.mkdir(this_output)

        #Creates the 'IDR' subfolder within the recently created folder to store the fasta sequences in an orderly manner
        this_output_IDRs = this_output+ '/' +'IDRs'
        os.mkdir(this_output_IDRs)
        
        #Creates Exceptions dictionary to keep track of exceptional proteins
        #So far, it consists of proteins whose size is smaller than the specified tolerance
        excepts = {}
        excepts['name'] = []
        excepts['ID'] = []

################################################################################################################################
################################################################################################################################
################################################################################################################################
#                                            Start of extract-iupred2A chunk                                                   #
################################################################################################################################
################################################################################################################################
################################################################################################################################

        #This chunk loads the IUPRED results from the txt file into a dictionaries 

        #This is important for the program to know when to stop
        c = open(f)
        lines1 = c.readlines()
        total_lines = len(lines1)
        #print("\n"+"iupred2A Results:"+"\n")
        #print("Running...")
        c.close()

        #Read the data
        f = open(f)

        #creates dictionary
        proteins_iupred2a = {}
        #print(proteins_iupred2a)


        #creates counters
        i = 1
        n = 0
        ndis = 0    


        #Starts the loop trhough the file. Line by line.
        while i <= total_lines:
            this_line = f.readline()
            i = i + 1


            #Checks if a new protein starts in this line. If it does, the program creates a new entry in the dictionary for this protein
            #Next, the program extracts the line of the protein name and adds it to the dictionary under 'name'
            if '>' in this_line:

                ################################## ################################################################################
                #If the 1st protein has already finished, the program estimates the %disorder for each protein, one by one, and saves it
                #But the final protein has to be proccessed outside of this loop
                if n != 0:
                    this_dis = 10000 * (ndis/this_resN)
                    this_dis = floor(this_dis)/100                          
                    #print("This protein is "+str(this_dis)+"% disordered")
                    proteins_iupred2a["protein"+str(n)]['dis'].append(this_dis)

                    ndis = 0
                    pass
                ##################################################################################################################

                n = n + 1
                proteins_iupred2a["protein"+str(n)] = {'name':[], 'ID':[], 'long_name':[], 'resN':[], 'resname':[], 'iupred':[], 'anchor':[], 'dis':[]}
                this_line= this_line.strip("\n")
                proteins_iupred2a["protein"+str(n)]['long_name'].append(this_line)

                #extracts Uniprot ID and name from protein's long name line
                #removes > and starts splitting the long name line
                this_line = this_line.strip(">")
                this_line = this_line.split()
                this_line_name = str(this_line[0])  
                #print(this_line_name)
                split_this_line_name = this_line_name.split("|")
                proteins_iupred2a["protein"+str(n)]['ID'].append(split_this_line_name[1])
                proteins_iupred2a["protein"+str(n)]['name'].append(split_this_line_name[2])
                #print(split_this_line_name[1])

                #print(proteins_iupred2a)

            #Checks if this line begins with # and discards it
            elif this_line.startswith('#'):
                pass


              #Checks if this line begins with "enter" and discards it. In other words, checks if it is an empty line
            elif this_line.startswith('\n'):
                pass

            #Checks if this line has information on a protein. If it does, the program splits the line into a list
            #Next, the program saves each item of the list into a temporary var and uses the vars to save them into the dictionary
            else:

                split_line = this_line.split()

                #For debug purposes
                #print(split_line)

                this_resN = int(split_line[0])
                this_resname = str(split_line[1])
                this_iupred = float(split_line[2])
                this_anchor = float(split_line[3])

                if this_iupred >= 0.5:
                    ndis = ndis + 1

                #print(this_resN)
                #print(this_resname)
                #print(this_line)

                proteins_iupred2a["protein"+str(n)]['resN'].append(this_resN)
                proteins_iupred2a["protein"+str(n)]['resname'].append(this_resname)
                proteins_iupred2a["protein"+str(n)]['iupred'].append(this_iupred)
                proteins_iupred2a["protein"+str(n)]['anchor'].append(this_anchor)


                #Plotting the data of current protein
                #plt.plot(proteins_iupred2a["protein"+str(n)]['resN'],proteins_iupred2a["protein"+str(n)]['iupred'], 'ro', label="IuPred Score")
                #plt.plot(proteins_iupred2a["protein"+str(n)]['resN'],proteins_iupred2a["protein"+str(n)]['anchor'], 'bo', label="Anchor Score")
                #anchor score won't be necessary

                # Decorating the plot
                #plt.xlabel('Residue Number')
                #plt.ylabel('Residue Score')


                # Legend
                #legend = plt.legend(bbox_to_anchor=(1, 1), loc=1, borderaxespad=1.)
                #legend.get_frame().set_edgecolor('grey')
                #legend.get_frame().set_linewidth(2.0)

        #The program estimates the %disorder for the last protein and saves it
        this_dis = 10000 * (ndis/this_resN)
        this_dis = floor(this_dis)/100                          
        #print("This protein is "+str(this_dis)+"% disordered")
        proteins_iupred2a["protein"+str(n)]['dis'].append(this_dis)

        #print(proteins_iupred2a["protein"+str(n)]["long_name"])
        #plt.show()
        print("Total number of lines in the file is: "+str(total_lines))
        print("Total number of proteins in the file is: "+str(n))


        #For debug purposes
        #print(proteins_iupred2a)        


        #f.close()        

################################################################################################################################
################################################################################################################################
################################################################################################################################
#                                                     end of iupred2A chunk                                                    #
################################################################################################################################
################################################################################################################################
################################################################################################################################

        #Creates the text file to store the compiled multi-FASTA results
        #f2 = open(str(this_output)+'/'+"IDRs_multi_FASTA.txt", "a")
        #f2.close()

        #Creates the text file to save the results table
        f3 = open(str(this_output)+'/'+"#Results_table.txt", "a")
        f3.write("ID"+' '+"Gene_Name"+' '+"FL_Protein_Length"+' '+"FL_Protein_Disorder_IUPred(%)"+' '+"IDR_Number"+' '+"IDRs Length Start End Confidence(%) Avg_Disorder(%)"+' '+"Long_Name"+" Tolerance:"+str(tolerance)+" Cut-off:"+str(min_IDR_length)+" Threshold:"+str(threshold)+"\n")
        f3.close()

        #Creates handy counters
        p = 0      #marks the current protein dictionary in the dictionaries
        n2 = 0      #marks the current residue in the dictionary
        I = 0      #marks the current IDR
        C = 0       #marks current IDR when applying cut_off
        I2 = 0     #marks number of IDRs that passed the cut_off test
        L = 0       #marks number of linkers in the current protein
        IT = 0      #marks total number of fasta files created for IDRs
        proteins_w_IDRs = 0 #marks number of proteins with IDRs
        avg_IDR_dis = 0 #marks average IDR confidence for current protein only considering the IDRs that pass the cut-off
        group_avg_IDR_dis = 0 #marks average IDR confidence for all the IDRs within the group, average of the average
        percentage_IDRs = 0  #marks percentage of proteins with at least one IDR
        
        #Creates decision making variables:
        previous = "none"  #information regarding previous residue to decide if we're in an IDR or in a folded region (fold)
        res = 0  #current residue's index

        #Switches the names of the Parameters established above to shorter names
        tol = tolerance
        cut_off = min_IDR_length
        T = threshold

        #creates master dictionary
        proteins_IDRs = {}



        #Starts to loop through the IUPRED dictionary protein by protein (the first protein is protein1 according to indexing)

        while p != n: #n is the total number of proteins from iupred2a chunk above
            #updates counter n2
            p = p + 1
            #resets residue counter n2, IDR counter, C, I2, and decision making variable previous
            n2 = 0
            I = 0
            previous = "none"
            C = 0    #resets handy cut_off counter for current IDR
            I2 = 0    #resets the handy counter for number of IDRs that passed the cut_off test
            L = 0     #resets the handy Linker counter
            avg_IDR_dis = 0 #resets average IDR confidence for current protein

            #Creates new dicionary entries for IDRs and transfers previous values
            proteins_IDRs["protein"+str(p)] = {'name':[], 'ID':[], 'long_name':[], 'resN':[], 'resname':[], 'iupred':[], 'anchor':[], 'dis':[], 'IDRs':{}, 'IDRs_amount':[]}

            proteins_IDRs["protein"+str(p)]['ID'] = proteins_iupred2a["protein"+str(p)]['ID']
            proteins_IDRs["protein"+str(p)]['name'] = proteins_iupred2a["protein"+str(p)]['name']
            proteins_IDRs["protein"+str(p)]['long_name'] = proteins_iupred2a["protein"+str(p)]['long_name']
            proteins_IDRs["protein"+str(p)]['iupred'] = proteins_iupred2a["protein"+str(p)]['iupred']
            proteins_IDRs["protein"+str(p)]['anchor'] = proteins_iupred2a["protein"+str(p)]['anchor']
            proteins_IDRs["protein"+str(p)]['resN'] = proteins_iupred2a["protein"+str(p)]['resN']
            proteins_IDRs["protein"+str(p)]['resname'] = proteins_iupred2a["protein"+str(p)]['resname']
            proteins_IDRs["protein"+str(p)]['dis'] = proteins_iupred2a["protein"+str(p)]['dis']

            #Starts to loop through the protein residue by residue (the first residue is res[0] according to indexing)
            #This part of the program has three distinct stages for proper decision making (begginning, middle, end)

            #Gets protein length to help with the loop and to define the distinct stages
            protein_length = len(proteins_IDRs["protein"+str(p)]['resN'])
            #print(proteins_IDRs["protein"+str(p)]['name'][0])

            #starts residue by residue loop
            while n2 <= protein_length - 1:

                #Beginning: (For the first residue)
                
                #Handling exception for when protein length is less than the specified tolerance
                if n2 == 0 and n2+tol > protein_length - 1:
                    #This protein is too short to provide an accurate analysis, so it will be skipped and a message we be printed out
                    print(" Sorry!\n "+str(proteins_IDRs["protein"+str(p)]['name'])+ " of Uniprot ID: "+str(proteins_IDRs["protein"+str(p)]['ID'])+" is too small to provide an accurate analysis.\n Protein size is smaller than the specified tolerance: "+ str(tol)+". Therefore, no IDRs will be extracted from this protein.")
                    #Saves the name and Uniprot ID of this exceptional protein into dictionary excepts
                    excepts['name'].append(proteins_IDRs["protein"+str(p)]['name'])
                    excepts['ID'].append(proteins_IDRs["protein"+str(p)]['ID'])
                    
                    #Updates n2 to escape the residue loop
                    n2 = protein_length
                    pass
                                
                #Back to the Beggining part of the program on a normal occasion: (For the first residue)
                if n2 == 0 and n2+tol <= protein_length - 1:
                    #Extracts iupred scores for the 2 residues taken into account for decision making in this iteration of the loop:
                    #current residue (res)
                    res = proteins_IDRs["protein"+str(p)]['iupred'][n2]
                    #incoming residue with respect to the established tolerance (res+tol)
                    res_tol = proteins_IDRs["protein"+str(p)]['iupred'][n2+tol]
                    #print(res)
                    #print(res_tol)
                    #print(' ')

                    #Starts decision making proccess to say if we're in an IDR or in a folded region (fold)
                    if res >= T and res_tol >= T:
                        previous = "IDR"   #for the next iteration, we are in an IDR
                        #Generates dictionary entries to store the first IDR information
                        I = I + 1
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)] = {'IDR_index':[], 'IDR_type':[], 'IDR_length':[], 'avg_IDR_score':[], 'IDR_seq':[], 'IDR_dis':[], 'IDR_start':[], 'IDR_end':[], 'IDR_resN':[], 'IDR_iupred':[]}
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_index'] = I
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] = proteins_IDRs["protein"+str(p)]['iupred'][n2]
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] = 1   #above threshold
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] = proteins_IDRs["protein"+str(p)]['resname'][n2]
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_type'] = "N_IDR"
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_start'] = proteins_IDRs["protein"+str(p)]['resN'][n2]

                        #print(proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'])
                        #print(previous)

                    if res >= T and res_tol < T:
                        previous = "IDR"   #for the next iteration, we are in an IDR
                        #Generates dictionary entries to store the first IDR information
                        I = I + 1
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)] = {'IDR_index':[], 'IDR_type':[], 'IDR_length':[], 'avg_IDR_score':[], 'IDR_seq':[], 'IDR_dis':[], 'IDR_start':[], 'IDR_end':[], 'IDR_resN':[], 'IDR_iupred':[]}
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_index'] = I
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] = proteins_IDRs["protein"+str(p)]['iupred'][n2]
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] = 1   #above threshold
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] = proteins_IDRs["protein"+str(p)]['resname'][n2]
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_type'] = "N_IDR"
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_start'] = proteins_IDRs["protein"+str(p)]['resN'][n2]

                        #print(proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'])
                        #print(previous)

                    if res < T and res_tol >= T:
                        previous = "IDR"   #for the next iteration, we are in an IDR
                        #Generates dictionary entries to store the first IDR information
                        I = I + 1
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)] = {'IDR_index':[], 'IDR_type':[], 'IDR_length':[], 'avg_IDR_score':[], 'IDR_seq':[], 'IDR_dis':[], 'IDR_start':[], 'IDR_end':[], 'IDR_resN':[], 'IDR_iupred':[]}
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_index'] = I
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] = proteins_IDRs["protein"+str(p)]['iupred'][n2]
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] = 0   #below threshold
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] = proteins_IDRs["protein"+str(p)]['resname'][n2]
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_type'] = "N_IDR"
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_start'] = proteins_IDRs["protein"+str(p)]['resN'][n2]

                        #print(proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'])
                        #print(previous)

                    if res < T and res_tol < T:
                        previous = "fold"  #for the next iteration, we are in a fold
                        #print(previous)

                    #Moves to the next residue
                    #print(n2)
                    n2 = n2 + 1


                #Middle: (For the central residues. The decision making will differ depending on whether we're in an IDR or in a fold)
                if n2 > 0 and n2+tol <= protein_length - 1:
                    #Extracts iupred scores for the 2 residues taken into account for decision making in this iteration of the loop:
                    #current residue (res)
                    res = proteins_IDRs["protein"+str(p)]['iupred'][n2]
                    #incoming residue with respect to the established tolerance (res+tol)
                    res_tol = proteins_IDRs["protein"+str(p)]['iupred'][n2+tol]

                    #print(proteins_IDRs["protein"+str(p)]['resN'][n2_tol])
                    #print(protein_length)
                    #print(res)
                    #print(res_tol)
                    #print(' ')

                    #Checks whether we're in an IDR:
                    if previous == "IDR":              

                        #Starts decision making proccess when we're in an IDR. This is biased towards writing the residue into an IDR sequence text
                        #This step uses the tolerance value to help minize mistakes due to outliers and "gaps" from the IUPRED results
                        #The cut_off stage at the end will help minize mistakes.
                        if res >= T and res_tol >= T: #we're very likely in a true IDR
                            previous = "IDR"   #for the next iteration, we are still in an IDR
                            #Adds IDR information to existing IDR entry
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] + proteins_IDRs["protein"+str(p)]['iupred'][n2]
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] + 1   #above threshold
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] + proteins_IDRs["protein"+str(p)]['resname'][n2]
                            #print(proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'])
                            #print(previous)

                        if res >= T and res_tol < T: #we're very likely in a true IDR, but approaching either a gap/outlier or the boundary
                            previous = "IDR"   #for the next iteration, we are still in an IDR
                            #Adds IDR information to existing IDR entry
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] + proteins_IDRs["protein"+str(p)]['iupred'][n2]
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] + 1   #above threshold
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] + proteins_IDRs["protein"+str(p)]['resname'][n2]
                            #print(previous)

                        if res < T and res_tol >= T: #we're very likely going through a gap/outlier of a true IDR, but approaching the gap/oulier's end
                            previous = "IDR"   #for the next iteration, we are still in an IDR
                            #Adds IDR information to existing IDR entry
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] + proteins_IDRs["protein"+str(p)]['iupred'][n2]
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] + 0   #below threshold
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] + proteins_IDRs["protein"+str(p)]['resname'][n2]
                            #print(previous)

                        if res < T and res_tol < T: #we're very likely in a fold now or going through a gap/outlier greater than the tolerance can handle
                            previous = "fold"  #for the next iteration, we are in a fold
                            #Wraps up this IDR of index I
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_length'] = len(proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'])
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_end'] = proteins_IDRs["protein"+str(p)]['resN'][n2-1]
                            #Finishing up calculation for average and confidence:
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] = floor(10000*proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score']/proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_length'])/100
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] = floor(10000*proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis']/proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_length'])/100
                            #print(previous)


                    #Checks whether we're in a fold:
                    if previous == "fold":              

                        #Starts decision making proccess when we're in a fold. This is biased towards not writing the residue into an IDR sequence text
                        #This step uses the tolerance value to help minize mistakes due to outliers and "gaps" from the IUPRED results
                        #The cut_off stage at the end will help minize mistakes.
                        if res >= T and res_tol >= T: #we're very likely no longer in a fold anymore, we very likely entered an IDR or are going through a gap/outlier greater than the tolerance can handle
                            previous = "IDR"   #for the next iteration, we are in an IDR
                            #Generates dictionary entries to store a new IDR's information
                            I = I + 1
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)] = {'IDR_index':[], 'IDR_type':[], 'IDR_length':[], 'avg_IDR_score':[], 'IDR_seq':[], 'IDR_dis':[], 'IDR_start':[], 'IDR_end':[], 'IDR_resN':[], 'IDR_iupred':[]}
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_index'] = I
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] = proteins_IDRs["protein"+str(p)]['iupred'][n2]
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] = 1   #above threshold
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] = proteins_IDRs["protein"+str(p)]['resname'][n2]
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_type'] = "Linker"
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_start'] = proteins_IDRs["protein"+str(p)]['resN'][n2]
                            #print(previous)

                        if res >= T and res_tol < T: #we're very likely going trhough a gap/outlier of a true fold, but approaching the gap/oulier's end
                            previous = "fold"   #for the next iteration, we are in a fold
                            #print(previous)

                        if res < T and res_tol >= T: #we're very likely in a true fold, but approaching either a gap/outlier or the boundary
                            previous = "fold"   #for the next iteration, we are in a fold
                            #print(previous)

                        if res < T and res_tol < T: #we're very likely in a true fold
                            previous = "fold"  #for the next iteration, we are in a fold
                            #print(previous)

                    #Moves to the next residue
                    #print(n2)
                    n2 = n2 + 1


                #End: (For the final residues. This region corresponds to the size of the tolerance)
                if n2 > 0 and n2+tol > protein_length - 1 and previous != "none":  
                # (and previous != "none") handles the exception of protein length being smaller than tolerance
                    #Extracts iupred scores for the current residue for decision making in this iteration of the loop:
                    #current residue (res)
                    res = proteins_IDRs["protein"+str(p)]['iupred'][n2]
                    #print(res)
                    #print(res_tol)
                    #print(' ')

                    #Starts decision making depending on whether we're in an IDR or in a fold. All the final residues will be handled in the same manner
                    if previous == "IDR":
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] + proteins_IDRs["protein"+str(p)]['iupred'][n2]
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_type'] = "C_IDR"   #the last IDR is a C-IDR
                        proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'] + proteins_IDRs["protein"+str(p)]['resname'][n2]

                        if res >= T:
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] + 1   #above threshold
                        else:
                            proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] = proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] + 0   #below threshold                        

                    if previous == "fold":    
                        #print(previous)
                        pass

                    #Moves to the next residue
                    #print(n2)
                    n2 = n2 + 1


            #Wrapping up this protein and its last IDR:
            
            #Handling exception for when protein length is less than the specified tolerance- this might not be necessary, but just in case
            if previous == "none":
                pass
                        
            #Back to wrapping up this protein and its last IDR:
            if previous == "IDR":            
                proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_length'] = len(proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_seq'])
                proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_end'] = proteins_IDRs["protein"+str(p)]['resN'][n2-1]
                #Finishing up calculation for average and confidence:
                proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score'] = floor(10000*proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['avg_IDR_score']/proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_length'])/100
                proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis'] = floor(10000*proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_dis']/proteins_IDRs["protein"+str(p)]["IDR"+str(I)]['IDR_length'])/100


            #Applying cut_off accoding to IDR minimum length parameter and exporting results into txt files and plot figures


            #Creating the general IUPRED plot
            fig, (ax1) = plt.subplots(1, sharex=True, figsize=(15,15))
            ax1.plot(proteins_IDRs["protein"+str(p)]['resN'],proteins_IDRs["protein"+str(p)]['iupred'], label="IUPred Score")

            while C != I:
                C = C + 1

                if proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_length'] >= cut_off:
                    I2 = I2 + 1

                    #Adds one more IDR to the count of average IDR confidence for current protein only considering IDRs that pass the cut-off
                    avg_IDR_dis = avg_IDR_dis + proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_dis']

                    #Counts the number of proteins with at least one IDR
                    if I2 == 1:
                        proteins_w_IDRs = proteins_w_IDRs + 1 #Adds one more protein to the count of proteins with IDRs

                    #Names the linkers properly according to number of linkers that passed the cut_off test
                    if proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_type'] == "Linker":
                        L = L + 1
                        proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_type'] = proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_type']+"_"+str(L)

                    #writes this protein's IDR sequences on a fasta.txt file within the subfolder IDRs
                    f2 = open(str(this_output_IDRs)+'/'+str(proteins_IDRs["protein"+str(p)]['name'][0])+"_"+str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_start'])+'_'+str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_end'])+'_'+str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_type'])+".txt", "a")
                    f2.write(">"+proteins_IDRs["protein"+str(p)]['ID'][0]+"|"+proteins_IDRs["protein"+str(p)]['name'][0]+"|"+"Confidence: "+str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_dis'])+'%'+"|"+str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_start'])+'_'+str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_end'])+"|"+proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_type']+'\n')
                    f2.write(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_seq']+'\n')
                    f2.write('\n')       
                    f2.close()

                    #writes this protein's IDR information on the summary results.txt file
                    f3 = open(str(this_output)+'/'+"#Results_table.txt", "a")
                    f3.write(str(proteins_IDRs["protein"+str(p)]['ID'][0])+' '+str(proteins_IDRs["protein"+str(p)]['name'][0])+' '+str(len(proteins_IDRs["protein"+str(p)]['resN']))+' '+str(proteins_IDRs["protein"+str(p)]['dis'][0])+' '+str(I2)+' '+str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_type'])+' '+ str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_length'])+' '+str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_start'])+' '+str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_end'])+' '+ str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_dis'])+' '+ str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['avg_IDR_score'])+' '+str(proteins_IDRs["protein"+str(p)]['long_name'][0])+"\n")
                    f3.close()

                    #Plotting data for each IDR
                    #fig, (ax1) = plt.subplots(1, sharex=True, figsize=(15,15))
                    ax1.axhline(y=T,alpha=0.7,c='r',linestyle='--', label="IUPred Threshold")
                    plt.axvline(x=proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_start'],alpha=0.7,c='g',linestyle='--', label="IUPred Threshold")
                    plt.axvline(x=proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_end'],alpha=0.7,c='g',linestyle='--', label="IUPred Threshold")
                    ax1.axvspan(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_start'], proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_end'], alpha=0.5, color='green')




                    #Saving the plot
                    #fig.savefig(this_output+'/'+str(proteins_IDRs["protein"+str(p)]['name'][0])+"_"+str(proteins_IDRs["protein"+str(p)]["IDR"+str(C)]['IDR_type'])+'.png')
                    #plt.close()

            #Wraps up the calculation of average IDR confidence for the current protein only considering IDRs that pass the cut-off
            if avg_IDR_dis != 0:
                avg_IDR_dis = floor(100*avg_IDR_dis/I2)/100
                #Adds this IDR to the calculation of the group average confidence only considering IDRs that pass the cut-off
                group_avg_IDR_dis = group_avg_IDR_dis + avg_IDR_dis
                #Creates plot title when IDRs were extracted
                ax1.set_title('IUPred2A Plot - '+str(proteins_IDRs["protein"+str(p)]['ID'][0])+' | '+ str(proteins_IDRs["protein"+str(p)]['name'][0])+' | Estimated Disorder: '+str(proteins_IDRs["protein"+str(p)]['dis']).strip(">'[]'")+'%'+' | #IDRs: '+str(I2)+'\n'+'Min. IDR Length: '+str(min_IDR_length)+' | Threshold: '+str(threshold)+' | Tolerance: '+str(tolerance)+' | Protein/IDRs Average Confidence: '+str(avg_IDR_dis)+'%', fontsize=18)
            if avg_IDR_dis == 0:
                avg_IDR_dis = "Not Apply"
                #Creates plot title when no IDRs were extracted
                ax1.set_title('IUPred2A Plot - '+str(proteins_IDRs["protein"+str(p)]['ID'][0])+' | '+ str(proteins_IDRs["protein"+str(p)]['name'][0])+' | Estimated Disorder: '+str(proteins_IDRs["protein"+str(p)]['dis']).strip(">'[]'")+'%'+' | #IDRs: '+str(I2)+'\n'+'Min. IDR Length: '+str(min_IDR_length)+' | Threshold: '+str(threshold)+' | Tolerance: '+str(tolerance)+' | Protein/IDRs Average Confidence: '+str(avg_IDR_dis), fontsize=18)

            #Saving and exporting the plot        
            ax1.set_ylabel('IUPred Score', fontsize=18)
            fig.savefig(this_output+'/'+str(proteins_IDRs["protein"+str(p)]['name'][0])+'.png')
            plt.close()

            #Saving the protein's number of IDRs into the protein's dictionary and into total count of fasta files created for IDRs
            proteins_IDRs["protein"+str(p)]['IDRs_amount'] = I2
            IT = IT + I2


        #Wraps up the calculation of average IDR confidence for the group only considering IDRs that pass the cut-off
        if group_avg_IDR_dis != 0:
            #This is an average of the averages type of calculation, thus has potential to be improved
            group_avg_IDR_dis = floor(100*group_avg_IDR_dis/proteins_w_IDRs)/100       


        #Calculates percentage of proteins with at least one IDR in the file
        if proteins_w_IDRs != 0:
            percentage_IDRs = floor(10000*proteins_w_IDRs/n)/100   

        #Saves final messages into Summary.txt to summarize run    
        f4 = open(str(summary_output)+"Summary.txt", "a")
        if group_avg_IDR_dis != 0:
            f4.write(str(this_name)+':'+'  min. IDR length = '+str(min_IDR_length)+'  threshold = '+str(threshold)+'  tolerance = '+str(tolerance)+"\n"+"Total number of lines in the file is: "+str(total_lines)+"\n"+"Total number of proteins in the file is: "+str(n)+"\n"+"Total number of IDR fasta files created is: "+str(IT)+"\n"+"Total number of proteins with at leats one IDR is: "+str(proteins_w_IDRs)+"\n"+"Percentage of proteins with at leats one IDR is: "+str(percentage_IDRs)+"%"+"\n"+"Average confidence for the entire group is: "+str(group_avg_IDR_dis)+"%"+"\n"+"Proteins whose length is less than specified tolerance: "+str(excepts['name'])+"\n"+"# --------------------------------------------- #"+"\n")
        else:
            f4.write(str(this_name)+':'+'  min. IDR length = '+str(min_IDR_length)+'  threshold = '+str(threshold)+'  tolerance = '+str(tolerance)+"\n"+"Total number of lines in the file is: "+str(total_lines)+"\n"+"Total number of proteins in the file is: "+str(n)+"\n"+"Total number of IDR fasta files created is: "+str(IT)+"\n"+"Total number of proteins with at leats one IDR is: "+str(proteins_w_IDRs)+"\n"+"# --------------------------------------------- #"+"\n")        
        f4.close()


        #Final messages to summarize run
        print("Plots were successfully created!")
        print("Results.txt table was successfully created!")
        print("Total number of IDR fasta files created is: "+str(IT))
        print("Total number of proteins with at leats one IDR is: "+str(proteins_w_IDRs))
        print("Percentage of proteins with at leats one IDR is: "+str(percentage_IDRs)+"%")
        if group_avg_IDR_dis != 0:
            print("Average confidence for the entire group is: "+str(group_avg_IDR_dis)+"%")
        print("Done"+"\n"+"\n"+"# --------------------------------------------- #")
        
    
    #Ending message
    print("\nExtract_IDRs is done!\n")
    print("# ------------------------------------------------------------------------------------------------------------------------ #")
    print("\n")    
    
    #This is what this function returns to be subsequently used in the List Generator Module/Function
    return summary_output

