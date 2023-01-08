#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
Created by Leo Dettori on 2021.07.30.
"""

def Load_Dett(data_path):

    import sys, os
    import numpy as np

    #Switches directory to the input location, which should contain the database .dett files    
    os.chdir(data_path)

    #Formats input location if need be
    if data_path.endswith("/"):
        pass
    else:
        data_path = data_path +"/"

    #Compiles list of files within the input directory
    files = next(os.walk(data_path))[2]
    #print(files)

    #Initial message
    print("Loading data..." + "\n")

    #Creates master dictionary that will organize data from all different .dett files
    data_dett = {}

    #Creates handy counters
    idett = 0

    #Starts to loop through each .dett file
    for file in files:
        if file.endswith(".dett"):
            #Updates .dett file counter to aid with sub-dictionary creation
            idett = idett + 1
            #Updates temporary path avriable to current .dett file's path
            current_data_path = data_path+file

            #This is important for the program to know when to stop
            c = open(current_data_path)
            lines1 = c.readlines()
            total_lines = len(lines1)
            #print("\n"+"UV-Vis Results:"+"\n")
            c.close()

            # Read the database file
            f = open(current_data_path)

            #creates/resets handy counters
            i = 0   #tracks which line of dett file we're in 
            iprot = 0  #tracks current protein number

            #Creates sub-dictionary for current .dett file
            data_dett['data_'+str(idett)] = {}

            #Extracts current file name for dictionary organization purposes
            file_name_split = current_data_path.split("/")
            file_name = file_name_split[-1].split(".dett")[0]

            #Saves current file name into its respective dictionary entry for labeling purposes
            data_dett['data_'+str(idett)]['name'] = file_name
            print(data_dett['data_'+str(idett)]['name'])

            #Starts the loop trhough the file. Line by line.
            while i <= total_lines:
                this_line = f.readline()
                i = i + 1

                #Reads and saves into lists the amino acids, dipeptides and tripeptides and their respective coded/clustered versions
                #Reads and saves into a list
                if this_line.startswith("pepname"):
                    pepname = this_line.split("<")[-1].strip("\n").split(",")        
                    #print(pepname)

                #Reads and saves into a list
                if this_line.startswith("codepepname"):
                    codepepname = this_line.split("<")[-1].strip("\n").split(",")        
                    #print(codepepname)

                #Reads and saves into a list
                if this_line.startswith("dipepname"):
                    dipepname = this_line.split("<")[-1].strip("\n").split(",")        
                    #print(dipepname)

                #Reads and saves into a list
                if this_line.startswith("codedipepname"):
                    codedipepname = this_line.split("<")[-1].strip("\n").split(",")        
                    #print(codedipepname)

                #Reads and saves into a list
                if this_line.startswith("tripepname"):
                    tripepname = this_line.split("<")[-1].strip("\n").split(",")        
                    #print(tripepname)

                #Reads and saves into a list
                if this_line.startswith("codetripepname"):
                    codetripepname = this_line.split("<")[-1].strip("\n").split(",")        
                    #print(codetripepname)


                #Starts to loop through each protein's data and assemble the dictionary
                #Reads current protein/truncation/slice name and adds it to respective dictionary entry
                if this_line.startswith("name"):  #we encountered a new protein
                    
                    #updates protein counter
                    iprot = iprot + 1

                    #Creates dictionary entry for current protein
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)] = {}

                    #Prepares the sub-dictionary for the current protein
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)] = {'name':[], 'seq':[], 'pep':{}, 'codepep':{}, 'dipep':{}, 'codedipep':{}, 'tripep':{}, 'codetripep':{}}                    
                    
                    this_name = this_line.split("<")[-1].strip("\n")
                    #Adds current protein name to its dictionary entry
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["name"] = this_name
                    #print(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["name"])

                #Reads current protein/truncation/slice sequence and adds it to respective dictionary entry
                if this_line.startswith("seq"):
                    this_seq = this_line.split("<")[-1].strip("\n")
                    #Adds current protein sequence to its dictionary entry
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["seq"] = this_seq
                    #print(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["seq"])

                #For the single peptide
                #Reads current protein/truncation/slice single amino acid/peptide scores and adds it to respective dictionary entry
                if this_line.startswith("pepfreq"):
                    this_pep = this_line.split("<")[-1].strip("\n").split(",")
                    this_pep_float = []
                    #converts the score from string to float
                    for j in this_pep:
                        this_pep_float.append(float(j))
                    #Converts float list into a numpy array to facilitate math operations later    
                    this_pep_float = np.array(this_pep_float)

                    #Adds current protein score to its dictionary entry
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["pep"]["freq"] = this_pep_float
                    #Adds the respective label to the sub-dictionary
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["pep"]["name"] = pepname
                    #print(len(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["pep"]["name"]))
                    #print(len(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["pep"]["freq"]))

                #Reads current protein/truncation/slice single amino acid/peptide code/cluster scores and adds it to respective dictionary entry
                if this_line.startswith("codepepfreq"):
                    this_codepep = this_line.split("<")[-1].strip("\n").split(",")
                    this_codepep_float = []
                    #converts the score from string to float
                    for j in this_codepep:
                        this_codepep_float.append(float(j))
                    #Converts float list into a numpy array to facilitate math operations later    
                    this_codepep_float = np.array(this_codepep_float)

                    #Adds current protein score to its dictionary entry
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codepep"]["freq"] = this_codepep_float
                    #Adds the respective label to the sub-dictionary
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codepep"]["name"] = codepepname
                    #print(len(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codepep"]["name"]))
                    #print(len(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codepep"]["freq"]))


                #For the dipeptide    
                #Reads current protein/truncation/slice double amino acid/peptide scores and adds it to respective dictionary entry
                if this_line.startswith("dipepfreq"):
                    this_dipep = this_line.split("<")[-1].strip("\n").split(",")
                    this_dipep_float = []
                    #converts the score from string to float
                    for j in this_dipep:
                        this_dipep_float.append(float(j))
                    #Converts float list into a numpy array to facilitate math operations later    
                    this_dipep_float = np.array(this_dipep_float)

                    #Adds current protein score to its dictionary entry
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["dipep"]["freq"] = this_dipep_float
                    #Adds the respective label to the sub-dictionary
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["dipep"]["name"] = dipepname
                    #print(len(data_1["protein_" + str(iprot)]["dipep"]["name"]))
                    #print(len(data_1["protein_" + str(iprot)]["dipep"]["freq"]))

                #Reads current protein/truncation/slice double amino acid/peptide code/cluster scores and adds it to respective dictionary entry
                if this_line.startswith("codedipepfreq"):
                    this_codedipep = this_line.split("<")[-1].strip("\n").split(",")
                    this_codedipep_float = []
                    #converts the score from string to float
                    for j in this_codedipep:
                        this_codedipep_float.append(float(j))
                    #Converts float list into a numpy array to facilitate math operations later    
                    this_codedipep_float = np.array(this_codedipep_float)

                    #Adds current protein score to its dictionary entry
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codedipep"]["freq"] = this_codedipep_float
                    #Adds the respective label to the sub-dictionary
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codedipep"]["name"] = codedipepname
                    #print(len(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codedipep"]["name"]))
                    #print(len(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codedipep"]["freq"]))


                #For the tripeptide
                #Reads current protein/truncation/slice triple amino acid/peptide scores and adds it to respective dictionary entry
                if this_line.startswith("tripepfreq"):
                    this_tripep = this_line.split("<")[-1].strip("\n").split(",")
                    this_tripep_float = []
                    #converts the score from string to float
                    for j in this_tripep:
                        this_tripep_float.append(float(j))
                    #Converts float list into a numpy array to facilitate math operations later    
                    this_tripep_float = np.array(this_tripep_float)

                    #Adds current protein score to its dictionary entry
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["tripep"]["freq"] = this_tripep_float
                    #Adds the respective label to the sub-dictionary
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["tripep"]["name"] = tripepname
                    #print(len(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["tripep"]["name"]))
                    #print(len(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["tripep"]["freq"]))

                #Reads current protein/truncation/slice triple amino acid/peptide code/cluster scores and adds it to respective dictionary entry
                if this_line.startswith("codetripepfreq"):
                    this_codetripep = this_line.split("<")[-1].strip("\n").split(",")
                    this_codetripep_float = []
                    #converts the score from string to float
                    for j in this_codetripep:
                        this_codetripep_float.append(float(j))
                    #Converts float list into a numpy array to facilitate math operations later    
                    this_codetripep_float = np.array(this_codetripep_float)

                    #Adds current protein score to its dictionary entry
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codetripep"]["freq"] = this_codetripep_float
                    #Adds the respective label to the sub-dictionary
                    data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codetripep"]["name"] = codetripepname
                    #print(len(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codetripep"]["name"]))
                    #print(len(data_dett['data_'+str(idett)]["protein_" + str(iprot)]["codetripep"]["freq"]))

            f.close()


    print("\n")
    print("Done!")
    print("# ------------------------------------------------------------------------------------------------------------------------ #")
    print("\n")
    
    #Returns database as a dictionary
    return data_dett

