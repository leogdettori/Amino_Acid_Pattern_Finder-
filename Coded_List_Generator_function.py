#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
Adapted on Mon Feb 23 2021

Script for preparing for the program that finds and plots the frequency of different dipeptide combinations
in selected protein sequences

@authors: Csaba Papp, Leo Dettori

"""

#Run this script to generatae a list file for Csabas Program

################################################################################################################################

#Specify Group Code and Groups:
#Make sure you run the Pre Set Tab in case you want to use Pre Set Groups

#group_code = Cluster_Code_5
#groups = Cluster_5

#All Amino Acids: ['I','V','L','F','C','M','A','G','T','W','S','Y','P','H','N','Q', 'E','D','K','R']
#Uncalled Amino Acids will remain as they are instead of being grouped



#Specify the path to the folde where the protein sequences are (either .fasta or .txt files):
#Don't delete the r in r'.........'

#seq_path = r'C:/Users/leogd/New_Readers_List_IDRs_40/IDRs'

#Specify the name of the output list file:

#list_name = "New_Readers_List_IDRs_40_Cluster_5"

#Specify the output address (standard is the destination where Csabas program runs):

#output_path = r'C:/Users/leogd/CsabaPapp_RotationProject/Sequence_Lists'

################################################################################################################################


import sys, os

def Coded_List_Generator(groups, group_code, seq_path, list_name, output_path):

    os.chdir(seq_path)
    files = [f for f in os.listdir('.') if f.endswith(".txt") | f.endswith(".fasta")]
    #print(files)

    cwd = os.getcwd()


    os.chdir(output_path)
    file = open(list_name+'.txt',"w") 
    counter = 0

    for f in files:
        counter = counter + 1
        name = f.split(".")[0]
        location = cwd + "\\" + f
        l = [name, location, str(group_code), str(groups)]
        for word in l:
            file.write(word + '\t')
        #for group in groups:
            #file1.write(str(group) + ", ")
        #    file1.write(', '.join(str(group)))
        file.write('\n')
        #file1.write(name + "\t" + location + "\t" + str(group_code) + "\t" + str(groups) + "\n")

    print(str(counter) + " sequences were successfully compiled!")    
    file.close()

