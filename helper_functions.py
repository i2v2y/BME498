#!/usr/bin/env python
# coding: utf-8

# In[36]:


def upload(filename):
    print(filename, " uploaded !")
    #TODO
    return

def download(email):
    print("link sent to ", email)
    #TODO
    return

def quality_control():
    print("Preprocessing: Executing QC...")
    return

def dimension_reduction():
    print("Preprocessing: Executing Dimension Reduction...")
    return

def alignment():
    print("Preprocessing: Executing alignment...")
    return

def normalization():
    print("Preprocessing: Executing normalization...")
    return

def clustering():
    print("Analyzation: Executing clustering...")
    return

def DE():
    print("Analyzation: Executing Differential Gene Analysis...")
    return

def TI():
    print("Analyzation: Executing Trajectory Inference...")
    return

def cell_annotation():
    print("Analyzation: Executing cell_annotation...")
    return

def visualize():
    print("Plotting")
    return

def preprocess(feature, userid ):
    if feature == "quality_control":
        quality_control()
    elif feature == "dimension_reduction": 
        dimension_reduction()
    elif feature == "alignment": 
        alignment() 
    elif feature == "normalization": 
        normalization() 
    else:
        print("Error, feature %s not supported" %(feature))
        
    visualize()
    
    return

def analyze(feature, userid ):
    if feature == "clustering":
        clustering()
    elif feature == "DE": 
        DE()
    elif feature == "TI": 
        TI() 
    elif feature == "cell_annotation": 
        cell_annotation() 
    else:
        print("Error, feature %s not supported" %(feature))
        
    visualize()
    return

def login(userid, password):
    print("user %s successfully logged in" %(userid))
    return


# In[40]:


def parse_request(data):
    for i in range (len(data)):
        userid = data[i][0]
        option = data[i][1]
        detail = data[i][2]
        #print(userid, option, detail)
        if option == "upload":
            upload(detail)
        elif option == "download":
            download(detail)
        elif option == "preprocess":
            preprocess(detail,userid)
        elif option == "analyze":
            analyze(detail,userid)
        elif option == "login":
            login(userid, detail)
        else:
            print("Unable to parse request: %s" %(data[i]))
    return
    
    
    


# In[43]:


import csv

with open('test.csv', newline='') as csvfile:
    data = list(csv.reader(csvfile))

#print(data)
parse_request(data)


# In[42]:


print(data)


# In[ ]:




