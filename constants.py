
# coding: utf-8

# In[5]:

import math
import numpy as np


# In[6]:

"""Constants"""
ts  = 0*10**8   # Time in yrs
te  = 4.5*10**9
dt  = 5*10**6

r_cmb = 3480e3

hl_Ur238 = 4.46e9
secingyr = 3.15569e16
secinyr  = math.pi*1e7


# In[7]:

def radius(r_l,r_t,N):
    r = np.linspace(r_l,r_t,num=N)
    return r


# In[ ]:



