#NN based implementation of acrobot swingup problem
#Alex Beyer - ENME743 - Final Project

#imports
#backend
import numpy as np
import time
from collections import deque
#environment
import gym
#ML setup
import torch
import torch.nn as nn
import torch.nn.functional as funct
import torch.optim as opt
from torch.distributions import Categorical

#init environment - seed rng for consistency across runs 
torch.manual_seed(0)
device = torch.device("cpu") #replace with "cuda" if running on cuda capable hardware
env = gym.make('Acrobot-v1')
env.seed(0) 

#build learning policy
class acroNN(nn.Module):
    def __init__(self, nIn = 6, nHide = 32, nOut = 3):