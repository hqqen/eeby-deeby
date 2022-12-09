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
import torch.nn.functional as func
import torch.optim as optim
from torch.distributions import Categorical

#build learning policy
class acroNN(nn.Module):
    def __init__(self, nIn = 6, nHide = 32, nOut = 3):#6, 32, 3
        super(acroNN, self).__init__()
        self.fc1 = nn.Linear(nIn, nHide)
        self.fc1a = nn.Linear(nHide, nHide)
        self.fc2 = nn.Linear(nHide,nOut)

    def forward(self, state):
        state = func.leaky_relu(self.fc1(state))
        state = func.leaky_relu(self.fc1a(state))
        state = self.fc2(state)
        return func.softmax(state, dim = 1)

    def act(self, state):
        if type(state) == tuple:
            state = torch.from_numpy(state[0]).float().unsqueeze(0).to(device)
        else:
            state = torch.from_numpy(state).float().unsqueeze(0).to(device)
        probs = self.forward(state).cpu()
        m = Categorical(probs)
        nextAction = m.sample()
        return nextAction.item() - 1, m.log_prob(nextAction)

#define reinforcement learning policy
def reinforce(network, nEps = 1500, tMax = 250, gamma = .1, nPrint = 100): #set tiem back to 100????
    scoresDEQ = deque(maxlen = 100)
    score = []
    for ep in range(1,nEps+1):
        logProbs = []
        rewards = []
        state = env.reset()
        for t in range(tMax):
            nextAction, logProb = network.act(state)
            logProbs.append(logProb)
            state, reward, truncated, terminated, _ = env.step(nextAction)
            rewards.append(reward)
            if truncated or terminated:
                break
        scoresDEQ.append(sum(rewards))
        score.append(sum(rewards))
        discounts = [gamma**i for i in range(len(rewards) + 1)]
        #R = sum([(disc, reward) for (disc, reward) in zip(discounts,rewards)])
        R = 0
        for i in range(np.amin([np.size(discounts), np.size(rewards)])):
            R += discounts[i]*rewards[i]
        loss = []
        for logProb in logProbs:
            loss.append(-logProb*R)
        lossSum = torch.cat(loss).sum()
        optimizer.zero_grad()
        lossSum.backward()
        optimizer.step()
    return loss, lossSum, score


#init environment - seed rng for consistency across runs 
torch.manual_seed(0)
device = torch.device("cpu") #replace with "cuda" if running on cuda capable hardware
env = gym.make('Acrobot-v1')
network = acroNN()
optimizer = optim.Adam(network.parameters(), lr=0.05)
loss, lossSum, score = reinforce(network)

1 == 1 