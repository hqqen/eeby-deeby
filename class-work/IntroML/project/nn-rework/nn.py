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
    def __init__(self, nIn = 6, nHide = 32, nOut = 3):
        super(acroNN, self).__init__()
        self.fc1 = nn.Linear(nIn, nHide)
        self.fc2 = nn.Linear(nHide,nOut)

    def forward(self, state):
        state = func.relu(self.fc1(state))
        state = self.fc2(state)
        return func.softmax(state, dim = 1)

    def act(self, state):
        state = torch.from_numpy(state).float().unsqueeze(0).to(device)
        probs = self.forward(state).cpu()
        m = Categorical(probs)
        nextAction = m.sample()
        return nextAction.item() - 1, m.log_prob(nextAction)

#define reinforcement learning policy
def reinforce(network, nEps = 5000, tMax = 100, gamma = 1.0, nPrint = 100):
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
        R = sum([disc, reward for disc, reward in zip(discounts,rewards)])
        loss = []
        for logProb in logProbs:
            loss.append(-logProb*R)
        loss = torch.cat(loss).sum()
        optimizer.zero_grad()
        loss.backwards()
        optimizer.step()


#init environment - seed rng for consistency across runs 
torch.manual_seed(0)
device = torch.device("cpu") #replace with "cuda" if running on cuda capable hardware
env = gym.make('Acrobot-v1')
network = acroNN()
optimizer = optim.Adam(network.parameters(), lr=0.001)
reinforce(network)