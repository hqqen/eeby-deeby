#NN based implementation of acrobot swingup problem
#Alex Beyer - ENME743 - Final Project

#imports
#backend
import numpy as np
import time
from collections import deque
import matplotlib.pyplot as plt
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
        self.fc1b = nn.Linear(nHide, nHide)
        self.fc1c = nn.Linear(nHide, nHide)
        self.fc2 = nn.Linear(nHide,nOut)

    def forward(self, state):
        state = func.softmax(self.fc1(state), dim = 1)
        state = func.softmax(self.fc1a(state), dim = 1)
        state = func.softmax(self.fc1b(state), dim = 1)
        state = func.softmax(self.fc1c(state), dim = 1)
        state = self.fc2(state)
        return func.softmax(state, dim = 1)

    def act(self, state):
        if type(state) == tuple:
            state = torch.from_numpy(state[0]).float().unsqueeze(0).to(device)
        else:
            state = torch.from_numpy(state).float().unsqueeze(0).to(device)
        probs = self.forward(state).cpu()
        m = Categorical(probs = probs)
        nextAction = m.sample()
        return nextAction.item() - 1, m.log_prob(nextAction) #nextAction.item() - 1

#define reinforcement learning policy
def reinforce(network, nEps = 1500, tMax = 1000, gamma = 1., nPrint = 100): #set tiem back to 100????
    env = gym.make('Acrobot-v1')
    scoresDEQ = deque(maxlen = 100)
    score = []
    opt= optim.Adam(network.parameters(), lr=3e-4)
    for ep in range(1,nEps+1):
        logProbs = []
        rewards = []
        state = env.reset()
        #env.state = np.array([0,0,0,0,0,0])
        #state = env.state
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
        R = sum([pair[0]*pair[1] for pair in zip(discounts,rewards)])
        #R = 0
        #for i in range(np.amin([np.size(discounts), np.size(rewards)])):
        #    R += discounts[i]*rewards[i]
        loss = []
        for logProb in logProbs:
            loss.append(-logProb*R)
        lossSum = torch.cat(loss).sum()

        #if score[-1] >= score[-2]:
        opt.zero_grad(set_to_none = True)
        lossSum.backward()
        #print("Pre Update Weights:")
        #print(list(network.parameters()))
        opt.step()
        #print("Post Update Weights:")
        #print(list(network.parameters()))
        1 == 1 #IF SCORE FOR THIS ROUND IS LOWER THAN THE LAST DONT UPDATE!!!!
    return loss, lossSum, score, nEps


#init environment - seed rng for consistency across runs 
torch.manual_seed(0)
device = torch.device("cpu") #replace with "cuda" if running on cuda capable hardware
env = gym.make('Acrobot-v1')
network = acroNN()
#optimizer = optim.Rprop(network.parameters(), lr=0.0025)
i = 0
#score = [-500]
#while score[-1] == -98989: 
loss, lossSum, score, nEps = reinforce(network)
i += 1
print(i)

print(f"Model took {i} runs to find a solution")


plt.plot(score)
plt.plot([0, nEps], [-100, -100], linestyle = 'dashed')
plt.title(r"Training Score Curve for Acrobot ($\alpha$ = .05)")
plt.xlabel('nEpoch #')
plt.ylabel('Average Reward')
plt.show()