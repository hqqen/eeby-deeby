#Alex Beyer - ILA5 - Linear policy Gradients
#implements linear policy gradients with a gaussian potential (i think?) via a NN using the Advatage Actor-Critic (A2C) Method

#necessary imports
import sklearn.preprocessing
import numpy as np
import time
import gym
import matplotlib.pyplot as plt

import torch.nn.functional as F
import torch.optim as optim
import torch.nn as nn
import torch

import time

#boilerplate setup
device = torch.device("cpu") #i ran this via CUDA but that throws errors for unsupported cards - change "cpu" to "cuda" to reenable

#plotter class with data storage for later
class Plotter():
    def __init__(self):
        self.data = []

#actor-critic agent
class ActorCritic(nn.Module):
    def __init__(self, nStates, nActions, nHidden = 16):
        super(ActorCritic, self).__init__()
        #build NN
        self.nActions = nActions
        self.layer1 = nn.Linear(nStates, nHidden)
        self.layer2 = nn.Linear(nHidden, nHidden)
        self.layer3 = nn.Linear(nHidden, nActions)
        self.value = nn.Linear(nHidden, 1)
        self.to(device)

    #forward pass - I'm trying a different way to build this than in the NN homework where I defined everything as part of the same network; here I have the activation functions separate and directly call each layer in the activation function
    #personally I think this makes more sense but that's just how I think about the problem
    #expirimented a bit and relu for all 3 layers gives some interesting properties; namely that it seems to be able to get a positive score - this score might be entirely artificial (I'm not sure the model can even score anything higher than 1 and even that might not be physically possible)
    # but it still seems to correspond to great model preformance
    def forward(self, x):
        x = F.relu(self.layer1(x))
        x = F.relu(self.layer2(x))
        mu = F.relu(self.layer3(x)) #try different activators
        sigma = F.softmax(self.layer3(x), dim = -1) + 1E-5
        dist = torch.distributions.Normal(mu.view(self.nActions).data, sigma.view(self.nActions).data)
        value = self.value(x)
        return dist, value

#build the A2C trainer as a class which calls the ActorCritic agent class in itself
class A2C:
    def __init__(self, envName, gamma = .5, learnRate = .05, nEps = 100, nSteps = 100, nEpsTest = 10):
        
        self.envName = envName
        self.env = gym.make(envName)
        self.model = ActorCritic(self.env.observation_space.shape[0], self.env.action_space.shape[0]).to(device)
        self.opt = optim.Adam(self.model.parameters(),learnRate)

        self.data = {"loss": []}
        self.startTime = None

        self.nEps = nEps
        self.nEpsTest = nEpsTest
        self.nSteps = nSteps
        self.gamma = gamma

    def initStateScaler(self):
        ssSample = np.array([self.env.observation_space.sample() for x in range(10000)])
        self.scaler = sklearn.preprocessing.StandardScaler()
        self.scaler.fit(ssSample)

    def scaleState(self, state):
        scaled = self.scaler.transform(np.array([state[0]]).reshape(1,-1))
        return scaled[0]

    def getNextAction(self, state):
        if type(state) is tuple:
            dist, value = self.model(torch.Tensor(state[0].reshape(1,-1))) 
        else:
            dist, value = self.model(torch.Tensor(state)) 
        action = dist.sample().numpy()
        nextProb = dist.log_prob(torch.FloatTensor(action))
        return action, nextProb, value

    def a2cUpdate(self, rewards, lProbs, values, state):

        qVals = []
        nextQ = 0
        pw = 0
        for reward in rewards[::-1]:
            nextQ += self.gamma ** pw * reward
            pw += 1
            qVals.append(nextQ)

        qVals = qVals[::-1]
        qVals = torch.tensor(qVals)
        qVals = (qVals - qVals.mean()) / (qVals.std() + 1e-5)

        loss = 0
        for nextProb, value, nextQ in zip(lProbs, values, qVals):

            advantage = nextQ - value.item()
            actorLoss = -nextProb * advantage
            criticLoss = F.smooth_l1_loss(value[0], nextQ)
            loss += criticLoss + actorLoss

        self.opt.zero_grad()
        loss.min().backward()  
        self.opt.step()

    # Main training loop.
    def train(self):

        score = 0.0
        rewards = []
        muRewards = []
        sigmaRewards = []
        
        self.initStateScaler()

        self.startTime = time.time()
        for e in range(self.nEps):
            state = self.env.reset()
            score = 0.0
            stepNum = 0

            rewards = []
            lProbs = []
            values = []

            for t in range(self.nSteps):
                
                stepNum += 1

                if type(state) is tuple:
                    state = self.scaleState(state[0].reshape(1,-1))
                else:
                    state = self.scaleState(state.reshape(1,-1))
                action, nextProb, value = self.getNextAction(state)
                state, reward, truncated, terminated, _ = self.env.step(action)
                score += reward
                rewards.append(reward)
                values.append(value)
                lProbs.append(nextProb)
                if truncated or terminated:
                    break

            rewards.append(score)

             # Update Actor - Critic 
            self.a2cUpdate(rewards, lProbs, values, state)

            if (e+1) % 2 == 0:
                print("ep {} got reward: {} in {} steps ".format(e+1, rewards[e], stepNum))

            if (e+1) % 10 == 0:
                nextMuReward, nextSigmaReward = self.test(self.nEpsTest,e)   
                print('ave reward: {} \n w/stdev: {}'.format(nextMuReward, nextSigmaReward))

            nextMuReward, nextSigmaReward = self.test(self.nEpsTest,e)
            muRewards.append(nextMuReward)
            sigmaRewards.append(nextSigmaReward)


        self.env.close()
        return rewards, lProbs, values, muRewards, sigmaRewards
 
    def test(self, nEps, tEp):
        testReward = []
        for e in range(self.nEpsTest):
            state = self.env.reset()
            nextReward = []
            for t in range(self.nSteps):
                action, _, _ = self.getNextAction(state)
                _, reward, truncated, terminated, _ = self.env.step(action)
                nextReward.append(reward)
                if truncated or terminated:
                    break
            testReward.append(sum(nextReward))
        return np.mean(testReward), np.std(testReward)

if __name__ == "__main__":
    A2C = A2C("MountainCarContinuous-v0")
    rewards, lProbs, values, muRewards, sigmaRewards = A2C.train()
    plt.plot(muRewards) #final element is the mean reward and shouldn't be plotted with the others
    plt.xlabel("Episode")
    plt.ylabel("Average Reward")
    plt.title("Average Reward by Episode")
    plt.show()