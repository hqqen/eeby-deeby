#basic imports
import random
import numpy as np
import gym
import matplotlib.pyplot as plt
import glob
import io
import base64
import os
from collections import deque,namedtuple #will be used for datastacks later

#code uses pytorch DNNs
import torch
from torch import nn

#get rid of tqdm
from tqdm.notebook import tqdm

#imports for getting & storing episode videos
from IPython.display import HTML
from IPython import display as ipythondisplay
from pyvirtualdisplay import Display
from gym.wrappers import RecordEpisodeStatistics
from gym.wrappers import RecordVideo 

class DQN(nn.Module):
    #deep Q-Network Class to call later
    #we pass in the dimensionality of the state and action spaces (6 and 3, respectively) and returns the optimal law Q
    def __init__(self, state_space_dim, action_space_dim):
        super().__init__()

        #define the NN
        self.linear = nn.Sequential(
                  nn.Linear(state_space_dim,64),
                  nn.ReLU(),
                  nn.Linear(64,64*2),
                  nn.ReLU(),
                  nn.Linear(64*2,action_space_dim)
                )

    def forward(self, x):
        x = x.to()#device)
        return self.linear(x)

class ReplayMemory(object):
    #function to store the model memory
    def __init__(self, memLen):
        #want a double-ended queue, not a normal stack so we can append to/pop from either end
        self.memory = deque(maxlen = memLen)

    def push(self, s, a, sPrime, r):
        # append a tuple (state, action, next state, reward) to the memory
        self.memory.append((s, a, sPrime, r))

    def popRand(self, batch_size):
        #given a batch size return a normally distributed sample of data points (pop a random number of points off the stack)
        #need this to be from a distribution to limit the effects of data correlation; gives poor model preformance
        batch_size = min(batch_size, len(self)) 
        return random.sample(self.memory,batch_size)

    def __len__(self):
        return len(self.memory) 

def show_videos():
    #function to load and play the set of mp4s corresponding with a test run
    # I'm not going to lie, a lot of this function is taken from code online 

    #get a list of all mp4s in the folder and srot by record order
    mp4list = glob.glob('video/*.mp4')
    mp4list.sort()

    #play each video in order, showing the name of each
    for mp4 in mp4list:
        print(f"\nSHOWING VIDEO {mp4}")
        video = io.open(mp4, 'r+b').read()
        encoded = base64.b64encode(video)
        ipythondisplay.display(HTML(data='''<video alt="test" autoplay 
            loop controls style="height: 400px;">
            <source src="data:video/mp4;base64,{0}" type="video/mp4" />
        </video>'''.format(encoded.decode('ascii'))))
    
def wrap_env(env, video_callable=None):
    #basic environment wrapper so we can record videos later
    env = RecordVideo(env, './video', force=True, video_callable=video_callable)
    return env

def choose_action_epsilon_greedy(net, state, epsilon):
    
    if epsilon > 1 or epsilon < 0:
        raise Exception('The epsilon value must be between 0 and 1')
                
    # Evaluate the network output from the current state
    with torch.no_grad():
        net.eval()
        state = torch.tensor(state, dtype=torch.float32) # Convert the state to tensor
        net_out = net(state)

    # Get the best action (argmax of the network output)
    best_action = int(net_out.argmax())
    # Get the number of possible actions
    action_space_dim = net_out.shape[-1]

    # Select a non optimal action with probability epsilon, otherwise choose the best action
    if random.random() < epsilon:
        # List of non-optimal actions (this list includes all the actions but the optimal one)
        non_optimal_actions = [a for a in range(action_space_dim) if a != best_action]
        # Select randomly from non_optimal_actions
        action = random.choice(non_optimal_actions)
    else:
        # Select best action
        action = best_action
        
    return action, net_out.cpu().numpy()

def choose_action_softmax(net, state, temperature):
    
    if temperature < 0:
        raise Exception('The temperature value must be greater than or equal to 0 ')
        
    # If the temperature is 0, just select the best action using the eps-greedy policy with epsilon = 0
    if temperature == 0:
        return choose_action_epsilon_greedy(net, state, 0)
    
    # Evaluate the network output from the current state
    with torch.no_grad():
        net.eval()
        state = torch.tensor(state, dtype=torch.float32)
        net_out = net(state)

    # Apply softmax with temp
    temperature = max(temperature, 1e-8) # set a minimum to the temperature for numerical stability
    softmax_out = nn.functional.softmax(net_out/temperature, dim=0).cpu().numpy()
                
    # Sample the action using softmax output as mass pdf
    all_possible_actions = np.arange(0, softmax_out.shape[-1])
    # this samples a random element from "all_possible_actions" with the probability distribution p (softmax_out in this case)
    action = np.random.choice(all_possible_actions,p=softmax_out)
    
    return action, net_out.cpu().numpy()

env = gym.make('Acrobot-v1', render_mode="rgb_array")
#env.seed(0)
#env = wrap_env(env, video_callable=lambda episode_id: True)

for num_episode in range(10): 
    state = env.reset()
    score = 0
    done = False
    # Go on until the pole falls off or the score reach -500
    while not done and score > -500:
      # Choose a random action
      action = random.choice([0, 1, 2])
      next_state, reward, truncated, terminated, info = env.step(action)
      done = truncated or terminated
      # Visually render the environment 
      env.render()
      # Update the final score (-1 for each step)
      score += reward 
      state = next_state
      # Check if the episode ended (the pole fell down)
    print(f"EPISODE {num_episode + 1} - FINAL SCORE: {score}") 


display = Display(visible=0, size=(1400, 900))
display.start()
#show_videos()