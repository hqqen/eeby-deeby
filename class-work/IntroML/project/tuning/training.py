#library imports
#get openAIgym 
import gym 
#utils
import numpy as np
import matplotlib.pyplot as plt

def main(): 
 #initialize environment, train is and close it when done
  bot = gym.make('Acrobot-v1')
  policyGrad(bot)#, 100, nSteps=200)
  bot.close()



def policyGrad(env, nEpoch=100, alpha=0.01, nSamples=5,\
                    nSteps=500, drawPlots=True):
  ''' 
  given:
  an environment "env" to train in
  a number of epohcs "nEpoch"
  gradient descent rate "alpha"
  number of trajecories to sample "nSamples"
  plotting flag "drawPlots"
  
  preform stochastic swingup training on the acrobot model, returning an ndarray of softmax weights and plot, depending on relevant flags
  ''' 

  #pull in state/accion spaces
  nStates = np.prod(env.observation_space.shape)
  nActions = env.action_space.n

  #xavier init (0 mean, variance of (1/nStates^(1/2))) the initial weights for softmax policy --CHECK THIS
  weights = np.random.normal(0, \
    1 / np.sqrt(nStates), \
    (nActions, nStates))

  #call plotter
  if drawPlots:
    epochs = [i for i in range(nEpoch)]
    rewards = [0] * nEpoch 

  #train
  for ep in range(nEpoch):
    currReward = 0 
    grad = np.zeros_like(weights)
    for tau in range(nSamples):
        #at the start of the loop reset the agent to restart training
      state = env.reset()
      
      for t in range(nSteps):
        #softmax the current state
        if isinstance(state,tuple):
            candidateAction = _softMax(np.matmul(weights, state[0][:].reshape(nStates, 1)))
        else:
            candidateAction = _softMax(np.matmul(weights, state[:].reshape(nStates, 1)))
        #randomly choose a next step from the action space
        nextAction = np.random.choice(nActions, p=candidateAction.flatten())
        #run the chosen action
        state, reward, truncated, terminated, _ = env.step(nextAction)

        if truncated or terminated: 
          break

        currReward += reward

        #one-hot encoding because thats the easiest way to do this (a discrete action space behaves like a categorical space fo my purposes)
        candidateAction = -candidateAction
        candidateAction[nextAction] += 1
        #find the local gradient
        grad += np.matmul( \
          (candidateAction * reward).reshape(nActions, 1), \
          state.reshape(1, nStates))

    #update weights and calculate current reward
    weights -= alpha * grad / nSamples
    rewards[ep] = currReward / nSamples

  #draw plots if flag is set
  if drawPlots:
    makePlots(env, alpha, nSamples, nSteps, epochs, rewards, weights)


def makePlots(env, alpha, nSamples, nSteps, epochs, rewards, weights):
  #given the same system parameters as the trainer, plot training results

  # Plot the rewards over epochs
  plt.plot(epochs, rewards)
  plt.xlabel('nEpoch #')
  plt.ylabel('Average Reward')
  plt.title(f'Policy Gradient Training Performance w/α = {alpha},'\
    f' τ = {nSamples} ')
  plt.show()

  #run som number of training episodes
  for ep in range(1, 4): 
    state = env.reset()
    # Step of trajectory 
    for t in range(nSteps):
        env.render()
        #softmax the state, taking into account datatypes (iteration 1 is fine, 2 and on is not)
        #python is deeply upsetting
        if isinstance(state, tuple):
            candidateAction = _softMax(np.matmul(weights, state[0][:]))
        else:
            candidateAction = _softMax(np.matmul(weights, state[:]))
        #use the gradient probabilities we found to probabilistically choose the next action and run the environment
        nextAction = np.random.choice(env.action_space.n, p=candidateAction)
        state, reward, truncated, terminated, _ = env.step(nextAction)
        if truncated or terminated:
            print(f'Episode {ep} lasted {t} iterations')
            break


def _softMax(x): 
#apply softmax policy to input x
  return np.exp(x - np.amax(x)) / np.sum(np.exp(x - np.amax(x)))

if __name__ == "__main__":
  main()
