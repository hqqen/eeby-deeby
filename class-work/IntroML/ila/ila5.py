#library imports
#get openAIgym 
import gym 
#utils
import numpy as np
import matplotlib.pyplot as plt

def main(): 
    nTrain = 10
    nEpoch = 100
    tau = 5
    alpha = .01
    
    #initialize environment, train is and close it when done
    bot = gym.make('MountainCarContinuous-v0')
    #we want to average training preformance over multiple runs - initialize storage
    ePlot = np.array(np.zeros(nEpoch))
    rPlot = np.array(np.zeros(nEpoch))
    for i in range(nTrain):
        e, r = policyGrad(bot, nEpoch = nEpoch, alpha = alpha, nSamples = tau)
        ePlot = ePlot + e
        rPlot = rPlot + r
        print(i+1)
    #averge over number of training runs
    ePlot = windowedAve(np.array(ePlot)/nTrain,1)
    rPlot = windowedAve(np.array(rPlot)/nTrain,1)
    
    makePlots(alpha, tau, ePlot, rPlot)
    
    bot.close()



def policyGrad(env, nEpoch=100, alpha=0.01, nSamples=5,\
                    nSteps=500, drawPlots=False):
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
    nActions = env.action_space

    #xavier init (0 mean, variance of (1/nStates^(1/2))) the initial weights for softmax policy --CHECK THIS
    weights = np.random.normal(0, \
        1 / np.sqrt(nStates), \
        (nActions, nStates))

    #call plotter
    #if drawPlots:
    epochs = [i for i in range(nEpoch)]
    rewards = [0] * nEpoch 

    #train
    for ep in range(nEpoch):
        currReward = 0 
        grad = np.zeros_like(weights)
        for tau in range(nSamples):
            #at the start of the loop reset the agent to restart training
            state = env.reset()
            #run for nSteps iterations to allow for swingup
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
                #check sim flags for having completed swingup, breaking run if succesful
                if truncated or terminated: 
                    break
                #if we havent broken, add the rewards
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
        makePlots(alpha, nSamples, epochs, rewards)
    return epochs, rewards


def makePlots(alpha, nSamples, epochs, rewards):
    #given the same system parameters as the trainer, plot training results

    # Plot the rewards over epochs
    plt.plot(epochs, rewards)
    plt.plot([epochs[0],epochs[-1]],[-100, -100],linestyle='dashed')
    plt.xlabel('nEpoch #')
    plt.ylabel('Average Reward')
    plt.title(f'Window Averaged (windowSize = 1) Policy Gradient Training Performance w/α = {alpha},'\
        f' τ = {nSamples} ')
    plt.show()

    #run som number of training episodes
    '''
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
    '''


def _softMax(x): 
    #apply softmax policy to input x
    return np.exp(x - np.amax(x)) / np.sum(np.exp(x - np.amax(x)))

def windowedAve(x,n):
    #given:
    #1d dataset "x"
    #window size "n"
    #return:
    #trailing window averaged vector "y"
    #the window will move through "x" with the final element in the window being the current value of "x"; 
    # values on the left of the dataset (i.e. elements with ordinality < n) will have a truncated window
    y = 0*x
    for i in range(np.size(np.array(x))):
        nWindow = n
        windowSum = 0
        if i < n:
            nWindow = i + 1
        for j in range(0,nWindow): #zero is included, zero indexing means we dont need to add 1
            windowSum += x[i-j] 
        y[i] = windowSum/nWindow
    return y

if __name__ == "__main__":
    main()
