'''
Alex Beyer
ENME743 - ILA5 - Linear Policy Gradients for Mountaincar-v1

Done in tensorflow to try and learn those tools better
The goal is to use policy gradients to solve the moutaincar problem in uder 90 timesteps (i.e. a total score above -90)
I also read into recurrent neural networks as I was figuring this problem out and I'm going to implement them here as practice; 
I think they'd be a useful tool for my research and I have not expirimented with them enough yet
To this end I'm using a fully connected multilayer perceptron

Note: The tutorial I followed to help understand some of this used TensorFlow 1.15.0, so I'm using some calls to the TF1 compaibility API.
'''

#imports
import numpy as np
import sklearn
import sklearn.preprocessing
import tensorflow.compat.v1 as tf
#from tensorflow.contrib import rnn #<---- cant use any rnn calls!!!!!!
import matplotlib as plt
import gym
        
#disable TF v2 features; the tutorial I used is so old that it not only doesn't use v2 features but the placeholder method it uses is no longer supported
tf.disable_v2_behavior()

#setup agent class
class pGradAgent:
    def __init__(self, env, num_input, init_learning_rate=5e-6, min_learning_rate=1e-9, learning_rate_N_max=2000,
                 shuffle=True, batch_size=1, sigma=None):
        self._env = env
        self._sess = tf.Session()
        self._states = tf.placeholder(tf.float32, (None, nStates), name="states")

        self._init_learning_rate = init_learning_rate
        self._min_learning_rate = min_learning_rate
        self._learning_rate_N_max = learning_rate_N_max
        self._learning_rate = tf.placeholder(tf.float32, shape=[])

        self._phi_hidden = 128
        self._sigma_hidden = 32

        # policy parameters
        self._mu_theta = tf.get_variable("mu_theta", [self._phi_hidden, 1],
                                         initializer=tf.zeros_initializer())
        if sigma is None:
            self._sigma_theta = tf.get_variable("sigma_theta", [self._sigma_hidden],
                                                initializer=tf.zeros_initializer())

        # neural featurizer parameters
        self._W1 = tf.get_variable("W1", [nStates, self._phi_hidden],
                                   initializer=tf.random_normal_initializer())
        self._b1 = tf.get_variable("b1", [self._phi_hidden],
                                   initializer=tf.constant_initializer(0))
        self._h1 = tf.nn.tanh(tf.matmul(self._states, self._W1) + self._b1)
        self._W2 = tf.get_variable("W2", [self._phi_hidden, self._phi_hidden],
                                   initializer=tf.random_normal_initializer())
        self._b2 = tf.get_variable("b2", [self._phi_hidden],
                                   initializer=tf.constant_initializer(0))
        self._phi = tf.nn.tanh(tf.matmul(self._h1, self._W2) + self._b2)

        self._mu = tf.matmul(self._phi, self._mu_theta)
        if sigma is None:
            self._sigma = tf.reduce_sum(self._sigma_theta)
            self._sigma = tf.exp(self._sigma)
        else:
            self._sigma = tf.constant(sigma, dtype=tf.float32)

        self._optimizer = tf.train.GradientDescentOptimizer(learning_rate=self._learning_rate)

        self._discounted_rewards = tf.placeholder(tf.float32, (None, 1), name="discounted_rewards")
        self._taken_actions = tf.placeholder(tf.float32, (None, 1), name="taken_actions")

        # we'll get the policy gradient by using -log(pdf), where pdf is the PDF of the Normal distribution
        self._loss = -tf.log(tf.sqrt(1/(2 * np.pi * self._sigma**2)) * tf.exp(-(self._taken_actions - self._mu)**2/(2 * self._sigma**2))) * self._discounted_rewards

        self._train_op = self._optimizer.minimize(self._loss)

        self._sess.run(tf.global_variables_initializer())

        self._num_input = nStates
        self._shuffle = shuffle
        self._batch_size = batch_size
        # rollout buffer
        self._state_buffer  = []
        self._reward_buffer = []
        self._action_buffer = []
        # record reward history for normalization
        self._all_rewards = []
        self._max_reward_length = 1000000
        self._discount_factor = 0.99

        observation_examples = np.array([env.observation_space.sample() for x in range(10000)])
        self._scaler = sklearn.preprocessing.StandardScaler()
        self._scaler.fit(observation_examples)

    def sample_action(self, system_state):
        1 == 1
        if np.size(system_state[0]) == 1:
            system_state = self._scaler.fit_transform(system_state[0].reshape(1, -1))    
        else:
            system_state = self._scaler.transform(system_state[0].reshape(1, -1))
        # Gaussian policy
        mu, sigma = self._sess.run([self._mu, self._sigma], feed_dict={
            self._states: system_state
        })
        action = np.random.normal(mu, sigma)
        action = np.clip(action, self._env.action_space.low[0], self._env.action_space.high[0])
        return action[0], sigma

    def store_rollout(self, state, action, reward):
        self._action_buffer.append(action)
        self._reward_buffer.append(reward)
        self._state_buffer.append(state)

    def update_model(self, iteration):
        N = len(self._reward_buffer)
        r = 0  # use discounted reward to approximate Q value

        discounted_rewards = np.zeros(N)
        for t in reversed(range(N)):
            r = self._reward_buffer[t] + self._discount_factor * r
            discounted_rewards[t] = r

        # reduce gradient variance by normalization
        self._all_rewards += discounted_rewards.tolist()
        self._all_rewards = self._all_rewards[:self._max_reward_length]
        discounted_rewards -= np.mean(self._all_rewards)
        discounted_rewards /= np.std(self._all_rewards)

        learning_rate = self._gen_learning_rate(iteration, l_max=self._init_learning_rate,
                                                l_min=self._min_learning_rate, N_max=self._learning_rate_N_max)

        all_samples = []
        for t in range(N-1):
            state  = self._state_buffer[t]
            action = self._action_buffer[t]
            reward = [discounted_rewards[t]]
            sample = [state, action, reward]
            all_samples.append(sample)
        if self._shuffle:
            np.random.shuffle(all_samples)

        batches = list(self._minibatches(all_samples, batch_size=self._batch_size))

        for b in range(len(batches)):
            batch = batches[b]
            states = [row[0] for row in batch]
            actions = [row[1] for row in batch]
            rewards = [row[2] for row in batch]

            self._sess.run([self._train_op], {
                self._states:             states,
                self._taken_actions:      actions,
                self._discounted_rewards: rewards,
                self._learning_rate:      learning_rate
            })

        self._clean_up()

    def _minibatches(self, samples, batch_size):
        for i in range(0, len(samples), batch_size):
            yield samples[i:i + batch_size]

    def _gen_learning_rate(self, iteration, l_max, l_min, N_max):
        if iteration > N_max:
            return l_min
        alpha = 2 * l_max
        beta = np.log((alpha / l_min - 1)) / N_max
        return alpha / (1 + np.exp(beta * iteration))

    def _clean_up(self):
        self._state_buffer  = []
        self._reward_buffer = []
        self._action_buffer = []


#main logic loop
if __name__ == "__main__":
    #set up simulation and associated variables for later
    env = gym.envs.make("MountainCarContinuous-v0")
    nStates = env.observation_space.shape[0]
    nEpisodes = 200
    nSteps = 1000
    render = False
    agent = pGradAgent(env,nStates)
    
    #begin training loop by setting up per-epoch logic
    for ep in range(nEpisodes):
        state = env.reset()
        totalReward = 0
        stdevs = [] #gaussian stdevs from policy gradient
        mus = []
        done = False
        
        #per-episode logic
        for step in range(nSteps):
            #debug - show the environment to watch whats going on
            if render:
                env.render()
            action, sigma = agent.sample_action(state) #<---------------- this needs to also return the mean mu for storage
            stateFuture, reward, truncated, terminated, _ = env.step(action)
            totalReward += reward
            stdevs.append(sigma)
            agent.store_rollout(state,action,reward) #<---------- rename, stores states within agent for later
            state = stateFuture
            
            if truncated or terminated:
                break
            
        agent.update_model(ep)
        
        print("Epoch " + ep + " done with total reward: " + totalReward)