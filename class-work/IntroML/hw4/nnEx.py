#library imports
from sklearn.model_selection import train_test_split
import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#environment setup w/rng seeding for consistency
sns.set_context('poster')
np.random.seed(1)

#define NN models
#im being a bad programmer and explicitly defining these instead of a more extensible method
class Net1(nn.Module):
    #PS4.5.1.1 NN model
    def __init__(self):
        super(Net1, self).__init__()  
        # Fully-Connected Layer: 1 (input data) -> 4 (next layer width)
        self.input = nn.Linear(1, 4)
        #hidden layer one 4 -> 4
        self.fc1 = nn.Linear(4, 4)  
        
        # No
        self.sigmoid = nn.Sigmoid()
        #self.sigmoid2 = nn.Sigmoid() #model.eval() seems to indicate I can only use each activation layer once
        # You can try other kinds as well
        # self.relu = nn.ReLU()
        # self.elu = nn.ELU()
        
        
        # Fully-Connected Layer: 4 (prev layer) -> 4 (next layer)
        self.fc2 = nn.Linear(4, 4) 
        # Output layer: 4 (prev layer) -> 1 (output)
        self.fc3 = nn.Linear(4, 1) 

        #our model is not going to run into numerical issues, disable float support by enabling double precision support
        self.double()
    
    # Forward pass builds the model prediction from the inputs
    def forward(self, x):  
        out = self.input(x)
        out = self.sigmoid(out)                            
        out = self.fc1(out)
        out = self.sigmoid(out)
        out = self.fc2(out)
        out = self.sigmoid(out)
        out = self.fc3(out)
        return out
    
class Net2(nn.Module):
    #PS4.5.1.2 NN model
    def __init__(self):
        super(Net, self).__init__()  
        # Fully-Connected Layer: 1 (input data) -> 4 (next layer width)
        self.input = nn.Linear(1, 5)
        #hidden layer one 4 -> 4
        self.fc1 = nn.Linear(5, 5)  
        
        # Non-Linear Layer
        self.relu = nn.ReLU()
        # self.elu = nn.ELU()
        
        
        # Fully-Connected Layer: 4 (prev layer) -> 4 (next layer)
        self.fc2 = nn.Linear(5, 1) 
        # Output layer: 4 (prev layer) -> 1 (output)

        #our model is not going to run into numerical issues, disable float support by enabling double precision support
        self.double()
    
    # Forward pass builds the model prediction from the inputs
    def forward(self, x):  
        out = self.input(x)
        out = self.relu(out)                            
        out = self.fc1(out)
        out = self.relu(out)
        out = self.fc2(out)
        return out

class Net3(nn.Module):
    #PS4.5.1.2 NN model
    def __init__(self):
        super(Net, self).__init__()  
        # Fully-Connected Layer: 1 (input data) -> 4 (next layer width)
        self.input = nn.Linear(1, 3)
        #hidden layer one 4 -> 4
        self.fc1 = nn.Linear(3, 4)  
        
        # Non-Linear Layer
        #self.sigmoid2 = nn.Sigmoid() #model.eval() seems to indicate I can only use each activation layer once
        # You can try other kinds as well
        # self.relu = nn.ReLU()
        self.elu = nn.ELU()
        
        
        # Fully-Connected Layer: 4 (prev layer) -> 4 (next layer)
        self.fc2 = nn.Linear(4, 5) 
        # Output layer: 4 (prev layer) -> 1 (output)
        self.fc3 = nn.Linear(5, 6)
        self.fc4 = nn.Linear(6, 1) 


        #our model is not going to run into numerical issues, disable float support by enabling double precision support
        self.double()
    
    # Forward pass builds the model prediction from the inputs
    def forward(self, x):  
        out = self.input(x)
        out = self.elu(out)                            
        out = self.fc1(out)
        out = self.elu(out)
        out = self.fc2(out)
        out = self.elu(out)
        out = self.fc3(out)
        out = self.elu(out)
        out = self.fc4(out)
        return out

#extensible NN definition - python throws a fit if I have a class be an argument for __init__ in another class so I'm explicitly coding in the different activation functions :(
class extNetSigmoid(nn.Module):
    #arbitrarily deep/wide pytorch NN generator -- sigmoid activator
    #inits 1 layer deep, 2 neuron wide net by default 
    def __init__(self, numHide = 1, numNodes = 2):
        
        super(extNetSigmoid, self).__init__()
        
        #store inputs as class properties
        self.numHide = numHide
        self.numNodes = numNodes

        #hardcode activator to avoid errors from passing classes to classes (I don't know but at least this works)
        self.activation = nn.Sigmoid()  

        #make a list of layers to iterate over
        self.layers = nn.ModuleList()
        #define the first layer to be of numNodes depth
        self.layers.append(nn.Linear(1,self.numNodes))
        #define numHide hidden layers to be of numNodes depth
        for i in range(self.numHide):
            self.layers.append(nn.Linear(self.numNodes,self.numNodes))
        #define output layer to be of numNodes depth
        self.layers.append(nn.Linear(self.numNodes,1))

        #our model is not going to run into numerical issues, disable float support by enabling double precision support
        self.double()
    
    # Forward pass builds the model prediction from the inputs
    def forward(self, x):  
        #step into loop
        out = x
        for layer in self.layers: #iterate over all layers in the list from before
            out = layer(out)
            #dont apply activaiton function to first or last layer
            if layer.out_features != 1 and layer.in_features != 1:
                out = self.activation(out)
        return out



def getLoss(model, inputs, targets, learnRate = .1, numEpochs = 10000):
    lossPlot = []
    for epoch in range(numEpochs):

        ## Do Forward pass
        # Make predictions
        optimizer = torch.optim.SGD(model.parameters(), lr = learnRate, weight_decay = 0)
        outputs = model(inputs.reshape(18,1).double())
        # Compute the loss function
        loss = criterion(outputs, targets)
        lossPlot.append(loss.cpu().detach().numpy())
        ## Update the model
        # Reset the optimizer gradients
        optimizer.zero_grad()
        # Compute the gradient of the loss function
        loss.backward()
        # Do an optimization step
        optimizer.step()
        
        # Print the loss
        if (epoch+1) % 200 == 0:
            print ('Epoch [{:4}/{}], Loss: {:.4f}'.format(epoch+1, numEpochs, loss.item()))
    return lossPlot

n_samples = 30

# What Loss function should we use? MSE!
criterion = nn.MSELoss()

# True Function we want to estimate
true_fun = lambda X: np.cos(1.5 * np.pi * X)

# Noisy Samples from the true function
X = np.sort(2*np.random.rand(n_samples)-1)

y = true_fun(X) + np.random.randn(n_samples) * 0.1

X_train, X_test, y_train, y_test = train_test_split(
     X, y, test_size=0.4, random_state=0)
plt.figure(figsize=(7,7))

X_train = X_train.reshape(-1,1)
y_train = y_train.reshape(-1,1)     
# Plot the data samples
plt.scatter(X_train,y_train, label="Train", c='Blue', s=20, edgecolors='none')
plt.scatter(X_test,y_test, label="Test", c='Red', s=50, edgecolors='none')
#plt.plot(xPlot, true_fun(xPlot), 'g--',label="True function")
plt.legend(loc="best")
sns.despine()
plt.ion()
  
# Convert numpy arrays to torch tensors
inputs = torch.from_numpy(X_train)
targets = torch.from_numpy(y_train)

model1 = extNetSigmoid(3,4)#extNetSigmoid(2, 4)
#lossPlot = getLoss(model1, inputs, targets)
lossPlot = []
numEpochs = 10000
learnRate = .05
optimizer = torch.optim.SGD(model1.parameters(), lr = learnRate, weight_decay = 0)
for epoch in range(numEpochs):

    ## Do Forward pass
    # Make predictions
    outputs = model1(inputs.reshape(18,1).double())
    # Compute the loss function
    loss = criterion(outputs, targets)
    lossPlot.append(loss.cpu().detach().numpy())
    ## Update the model
    # Reset the optimizer gradients
    optimizer.zero_grad()
    # Compute the gradient of the loss function
    loss.backward()
    # Do an optimization step
    optimizer.step()
    
    # Print the loss
    if (epoch+1) % 200 == 0:
        print ('Epoch [{:4}/{}], Loss: {:.4f}'.format(epoch+1, numEpochs, loss.item()))

testOutput = []
xPlot = np.linspace(-1.5,1.5)
#for i in range(len(X_test)):        
#    testOutput.append( model(torch.from_numpy(np.array(X_test[i]))) )#np.linspace(-1.5,1.5)))
testOutput = model1(torch.from_numpy(np.array(X_test.reshape(12,1))))
plt.figure()
plt.title("Model Function Estimate vs Real Function")
plt.plot(X_test, testOutput.cpu().detach().numpy())
plt.plot(xPlot, true_fun(xPlot))


plt.figure()
plt.title("Trainging Loss Over Epochs")
plt.plot(np.array(lossPlot))


plt.figure()
plt.title("Predicted Values")

plt.plot(X_train, y_train, 'ro', label='Data')
predicted = model1(torch.from_numpy(xPlot.reshape(50,1))).detach().numpy()
plt.plot(xPlot, predicted, 'b', label='Prediction')
plt.ylim([1.5,-1.5])
plt.legend()
plt.show(block=True)