from autograd import grad

def bounce(g ,e, v0, h0):
    num_timesteps = 4
    v = [0]*num_timesteps; h = [0]*num_timesteps
    v[0] = v0 # Meters/sec
    h[0] = h0 # Meters
    dt = 0.5  # seconds
    for i in range(1,num_timesteps):
        h[i] = h[i-1] + v[i-1]*dt # Calculate new position
        v[i] = v[i-1] - g*dt # Calculate new velocity
        if h[i] <= 0 :   # Check if bounced
            h[i] = -h[i]    # Assume we regain the negative height
            v[i] = -e*v[i]  # flip the velocity, discounting by the CoR
            #print('i: {}, v: {}, h: {}'.format(i,v[i],h[i]))
    h_final = h[i]
    v_final = v[i]
    return v_final

bGrad = grad(bounce,(0,1,2,3)) #(g=10,e=0.9,v0=0.0,h0=5.0)

print(bGrad(10.,.9,0.,5.))
