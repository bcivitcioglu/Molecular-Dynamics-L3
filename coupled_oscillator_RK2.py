import numpy as np
import matplotlib.pyplot as plt
global num_steps
num_steps = 1000
def coupled_pendula_RK2(s=0.02,k=0.5,m=1,L=1):
    g = 9.8
    omega_0 = (g/L)**(0.5)
    omega_0_2 = (2*k + (g/L))**(0.5)

    # Let's start by defining our variables in a discrete manner
    t = np.linspace(0, num_steps*s, num_steps+1)
    alpha = np.zeros(num_steps+1)
    v_alpha = np.zeros(num_steps+1)
    beta = np.zeros(num_steps+1)
    v_beta = np.zeros(num_steps+1)

    # Initial conditions
    # If theta[0] = 10 and beta[0]=10 then alpha[0] = 10, beta[0] = 0
    alpha[0] = 0
    v_alpha[0] = 0
    beta[0] = 10
    v_beta[0] = 0

    # The RK2 Method
    for n in np.arange(num_steps):
        k_1 = s*v_alpha[n]
        l_1 = -s*(omega_0**2)*alpha[n]
        m_1 = s*v_beta[n]
        o_1 = -s*(omega_0_2**2)*beta[n]

        k_2 = k_1 + s*alpha[n]*(s/2)*(-omega_0**2)
        l_2 = l_1 + v_alpha[n]*s*(-omega_0**2)*(s/2)
        m_2 = m_1 + s*beta[n]*(s/2)*(-omega_0_2**2)
        o_2 = o_1 + v_beta[n]*s*(-omega_0_2**2)*(s/2)

        alpha[n+1] = alpha[n] + k_2
        v_alpha[n+1] = v_alpha[n] + l_2
        beta[n+1] = beta[n] + m_2
        v_beta[n+1] = v_beta[n] + o_2
        
    theta_1 = (1/2)*(alpha+beta)
    theta_2 = (1/2)*(alpha-beta)
    # Plotting the result and comparing with the exact solution
    fig = plt.figure()
    line_1, line_2 = plt.plot(t, theta_1, 'b-', t, theta_2, 'r--')
    fig.legend((line_1, line_2), ('$\theta_1$', '$\theta_2$'), 'upper left')
    plt.xlabel('t')
    plt.ylabel('$\theta$')
    plt.show()

coupled_pendula_RK2()
