

import numpy as np
from sympy import *
from scipy.integrate import odeint
from scipy.linalg import solve_continuous_are
import control
import matplotlib.pyplot as plt

class Acrobot(object):
	def __init__(self):
		
		self.init_state_ = np.array([np.pi,0.0,0.0,0.0])
		self.x0 = np.array([np.pi-0.001,0.0,0.0,0.0])
		self.tspan_ = np.arange(0.0,10.0,0.01)
		self.mass_ = np.array([8.0,8.0])
		self.length_ = np.array([0.5,1.0])
		self.g_ = 9.81
		self.inertia_ = self.mass_*pow(self.length_,2)
		self.q1, self.q2,self.q1_dot,self.q2_dot= symbols('q1 q2 q1_dot q2_dot')

	def simulateAcrobot(self):

		m1 = self.inertia_[0]+self.inertia_[1]+self.mass_[0]*pow(self.length_[0],2) + 2*self.mass_[1]*self.length_[0]*self.length_[1]*cos(self.q2)
		m2 = self.inertia_[1] + self.mass_[1]*self.length_[0]*self.length_[1]*cos(self.q2)
		m3 = self.inertia_[1] + self.mass_[1]*self.length_[0]*self.length_[1]*cos(self.q2)
		m4 = self.inertia_[1]
		# M = np.array([[m1,m2],[m3,m4]])
		M = Matrix([[m1,m2],[m3,m4]])
		# print("M = " , M)
		# print("******************************" )

		c1 = -2*self.mass_[1]*self.length_[0]*self.length_[1]*sin(self.q2)*self.q2_dot
		c2 = -self.mass_[1]*self.length_[0]*self.length_[1]*sin(self.q2)*self.q2_dot
		c3 = self.mass_[1]*self.length_[0]*self.length_[1]*sin(self.q2)*self.q2_dot
		c4 = 0
		# C = np.array([[c1,c2],[c3,c4]])
		C = Matrix([[c1,c2],[c3,c4]])
		# print("C = " , C)
		# print("******************************" )

		t1 = -self.mass_[0]*self.g_*self.length_[0]*sin(self.q1) - self.mass_[1]*self.g_*(self.length_[0]*sin(self.q1) + self.length_[1]*sin(self.q1 + self.q2))
		t2 = -self.mass_[1]*self.g_*self.length_[1]*sin(self.q1 + self.q2)
		# T = np.array([[t1],[t2]])
		T = Matrix([[t1],[t2]])
		print("T = " , T)
		print("******************************" )

		# f1 = np.array([[self.q1],[self.q2]])
		f1 = Matrix([[self.q1_dot],[self.q2_dot]])

		# print("f1 = " , f1)
		# print("******************************" )

		M_inv = M.inv()

		# cf1 = C*f1
		print(sin(np.pi))
		print("T= ", T.subs([(self.q1,self.init_state_[0]),(self.q2,self.init_state_[1])]))
		# print("T= ", M_inv)
		# print("C_= ", M_inv)


		f2 = M_inv*(T - C*f1)
		# print("f2 = " , f2)
		# print("******************************" )
		# f2_init = f2.subs([(self.q1,self.init_state_[0]),(self.q2,self.init_state_[1])])
		# print(f2_init)
		

		f = Matrix([[f1],[f2]])

		# print("f = " , f)
		# print("******************************" )

		f0 = f.subs([(self.q1,self.init_state_[0]),(self.q2,self.init_state_[1]),(self.q1_dot,self.init_state_[2]),(self.q2_dot,self.init_state_[3])])


		diff_f = f.jacobian([self.q1,self.q2,self.q1_dot,self.q2_dot])

		A = diff_f.subs([(self.q1,self.init_state_[0]),(self.q2,self.init_state_[1]),(self.q1_dot,self.init_state_[2]),(self.q2_dot,self.init_state_[3])])

		delf_delu = Matrix([[0],[0],[M_inv*Matrix([[0],[1]])]])

		B = delf_delu.subs([(self.q1,self.init_state_[0]),(self.q2,self.init_state_[1]),(self.q1_dot,self.init_state_[2]),(self.q2_dot,self.init_state_[3])])

		A = np.array(A).astype(np.float64)
		B = np.array(B).astype(np.float64)
		f0 = np.array(f0).astype(np.float64)


		xout = odeint(rhs,self.x0,self.tspan_,args = (A,B,f0,self.x0))
		
		q1_plot = xout[:,0]
		q2_plot = xout[:,1]
		q1dot_plot = xout[:,2]
		q2dot_plot = xout[:,3]
		print(q1_plot)
		# print(q2_plot)

		plt.plot(np.arange(0,10,0.01),q2dot_plot)
		plt.title("Angular velocity with time steps")
		plt.ylabel("Derivative of theta2 with time ")
		plt.xlabel("Time")

		# plt.plot(np.arange(0,10,0.01),q2_plot)

		plt.show()

def rhs(x,t,A,B,f0,x0):
	acrobot_obj = Acrobot()
	x0 = x0.reshape((4,1))
	x = x.reshape((4,1))

	Q = np.array(100*eye(4)).astype(np.float64)
	R = 100
	X = solve_continuous_are(A,B,Q,R)
	K = (np.matmul(B.T,X))/R
	K = K.reshape((1,4))
	print("k=",K)
	print("x=",x-x0 )
	# print("x0=",x0 )
	

	u = np.matmul(-K,(x-x0))
	state_change = np.matmul(A,(x-x0))
	input_change = np.matmul(B,u)


	dxdt = state_change + input_change + f0
	dxdt = dxdt.reshape(4)

	return dxdt


if __name__ == "__main__":
	acrobot_obj = Acrobot()
	acrobot_obj.simulateAcrobot()

