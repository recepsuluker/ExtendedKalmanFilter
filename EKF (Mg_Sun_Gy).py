# Purpose:  Integrated test showing how to setup and run a Python BSK module 
# Author:   Recep Suluker
# Creation and Update Date: 5th of Aug, 2022
# source : https://ahrs.readthedocs.io/en/latest/filters/ekf.html


import math
import numpy as np 
import matplotlib.pyplot as plt  
from operator import matmul
from pandas import array


  #The Kalman Gain can be thought of as a blending value within the range  [0.0,1.0]
def KalmanGain (qHat_t, V_t):                   
    K = [0.0, 1.0]
    if K == 0:
        q_t = qHat_t 
    else: 
        K == 1
        q_t = qHat_t + V_t
    return K,q_t


    #To rotate any vector x∈R3 through a quaternion "q"   
def CMatrix (qHat_w, qHat_x, qHat_y, qHat_z):           
    C00 =  1 - 2*(qHat_y**2 + qHat_z**2) 
    C01 =  2 * (qHat_x*qHat_y - qHat_w*qHat_z)
    C02 =  2 * (qHat_x*qHat_z + qHat_w*qHat_y)
     
    C10 =  2 * (qHat_x*qHat_y + qHat_w*qHat_z)
    C11 =  1 - 2*(qHat_x**2 + qHat_z**2) 
    C12 =  2 * (qHat_y*qHat_z + qHat_w*qHat_x)
        
    C20 =  2 * (qHat_x*qHat_z - qHat_w*qHat_y)
    C21 =  2 * (qHat_w*qHat_x + qHat_w*qHat_x)
    C22 =  1 - 2*(qHat_x**2 + qHat_z**2)

    C =  np.array([ [C00, C01, C02],
                    [C10, C11, C12],
                    [C20, C21, C22]])
    return C

def Rotation_Matrix (x1, x2, x3):
    x = np.array([x1, x2, x3]).transpose()
    xdot = np.matmul(CMatrix , x)
    return xdot


class System:
    def __init__(self,q_x,q_y,q_z,q_w, w_x, w_y, w_z):
        
        quaternion = np.array([q_w, q_x, q_y, q_z])     # Initial estimate of the quaternion
        bias = np.array([w_x, w_y, w_z])                # Initial estimate of the gyro bias
        self.xHat_t = np.concatenate((quaternion, bias)).transpose()

        self.P = np.identity(7)
        self.Q_t = np.identity(7)
        self.R = np.zeros(6,6)
        self.C = np.array(3,3)
        self.xdot = np.array(3,1)
        self.K = None
        self.S = None
        self.xHat = None
        self.PHat_t = None
        self.xHatPrev = None
        self.PHatPrev_t = None
        self.jacobianMatrix = np.array(6,4)
        self.zHat_t = np.array(6,1)
        self.V_t = np.array(6,1)
        self.FMatrix = np.array(4,4)
        self.E_w = np.array(3,3)
        self.W_t = np.array(4,3)
        self.qHat_t = np.array(4,1)
        self.gamma_ω = np.array(3,1)
        self.W_t = np.array(4,1)
        self.hqHat_t = np.array(6,1)

    


       
       #xhat_t is state matrix in prediction
    def State_Matrix (qHat_x,qHat_y, qHat_z, dt, w_x, w_y, w_z, q_w, q_x, q_y, q_z): 
        qHat_w = q_w - (dt/2 * w_x *q_x) - (dt/2 * w_y*q_y) - (dt/2 * w_z * q_z)
        qHat_x = q_x + (dt/2 * w_x *q_w) - (dt/2 *w_y *q_z) + (dt/2 * w_z *q_y)
        qHat_y = q_y + (dt/2 * w_x *q_z) + (dt/2* w_y *q_w) - (dt/2 * w_z *q_x )
        qHat_z = q_z - (dt/2 * w_x*q_y) + (dt/2 * w_y * q_x) + (dt/2 * w_z * q_w)

        qHat_t = np.array([[qHat_w] 
                           [qHat_x]
                           [qHat_y] 
                           [qHat_z]])
        return qHat_t


        #xHat_t = F(qt_1, w_t)
    def FMatrix(self, dt, w_x,w_y, w_z):   
        F00 =  1 
        F01 = -dt/2 * w_x 
        F02 = -dt/2 * w_y 
        F03 = -dt/2 * w_z
        F10 =  dt/2 * w_x
        F11 =  1
        F12 =  dt/2 * w_z
        F13 = -dt/2 * w_y
        F20 =  dt/2 * w_y
        F21 = -dt/2 * w_z
        F22 =  1
        F23 =  dt/2 * w_x
        F30 =  dt/2 * w_z
        F31 =  dt/2 * w_y
        F32 = -dt/2 * w_x
        F33 =  1
        A =  np.array([ [F00, F01, F02, F03],
                        [F10, F11, F12, F13],
                        [F20, F21, F22, F23],
                        [F30, F31, F32,F33]])
        return 



        #the spectral density of the gyroscopic noises on each axis, whose standard deviation, 
        # σω, is specified as a scalar in rad/s.
    def Standard_Deviation(gamma_ωx,gamma_ωy, gamma_ωz ):
        gamma_ω = math.sqrt(np.array([gamma_ωx**2 + gamma_ωy**2 + gamma_ωz**2 ])).transpose() 
        return gamma_ω


        #spectral noise covariance matrix, is specified as a scalar in rad/s.
    def Spec_Noise_Cov_Matrix (gamma_ωx,gamma_ωy, gamma_ωz ):                                
        E00 =  gamma_ωx**2
        E01 =  0
        E02 =  0
     
        E10 =  0
        E11 =  gamma_ωy**2
        E12 =  0
        
        E20 =  0
        E21 =  0
        E22 =  gamma_ωz**2

        E_w =  np.array([ [E00, E01, E02],
                          [E10, E11, E12],
                          [E20, E21, E22]])
        return E_w



        #the Jacobian of the prediction, but with respect to the angular rate
    def Jacob_Predic (self, dt, q_w, q_x, q_y, q_z):

        W_t = (dt/2) * np.array ([[-q_x, -q_y, -q_z],
                                  [ q_w, -q_z,  q_y],
                                  [ q_z,  q_w, -q_x],
                                  [-q_y,  q_x,  q_w]])
        return W_t



        #sun and magnometer sensor x,y,z axis     
    def CurrentMeasurement_vector (Sun_x,Sun_y,Sun_z, Mag_x, Mag_y,Mag_z):
            zHat_t = np.array([Sun_x, Sun_y, Sun_z, Mag_x, Mag_y, Mag_z]).transpose()  
            return zHat_t 




     #Global References
      #NED defines the X-, Y-, and Z-axis colinear to the geographical North, East, and Down directions, respectively.
      #ENU defines the X-, Y-, and Z-axis colinear to the geographical East, North, and Up directions, respectively.
    def SunDirection (theta, NED, ENU):
        Sun = test
        if test == NED :
            test =  (np.array([0, math.cos(theta), - math.sin(theta)]).transpose()) / (math.sqrt(math.cos(theta)**2 + math.sin(theta)**2))
        else : 
            test == ENU
            test = (np.array([math.cos(theta), 0, math.sin(theta)]).transpose()) / (math.sqrt(math.cos(theta)**2 + math.sin(theta)**2))
        return test, Sun
     #
    def MagDirection (theta, NED, ENU):
        Mag = test
        if test == NED :
            test =  (np.array([math.cos(theta), 0, math.sin(theta)]).transpose()) / (math.sqrt(math.cos(theta)**2 + math.sin(theta)**2))
        else : 
            test == ENU
            test = (np.array([0, math.cos(theta), - math.sin(theta)]).transpose()) / (math.sqrt(math.cos(theta)**2 + math.sin(theta)**2))
        return test,Mag


     #Normalize the Sun sensor’ measurements
    def Normalize_Sun (Sun_x, Sun_y, Sun_z):
        NormSun = (np.array([Sun_x, Sun_y, Sun_z]).transpose()) / (math.sqrt(Sun_x**2 + Sun_y**2 + Sun_z**2))
        return 
     #Normalize the Mag. sensor’ measurements
    def Normalize_Mag (Mag_x, Mag_y, Mag_z):
        NormSun = (np.array([Mag_x, Mag_y, Mag_z]).transpose()) / (math.sqrt(Mag_x**2 + Mag_y**2 + Mag_z**2))
        return 


        #h(qHat_t)
    def Measurement_Model(self, Sun_x,Sun_y,Sun_z,Mag_x,Mag_y,Mag_z,q_w, q_x, q_y, q_z ):    
        hqHat_t = 2 * np.array([[ h00 ],
                                [ h10 ],
                                [ h20 ],
                                [ h30 ],
                                [ h40 ],
                                [ h50 ]])

        h00 = Sun_x*(0.5-q_y**2-q_z**2) + Sun_y*(q_w*q_z + q_x*q_y) + Sun_z*(q_x*q_z-q_w*q_y) 
        h10 = Sun_x*(q_x*q_y - q_w*q_z) + Sun_y*(0.5-q_x**2-q_z**2) + Sun_z*(q_w*q_x+q_y*q_z) 
        h20 = Sun_x*(q_w*q_y + q_x*q_z) + Sun_y*(q_y*q_z - q_w*q_x) + Sun_z*(0.5-q_x**2-q_y**2)
        h30 = Mag_x*(0.5-q_y**2-q_z**2) + Mag_y*(q_w*q_z+q_x*q_y) + Mag_z*(q_x*q_z-q_w*q_y)
        h40 = Mag_x*(q_x*q_y-q_w*q_z) + Mag_y*(0.5-q_x**2-q_z**2) + Mag_z*(q_w*q_x+q_y*q_z)
        h50 = Mag_x*(q_w*q_y+q_x*q_z) + Mag_y*(q_y*q_z-q_w*q_x) + Mag_z*(0.5-q_x**2-q_y**2)
        return hqHat_t



        # H(qHat_t)
    def getJacobianMatrix(self, Sun_x,Sun_y,Sun_z,Mag_x,Mag_y,Mag_z,q_w, q_x, q_y, q_z ): 
        e00 = Sun_y*q_z - Sun_z*q_y
        e01 = Sun_y*q_y + Sun_z*q_z
        e02 = -2*Sun_x*q_y + Sun_y*q_x - Sun_z*q_w
        e03 = -2*Sun_x*q_z + Sun_y*q_w + Sun_z

        e10 = -Sun_x*q_z + Sun_z*q_x
        e11 = Sun_x*q_y - 2*Sun_y*q_x + Sun_z*q_w
        e12 = Sun_x*q_x + Sun_z*q_z
        e13 = -Sun_x*q_w - 2* Sun_y*q_z + Sun_z

        e20 = Sun_x*q_y - Sun_y*q_x
        e21 = Sun_x*q_z - Sun_y*q_w - 2*Sun_z*q_x
        e22 = Sun_x*q_w + Sun_y*q_z - 2*Sun_z*q_y
        e23 = Sun_x*q_x + Sun_y*q_y
        
        e30 = Mag_y*q_z - Mag_z*q_z 
        e31 = Mag_y*q_y + Mag_z*q_z
        e32 = -2*Mag_x*q_y + Mag_y*q_x - Mag_z*q_w
        e33 = -2*Mag_x*q_z + Mag_y*q_w + Mag_z


        e40 = -Mag_x*q_z + Mag_z*q_x
        e41 = Mag_x*q_y - 2*Mag_y*q_x + Mag_z*q_w
        e42 = Mag_x*q_x + Mag_z*q_z
        e43 = -Mag_x*q_w -2*Mag_y*q_z + Mag_z

        e50 = Mag_x*q_y - Mag_y*q_x
        e51 = Mag_x*q_z - Mag_y*q_w -2*Mag_z*q_x
        e52 = Mag_x*q_w + Mag_y*q_z -2*Mag_z*q_y
        e53 = Mag_x*q_x + Mag_y*q_y

        jacobianMatrix = 2 * np.array([[e00, e01, e02, e03],
                                       [e10, e11, e12, e13],
                                       [e20, e21, e22, e23],
                                       [e30, e31, e32, e33],
                                       [e40, e41, e42, e43],
                                       [e50, e51, e52, e53]])
        return jacobianMatrix



     #Prediction of EKF 
    def Prediction(self):
        q = self.xHat[0:4]
        xHat_t = self.qhat_t
        self.xHatPrev = self.xHat

        Q_t = np.matmul(np.matmul(self.gamma_ω, self.W_t), self.W_t.transpose())
        PHat_t = np.matmul(np.matmul(self.FMatrix, self.PHatPrev_t),self.FMatrix.transpose()) + self.Q_t

        return q, xHat_t ,Q_t, PHat_t
      
      #Correction of EKF 
    def Correction (self, sigmaSun, sigmaMag, zHat_t, hqHat_t):
        tmp1 = np.concatenate(((np.identity(3) * sigmaSun**2 ) , np.identity(3)), axis=1)       
        tmp2 = np.concatenate((np.zeros((3)), ((np.identity(3) * sigmaMag**2 ))), axis=1)         
        self.R = np.concatenate((tmp1, tmp2), axis=0)

        self.V_t = np.subtract(zHat_t, hqHat_t)
        self.S = np.linalg.inv(np.matmul(np.matmul(self.jacobianMatrix, self.PHat_t), self.jacobianMatrix.transpose()) + self.R)
        self.K = np.matmul (np.matmul (self.PHat_t), self.jacobianMatrix.transpose(), self.S)

        self.qHat = self.qHat_t + np.matmul(self.K, self.V_t)
        self.P = np.matmul(np.identity(4) - np.matmul(self.K, self.jacobianMatrix), self.PHat_t)

        return
