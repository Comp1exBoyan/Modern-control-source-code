# backstepping 设计的 非线性控制器仿真
import numpy as np
import matplotlib.pyplot as plt
import cvxpy
import math
import datetime
import time


# 车辆固有参数
L = 2 #车辆轴距
v0 = 2 #车辆初始速度
x0 = 0 #车辆初始x轴位置
y0 = 0 # 车辆初始y轴位置
psi0= 0.5 #车辆初始位姿
k1 = -2.4
k2 = -1.8
k3 = -0.96
dt = 0.1

NX = 3  # x = x, y, yaw
NU = 2  # u = [v,delta]
T = 6

class vehicle_model:
    def __init__(self, x, y, psi, v, L, dt) -> None:
        self.x = x
        self.y = y
        self.psi = psi
        self.v = v
        self.L = L
        # 实现是离散的模型
        self.dt = dt

        # 初始化车辆控制模型

    def update(self, u1, u2) -> None:
        self.x = self.x+self.v*math.cos(self.psi)*self.dt
        self.y = self.y+self.v*math.sin(self.psi)*self.dt
        self.psi = self.psi+self.v/self.L*math.tan(u2)*self.dt
        self.v = u1
        
        # u1 -- 控制器算出来的速度  u2 -- 控制算出来的前轮转角

    def get_state(self) -> None:
        return self.x, self.y, self.psi, self.v

        # 便于记录中途车辆的状态


    ## 模型部分结束，该控制器是基于非线性模型设计，所以只需要做到这一步



    # 轨迹部分 轨迹的坐标点就是期望位置

class MyReferencePath:
    def __init__(self):
        # set reference trajectory
        # refer_path包括4维：位置x, 位置y， 轨迹点的切线方向, 曲率k
        self.refer_path = np.zeros((1000, 4))
        self.refer_path[:, 0] = np.linspace(0, 100, 1000)  # x
        self.refer_path[:, 1] =  np.sin(self.refer_path[:, 0]/20) + 0.5*np.cos(self.refer_path[:, 0]/8)  # y
        # 使用差分的方式计算路径点的一阶导和二阶导，从而得到切线方向和曲率
        for i in range(len(self.refer_path)):
            if i == 0:
                dx = self.refer_path[i+1, 0] - self.refer_path[i, 0]
                dy = self.refer_path[i+1, 1] - self.refer_path[i, 1]
                ddx = self.refer_path[2, 0] + \
                    self.refer_path[0, 0] - 2*self.refer_path[1, 0]
                ddy = self.refer_path[2, 1] + \
                    self.refer_path[0, 1] - 2*self.refer_path[1, 1]
            elif i == (len(self.refer_path)-1):
                dx = self.refer_path[i, 0] - self.refer_path[i-1, 0]
                dy = self.refer_path[i, 1] - self.refer_path[i-1, 1]
                ddx = self.refer_path[i, 0] + \
                    self.refer_path[i-2, 0] - 2*self.refer_path[i-1, 0]
                ddy = self.refer_path[i, 1] + \
                    self.refer_path[i-2, 1] - 2*self.refer_path[i-1, 1]
            else:
                dx = self.refer_path[i+1, 0] - self.refer_path[i, 0]
                dy = self.refer_path[i+1, 1] - self.refer_path[i, 1]
                ddx = self.refer_path[i+1, 0] + \
                    self.refer_path[i-1, 0] - 2*self.refer_path[i, 0]
                ddy = self.refer_path[i+1, 1] + \
                    self.refer_path[i-1, 1] - 2*self.refer_path[i, 1]
            self.refer_path[i, 2] = math.atan2(dy, dx)  # yaw
            # 计算曲率:设曲线r(t) =(x(t),y(t)),则曲率k=(x'y" - x"y')/((x')^2 + (y')^2)^(3/2).
            # 参考：https://blog.csdn.net/weixin_46627433/article/details/123403726
            self.refer_path[i, 3] = (
                ddy * dx - ddx * dy) / ((dx ** 2 + dy ** 2)**(3 / 2))  # 曲率k计算

    def calc_track_error(self, x, y):
        """计算跟踪误差

        Args:
            x (_type_): 当前车辆的位置x
            y (_type_): 当前车辆的位置y

        Returns:
            _type_: _description_
        """
        # 寻找参考轨迹最近目标点
        d_x = [self.refer_path[i, 0]-x for i in range(len(self.refer_path))]
        d_y = [self.refer_path[i, 1]-y for i in range(len(self.refer_path))]
        d = [np.sqrt(d_x[i]**2+d_y[i]**2) for i in range(len(d_x))]
        s = np.argmin(d)  # 最近目标点索引

        yaw = self.refer_path[s, 2]
        k = self.refer_path[s, 3]
        angle = normalize_angle(yaw - math.atan2(d_y[s], d_x[s]))
        e = d[s]  # 误差
        if angle < 0:
            e *= -1

        return e, k, yaw, s

    def calc_ref_trajectory(self, robot_state, dl=1.0):
        """计算参考轨迹点，统一化变量数组，便于后面MPC优化使用
            参考自https://github.com/AtsushiSakai/PythonRobotics/blob/eb6d1cbe6fc90c7be9210bf153b3a04f177cc138/PathTracking/model_predictive_speed_and_steer_control/model_predictive_speed_and_steer_control.py
        Args:
            robot_state (_type_): 车辆的状态(x,y,yaw,v)
            dl (float, optional): _description_. Defaults to 1.0.

        Returns:
            _type_: _description_
        """
        e, k, ref_yaw, ind = self.calc_track_error(
            robot_state[0], robot_state[1])

        xref = np.zeros((NX, T + 1))
        dref = np.zeros((NU, T))
        ncourse = len(self.refer_path)

        xref[0, 0] = self.refer_path[ind, 0]
        xref[1, 0] = self.refer_path[ind, 1]
        xref[2, 0] = self.refer_path[ind, 2]

        # 参考控制量[v,delta]
        ref_delta = math.atan2(L*k, 1)
        dref[0, :] = robot_state[3]
        dref[1, :] = ref_delta

        travel = 0.0

        for i in range(T + 1):
            travel += abs(robot_state[3]) * dt
            dind = int(round(travel / dl))

            if (ind + dind) < ncourse:
                xref[0, i] = self.refer_path[ind + dind, 0]
                xref[1, i] = self.refer_path[ind + dind, 1]
                xref[2, i] = self.refer_path[ind + dind, 2]

            else:
                xref[0, i] = self.refer_path[ncourse - 1, 0]
                xref[1, i] = self.refer_path[ncourse - 1, 1]
                xref[2, i] = self.refer_path[ncourse - 1, 2]

        return xref, ind, dref

# 归一化角度
def normalize_angle(angle):
    """
    Normalize an angle to [-pi, pi].

    :param angle: (float)
    :return: (float) Angle in radian in [-pi, pi]
    copied from https://atsushisakai.github.io/PythonRobotics/modules/path_tracking/stanley_control/stanley_control.html
    """
    while angle > np.pi:
        angle -= 2.0 * np.pi

    while angle < -np.pi:
        angle += 2.0 * np.pi

    return angle


# backstepping controller

def backstepping_ctrl(xref, delta_ref, veh, xref_reg):
    e1 =  veh.x - xref[0, 0] # x 误差
    e2 =  veh.y - xref[1, 0] # y 误差
    e3 =  veh.psi - xref[2, 0] # yaw 角度误差
    
    e = [e1,e2,e3]

    item1 = k1* e1 + (xref[0, 0] - xref_reg[0, 0])/dt
    item2 = k2* e2 + (xref[1, 0] - xref_reg[1, 0])/dt
    item3 = k3* e3 + (xref[2, 0] - xref_reg[2, 0])/dt
    u1 = np.sqrt(pow(item1,2) + pow(item2, 2))
    u2=  np.arctan(L/u1 * item3)
    return u1,u2,e 

itertions = 200
def main():
    reference_path = MyReferencePath()
    goal = reference_path.refer_path[-1,0:2]

    veh  = vehicle_model(x0, y0, psi0, v0, L, dt)
    x_ = []
    y_ = []
    fig = plt.figure(1)
    robot_state = np.zeros(4)
    robot_state[0] = x0
    robot_state[1] = y0
    robot_state[2]=  psi0
    robot_state[3]=  v0
    xref, target_ind, dref = reference_path.calc_ref_trajectory(robot_state)
    xref_reg = xref
    x_e =[]
    y_e =[]
    yaw_e = []
    start_dt = datetime.datetime.now()
    for i in range(itertions):
        robot_state[0] = veh.x
        robot_state[1] = veh.y
        robot_state[2]=  veh.psi
        robot_state[3]=  veh.v
        state_c = robot_state[0:3] # current state
        xref, target_ind, dref = reference_path.calc_ref_trajectory(robot_state)
        u1,u2,e = backstepping_ctrl(xref, dref, veh, xref_reg)
        xref_reg = xref #保存上一时刻的轨迹值
        x_e.append(e[0])
        y_e.append(e[1])
        yaw_e.append(e[2])
        if(u1 > 5):
            u1 = 5
        elif(u1 < -5):
            u1 = -5
        if(u2 > 0.785):
            u2 = 0.785
        elif(u2 < -0.785):
            u2 = -0.785

        veh.update(u1,u2)
        x_.append(veh.x * 0.1)
        y_.append(veh.y * 10)
        plt.cla()
        # plt.title("Simulation Experiment")
        plt.plot(reference_path.refer_path[:, 0]*0.1, reference_path.refer_path[:,
                 1]*10, "-r",  linewidth=1.0, label="reference")
        plt.plot(x_, y_, "-.b", label="trajectory")
        plt.plot(reference_path.refer_path[target_ind, 0]*0.1,
                 reference_path.refer_path[target_ind, 1]*10, "go", label="target")
        plt.xlabel("x-axis position (m)")
        plt.ylabel("y-axis position (m)")
        # plt.axis("equal")
        plt.legend()
        plt.grid(True)
        plt.pause(0.001)

        # camera.snap()
        # 判断是否到达最后一个点
        if np.linalg.norm(robot_state[0:2]-goal)<=0.1:
            print("reach goal")
            break
    end_dt = datetime.datetime.now()
    print("time consume", end_dt - start_dt)

    np.savetxt('pid_x.txt', x_e)
    # # with open('nonlinear_x.txt') as f:
    # #         for line in f:
    # #             print(line, end='')

    # np.savetxt('nonlinear_y.txt', y_e)


    # np.savetxt('nonlinear_yaw.txt', yaw_e)  
                        
    # plt.figure()
    plt.show()
    index = np.linspace(0,100,itertions)
    plt.title("Tracking Errors")
    plt.xlabel("x-axis position")
    plt.ylabel("y-axis position")
    # plt.plot(index, x_e, label = 'x-axis error')
    # plt.plot(index, y_e, label = 'y-axis error')
    plt.plot(index, yaw_e, label = 'yaw error')
    plt.legend()
    plt.grid()
    plt.show()

    

if __name__=='__main__':
    main()