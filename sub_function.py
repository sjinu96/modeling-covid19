import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import pandas as pd

import matplotlib as mpl
import random
import pandas as pd
import datetime
from itertools import product


#그리드 함수랑 같이 쓰인다. or 개별 실험때
# omega x 3, beta , c의 parameter 값을 설정함. 
#dealy : 며칠을 미룰 지.
# omega_A[10] : 11일차의 parameter, 즉 11-?12의 parameter임. 
# 그러면, lamda[10]은 11일차의 parameter이니까 ,하루 미루어야 하잖아..
def forgridsearch_parameter(start,end,*parameter):#/// 7개가 들어온다 ()
    #parameter의 수에 따라 달라질 수 있음.[omega_G,omega_A,omega_B,beta,c] 
    #2번재는 [og,ob,beta,c,upper,lower]순. xxx
    #parameter=(0.3,0.1,0.995,0.05) 등.
    global omega_G, omega_A,omega_B, beta, c,p,fever # 없어도 됨.
    global delay, prev_parameter # prev parameter : 마지막에서 두번째 값. array로서 정의
    # array=np.array([omega_G,omega_B,theta,beta,c]) array로 반환해버리면 global 유지 x
    if type(prev_parameter)==int: #윈도우 별 첫 실행
        prev_parameter=parameter # 단순 trick, 후에 수정해야함 ### (1-14)
    omega_G[start-1:start+delay-1]=np.linspace(prev_parameter[0],parameter[0],delay+2)[1:-1] #천천히 선형증가.
    omega_G[start+delay-1:end]=parameter[0]
    omega_A[start-1:start+delay-1]=np.linspace(prev_parameter[1],parameter[1],delay+2)[1:-1] #천천히 선형증가.
    omega_A[start+delay-1:end]=parameter[1]
    omega_B[start-1:start+delay-1]=np.linspace(prev_parameter[2],parameter[2],delay+2)[1:-1] #천천히 선형증가.
    omega_B[start+delay-1:end]=parameter[2]
    beta[start-1:end]=parameter[3]
    c[start-1:start+delay-1]=np.linspace(prev_parameter[4],parameter[4],delay+2)[1:-1] #천천히 선형증가.
    c[start+delay-1:end]=parameter[4]
    
    p[start-1:start+delay-1]=np.linspace(prev_parameter[5],parameter[5],delay+2)[1:-1] #천천히 선형증가.
    p[start+delay-1:end]=parameter[5]
    fever[start-1:start+delay-1]=np.linspace(prev_parameter[6],parameter[6],delay+2)[1:-1] #천천히 선형증가.
    fever[start+delay-1:end]=parameter[6]
    
    return

# Full grid search
# 출력 : (신규 / 누적 확진자 기준) 각각 2개씩, 총 6개
# 1. 최적 파라미터 2. 그 때의 예측 정보 3.총 손실 score 
# cont : init_val을 줄 시, (X,(Q_new,Q_total)) 로 인자를 준다. 기본은 0.
def searching_parameter(start_days,end_days,para,para_range,iternum,MSE_weight,cont=0):
    total_predict_new_array=[]
    total_predict_total_array=[]
    total_mse_score_array_new=[]
    total_mse_score_array_total=[]
    real_sol_new=confirmed_new_korea[start_days-1:end_days] #실제 신규 확진자
    real_sol_total=confirmed_korea[start_days-1:end_days] #실제 누적 확진자
    full_parameter=griding(para,para_range,iternum)
    for i in range(len(full_parameter)):
        #1. 전역 parameter초기화 
        forgridsearch_parameter(start_days,end_days,full_parameter[i])
        print(i,"번째, parameter : ", full_parameter[i])
        #2. simulation ()
        predict_array_new,predict_array_total=seir_simulation(start_days,end_days,cont)[1]
        # 예측 신규/누적 확진자 저장.
        #3. 출력 중 신규 확진자 수의 array 저장.
        total_predict_new_array.append(predict_array_new)
        total_predict_total_array.append(predict_array_total)
        #4. mse최소값 생성. (누적확진자, 신규확진자 둘다 비교하자.
        total_mse_score_array_new.append(MSE(predict_array_new,real_sol_new,MSE_weight))
        total_mse_score_array_total.append(MSE(predict_array_total,real_sol_total,MSE_weight))
    idx_new=np.argmin(total_mse_score_array_new)
    idx_total=np.argmin(total_mse_score_array_total)
    #만약 모든 array가 난수라면 오류가 뜰 것.(우선 history 보기 위해 비활)
    
    opt_para_new=full_parameter[idx_new]
    opt_para_total=full_parameter[idx_total]
    opt_pred_new_new=total_predict_new_array[idx_new] # 신규 기준 최적 파라미터의 신규 확진자
    opt_pred_new_total=total_predict_total_array[idx_new] # 신규 기준 최적 파라미터의 누적 확진자
    opt_pred_total_new=total_predict_new_array[idx_total] # 누적 기준 최적 파라미터의 신규 확진자
    opt_pred_total_total=total_predict_total_array[idx_total] # 누적 기준 최적 파라미터의 누적 확진
    print("신규 확진자 기준최적 파라미터 : ", opt_para_new)
    print("누적 확진자 기준 최적 파라미터 : ", opt_para_total)
    return ([opt_para_new,opt_para_total],[opt_pred_new_new,opt_pred_new_total],
            [opt_pred_total_new,opt_pred_total_total],
[total_mse_score_array_new,total_mse_score_array_total]
,[total_predict_new_array, total_predict_total_array])
#,[total_predict_array_new, total_predict_total_array]) 필요하면.
#출력은 차례대로 [ (신규,누적)(기준)의 파라미터,[신규기준 최적 파라미터의 plotting]
#              (누적 기준 최적 파라미터의 plotting)
#              (신규,누적)(기준) 총 손실 score
#              (full_parameter : 최적 찾기) - 물론, 없어도 gridding 따로 해보면 알긴 암.

# plotting (실제 확진자 <-> 예측 확진자)
# 역시 parameter 순서는 [om_G,om_B,theta,beta]
# predict : [신규확진자,누적확진자] 
# opt 0 : total //opt 1 : new// opt 2 : all
# input : predict - (신규,누적) <-- 2개의 시계열.
# parameter 이상함 , 한 번 simul 할 때 para가 변하면 안됨. 추후 수동 시뮬 단계에서 잘 나누기.
def comparing_real_pred(start_days,end_days,predict,parameter,opt=0):
    "우선 사용하지 않는다."
    para=["G","A","B","beta","c"]# theta는 time-variable로 변하면서 쓰징 낳음.
    real_total=korea_total[start_days-1:end_days] # 이미 평탄화 되어있어야
    real_new=korea_new[start_days-1:end_days]
    time=np.linspace(start_days,end_days,end_days-start_days+1)
    parameter=list(parameter)
    
    for i in range(len(parameter)):
            parameter[i]=round(parameter[i],2)
    # total 출력
    if (opt==0 and opt!=1)or (opt==2):
        
        plt.plot(time,predict[1],label="{}:{}, {}:{}, {}:{}, {}:{},{}:{}".format(para[0],parameter[0],
                                                                        para[1],parameter[1],
                                                                        para[2],parameter[2],
                                                                        para[3],parameter[3],
                                                                                  para[4],parameter[4]))
        plt.plot(time,real_total,'wo',ms=2,mec='k',label="real_total")
    # opt=0으로 하면 total만 출력함.
    if (opt==1 and opt!=0) or(opt==2):
        plt.plot(time,predict[0],label="{}:{}, {}:{}, {}:{}, {}:{},{}:{}".format(para[0],parameter[0],
                                                                    para[1],parameter[1],
                                                                    para[2],parameter[2],
                                                                    para[3],parameter[3],
                                                                    para[4],parameter[4]))
        plt.plot(time,real_new,'wo',ms=2,mec='b',label="real_new")
    plt.legend(loc="best")
    
    
    
    

# 신규, 누적 확진자 자르기.
def slicing(start_days,end_days):
    a=korea_new[start_days-1:end_days] # 이미 smoothing이 끝난 상태이다.
    b=korea_total[start_days-1:end_days]
    

    
    return a,b


# Gridding (5개용)(수정필요 - Full grid search아니면 잘 안씀)
from itertools import product
def griding(para_array,para_range,grid_number): # 리스트와 튜플로 받자.  
    #우선, para_array : [omega_G,omega_B,theta,beta,c]
    # para,range : [(0,1), (0.2,0.5), (0.7,0.8) ..]
    # grid_number : [3,4,5] : 총 60개.
    # 결과값은 1줄로 출력하자. (요소는 tuple)
    parameter_grid=[]
    for i in range(len(para_range)):
        parameter_grid.append(np.linspace(para_range[i][0],
                                                  para_range[i][1],grid_number[i]))
    full=list(product(parameter_grid[0],parameter_grid[1],parameter_grid[2],
                      parameter_grid[3], parameter_grid[4]))
                                             #수정해야함, *args로 list를 받을 수 있게끔
    return full                                 # Class 상속귀찮으니 수동수정 ㄱ
        


# Loss Function + 가중치 포함.
# weight : 0~1, 1이면 가중치 x, 0이면 가장 최근 것만 반영.
# 즉, (0.95^30, 0.95^29 ... 1) 
# norm2
def MSE(X,t,weight=1): #둘다 array가 될 것, time 잘 맞추기.
    X=np.array(X)
    t=np.array(t)
    weight_array=np.array([(weight)**(len(X)-x) for x in range(len(X))])
    #return (1/len(X)*np.sum(np.abs((X-t))))
    return (1/len(X)*np.sum((X-t)**2*weight_array))



# Data Smoothing
# a : 평탄화 할 array
# iter : 반복횟수.
def smoothing(a, iter): 
    
    for j in range(iter):
        a2=[]
        for i in range(len(a)):
            if i==0:
                a2.append((a[i]+a[i+1]+a[i+2])/3)
            elif 0<i<len(a)-1:
                a2.append((a[i-1]+a[i]+a[i+1])/3)  
            else :
                a2.append((a[i-2]+a[i-1]+a[i])/3)
        a=a2
    
    return a2

def smoothing2(a, iter):  # 전방 5일 이평.
    
    for j in range(iter):
        a2=[]
        for i in range(len(a)):
            if i>(3):
                a2.append((a[i]+a[i-1]+a[i-2]+a[i-3]+a[i-4])/5)
            elif i==3:
                a2.append((a[i]+a[i-1]+a[i-2]+a[i-3]+a[i+1])/5)
            elif i==2:
                a2.append((a[i]+a[i-1]+a[i-2]+a[i+1]+a[i+2])/5)    
            elif i==1:
                a2.append((a[i]+a[i-1]+a[i+1]+a[i+2]+a[i+3])/5)
            else:
                a2.append((a[i]+a[i+1]+a[i+2]+a[i+3]+a[i+4])/5)
                              
        a=a2
    
    return a2

def smoothing3(a, iter):  # 전방 5일 이평.
    
    for j in range(iter):
        a2=[]
        for i in range(len(a)):
            if i>(3):
                a2.append((a[i]+a[i-1]+a[i-2]+a[i-3]+a[i-4])/5)
            elif i==3:
                a2.append((a[i]+a[i-1]+a[i-2]+a[i-3])/4)
            elif i==2:
                a2.append((a[i]+a[i-1]+a[i-2])/3)    
            elif i==1:
                a2.append((a[i]+a[i-1])/2)
            else:
                a2.append((a[i])/1)
                              
        a=a2
    
    return a2



def create_sigmoid(start,end,upper,lower,dwell,Q): #dwell : sigmoid의 모양 결정.
    # dwell이 0.8이라면 본 시그모이드 모양에서 일부만 따온다. (1이여야 많이 따옴)

    var=-1/((start-end)/2)*np.log((1-dwell)/dwell)
    #x=np.linspace(start,end,end-start+1) #이렇게하면 안댐
    x=np.linspace(0,end-start,end-start+1)
   
    y=Q*np.exp(var*(x-np.full(end-start+1,(end-start)/2)))/(1+Q*np.exp(var*(x-np.full(end-start+1,(end-start)/2))))*(upper-lower)+lower

    return y




# not sigmoid. only by-polar value
# plt에서 startday-enddays에 따른 ticks 조정하는 visualization 함수 하나 만들면 됨
# def theta_first_scaling(start,end,height): # (시작일,종료일,,끝높이)
#     global theta
#     scaled_data=naver_scaled.values[start-1:end].flatten()*(height/naver_scaled.values[start-1:end].max())
#     theta[start-1:end]=scaled_data
    
#     return scaled_data

# opt 0 : theta_G
# opt 1 : theta_B
# 무조건 정책시행 전 (delay)일 전을 시작날짜로 해야함.
# theta_g_prev 등은 (delay)일 전으로 가져와야함 (시작보다 하루 전)
# opt : start --> 각 time-window 별 첫 시행에서는 선형증가 없애자.
# 물론 parameter에 비해 +1 된 경향이 있음(즉, 150일치의 simulation에는 149일만 필요한데..)
def theta_after_scaling(start,end,upper_G,lower_G,dwell_G,Q_G,lower_B,upper_B,dwell_B,Q_B,delay=3): #순서 바꾸면 그대로 sigmoid처럼 됨.
    global theta_G, theta_A, theta_B # 안 해도 됨.
    global theta_G_prev,theta_B_prev # 해야 함. 정책시행 이전 5일
#     if opt=='start': 일단버리자. 
#         y_g=create_sigmoid(start,end,upper_G,lower_G,dwell_G)
#         y_b=create_sigmoid(start,end,lower_B,upper_B,dwell_B)
#         theta_G[start-1:end]=y_g     
#         theta_B[start-1:end]=y_b
#    if opt!='start':
    y_g=create_sigmoid(start,end-delay,upper_G,lower_G,dwell_G,Q_G)
    y_b=create_sigmoid(start,end-delay,lower_B,upper_B,dwell_B,Q_B)
    theta_G[start-1:start+delay-1]=np.linspace(theta_G_prev,y_g[0],delay+2)[1:-1] #천천히 선형증가.
    theta_G[start+delay-1:end]=y_g
    theta_B[start-1:start+delay-1]=np.linspace(theta_B_prev,y_b[0],delay+2)[1:-1] #천천히 선형증가.
    theta_B[start+delay-1:end]=y_b
    theta_A[start-1:end]=np.full(end-start+1,1)-(theta_G[start-1:end]+
                                                 theta_B[start-1:end])
                                                
from itertools import product
# 아래는 sigmoid fitting을 하는지, 안 하는지에 따라 나뉜다. 
def sa_gridding(alpha=0,eta=0,bol=0, move=0):
    iternum=(len(alpha)*len(eta)*len(bol)*len(move))
    gridding=list(product(alpha,eta,bol,move))
    return gridding

def sa_visual(start,end,full_dict,maxnum,opt): # 상위 maxnum개 출력, 2번째 이상 기점.
    full_df=pd.DataFrame(full_dict)
    best_result=full_df.sort_values(['best_mse'])[0:maxnum].reset_index() #인덱스초기화
    best_params=best_result['best_x'].values #최적 파라미터
    # best_parmas : (num*dim)개의 matrix
    hyper_para=best_result['hyperpara'].values
    dwell=[x[6] for x in best_result['best_x']]
    num=0
    best_mse=best_result['best_mse'].values
    plt.figure(figsize=(24,maxnum*4))
    for i in range(maxnum):
        best_para=best_params[i]
        forgridsearch_parameter(start,end,*(np.append(best_para[:5],best_para[-2:]))) #변수설정.
        global delay
        theta_after_scaling(start,end,best_para[5],best_para[6],best_para[7],best_para[8],best_para[9],best_para[10],best_para[11],best_para[12],delay=delay) #upper.lower,dwell
        globalpara(best_para[13],best_para[14],best_para[15],best_para[16],best_para[17])
        full,pred=seir_simulation(start,end,opt) #첫실행엔 opt에 0
        num+=1
        plt.subplot(maxnum,4,num)
        comparing_real_pred(start,end,predict=pred,parameter=best_para,opt=0)
        if (num<5):
            plt.title('total_case',loc='left',fontsize=12,fontweight=0)
        num+=1
        plt.subplot(maxnum,4,num)
        comparing_real_pred(start,end,predict=pred,parameter=best_para,opt=1)
        if (num<5):
            plt.title('new_case',loc='left',fontsize=12,fontweight=0)
        num+=1
        plt.subplot(maxnum,4,num)
        plt.plot(theta_G[start-1:end], label='g')
        plt.plot(theta_A[start-1:end], label='a')
        plt.plot(theta_B[start-1:end], label='b')
        plt.legend()
        if (num<5):
            plt.title('theta : {}'.format(best_para[6]),loc='left',fontsize=12,fontweight=0)
        num+=1
        plt.subplot(maxnum,4,num)
        plt.plot(best_result['full_mse'][i])
        plt.title('mse_by_iter({})'.format(best_mse[i]),loc='left',fontsize=15,fontweight=0)
    plt.suptitle("Best {} results({}-{})".format(maxnum,start,end),fontsize=30,color='black')
    #plt.text(0.5, 0.02, 'Time', ha='center', va='center')
    #plt.text(0.06, 0.5, 'value', ha='center', va='center', rotation='vertical')
    
    return



#각 윈도우별 첫째.
def sa_visual_first(start,end,full_dict,maxnum,opt): # 상위 maxnum개 출력, 2번째 이상 기점
    global theta_G_prev,theta_B_prev
    full_df=pd.DataFrame(full_dict)
    best_result=full_df.sort_values(['best_mse'])[0:maxnum].reset_index() #인덱스초기화
    best_params=best_result['best_x'].values #최적 파라미터
    # best_parmas : (num*dim)개의 matrix
    hyper_para=best_result['hyperpara'].values
    dwell=[x[6] for x in best_result['best_x']]
    num=0
    best_mse=best_result['best_mse'].values
    plt.figure(figsize=(24,maxnum*4))
    for i in range(maxnum):
        best_para=best_params[i]
        forgridsearch_parameter(start,end,*(np.append(best_para[:5],best_para[-2:]))) #변수설정.
        
        theta_G_prev=best_para[5] #선형성 없애기 위하여. 
        theta_B_prev=best_para[9]
        global delay
        theta_after_scaling(start,end,best_para[5],best_para[6],best_para[7],best_para[8],best_para[9],best_para[10],best_para[11],best_para[12],delay=delay) #upper.lower,dwell
        globalpara(best_para[13],best_para[14],best_para[15],best_para[16],best_para[17])
        full,pred=seir_simulation(start,end,opt) #첫실행엔 opt에 0
        num+=1
        plt.subplot(maxnum,4,num)
        comparing_real_pred(start,end,predict=pred,parameter=best_para,opt=0)
        if (num<5):
            plt.title('total_case',loc='left',fontsize=12,fontweight=0)
        num+=1
        plt.subplot(maxnum,4,num)
        comparing_real_pred(start,end,predict=pred,parameter=best_para,opt=1)
        if (num<5):
            plt.title('new_case',loc='left',fontsize=12,fontweight=0)
        num+=1
        plt.subplot(maxnum,4,num)
        plt.plot(theta_G[start-1:end], label='g')
        plt.plot(theta_A[start-1:end], label='a')
        plt.plot(theta_B[start-1:end], label='b')
        plt.legend()
        if (num<5):
            plt.title('theta : {}'.format(best_para[6]),loc='left',fontsize=12,fontweight=0)
        num+=1
        plt.subplot(maxnum,4,num)
        plt.plot(best_result['full_mse'][i])
        plt.title('mse_by_iter({})'.format(best_mse[i]),loc='left',fontsize=15,fontweight=0)
    plt.suptitle("Best {} results({}-{})".format(maxnum,start,end),fontsize=30,color='black')
    #plt.text(0.5, 0.02, 'Time', ha='center', va='center')
    #plt.text(0.06, 0.5, 'value', ha='center', va='center', rotation='vertical')

    return 



def visualize_by_bestparams(start,end,best_params,inputs):
    plt.figure(figsize=(24,maxnum*4))
    forgridsearch_parameter(start,end,*np.append(best_para[:5],best_para[-2:])) #변수설정.t
    global delay
    theta_after_scaling(start,end,best_para[5],best_para[6],best_para[7],best_para[8],best_para[9],best_para[10],best_para[11],best_para[12],delay) #upper.lower,dwell
    globalpara(best_para[13],best_para[14],best_para[15],best_para[16],best_para[17])
    full,pred=seir_simulation(start,end,inputs)
    plt.title('Individual prediction by best parameter')
    # 이후 추가(1.1)
    
    return full,pred # full : S~E~A~I~Q...

def check_error(full_dict): #first_full_dict. (sa_pipe 이후. )
    tau1=full_dict['beta']*full_dict['c']*full_dict['S']*full_dict['ep_A']*full_dict['ep']*full_dict['sigma']*(full_dict['omega_G']*full_dict['E_G']+full_dict['omega_A']*full_dict['E_A']+full_dict['omega_B']*full_dict['E_B'])/51780000
    tau2=full_dict['beta']*full_dict['c']*full_dict['S']*full_dict['fever']*(full_dict['omega_G']*full_dict['I_G']+full_dict['omega_A']*full_dict['I_A']+full_dict['omega_B']*full_dict['I_B'])/51780000
    tau3=full_dict['beta']*full_dict['c']*full_dict['S']*full_dict['ep_A']*(full_dict['omega_G']*full_dict['A_G']+full_dict['omega_A']*full_dict['A_A']+full_dict['omega_B']*full_dict['A_B'])/51780000
    

    plt.figure(figsize=(10,20))
    plt.subplot(5,1,1)
    plt.title("Which group will affect how much?")
#     (full_dict['Q_total']*tau1/(tau1+tau2+tau3)).plot(kind='line',label='E')
#     (full_dict['Q_total']*tau2/(tau1+tau2+tau3)).plot(kind='line',label='I')
#     (full_dict['Q_total']*tau3/(tau1+tau2+tau3)).plot(kind='line',label='A')
    tau1.plot(kind='line',label='E')
    tau2.plot(kind='line',label='I')
    tau3.plot(kind='line',label='A')
    plt.legend()
    plt.subplot(5,1,2)
    plt.title("plotting theta")
    full_dict['theta_G'].plot(kind='line', label='theta_G')
    full_dict['theta_A'].plot(kind='line',label='theta_A')
    full_dict['theta_B'].plot(kind='line',label='theta_B')
    plt.legend()
    plt.subplot(5,1,3)
    plt.title("comparing with real new confirmed")
    full_dict['Q_new'].reset_index(drop=True).plot(kind='line',label='predict_new')
    korea_new[full_dict.index[0]-1:full_dict.index[-1]].reset_index(drop=True).plot(kind='line', label='real_new')
    plt.legend()
    plt.subplot(5,1,4)
    plt.title("Unchaning parameter")
    full_dict['beta'].plot(kind='line', label='beta({})'.format(full_dict['beta'].iloc[-1]))
    full_dict['ep_A'].plot(kind='line',label='ep_A({})'.format(full_dict['ep_A'].iloc[-1]))
    full_dict['sigma'].plot(kind='line',label='sigma({})'.format(full_dict['sigma'].iloc[-1]))
    full_dict['ep'].plot(kind='line',label='ep({})'.format(full_dict['ep'].iloc[-1]))
    full_dict['row'].plot(kind='line',label='row({})'.format(full_dict['row'].iloc[-1]))
    full_dict['fever'].plot(kind='line',label='fever({})'.format(full_dict['fever'].iloc[-1]))
    full_dict['p'].plot(kind='line',label='p({})'.format(full_dict['p'].iloc[-1]))
    plt.legend()
    plt.subplot(5,1,5)
    plt.title("chaning parameter")
    full_dict['c'].plot(kind='line', label='c')
    full_dict['omega_G'].plot(kind='line', label='omega_G')
    full_dict['omega_A'].plot(kind='line', label='omega_A')
    full_dict['omega_B'].plot(kind='line', label='omega_B')
    full_dict['delta'].plot(kind='line',label='delta')
    plt.legend()
        

    #return tau1,tau2,tau3