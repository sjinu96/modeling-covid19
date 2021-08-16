# 위처럼, t-days -> (t+1) days로의 방정식 구현
# parameter의 원본을 유지하기 위해 변수이름 뒤에 2를 붙혔다.
# 여기서의 변수 이름은, 후에 seir_simulation 함수와 맞춰 줘야 한다.
# var_para1 : 1개의 변수값들 (beta,c2,lamda 등)
# var_para2 : 1개의 변수값들 (omega,theta)

def seir_onestep(X,var_para1,var_para2): 
    '''
    이는, 이하 seir_simulation과 이어진다.
    X[t]=[S,E_G,E_A,E_B,I_G,I_A,I_B,A_G,A_A,A_B,Q,R,R_I,R_A,Q_total] , for days = t
    input : X[t]
    output : X[t+1]
    var_para1 : [beta,lamda,c]
    var_para2 : [omega x 3 , theta x 3 , p,fever2] : 8개, 모두 상수.
    '''
    global ep
    N=sum(list(X))
    beta2,lamda2,c2=var_para1
    omega_G2,omega_A2,omega_B2,theta_G2,theta_A2,theta_B2,p2,fever2=var_para2
    S,E_G,E_A,E_B,I_G,I_A,I_B,A_G,A_A,A_B,Q,R,R_I,R_A,Q_total = list(X)
    I=I_G+I_A+I_B
    A=A_G+A_A+A_B
    E=E_G+E_A+E_B
    
    Gamma=beta2*c2*(fever2*omega_G2*I_G+fever2*omega_A2*I_A+fever2*omega_B2*I_B*+ep_A*omega_G2*A_G
                 +ep_A*omega_A2*A_A+ep_A*omega_B2*A_B+omega_G2*ep_A*ep*sigma*E_G+omega_A2*ep_A*ep*sigma*E_A+
                   omega_B2*ep_A*ep*sigma*E_B)/N
    
    gamma=1/(9-ep) # 9일은 가정하자. 
    
    tau=delta*p2/(1-p2) # p : 사후 무증상 확진율
    S_next=S-Gamma*S-lamda2
    E_G_next=E_G+theta_G2*Gamma*S-sigma*E_G-tau*I_G #1-10
    E_A_next=E_A+theta_A2*Gamma*S-sigma*E_A-tau*I_A ###
    E_B_next=E_B+theta_B2*Gamma*S-sigma*E_B-tau*I_B
    I_G_next=I_G+(1-row)*sigma*E_G-delta*I_G
    I_A_next=I_A+(1-row)*sigma*E_A-delta*I_A###
    I_B_next=I_B+(1-row)*sigma*E_B-delta*I_B
    A_G_next=A_G+row*sigma*E_G-gamma*A_G
    A_A_next=A_A+row*sigma*E_A-gamma*A_A ###
    A_B_next=A_B+row*sigma*E_B-gamma*A_B
    Q_next=Q+delta*I-mu*Q+lamda2+tau*I
    R_I_next=R_I+mu*Q
    R_A_next=R_A+gamma*A
    R_next=R_I_next+R_A_next
    
    
    # 이하, tdays -> (t+1)days로 갈 때 신규 확진자 나타내기 위함
    # Q만 구하면 된다. # 나머지는, 각 구획에 머무르는 시간으로 나눠주면 된다.  
    
    Q_new=lamda2+delta*I+tau*I # 1=10
    Q_total_next=Q_total+Q_new
    return ([S_next,E_G_next,E_A_next,E_B_next,I_G_next,I_A_next,I_B_next,A_G_next,
             A_A_next,A_B_next,Q_next,R_next,R_I_next,R_A_next,Q_total_next],[Q_new,Q_total_next])

# 단, 이 코드는 단독으로 실행할 수 없고 , 하단의 seir_simulation 함수 내부에서만 작동한다.
# parameter 때문
# 여기서는, seir_simulation 함수 내에서 매 t번째 날의 parameter들을 "상수"로 고정해야 하기 때문에
# parameter 뒤에 2를 붙혔다.

# 함수 1 : 각 집단별 초기값 set 출력
# input : 시작 할 날짜 (우선 days, 후에 datetime 사용하여 수정)
# output : X0=[S,E_G,E_B.....] 단, 해당 일차부터 
# 개선 : days <=15인 경우는 예외로 처리.
# opt : init_val을 임의로 조정하는 부분이 필요할 때
def init_val(days,opt=0):
    # 아래 주석은 모두 30일차 기준 예시. 
    global init_I
    
    if opt!=0 :# opt에 인자가 주어 질 경우, 그 인자로 바로 반환한다.
        return opt # 단, 모양 자체가 (x,[a,b]) 꼴이여야 한다.
    if days<=20: # init_I : 초기 확진자 수.
        I=init_I
        A=row/(1-row)*I #항체검사 기반.
        E=I/((1-p)*sigma) # 무증상확진자 고려해서 배수.
        Q=0
        R=0
        S=50000000-E-I-A-Q-R
        #1일차만 고려하면 되므로. 
        Q_total_init=sum(korea_new[:days])
        Q_init=sum(korea_new[days-1:days+3])*(1/4) #조심해야함.. 뭐 상관은 없다 크게,, 
      
    N=51780000
    
    gamma=1/(9-ep) 

   
    Q_t_2=korea_qua.iloc[days-3]


    R_I_t_2=korea_total[days-3]-Q_t_2
   
    I_t_2=I_t_1
    
    E_t_2=delta/(sigma*(1-row))*I_t_2
    A_t_2=row*sigma/gamma*E_t_2
    R_A_t_2=((1-p[days-1])*(Q_t_2+R_I_t_2-lamda_total[days-3])+I_t_2)*(row/(1-row))-A_t_2
    R_t_2=R_I_t_2+R_A_t_2
   
    S_t_2=N-E_t_2-I_t_2-Q_t_2-R_I_t_2-A_t_2-R_A_t_2

    A_t_1=A_t_2+row*sigma*E_t_2-gamma*A_t_2
   
    OMEGA=beta[days-1]*c[days-1]*(omega_G[days-1]*theta_G[days-1]+omega_A[days-1]*theta_A[days-1]+
                                          omega_B[days-1]*theta_B[days-1])
    GAMMA_t_2=OMEGA*(fever[days-1]*I_t_2+ep_A*A_t_2+ep_A*ep*sigma*E_t_2)/N
    E_t_1=E_t_2+GAMMA_t_2*S_t_2-sigma*E_t_2-p[days-1]*delta/(1-p[days-1])*I_t_2
    S_t_1=S_t_2-GAMMA_t_2*S_t_2-lamda[days-2]
    R_I_t_1=R_I_t_2+mu*Q_t_2
    R_A_t_1=R_A_t_2+gamma*A_t_2
    R_t_1=R_I_t_1+R_A_t_1
    Q_t_1=Q_t_2+delta*I_t_2-mu*Q_t_2+lamda[days-2]+p[days-1]*delta/(1-p[days-1])*I_t_2
    
    
    
    GAMMA_t_1=OMEGA*(fever[days-1]*I_t_1+ep_A*A_t_1+ep_A*ep*sigma*E_t_1)/N
    
 
    tau=delta*p[days-1]/(1-p[days-1]) # p : 사후 무증상 확진율
    S=S_t_1-GAMMA_t_1*S_t_1-lamda[days-1]
    E=E_t_1+GAMMA_t_1*S_t_1-sigma*E_t_1-tau*I_t_1#1-10
    I=I_t_1+(1-row)*sigma*E_t_1-delta*I_t_1
    A=A_t_1+row*sigma*E_t_1-gamma*A_t_1
    Q=Q_t_1+delta*I_t_1-mu*Q_t_1+lamda[days-1]+tau*I_t_1
    R_I=R_I_t_1+mu*Q_t_1
    R_A=R_A_t_1+gamma*A_t_1
    R=R_I+R_A

    #R은 집계가 되지 않을 사람들까지 포함 된 수이다.
    E_G=E*theta_G[days-1]
    E_A=E*theta_A[days-1]
    E_B=E*theta_B[days-1]
    I_G=I*theta_G[days-1]
    I_A=I*theta_A[days-1]
    I_B=I*theta_B[days-1]
    A_G=A*theta_G[days-1]
    A_A=A*theta_A[days-1]
    A_B=A*theta_B[days-1]

    #얘네는 초기라 어쩔 수가 없음. 
    Q_init=tau*I_t_1+lamda[days-1]+delta*I_t_1
    Q_total_init=Q_t_1+R_I_t_1+Q_init
    for idx, x in enumerate([A_t_1,A_t_2,A,E_t_1,E_t_2,E,I_t_1,I_t_2,I,R,R_I,R_A]):
        if x < 0:
            print('초기값(DAYS-3~day-1 에서 음수가 있습니다)', idx,[A_t_1,A_t_2,A,E_t_1,E_t_2,E,I_t_1,I_t_2,I,R,R_I,R_A])
            #raise NotImplementedError
    
    return ([S,E_G,E_A,E_B,I_G,I_A,I_B,A_G,A_A,A_B,Q,R,R_I,R_A,Q_total_init]
            ,[Q_init,Q_total_init])





# 함수 2 : 각 집단별 초기 parameter set 출력 , para1은 고정 상수라 필요x
# input : 시작할 날짜
# output : para2=[beta,lamda,c], para3=[omega 3개, theta 3개.] 단, 해당 일차부터.
def init_para(start_days, end_days):
    para2_mod=np.array(para2)
    para2_from_day_t=para2_mod[:,start_days-1:end_days-1]
    para2_from_day_t[1]=lamda[start_days:end_days] #02-01 람다맞추기
    para3_mod=np.array(para3)
    para3_from_day_t=para3_mod[:,start_days-1:end_days-1]
    return (para2_from_day_t,para3_from_day_t)

# Optimization을 실행 해야하는 변수 : GOOD,BAD 가중치와 THETA 
# X0=[S,E_G,E_A,E_B,I_G,I_A,I_B,A_G,A_A,A_B,Q,R,Q_total] # 시뮬레이션 시작 시 초기값.
# para1 = [alpha,row,sigma,delta,mu,gamma,epsilon_A] # 상수로 고정 된 parameter
# para2 = [beta,lamda,c] # time-table로 표현되는 parameter
# para3 = [omega 3개, theta 3개, p,fever)optimization 필요    
# opt : 임의로 init_val을 실행할 시. 인자로 opt(X,[q_new,q_total])이 온다. (1.1 쭉쓴다)

def globalpara(sigma_, delta_, row_,ep_A_,ep_):

    global sigma, delta, row, ep_A,ep,p,fever
    global prev_pa
    sigma=1/sigma_ # 제발 역수로 해야함
    delta=1/delta_
    row=row_
    ep_A=ep_A_
    ep=ep_
    
    #     p=p_
    #     fever=fever_
    #dealy 적용하기.
    return


def seir_simulation(start,end,opt=0): # (start)일부터 (end)일까지 simulating.
    # input : 시작일, 종료일
    # output : X=(S,E,A,I.....R)의 Time-Series
    para2_cut,para3_cut=init_para(start,end) # parameter limitation
    X0=init_val(start,opt)[0] # 시작일의 X[t]=[S[t],E[t],I[t],,,Q_total[t]]
    Q_init,Q_total=init_val(start,opt)[1] # Q신규와, Q누적의 초기 집단
    Total_Timeseries=np.zeros([len(X0),end-start+1]) # Total array i정의.
    Total_Timeseries[:,0]=X0 # 초기값 대입. 
    Q_Count=[] # 신규, 누적확진자 집계용, 실질적인 mse에 도입할 값.
    Q_Count.append([Q_init,Q_total])
    for i in range(end-start):
        #@##여기부턴.
        beta2,lamda2,c2=para2_cut[0][i],para2_cut[1][i],para2_cut[2][i]
        omega_G2,omega_A2,omega_B2=para3_cut[0][i],para3_cut[1][i],para3_cut[2][i] # 변수이동.
        theta_G2,theta_A2,theta_B2=para3_cut[3][i],para3_cut[4][i],para3_cut[5][i]
        p2,fever2=para3_cut[6][i],para3_cut[7][i]
        var_para1=[beta2,lamda2,c2]
        var_para2=[omega_G2,omega_A2,omega_B2,theta_G2,theta_A2,theta_B2,p2,fever2]
        Total_Timeseries[:,i+1],Q_counting=seir_onestep(Total_Timeseries[:,i],var_para1,var_para2) # <--함수에 쓰임
        Q_Count.append(Q_counting)
    Q_new=np.array(Q_Count)[:,0]
    Q_total=np.array(Q_Count)[:,1]
    if Total_Timeseries.min()<0:
        print('음수가 있습니다')
        #raise NotImplementedError
        
    return ([Total_Timeseries,[Q_new,Q_total]]) # 이거 Fullgrid-search땐 벗겨야함. 




##########################
###
### 시뮬레이티드 어닐링 + SEIR(실제론 SEHIQR)
###
##########################
def sa_seir_main(start,end,para_space,opt,init_state=0,**hyperpara): 
    default = {
            'tol_num':50, # n step간 해 개선이 없으면 종료
            'tol_val':1, # n step간 해 개선이 적으면 종료
            'schedule': 'Boltzmann', # Boltzmann 스케쥴, 좋지 못함, 현재 실행 x
            'bol' : 1, # Acceptance probability를 높임 (더 잘 이동함)
            'eta' : 1, #학습률, 이동을 더 세밀하게 함
            'alpha' : 0.985,
        'move' : 1,# Annealing schedule, 지수감쇄
    'stop':0.01} #온도가 일정 이하 내려가면 종료함.
    
    default.update(hyperpara)
    n_iter=default['n_iter']
    tol_num=default['tol_num']
    tol_val=default['tol_val']
    schedule=default['schedule']
    bol=default['bol']
    eta=default['eta']
    alpha=default['alpha']
    stop=default['stop']
    move=default['move']
    # 이거 꼭 필요한가 ? 
    
    step=0
    
    dim=len(para_space) #parameter의 차원
    best={'loss':1e+10,'state':'변수9~10개'} #[최적손실, 최적상태.]
    initial={'state1':'og,ob,beta,c 외 5~6개'} # parameter의 초기값 
    full_mse=[] # 손실값 추이 
    full_s=[] # 상태 추이
    full_temp=[] # 온도 추이 (개별실행 빼면 필요 없을듯)
    full_move_count=[]
    # IF 초기값이 주어진다면, (논문 투고 전 정밀하게 추정 위해서) 
    if(init_state):
        print('초기값이 추가로 입력되었습니다.', init_state)
        s=init_state
        init_eval=sa_eval(s,start,end,opt)
        temp=1.2*init_eval #초기 온도
        best['loss']=init_eval
        best['state']=s
        initial['state']=s
        full_s.append(s)
        full_mse.append(init_eval)
    else: # init(초기치가 입력되지 않았을 때, 초기값 정하기)
        while True: #초기값, 온도 설정
            s=sa_initialize(para_space)
            init_eval=sa_eval(s,start,end,opt)
            
            if (np.isnan(init_eval)):
                print('nan입니다.(초기값)')
                continue # nan이여도 초기값 다시 할당.
            if init_eval>1e7:
                continue # 초기값이 너무 크면 다시 할당
            temp=1.2*init_eval # 초기온도 지정
            best['loss']=init_eval
            best['state']=s
            initial['state']=s #초기값 저장.
            
            break #아니라면 while문 탈출
    #loss_s=sa_eval(s,start,end,opt) # 초기 loss 설정.
    loss_s=init_eval
    print('loss_s', loss_s)
    loss_for_count=0 #50번이상 정체되면 스탑.
    local=0 # local에 갇혔으면 스탑
    count=0 # 조기 종료 해 미개선 count용
    while(temp>stop): # stop 이하가 되면 stop. 
        move_count=0 # 내부 루프에서 이동 횟수(더 안 좋은 곳으로)
        for iter_ in range(int(n_iter)): # 같은 온도 내 반복(우선 static)
            s_new=sa_move(s,alpha,step,temp,eta,move) # 함수 내부 특성 상 temp,eta를 인자로 주는게 좋음.
            loss_s_new=sa_eval(s_new,start,end,opt)
            if(loss_s_new<best['loss']): #베스트용.
                best['loss']=loss_s_new
                best['state']=s_new
            
            if(loss_s_new<loss_s): # 더 좋은 방향
                s=s_new
                loss_s=loss_s_new
            else: # 더 안 좋은 방향 --> 확률로 이동 bol이 작으면 더 보수적으로 움직임.
                if temp<1 :
                    bol=0.3
                if (temp<0.1):
                    bol=0.01 # 그래도 마지막엔 경사하강만 진행하게끔.
                if (np.exp(-(loss_s_new-loss_s)/(bol*temp))>np.random.rand()): #우선 상수 x.---
                    move_count+=1
                    s=s_new 
                    loss_s=loss_s_new #일정확률에 따라 안 좋은 방향으로도 이동한다
        #print("worst_accept : ",move_count) #안 좋은 곳으로 이동한 횟수(어느 정도 있어야 함)
        full_move_count.append(move_count)
        if(loss_for_count==loss_s): #만약 이전 것과 지금 것의 loss가 같다면 ? 
            count+=1
        else:
            count=0
        loss_for_count=loss_s # 해당 loss, loss가 이전과 같다는 것은, 이동을 아예 못 한다는 것.
        #이는, acceptance 확률이 낮을 수도, 혹은 이동 보폭이 클 수도. 
        if(count==tol_num):
                    print("\nlocal minimum 의심, 시뮬레이션 종료")
                    break # tol_num번동안 개선이 없으면 stop.
        ''' 일단패스
        #local에 갇힌 것으로 의심될 때.
        if (step>200): # 200번 이상 돌렸을 때 best에 심각하게 접근하지 못 한다면, 
                if(loss_s>best['loss']):
                    local+=1
                else:
                    local=0
        if local==200: # 200번동안 개선이 없다면 out.
            print('best_loss에 접근하기 힘들다고 판단됨, best_loss:{}, 현재:{}'.format(best['loss'],loss_s))
            break
        # stop. 일단 패스
        '''
        step+=1 # step 증가.
        if(n_iter<=70):
            n_iter+=0.2 # 3번돌릴 때 마다 내부 step수를 늘린다.
        full_mse.append(loss_s)
        full_s.append(s)
        full_temp.append(temp)
        #온도가 클땐 빠르게 감소.
        if(temp>500):
            temp=temp/np.log(1+step/3)
        elif(temp>10):
            temp=temp*0.97 # 조금 빠르게 감소하되, 충분히 탐색할 시간을 준다. 
        else:
            temp=temp*alpha #그 이하는 지수감쇄.
        
        #print("step:{}, loss={}\n".format(step,loss_s))
        #print("temp : ", temp)
    full_dict=dict(mse=full_mse,x=full_s,temp=full_temp,
                   best=[best['state'],best['loss']],
                   initial=initial['state'],move_count=full_move_count)   
    return full_dict
  

#이동함수 - T와 s에 의존
#시간에 따라 보폭이 짧아야 함
#k는 운동상수. (클수록 잘 움직임 : 기본은 1)
def sa_move(s,alpha,step,temp,eta,move):
    s=np.array(s) 
    
    upper=np.array([x[1] for x in para_space]) #상한리스트
    lower=np.array([x[0] for x in para_space])
    lim=(upper-lower)/2
    dim=len(s)
    '''
    if(schedule=='Boltzmann'):
        std=np.amin([np.sqrt(temp)*np.ones(dim),(upper-lower)/(5)],axis=0) 
        #가능한가 ?  
        y=np.random.normal(0,std,size=dim)
    print(y)'''
    
    move_method=(1/(1+(move/temp))) #산술이동 --> 빠르게 감소해서 후에 이동은 잘하지만(정밀) temp클 때 구림.
    #30은, 온도가 30일 때 부터 감소가 유의미하게끔. (1로 설정하면 온도가 1이상일 땐 아무것도 못함)
    #move_method=(1/(1+(np.log(1+temp)))) #로그 이동 --> 너무 느리게 감소함, 그래서 이동을 못함.
    out_count=0
    while(1):
       
        y=np.random.uniform(-lim,lim,size=dim)*(move_method) # 이 함수를 어떻게 할 것인가.
        if out_count>500:
            move_method=(1/(1+(move/temp*1.1)))
            out_count=0
            #print('500번이 넘었다.. 뭐냐..')
        # ****이동을 잘 못한다 --> eta를 줄여서 세밀하게
        # ****이동을 너무 잘한다 --? eta를 늘려서 보폭을 크게(이동하기 어려움)-temp가 줄어들면서 처리됨
        s_new=s+eta*y
        if s_new[5]<s_new[6]:
            out_count+=1
            #print('s[5]<s[6]')
            
            continue
        elif s_new[9]>s_new[10]:
            out_count+=1
            #print('s[9]>s[10]')
            continue
        if (s_new[5]+s_new[9])>1:
            out_count+=1
            #print('theta 합 > 1')
            continue
        elif (s_new[6]+s_new[10])>1:
            out_count+=1
            #print('theta 합 > 1')
            continue
        else:
            break
    
    #s_new=s+eta*y
    for x in range(len(s_new)):
        if(s_new[x]<para_space[x][0]): #범위가 좌측으로 벗어남.
            s_new[x]=s[x]-eta*y[x]
        elif(s_new[x]>para_space[x][1]): #범위가 우측으로 벗어남 
            s_new[x]=s[x]-eta*y[x]
        else:
            pass
        
    return s_new

import random
#input : parameter spaced와 초기 Temperature 설정방법
#output : 초기의 상태.
def sa_initialize(parameter_space):
    s=[]
    for (a,b) in parameter_space:
        s.append(random.uniform(a,b))
    
    return s

# para_space  1: ["omega_G", +om_A , omega_B", "beta", "c",]
#para_space 2 : [(upper,lower,dwell)*2,sigma.delta,row,ep_A]

def sa_eval(state,start,end,opt,loss='new'): #start,end *arg로 합쳐야함.
    global init_I
    if loss=='new':
        real_new=slicing(start,end)[0] # 실제 확진자수를 불러옴(누적으로)
    if loss=='total':
        real_total=slicing(start,end)[1] #
    forgridsearch_parameter(start,end,*(np.append(state[:5],state[-2:]))) # og,oa,ob,beta,c
    #if (opt==1) : #state 1 (1~39)일 할 때. 
    #    theta_first_scaling(start,end,state[4]) 
    #    globalpara(sigma_=state[5],delta_=state[6],row_=state[7],ep_A_=state[8]) # sigam,delta,row,ep_A
    # 위에는 사용하지 않음.  
    global delay
    
    theta_after_scaling(start,end,state[5],state[6],state[7],state[8],state[9],state[10],state[11],state[12],delay=delay) # sigmoid 의 upper.lower,dwel
    globalpara(sigma_=state[13],delta_=state[14],row_=state[15],ep_A_=state[16],ep_=state[17]) # sigam,delta,row,ep_A
    if (opt!=1): # 중간값을 줘야 하는 경우.
        if loss=='new':
            pred=seir_simulation(start,end,opt)[1][0] #이거 가볍게 해야함.
            RMSE=np.sqrt(MSE(pred,real_new))
        else :
            pred=seir_simulation(start,end,opt)[1][1] #이거 가볍게 해야함.
            RMSE=np.sqrt(MSE(pred,real_total))
            
    else: #첫번째 일 땐 이어서 하지 않기 때문에 opt가 들어가지 않는다. - 윈도우 별 첫 시도.
        if loss=='new':
            pred=seir_simulation(start,end)[1][0] #이거 가볍게 해야함.
            RMSE=np.sqrt(MSE(pred,real_new))
        else : #누적으로 fitting할 때 ? (변수 2개)
            pred=seir_simulation(start,end)[1][1] #이거 가볍게 해야함.
            RMSE=np.sqrt(MSE(pred,real_total))
    return RMSE


#각 Time-window 별 첫실행시!! 
# real_first=0이 기본, 0이 아니라면 upper와 lower는 모두 같음.(즉, Window 1의 첫번째 경우.)
def sa_seir_main_first(start,end,para_space,opt,init_state=0,real_first=0,**hyperpara): 
    default = {
            'tol_num':50, # n step간 해 개선이 없으면 종료
            'tol_val':1, # n step간 해 개선이 적으면 종료
            'schedule': 'Boltzmann', # Boltzmann 스케쥴, 좋지 못함, 현재 실행 x
            'bol' : 1, # Acceptance probability를 높임 (더 잘 이동함)
            'eta' : 1, #학습률, 이동을 더 세밀하게 함
            'alpha' : 0.985,
        'move':1,# Annealing schedule, 지수감쇄
    'stop':0.01} #온도가 일정 이하 내려가면 종료함.
    
    default.update(hyperpara)
    n_iter=default['n_iter']
    tol_num=default['tol_num']
    tol_val=default['tol_val']
    schedule=default['schedule']
    bol=default['bol']
    eta=default['eta']
    alpha=default['alpha']
    stop=default['stop']
    move=default['move']
    # 이거 꼭 필요한가 ? 
    
    step=0
    
    dim=len(para_space) #parameter의 차원
    best={'loss':1e+10,'state':'변수9~10개'} #[최적손실, 최적상태.]
    initial={'state1':'og,ob,beta,c 외 5~6개'} # parameter의 초기값 
    full_mse=[] # 손실값 추이 
    full_s=[] # 상태 추이
    full_temp=[] # 온도 추이 (개별실행 빼면 필요 없을듯)
    full_move_count=[]
    # IF 초기값이 주어진다면, (논문 투고 전 정밀하게 추정 위해서) 
    if(init_state):
        print('초기값이 추가로 입력되었습니다.', init_state)
        s=init_state
        init_eval=sa_eval(s,start,end,opt)
        temp=1.2*init_eval #초기 온도
        best['loss']=init_eval
        best['state']=s
        initial['state']=s
        full_s.append(s)
        full_mse.append(init_eval)
    else: # init(초기치가 입력되지 않았을 때, 초기값 정하기)
        while True: #초기값, 온도 설정
            s=sa_initialize(para_space)
            #_first 함수의 차별점(이하5줄)
            
            global theta_G_prev, theta_B_prev
            theta_G_prev=s[5] #선형성 없애기 위하여. 
            theta_B_prev=s[9]
            
            if real_first!=0: # Window 1의 첫 시행 때만 상수로 가정한다. 
                s[6]=s[5]
                s[10]=s[9]

            init_eval=sa_eval(s,start,end,opt)
            if (np.isnan(init_eval)):
                print('nan입니다.(초기값)')
                continue # nan이여도 초기값 다시 할당.
            if init_eval>1e7:
                continue # 초기값이 너무 크면 다시 할당
            temp=1.2*init_eval # 초기온도 지정
            best['loss']=init_eval
            best['state']=s
            initial['state']=s #초기값 저장.
            
            break #아니라면 while문 탈출
    #loss_s=sa_eval(s,start,end,opt) # 초기 loss 설정.
    loss_s=init_eval
    print('loss_s', loss_s)
    loss_for_count=0 #50번이상 정체되면 스탑.
    local=0 # local에 갇혔으면 스탑
    count=0 # 조기 종료 해 미개선 count용
    while(temp>stop): # stop 이하가 되면 stop. 
        move_count=0 # 내부 루프에서 이동 횟수(더 안 좋은 곳으로)
        for iter_ in range(int(n_iter)): # 같은 온도 내 반복(우선 static)
            s_new=sa_move(s,alpha,step,temp,eta,move) # 함수 내부 특성 상 temp,eta를 인자로 주는게 좋음.
            
            # 01-01 추가. eval전에 진행해야함.
            theta_G_prev=s_new[5] #선형성 없애기 위하여. _
            theta_B_prev=s_new[9]
            if real_first!=0: # 맨 최초 시작 때는 선형 가정.
                s_new[6]=s_new[5]
                s_new[10]=s_new[9]
                
            loss_s_new=sa_eval(s_new,start,end,opt)
            if(loss_s_new<best['loss']): #베스트용.
                best['loss']=loss_s_new
                best['state']=s_new
            
            if(loss_s_new<loss_s): # 더 좋은 방향
                s=s_new
                loss_s=loss_s_new
            else: # 더 안 좋은 방향 --> 확률로 이동 bol이 작으면 더 보수적으로 움직임.
                if temp<1 :
                    bol=0.3
                if (temp<0.1):
                    bol=0.01 # 그래도 마지막엔 경사하강만 진행하게끔.
                if (np.exp(-(loss_s_new-loss_s)/(bol*temp))>np.random.rand()): #우선 상수 x.---
                    move_count+=1
                    s=s_new 
                    loss_s=loss_s_new #일정확률에 따라 안 좋은 방향으로도 이동한다
        #print("worst_accept : ",move_count) #안 좋은 곳으로 이동한 횟수(어느 정도 있어야 함)
        full_move_count.append(move_count)
        if(loss_for_count==loss_s): #만약 이전 것과 지금 것의 loss가 같다면 ? 
            count+=1
        else:
            count=0
        loss_for_count=loss_s # 해당 loss, loss가 이전과 같다는 것은, 이동을 아예 못 한다는 것.
        #이는, acceptance 확률이 낮을 수도, 혹은 이동 보폭이 클 수도. 
        if(count==tol_num):
                    print("\nlocal minimum 의심, 시뮬레이션 종료")
                    break # tol_num번동안 개선이 없으면 stop.
        ''' 일단패스
        #local에 갇힌 것으로 의심될 때.
        if (step>200): # 200번 이상 돌렸을 때 best에 심각하게 접근하지 못 한다면, 
                if(loss_s>best['loss']):
                    local+=1
                else:
                    local=0
        if local==200: # 200번동안 개선이 없다면 out.
            print('best_loss에 접근하기 힘들다고 판단됨, best_loss:{}, 현재:{}'.format(best['loss'],loss_s))
            break
        # stop. 일단 패스
        '''
        step+=1 # step 증가.
        if(n_iter<=70):
            n_iter+=0.2 # 3번돌릴 때 마다 내부 step수를 늘린다.
        full_mse.append(loss_s)
        full_s.append(s)
        full_temp.append(temp)
        #온도가 클땐 빠르게 감소.
        if(temp>500):
            temp=temp/np.log(1+step/3)
        elif(temp>10):
            temp=temp*0.97 # 조금 빠르게 감소하되, 충분히 탐색할 시간을 준다. 
        else:
            temp=temp*alpha #그 이하는 지수감쇄.
        
        #print("step:{}, loss={}\n".format(step,loss_s))
        #print("temp : ", temp)
    full_dict=dict(mse=full_mse,x=full_s,temp=full_temp,
                   best=[best['state'],best['loss']],
                   initial=initial['state'],move_count=full_move_count)   
    return full_dict


def sa_main_iter(start,end,opt,para_space,gridding,init_state=0,**hyperpara): #hyperparameter 설정해야함. 
    #gridding
    #global init_I
    #init_I=1
    
    iternum=len(gridding)
    full_dict={}
    best_mse=[] # (internum)차원
    best_x=[] # (iternum)차원
    period_array=[]
    for i in range(iternum):
        period_array.append([start,end])
    initial_point=[]
    full_mse=[] # full mse array, (iternum)*(step)
    full_temp=[] # full temp array (iternum)*step
    final_mse=[] # 수렴을 잘 했는지 # (iternum) 차원
    final_x=[] # 수렴의 마지막 x.
    full_move=[] #얼마나 worst move를 했는지. 
    for i in range(iternum):
        alpha,eta,bol,move=gridding[i]
        print("\n\nSimulation :{}, alpha:{}, eta:{}, bol:{}, move:{}\n\n".format(i+1,alpha,eta,bol,move))
        #패러미터 다섯개임.
        alpha,eta,bol,move=gridding[i]
        default ={'n_iter':50, 'tol_num':80,'eta':eta 
                ,'schedule':'Boltzmann','alpha':alpha,
                'bol':bol,'move':move,'stop':0.001}

        default.update(hyperpara)
        hyperpara=default
        sa_pred=sa_seir_main(start,end,para_space,opt,init_state=init_state,**hyperpara)
        best_mse.append(sa_pred['best'][1])
        best_x.append(sa_pred['best'][0])
        
        full_mse.append(sa_pred['mse']) # full mse
        full_temp.append(sa_pred['temp'])
        
        final_mse.append(sa_pred['mse'][-1]) # 마지막 mse
        final_x.append(sa_pred['x'][-1]) # 마지막 parameter
        
        initial_point.append(sa_pred['initial'])
        full_move.append(sa_pred['move_count'])
        print('결과 : ', np.round(sa_pred['best'][0],3), '-->' ,sa_pred['best'][1])
        
    full_dict=dict(period=period_array,best_mse=best_mse,best_x=best_x,hyperpara=gridding
                   ,final_mse=final_mse,final_x=final_x,
                   initial=initial_point,full_move=full_move
                   ,full_mse=full_mse,full_temp=full_temp
                  )
        
    

    return full_dict

# 각 Time-window 별 첫 실행

def sa_main_iter_first(start,end,opt,para_space,gridding,init_state=0,real_first=0,**hyperpara): #hyperparameter 설정해야함. 
    #gridding
    #global init_I
    #init_I=1
    
    iternum=len(gridding)
    full_dict={}
    best_mse=[] # (internum)차원
    best_x=[] # (iternum)차원
    period_array=[]
    for i in range(iternum):
        period_array.append([start,end])
    initial_point=[]
    full_mse=[] # full mse array, (iternum)*(step)
    full_temp=[] # full temp array (iternum)*step
    final_mse=[] # 수렴을 잘 했는지 # (iternum) 차원
    final_x=[] # 수렴의 마지막 x.
    full_move=[] #얼마나 worst move를 했는지. 
    for i in range(iternum):
        alpha,eta,bol,move=gridding[i]
        print("\n\nSimulation :{}, alpha:{}, eta:{}, bol:{}, move:{}\n\n".format(i+1,alpha,eta,bol,move))
        #패러미터 다섯개임.
        alpha,eta,bol,move=gridding[i]
        default ={'n_iter':50, 'tol_num':80,'eta':eta 
                ,'schedule':'Boltzmann','alpha':alpha,
                'bol':bol,'move':move,'stop':0.001}

        default.update(hyperpara)
        hyperpara=default
        sa_pred=sa_seir_main_first(start,end,para_space,opt,init_state=init_state,real_first=real_first,**hyperpara)
        best_mse.append(sa_pred['best'][1])
        best_x.append(sa_pred['best'][0])
        
        full_mse.append(sa_pred['mse']) # full mse
        full_temp.append(sa_pred['temp'])
        
        final_mse.append(sa_pred['mse'][-1]) # 마지막 mse
        final_x.append(sa_pred['x'][-1]) # 마지막 parameter
        
        initial_point.append(sa_pred['initial'])
        full_move.append(sa_pred['move_count'])
        print('결과 : ', sa_pred['best'][0], '-->' ,sa_pred['best'][1])
        
    full_dict=dict(period=period_array,best_mse=best_mse,best_x=best_x,hyperpara=gridding
                   ,final_mse=final_mse,final_x=final_x,
                   initial=initial_point,full_move=full_move
                   ,full_mse=full_mse,full_temp=full_temp
                  )
        
    

    return full_dict


#### 파이프라인
def sa_pipe(start,end,theta_start,theta_end,best_parameter,opt=0,full_array=0,prev_dict=0): # opt : output, 
    global init_I
    if (full_array==0): # 첫번째 회차.
        full_array=[]
    forgridsearch_parameter(start,end,*np.append(best_parameter[:5],best_parameter[-2:])) #변수설정
    if (opt==0 or full_array==[]):# 첫실행
        #잠시만, 이거 안 씀.
        #theta_first_scaling(start,end,0.980)
        print("Window 별 첫 실행")
        global theta_G_prev, theta_B_prev
        theta_G_prev=best_parameter[5] #선형성 없애기 위하여. 
        theta_B_prev=best_parameter[9]
        global delay
        theta_after_scaling(theta_start,theta_end,best_parameter[5],best_parameter[6],best_parameter[7],
                           best_parameter[8],best_parameter[9],best_parameter[10],best_parameter[11],best_parameter[12],delay=delay)
        globalpara(best_parameter[13],best_parameter[14],best_parameter[15],best_parameter[16],best_parameter[17]) #sigma~epsilon~ep~p
        full,pred=seir_simulation(start,end) 
        full_store_collection=full #[:,start-1:end] #변형되면 안 되는 배열
        full_store_counting=np.array(pred) #[:,start-1:end] 
    elif (opt!=0 and full_array!=0):
        print("Window 별 첫 실행 이후")
        theta_after_scaling(theta_start,theta_end,best_parameter[5],best_parameter[6],best_parameter[7],
                           best_parameter[8],best_parameter[9],best_parameter[10],best_parameter[11],best_parameter[12],delay=delay)
        globalpara(best_parameter[13],best_parameter[14],best_parameter[15],best_parameter[16],
                  best_parameter[17])
        full,pred=seir_simulation(start,end,opt)
        full_store_collection=np.concatenate((full_array[0][:,:-1],full),axis=1)
        full_store_counting=np.concatenate((full_array[1][:,:-1],pred),axis=1)
    else :
        print("매개인자 입력 제대로 하기.")
    X_cont=full_store_collection[:,-1] # 마지막 집단 그대로 input으로 사용.
    Q_cont=full_store_counting[:,-1] #Q_new, Q_total
    full_array=[full_store_collection,full_store_counting]
    next_input=[X_cont,Q_cont]
    # history 저장
    parameter_name=['omega_G','omega_A','omega_B','beta','c','theta_G',
                    'theta_A','theta_B','sigma','delta','row','ep_A','ep','p','fever'] # 1-10
    parameter_variable=[omega_G,omega_A,omega_B,beta,c,theta_G,theta_A,theta_B,sigma,delta,row,ep_A,ep,p,fever] # 꼭 이렇게 해야하나.
    variable_name=['S','E_G','E_A','E_B','I_G','I_A','I_B','A_G','A_A','A_B','Q','R'] # 하나 빠짐
    variable2_name=['Q_new','Q_total']
    
    
    # 총 Tiemseries.
    full_dict={} # 모든 변수.
    for idx, x in enumerate(parameter_name):
        if (type(parameter_variable[idx])==float) or (type(parameter_variable[idx])==np.float64): # sigma,delta,row 등 float으로 정의하였을 때.
            full_dict[x]=np.full(end-start+1,parameter_variable[idx])
        else:
            full_dict[x]=parameter_variable[idx][start-1:end] # 마지막은 사실상 안 씀..(init_para에서 1개 빼고 정의했으니깐.) - 어차피 후에 사라짐.

    for idx,x in enumerate(variable_name):
        full_dict[x]=full[idx] #[start-1:end]
    for idx,x in enumerate(variable2_name):
        full_dict[x]=pred[idx] #[start-1:end]
    full_dict['abroad']=lamda[start-1:end] # 위에 식도 변형했으므로, 간단히 맞춰주자.
    if opt==0:
        full_dict['input']=[init_I]+list(np.full(end-start,0)) # [1,0,0,0...으로 패딩]
    else:
        full_dict['input']=[opt]+list(np.full(end-start,0))
    full_dict['output']=list(np.full(end-start,0))+[next_input]
    full_dict['best_parameter']=[best_parameter]+list(np.full(end-start,0))
    full_dict['full_array']=list(np.full(end-start,0))+[full_array]
    full_dict['dwell_g']=[best_parameter[7]]+list(np.full(end-start,0))
    full_dict['Q_g']=[best_parameter[8]]+list(np.full(end-start,0))
    full_dict['dwell_b']=[best_parameter[11]]+list(np.full(end-start,0))
    full_dict['Q_b']=[best_parameter[12]]+list(np.full(end-start,0))
    next_theta_prev=(full_dict['theta_G'][-2], full_dict['theta_B'][-2])
    full_dict=pd.DataFrame(full_dict,index=list(range(start,end+1)))
    # prev dict의 output,full array도 저장하자.
    if type(prev_dict)!=int : # 이후 실행 때. 
        prev_output=prev_dict['output'].iloc[-1]
        prev_full_array=prev_dict['full_array'].iloc[-1]
        full_dict=pd.DataFrame(full_dict,index=list(range(start,end+1))) # dataframe.
        full_dict['output']=[prev_output]+list(np.full(end-start-1,0))+[next_input] # output조정.
        full_dict['full_array']=[prev_full_array]+list(np.full(end-start-1,0))+[full_array]
        # 그 다음 실행 시.
        full_dict=pd.concat([prev_dict[:-1],full_dict]) #한칸 떼고 합친다. 

    return next_input, full_array, full_dict, next_theta_prev   # return이 곧바로 다음 input이다. 