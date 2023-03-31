import pandas as pd
import numpy as np
import streamlit as st
import pandas as pd
import fluids
import matplotlib.pyplot as plt
import neqsim
from neqsim.thermo.thermoTools import *
from neqsim.process.processTools import *
from neqsim import methods
from neqsim.thermo import fluid, TPflash, createfluid2
from neqsim.process import pipe, pipeline, clearProcess, stream, runProcess
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import jpype 

url ='http://raw.githubusercontent.com/Ahmedhassan676/pressure_drop/main/table.csv'
df_gas = pd.read_csv(url, index_col=[0])
url_1 = 'http://raw.githubusercontent.com/Ahmedhassan676/pressure_drop/main/comp.csv'
df_comp_table = pd.read_csv(url_1)
url_2 ='http://raw.githubusercontent.com/Ahmedhassan676/pressure_drop/main/summary.csv'
df_summary = pd.read_csv(url_2, index_col=[0])
url_3 ='http://raw.githubusercontent.com/Ahmedhassan676/pressure_drop/main/table_liq.csv'
df_liq = pd.read_csv(url_3, index_col=[0])

def convert_data(df):
     csv = df.to_csv(index=False).encode('utf-8')
     return csv
def find_nearest(D):
    array = np.array([0.5,0.75,1,1.5,2,3,4,5,6,8,10,12,14,16,18,20,24])
    idx = (np.abs(array - D)).argmin()
    return array[idx]
def k_calculations(df,df_comp_table,suc_t,disch_t):
        
        temperatures = np.array([suc_t,disch_t])*1.8+ 32 
        df['y_MCp_suc']=[np.interp(temperatures[0],df_comp_table['Gas'][3:],df_comp_table['{}'.format(compound)][3:]) for compound in df.index]
        suc_MCp = (df['mol%']*df['y_MCp_suc']).sum()*0.01
        k_suc = suc_MCp/(suc_MCp-1.986)
        df['y_MCp_disch']=[np.interp(temperatures[1],df_comp_table['Gas'][3:],df_comp_table['{}'.format(compound)][3:]) for compound in df.index]
        disch_MCp = (df['mol%']*df['y_MCp_disch']).sum()*0.01
        k_disch = disch_MCp/(disch_MCp-1.986)
        k = (k_suc + k_disch)/2

        
        return k
def Summary_calculations(Q_std,D,G,mu,f_E,p1,p2,t,m_wt,k,rho2,L,z,dp_100m):
        tb = 273.15+15.55556
        Pb = 101.325
        D_mm = D * 1000
        A = np.pi*0.25*((D_mm/1000)**2)
        Q_normal = Q_std*(273.15/(273.15+15))
        Q_actual = Q_std*(1.033023/p2)*((t+273.15)/tb)
        v = Q_actual/(A*3600) #Velocity (m/s)
        sonic_velocity =((9.81*k*847.9*(t+273.15))/m_wt)**0.5
        Re = 0.5134*(Pb/tb)*((G*Q_std*24)/(mu*0.01*D_mm))
        mach = v/sonic_velocity
        dp_percent = ((p1-p2)/p1)*100
        rho_v_2 = rho2*(v**2)
        summary_list = [p1,p2,t,L,D,Q_std,Q_normal,Q_actual,dp_percent,dp_100m,f_E,Re,v,sonic_velocity,mach,rho_v_2,m_wt,z,k,1/np.sqrt(k),mu]
        return summary_list
def Z_calculations(df,t_suc,p_suc):
        pc = np.sum(df['mol%']*df['Pc']) * 0.01
        tc = np.sum(df['mol%']*(df['Tc']+460)) * 0.01  
        m_wt = np.sum(df['mol%']*df['m.wt'])*0.01
        Tr = (t_suc*1.8 + 32+460)/tc
        Pr = (p_suc +1.03323)*14.2233/pc
        A = [1,0.31506237,-1.04670990,-0.57832729,0.53530771,-0.61232032,-0.10488813,0.68157001,0.68446549]
        rho = (0.27*Pr)/Tr
        Z = 1
        error = 10
        y = 1/0.9
        while error > 0.001:
            part_1 = (A[1]+(A[2]/Tr)+(A[3]/(Tr**3)))*rho*y
            part_2 = (A[4]+(A[5]/Tr))*(rho**2)*(y**2)
            part_3 = ((A[5]*A[6]*(rho**5)*(y**5))/Tr)
            part_4 = (A[7]*(rho**3)*(y**3))/((Tr**3)*(1+(A[8]*(rho**2)*(y**2)))*np.exp(-A[8]*(rho**2)*(y**2)))
            Z = A[0]+part_1+part_2+ part_3 + part_4 
            error = abs(Z-(1/y))
            y = 1/Z 
        
        return Z, m_wt
def choose_composition():
            
            df = pd.DataFrame({'Composition/property':df_comp_table.columns[1:],'mol%':np.zeros(len(df_comp_table.columns)-1), 'm.wt':df_comp_table.iloc[0,1:],'Pc':df_comp_table.iloc[1,1:],'Tc':df_comp_table.iloc[2,1:]})
            try:
                sum_of_comp = 0 
                c1,c2,c3,c4,c5,c6,c7,c8,c9,c16,c10,c11,c13,h2o = 0,0,0,0,0,0,0,0,0,0,0,0,0,0
                options_list = [df_comp_table.columns[i] for i in [22,1,4,6,11,10,13,15,20,21,24,25,23,19]]
                while sum_of_comp != 100:
                    options = st.multiselect(
                    'Select your components', options_list
                    )
                    if df_comp_table.columns[22] in options:
                        c1 = st.number_input('hydrogen%', key = 'c1')
                    if df_comp_table.columns[1] in options:
                        c2 = st.number_input('methane%', key = 'c2')
                    if df_comp_table.columns[4] in options:
                        c3 = st.number_input('ethane%', key = 'c3')
                   # if df_comp_table.columns[3] in options:
                    #    c15 = st.number_input('ethylene%', key = 'c15')
                    if df_comp_table.columns[6] in options:
                        c4 = st.number_input('propane%', key = 'c4')
                    if df_comp_table.columns[11] in options:
                        c5 = st.number_input('nbutane%', key = 'c5')
                    if df_comp_table.columns[10] in options:
                        c6 = st.number_input('isobutane%', key = 'c6')
                    if df_comp_table.columns[13] in options:
                        c7 = st.number_input('pentane%', key = 'c7')
                    if df_comp_table.columns[15] in options:
                        c8 = st.number_input('hexane%', key = 'c8')
                    if df_comp_table.columns[20] in options:
                        c9 = st.number_input('Oxygen%', key = 'c9')
                    if df_comp_table.columns[21] in options:
                        c16 = st.number_input('nitrogen%', key = 'c16')
                    if df_comp_table.columns[24] in options:
                        c10 = st.number_input('carbon monoxide%', key = 'c10')
                    if df_comp_table.columns[25] in options:
                        c11 = st.number_input('carbon dioxide%', key = 'c11')
                    if df_comp_table.columns[23] in options:
                        c13 = st.number_input('hydrogen sulphide%', key = 'c13')
                    #if df_comp_table.columns[17] in options:
                     #   nh3 = st.number_input('Ammonia%', key = 'nh3')
                    if df_comp_table.columns[19] in options:
                        h2o = st.number_input('Water vapor%', key = 'h2o')
                    if c1 or c2 or c3 or c4 or c5 or c6 or c7 or c8 or c9 or c16 or c10 or c11 or c13 or h2o:
                        c = []
                        for i in (c1,c2,c3,c4,c5,c6,c7,c8,c9,c16,c10,c11,c13,h2o):
                            c.append(i)
                        
                        for (i,j) in zip(options_list,c):
                            if j != None:
                                    df.loc[i,'mol%'] = j
                        
                        sum_of_comp = np.sum(df['mol%'])
                        
                st.success('Composition in Mol. percent completed!', icon="✅")
                
                return df[df['mol%'] != 0] ,True

            except (ValueError, st.errors.DuplicateWidgetID): pass
            except (TypeError, KeyError, ZeroDivisionError):st.write('Please Check your data')
def NeqSim_calculations(q,D,df_comp,t,p1,p2,L,type):
        if type == 'estimate quantity': 
            p1 = p1*98.066*0.01
            p2 = p2*98.066*0.01
            q = q*24/(1000000)
            #Creating inlet fluid using SRK-EoS   
            names =  list(df_comp.index.str.lower())
            molefractions = list(df_comp['mol%']*0.01)
            fluid1 = createfluid2(names, molefractions)
            fluid1.setMixingRule('classic')
            fluid1.setTemperature(t, "C")
            fluid1.setPressure(p1, "bara")
            fluid1.setTotalFlowRate(q, "MSm3/day")
            
            TPflash(fluid1)

            

            deltaElevation = 0.0
            pipeLength = L
            roughness= 0.00005
            diameter = D
            error = 10
            method = "friction theory"
            fluid1.getPhase('gas').getPhysicalProperties().setViscosityModel(method)
            fluid1.initProperties()
            while error > 0.01:
                fluid1.setTotalFlowRate(q, "MSm3/day")
                clearProcess()
                
                stream1 = stream(fluid1)
                
                pipeSimple = pipe(stream1, pipeLength, deltaElevation, diameter, roughness)
                runProcess()
                
                
                error = abs(p2 - pipeSimple.getOutStream().getFluid().getPressure('bara'))
                if p2 - pipeSimple.getOutStream().getFluid().getPressure('bara') > 0.01:
                    q = q - q*0.01
                elif p2 - pipeSimple.getOutStream().getFluid().getPressure('bara') < -0.01:
                    q = q + q*0.01
                else: error = error
                result = q*1000000/24
        if type == "estimate upstream pressure":
            p1 = p1*98.066*0.01
            p2 = p2*98.066*0.01
            q = q*24/(1000000)
            #Creating inlet fluid using SRK-EoS   
            names =  list(df_comp.index.str.lower())
            molefractions = list(df_comp['mol%']*0.01)
            fluid1 = createfluid2(names, molefractions)
            fluid1.setMixingRule('classic')
            fluid1.setTemperature(t, "C")
            fluid1.setPressure(p1, "bara")
            fluid1.setTotalFlowRate(q, "MSm3/day")
            
            TPflash(fluid1)

            

            deltaElevation = 0.0
            pipeLength = L
            roughness= 0.00005
            diameter = D
            error = 10
            method = "friction theory"
            fluid1.getPhase('gas').getPhysicalProperties().setViscosityModel(method)
            fluid1.initProperties()
            while error > 0.01:
                fluid1.setPressure(p1, "bara")
                clearProcess()
                
                stream1 = stream(fluid1)
                
                pipeSimple = pipe(stream1, pipeLength, deltaElevation, diameter, roughness)
                runProcess()  
                error = abs(p2 - pipeSimple.getOutStream().getFluid().getPressure('bara'))
                if p2 - pipeSimple.getOutStream().getFluid().getPressure('bara') > 0.01:
                    p1 = p1 + p1*0.01
                elif p2 - pipeSimple.getOutStream().getFluid().getPressure('bara') < -0.01:
                    p1 = p1 - p1*0.01
                else: error = error
            result = pipeSimple.getInStream().getFluid().getPressure('bara')/0.9806
        if type == "estimate downstream pressure":
            p1 = p1*98.066*0.01
            p2 = p2*98.066*0.01
            q = q*24/(1000000)
            #Creating inlet fluid using SRK-EoS   
            names =  list(df_comp.index.str.lower())
            molefractions = list(df_comp['mol%']*0.01)
            fluid1 = createfluid2(names, molefractions)
            fluid1.setMixingRule('classic')
            fluid1.setTemperature(t, "C")
            fluid1.setPressure(p1, "bara")
            fluid1.setTotalFlowRate(q, "MSm3/day")
            
            TPflash(fluid1)

            

            deltaElevation = 0.0
            pipeLength = L
            roughness= 0.00005
            diameter = D
            error = 10
            method = "friction theory"
            fluid1.getPhase('gas').getPhysicalProperties().setViscosityModel(method)
            fluid1.initProperties()
               
            clearProcess()
            
            stream1 = stream(fluid1)
            
            pipeSimple = pipe(stream1, pipeLength, deltaElevation, diameter, roughness)
            runProcess()  
            
            result = pipeSimple.getOutStream().getFluid().getPressure('bara')/0.9806
        return result
def Calculate_OutVelocity(D,Q_std,p2,t):
        tb = 273.15+15.55556
        
        D_mm = D * 1000
        A = np.pi*0.25*((D_mm/1000)**2)
        
        Q_actual = Q_std*(1.033023/p2)*((t+273.15)/tb)
        v = Q_actual/(A*3600) #Velocity (m/s) 
        return v
def graph_NeqSim(q,D,df_comp,t,p1,L):
            
            p2 = []
            p2.append(p1)
            L_list = []
            L_list.append(0)
            velocity = []
            velocity.append(Calculate_OutVelocity(D,q,p1,t))
            p1 = p1*98.066*0.01
            
            q = q*24/(1000000)
            #Creating inlet fluid using SRK-EoS   
            names =  list(df_comp.index.str.lower())
            molefractions = list(df_comp['mol%']*0.01)
            fluid1 = createfluid2(names, molefractions)
            fluid1.setMixingRule('classic')
            fluid1.setTemperature(t, "C")
            
            fluid1.setTotalFlowRate(q, "MSm3/day")
            
            TPflash(fluid1)

            

            deltaElevation = 0.0
            pipeLength = L*0.1
            roughness= 0.00005
            diameter = D
            
            method = "friction theory"
            fluid1.getPhase('gas').getPhysicalProperties().setViscosityModel(method)
            fluid1.initProperties()
            if 'hydrogen' in names:
                 for i in range(2):
                    pipeLength = L*0.5
                    fluid1.setPressure(p1, "bara") 
                    clearProcess()
                    
                    stream1 = stream(fluid1)
                    
                    pipeSimple = pipe(stream1, pipeLength, deltaElevation, diameter, roughness)
                    runProcess()
                    if p2[-1] != pipeSimple.getOutStream().getFluid().getPressure('bara')/0.9806:
                        p2.append(pipeSimple.getOutStream().getFluid().getPressure('bara')/0.9806)
                        velocity.append(Calculate_OutVelocity(D,q*1000000/24,p2[-1],t))
                        L_list.append((1/3)*len(p2)*L)
                        p1 = pipeSimple.getOutStream().getFluid().getPressure('bara')
            else:
                for i in range(9):
                    fluid1.setPressure(p1, "bara") 
                    clearProcess()
                    
                    stream1 = stream(fluid1)
                    
                    pipeSimple = pipe(stream1, pipeLength, deltaElevation, diameter, roughness)
                    runProcess()
                    if p2[-1] != pipeSimple.getOutStream().getFluid().getPressure('bara')/0.9806:
                        p2.append(pipeSimple.getOutStream().getFluid().getPressure('bara')/0.9806)
                        velocity.append(Calculate_OutVelocity(D,q*1000000/24,p2[-1],t))
                        L_list.append(0.1*len(p2)*L)
                        p1 = pipeSimple.getOutStream().getFluid().getPressure('bara')
            if L_list[-1] != L:
                 st.warning('Check your dp percent or volume flow rate input as Neqsim couldnt fully converge')
                      
            
            fig,axs = plt.subplots(1,2)
            fig.set_figheight(4)
            fig.set_figwidth(12)
            axs[0].plot(L_list,p2)
            axs[0].grid()
            axs[0].set_xlabel("Length (m)")
            axs[0].set_ylabel("Pressure (kg/cm2.a)")

            axs[1].plot(L_list,velocity)
            axs[1].grid()
            axs[1].set_xlabel("Length (m)")
            axs[1].set_ylabel("Velocity (m/s)")
            st.pyplot(fig)   
            
def get_viscosity(df_comp,p1,t):
        
        
        #Creating inlet fluid using SRK-EoS
        names =  list(df_comp.index.str.lower())
        molefractions = list(df_comp['mol%']*0.01)
        fluid1 = createfluid2(names, molefractions)
        fluid1.setMixingRule('classic')
        fluid1.setTemperature(t, "C")
        fluid1.setPressure(p1, "bara")
        fluid1.setTotalFlowRate(1, "MSm3/day")
        TPflash(fluid1)
        method = "friction theory"
        fluid1.getPhase('gas').getPhysicalProperties().setViscosityModel(method)
        fluid1.initProperties()
        mu = fluid1.getViscosity('cP')
        rho = fluid1.getDensity('kg/m3')
        return mu, rho

def general_gas_equation(q,p1,p2,D,G,z,L,t,mu,type):
    if type == 'estimate quantity':
        L = L*0.001
        f = 0.02
        e = 0.00005 #Roughness
        tb = 273+15.55556
        Pb = 101.325
        D = D *1000
        calc_port =G*(t+273.15)*L*z
        dp_2 = ((p1*98.066)**2)-((p2*98.066)**2)
        D_mm = D*1000
        error =10
        while error > 0.0001:
            Q_std = ((11.4946*(10**-4))*(1/(f**0.5))*(tb/Pb)*((dp_2/calc_port)**0.5)*((D)**2.5))/24 #*0.947942947917463/24
            Re = 0.5134*(Pb/tb)*((G*Q_std*24)/(mu*0.01*D_mm))
            f1 = (1/(-2*np.log((e/(3.7*D))+(2.51/(Re*(f**0.5))))))**2
            
            error = abs(f - f1)
            f = f1
        result1,result2 = Q_std,f
    if type == "estimate upstream pressure":   
        Q_std = q
        L = L*0.001
        f = 0.02
        e = 0.00005 #Roughness
        tb = 273+15.55556
        Pb = 101.325
        D = D *1000
        calc_port =G*(t+273.15)*L*z
        #dp_2 = ((p1*98.066)**2)-((p2*98.066)**2)
        D_mm = D*1000
        error =10
        Re = 0.5134*(Pb/tb)*((G*Q_std*24)/(mu*0.01*D_mm))
        while error > 0.0001:
            f1 = (1/(-2*np.log((e/(3.7*D))+(2.51/(Re*(f**0.5))))))**2
            error = abs(f - f1)
            f = f1
        p1 = np.sqrt(((((Q_std*24*(f**0.5))/((11.4946*(10**-4))*(tb/Pb)*((D)**2.5)))**2)*calc_port)+((p2*98.066)**2))/98.066
        result1,result2 = p1,f
    if type == "estimate downstream pressure":   
        Q_std = q
        L = L*0.001
        f = 0.02
        e = 0.00005 #Roughness
        tb = 273+15.55556
        Pb = 101.325
        D = D *1000
        calc_port =G*(t+273.15)*L*z
        #dp_2 = ((p1*98.066)**2)-((p2*98.066)**2)
        D_mm = D*1000
        error =10
        Re = 0.5134*(Pb/tb)*((G*Q_std*24)/(mu*0.01*D_mm))
        while error > 0.0001:
            f1 = (1/(-2*np.log((e/(3.7*D))+(2.51/(Re*(f**0.5))))))**2
            error = abs(f - f1)
            f = f1
        p2 = np.sqrt(-((((Q_std*24*(f**0.5))/((11.4946*(10**-4))*(tb/Pb)*((D)**2.5)))**2)*calc_port)+((p1*98.066)**2))/98.066
        result1,result2 = p2,f    
    return result1,result2

def gas_equations(q,p1,p2,D,G,L,t,z,mu,type):
    
    t = t +273.15
    if type=="estimate quantity":
        Q = [0,0,0,0,0,0,0,0,0]
        p1=p1*98.066*1000
        p2=p2*98.066*1000
        t = t+273.15 
        Q[0] = fluids.compressible.Panhandle_A(SG=G,Tavg=t,L=L,D=D,P1=p1,P2=p2,Q=None,Ts=288.7,Ps=101325.0,Zavg=z,E=0.95)

        Q[1]=fluids.compressible.Panhandle_B(SG=G,Tavg=t,L=L,D=D,P1=p1,P2=p2,Q=None,Ts=288.7,Ps=101325.0,Zavg=z,E=0.95)

        Q[2]=fluids.compressible.Weymouth(SG=G,Tavg=t,L=L,D=D,P1=p1,P2=p2,Q=None,Ts=288.7,Ps=101325.0,Zavg=z,E=0.95)


        result =np.array(Q)*3600
    if type=="estimate upstream pressure":
        q = q/3600
        P = [0,0,0,0,0,0,0,0,0]
        
        p2=p2*98.066*1000
        
        P[0] = fluids.compressible.Panhandle_A(SG=G,Tavg=t,L=L,D=D,P1=None,P2=p2,Q=q,Ts=288.7,Ps=101325.0,Zavg=z,E=0.95)

        P[1]=fluids.compressible.Panhandle_B(SG=G,Tavg=t,L=L,D=D,P1=None,P2=p2,Q=q,Ts=288.7,Ps=101325.0,Zavg=z,E=0.95)

        P[2]=fluids.compressible.Weymouth(SG=G,Tavg=t,L=L,D=D,P1=None,P2=p2,Q=q,Ts=288.7,Ps=101325.0,Zavg=z,E=0.95)

  
        result=np.array(P)/98066
    if type=="estimate downstream pressure":
        q = q/3600
        P = [0,0,0,0,0,0,0,0,0]
        p1=p1*98.066*1000
        
        P[0] = fluids.compressible.Panhandle_A(SG=G,Tavg=t,L=L,D=D,P1=p1,P2=None,Q=q,Ts=288.7,Ps=101325.0,Zavg=z,E=0.95)

        P[1]=fluids.compressible.Panhandle_B(SG=G,Tavg=t,L=L,D=D,P1=p1,P2=None,Q=q,Ts=288.7,Ps=101325.0,Zavg=z,E=0.95)

        P[2]=fluids.compressible.Weymouth(SG=G,Tavg=t,L=L,D=D,P1=p1,P2=None,Q=q,Ts=288.7,Ps=101325.0,Zavg=z,E=0.95)

        
        result=np.array(P)/98066
    
    return result

def main():
    html_temp="""
    <div style="background-color:lightblue;padding:16px">
    <h2 style="color:black"; text-align:center> Line sizing tool </h2>
    </div>
    
        """
    st.markdown(html_temp, unsafe_allow_html=True)
    s1 = st.selectbox('Chooose your line sizing?',('Gas - estimate Quantity (Std.m3/hr)','Gas - estimate Upstream pressure (kg/cm2.a)','Gas - estimate Downstream pressure (kg/cm2.a)','Liquid pressure drop/NPSHa','Estimate Equivalent Length','Check Standards for Line Sizing'), key = 'type')
    if s1 == 'Gas - estimate Quantity (Std.m3/hr)':
        st.write('## Estimation of Equivalent Length') 
        st.write("""When the piping layout is not available, the equivalent length (Le) of the piping will be estimated based on the straight length (Ls) as follows:\n 1. Process area: 3.0 times Ls\n 2. Common area: 1.5 times Ls\n 3. Offsite area: 1.3 times Ls""")
        st.write("""Note that Ls is the sum of XYZ coordinate length. \n For large size or high pressure piping, it is recommended to estimate the number of elbows tees and valves and evaluate the equivalent length, assuming piping layout.""")
        
        df_gas['input'] = 0.00
        edited_df = st.experimental_data_editor(df_gas)
        q=0
        p1,p2,t,L,D = edited_df.iloc[0,0],edited_df.iloc[1,0],edited_df.iloc[2,0],edited_df.iloc[3,0],edited_df.iloc[4,0]
        if D != 0 and p1 != 0 and p2 !=0 :
            try:
                D = fluids.nearest_pipe(NPS=find_nearest(D))[1]
            except ValueError: pass
            Q_std = [0,0,0,0,0,0,0,0,0,0,0]
            s_select = st.selectbox('Estimate M.wt, Cp/Cv and Z, Density and Viscosity?',('Yes','I already have these values'), key = 'k_calculations')
            if s_select == 'I already have these values':
                    m_wt= st.number_input('Molecular weight' , key = 'mwt')
                    z= st.number_input('Compressibility factor', key = 'z')
                    k= st.number_input('Cp/Cv', key = 'k')
                    rho2= st.number_input('Density (kg/m3)', key = 'rho')
                    mu= st.number_input('Viscosity (Cp)', key = 'vis')
                    G = m_wt/29
            else:
                    try:
                            df_comp, done = choose_composition()
                            if done == True:
                                z1, m_wt = Z_calculations(df_comp,t,p1)
                                z2, m_wt = Z_calculations(df_comp,t,p2)
                                G = m_wt/29
                                z = (z1+z2)*0.5
                                mu,rho1 = get_viscosity(df_comp,p1,t)
                                mu,rho2 = get_viscosity(df_comp,p2,t)
                                k = k_calculations(df_comp,df_comp_table,t,t)
                            
                    except (ValueError,TypeError, KeyError, ZeroDivisionError):st.write('your total mol. percent should add up to 100')
                    except UnboundLocalError: pass

        if st.button("Reveal Calculations", key = 'calculations_table22'):
            try:
                if p2 > p1:
                     st.warning('Downstream presure is higher than Upstream pressure!')
                else:
                    dp_100m = ((p1 - p2)*100)/L
                    Q_std[:3] = gas_equations(q,p1,p2,D,G,L,t,z,mu,'estimate quantity')
                    Q_std[10],f = general_gas_equation(Q_std[0],p1,p2,D,G,z,L,t,mu,'estimate quantity')
                    df_result = pd.DataFrame(df_summary)
                    if s_select != 'I already have these values':
                        Q_std[9] = NeqSim_calculations(Q_std[10],D,df_comp,t,p1,p2,L,'estimate quantity')
                        df_result['NeqSim Simulator'] = Summary_calculations(Q_std[9],D,G,mu,np.nan,p1,p2,t,m_wt,k,rho2,L,z,dp_100m)
                    else: del df_result['NeqSim Simulator'] 
                    
                    df_result['General Gas'] = Summary_calculations(Q_std[10],D,G,mu,f,p1,p2,t,m_wt,k,rho2,L,z,dp_100m)
                    
                    
                    df_result['Panhandle_A'] = Summary_calculations(Q_std[0],D,G,mu,0.95,p1,p2,t,m_wt,k,rho2,L,z,dp_100m)
                    df_result['Panhandle_B'] = Summary_calculations(Q_std[1],D,G,mu,0.95,p1,p2,t,m_wt,k,rho2,L,z,dp_100m)
                    df_result['Weymouth'] = Summary_calculations(Q_std[2],D,G,mu,0.95,p1,p2,t,m_wt,k,rho2,L,z,dp_100m)
                    

                    st.dataframe(df_result)
                    st.download_button("Click to download your calculations table!", convert_data(df_result.reset_index()),"pressure_drop_calculations.csv","text/csv", key = "download1")
                    if s_select != 'I already have these values':
                        graph_NeqSim(Q_std[9],D,df_comp,t,p1,L)
            except UnboundLocalError: st.write('Check your data input!')
            except jpype.JException: 
                dp_100m = ((p1 - p2)*100)/L
                Q_std[:3] = gas_equations(q,p1,p2,D,G,L,t,z,mu,'estimate quantity')
                Q_std[10],f = general_gas_equation(Q_std[0],p1,p2,D,G,z,L,t,mu,'estimate quantity')
                
                
                df_result = pd.DataFrame(df_summary)
                df_result['General Gas'] = Summary_calculations(Q_std[10],D,G,mu,f,p1,p2,t,m_wt,k,rho2,L,z,dp_100m)
                del df_result['NeqSim Simulator'] 
                
                df_result['Panhandle_A'] = Summary_calculations(Q_std[0],D,G,mu,0.95,p1,p2,t,m_wt,k,rho2,L,z,dp_100m)
                df_result['Panhandle_B'] = Summary_calculations(Q_std[1],D,G,mu,0.95,p1,p2,t,m_wt,k,rho2,L,z,dp_100m)
                df_result['Weymouth'] = Summary_calculations(Q_std[2],D,G,mu,0.95,p1,p2,t,m_wt,k,rho2,L,z,dp_100m)
                

                st.dataframe(df_result)
                st.download_button("Click to download your calculations table!", convert_data(df_result.reset_index()),"pressure_drop_calculations.csv","text/csv", key = "download1")
                st.warning('Neqsim simulator couldnt converge')
    elif s1 == "Gas - estimate Upstream pressure (kg/cm2.a)":
        st.write('## Estimation of Equivalent Length') 
        st.write("""When the piping layout is not available, the equivalent length (Le) of the piping will be estimated based on the straight length (Ls) as follows:\n 1. Process area: 3.0 times Ls\n 2. Common area: 1.5 times Ls\n 3. Offsite area: 1.3 times Ls""")
        st.write("""Note that Ls is the sum of XYZ coordinate length. \n For large size or high pressure piping, it is recommended to estimate the number of elbows tees and valves and evaluate the equivalent length, assuming piping layout.""")
        
        df_gas_modified = df_gas.copy()
        df_gas_modified['input'] = 0.00
        df_gas_modified.rename(index={'P1 (Kg/cm2)': 'Flow Rate (Std.m3/hr)'}, inplace=True)
        edited_df = st.experimental_data_editor(df_gas_modified)
        
        
        
        q,p2,t,L,D = edited_df.iloc[0,0],edited_df.iloc[1,0],edited_df.iloc[2,0],edited_df.iloc[3,0],edited_df.iloc[4,0]
        if D != 0 and p2 != 0:
            try:
                D = fluids.nearest_pipe(NPS=find_nearest(D))[1]
            except ValueError: pass
            P1 = [0,0,0,0,0,0,0,0,0,0,0]
            s_select = st.selectbox('Estimate M.wt, Cp/Cv and Z, Density and Viscosity?',('Yes','I already have these values'), key = 'k_calculations')
            if s_select == 'I already have these values':
                    m_wt= st.number_input('Molecular weight' , key = 'mwt')
                    z= st.number_input('Compressibility factor', key = 'z')
                    k= st.number_input('Cp/Cv', key = 'k')
                    rho2= st.number_input('Density (kg/m3)', key = 'rho')
                    mu= st.number_input('Viscosity (Cp)', key = 'vis')
                    G = m_wt/29
            else:
                    try:
                        done = False
                        df_comp, done = choose_composition()
                        if done == True:
                            z2, m_wt = Z_calculations(df_comp,t,p2)
                            z = z2
                            G = m_wt/29
                            mu,rho2 = get_viscosity(df_comp,p2,t)
                            k = k_calculations(df_comp,df_comp_table,t,t)
                        
                        
                        
                    except (ValueError,TypeError, KeyError, ZeroDivisionError):st.write('your total mol. percent should add up to 100')
                    except UnboundLocalError: pass

        if st.button("Reveal Calculations", key = 'calculations_table_P1'):
            try:
                P1[:3] = gas_equations(q,None,p2,D,G,L,t,z,mu,"estimate upstream pressure")
                P1[4],f = general_gas_equation(q,None,p2,D,G,z,L,t,mu,"estimate upstream pressure")
                df_result = pd.DataFrame(df_summary)
                if s_select != 'I already have these values':
                    P1[3] = NeqSim_calculations(q,D,df_comp,t,P1[0],p2,L,"estimate upstream pressure")
                    dp_100m = ((np.array(P1) -p2)*100)/L
                    df_result['NeqSim Simulator'] = Summary_calculations(q,D,G,mu,0,P1[3],p2,t,m_wt,k,rho2,L,z,dp_100m[3])
                else: 
                     dp_100m = ((np.array(P1) -p2)*100)/L
                     del df_result['NeqSim Simulator']
                
                
                df_result['General Gas'] = Summary_calculations(q,D,G,mu,f,P1[4],p2,t,m_wt,k,rho2,L,z,dp_100m[4])
                
                
                df_result['Panhandle_A'] = Summary_calculations(q,D,G,mu,0.95,P1[0],p2,t,m_wt,k,rho2,L,z,dp_100m[0])
                df_result['Panhandle_B'] = Summary_calculations(q,D,G,mu,0.95,P1[1],p2,t,m_wt,k,rho2,L,z,dp_100m[1])
                df_result['Weymouth'] = Summary_calculations(q,D,G,mu,0.95,P1[2],p2,t,m_wt,k,rho2,L,z,dp_100m[2])
                st.dataframe(df_result)
                st.download_button("Click to download your calculations table!", convert_data(df_result.reset_index()),"pressure_drop_calculations.csv","text/csv", key = "download2")
                if s_select != 'I already have these values':
                        graph_NeqSim(q,D,df_comp,t,P1[3],L)
            except UnboundLocalError: st.write('Check your data input!')
            except jpype.JException: 
                P1[:3] = gas_equations(q,None,p2,D,G,L,t,z,mu,"estimate upstream pressure")
                P1[4],f = general_gas_equation(q,None,p2,D,G,z,L,t,mu,"estimate upstream pressure")
                
                
                dp_100m = ((np.array(P1) -p2)*100)/L
                df_result = pd.DataFrame(df_summary)
                df_result['General Gas'] = Summary_calculations(q,D,G,mu,f,P1[4],p2,t,m_wt,k,rho2,L,z,dp_100m[4])
                del df_result['NeqSim Simulator'] 
                
                df_result['Panhandle_A'] = Summary_calculations(q,D,G,mu,0.95,P1[0],p2,t,m_wt,k,rho2,L,z,dp_100m[0])
                df_result['Panhandle_B'] = Summary_calculations(q,D,G,mu,0.95,P1[1],p2,t,m_wt,k,rho2,L,z,dp_100m[1])
                df_result['Weymouth'] = Summary_calculations(q,D,G,mu,0.95,P1[2],p2,t,m_wt,k,rho2,L,z,dp_100m[2])
                st.dataframe(df_result)
                st.download_button("Click to download your calculations table!", convert_data(df_result.reset_index()),"pressure_drop_calculations.csv","text/csv", key = "download2")
                st.warning('Neqsim simulator couldnt converge')
            
    elif s1 == "Gas - estimate Downstream pressure (kg/cm2.a)":
        st.write('## Estimation of Equivalent Length') 
        st.write("""When the piping layout is not available, the equivalent length (Le) of the piping will be estimated based on the straight length (Ls) as follows:\n 1. Process area: 3.0 times Ls\n 2. Common area: 1.5 times Ls\n 3. Offsite area: 1.3 times Ls""")
        st.write("""Note that Ls is the sum of XYZ coordinate length. \n For large size or high pressure piping, it is recommended to estimate the number of elbows tees and valves and evaluate the equivalent length, assuming piping layout.""")
        
        df_gas_modified = df_gas.copy()
        df_gas_modified['input'] = 0.00
        df_gas_modified.rename(index={'P1 (Kg/cm2)': 'Flow Rate (Std.m3/hr)','P2(Kg/cm2)':'P1 (Kg/cm2)'}, inplace=True)
        edited_df = st.experimental_data_editor(df_gas_modified)
        
        
        
        q,p1,t,L,D = edited_df.iloc[0,0],edited_df.iloc[1,0],edited_df.iloc[2,0],edited_df.iloc[3,0],edited_df.iloc[4,0]
        if D != 0 and p1 != 0:
            try:
                D = fluids.nearest_pipe(NPS=find_nearest(D))[1]
            except ValueError: pass
            P2 = [0,0,0,0,0,0,0,0,0,0,0]
            s_select = st.selectbox('Estimate M.wt, Cp/Cv and Z, Density and Viscosity?',('Yes','I already have these values'), key = 'k_calculations')
            if s_select == 'I already have these values':
                    m_wt= st.number_input('Molecular weight' , key = 'mwt')
                    z= st.number_input('Compressibility factor', key = 'z')
                    k= st.number_input('Cp/Cv', key = 'k')
                    rho1= st.number_input('Density (kg/m3)', key = 'rho')
                    mu= st.number_input('Viscosity (Cp)', key = 'vis')
                    G = m_wt/29
            else:
                    try:
                            df_comp, done = choose_composition()
                            if done == True:
                            
                                z1, m_wt = Z_calculations(df_comp,t,p1)
                                
                                z = z1
                                G = m_wt/29
                                mu,rho1 = get_viscosity(df_comp,p1,t)
                                k = k_calculations(df_comp,df_comp_table,t,t)
                            
                            
                            
                    except (ValueError,TypeError, KeyError, ZeroDivisionError):st.write('your total mol. percent should add up to 100')
                    except UnboundLocalError: pass

        if st.button("Reveal Calculations", key = 'calculations_table_P2'):
            try:
                P2[:3] = gas_equations(q,p1,None,D,G,L,t,z,mu,"estimate downstream pressure")
                P2[4],f = general_gas_equation(q,p1,None,D,G,z,L,t,mu,"estimate downstream pressure")
                df_result = pd.DataFrame(df_summary)
                if s_select != 'I already have these values':
                    P2[3] = NeqSim_calculations(q,D,df_comp,t,p1,P2[0],L,"estimate downstream pressure")
                    dp_100m = ((p1 - np.array(P2) )*100)/L
                    df_result['NeqSim Simulator'] = Summary_calculations(q,D,G,mu,0,p1,P2[3],t,m_wt,k,rho1,L,z,dp_100m[3])
                else: 
                     del df_result['NeqSim Simulator']
                     dp_100m = ((p1 - np.array(P2) )*100)/L
                df_result['General Gas'] = Summary_calculations(q,D,G,mu,f,p1,P2[4],t,m_wt,k,rho1,L,z,dp_100m[4])
                
                
                df_result['Panhandle_A'] = Summary_calculations(q,D,G,mu,0.95,p1,P2[0],t,m_wt,k,rho1,L,z,dp_100m[0])
                df_result['Panhandle_B'] = Summary_calculations(q,D,G,mu,0.95,p1,P2[1],t,m_wt,k,rho1,L,z,dp_100m[1])
                df_result['Weymouth'] = Summary_calculations(q,D,G,mu,0.95,p1,P2[2],t,m_wt,k,rho1,L,z,dp_100m[2])
                st.dataframe(df_result)
                st.download_button("Click to download your calculations table!", convert_data(df_result.reset_index()),"pressure_drop_calculations.csv","text/csv", key = "download3")
                if s_select != 'I already have these values':
                    graph_NeqSim(q,D,df_comp,t,p1,L)
            except UnboundLocalError: st.write('Check your data input!')
            except ValueError: st.warning('Check your assumptions: Empirical equations couldnt find a valid solution')
            except jpype.JException: 
                P2[:3] = gas_equations(q,p1,None,D,G,L,t,z,mu,"estimate downstream pressure")
                P2[4],f = general_gas_equation(q,p1,None,D,G,z,L,t,mu,"estimate downstream pressure")
                
                
                dp_100m = ((p1 - np.array(P2) )*100)/L
                df_result = pd.DataFrame(df_summary)
                df_result['General Gas'] = Summary_calculations(q,D,G,mu,f,p1,P2[4],t,m_wt,k,rho1,L,z,dp_100m[4])
                del df_result['NeqSim Simulator'] 
                
                df_result['Panhandle_A'] = Summary_calculations(q,D,G,mu,0.95,p1,P2[0],t,m_wt,k,rho1,L,z,dp_100m[0])
                df_result['Panhandle_B'] = Summary_calculations(q,D,G,mu,0.95,p1,P2[1],t,m_wt,k,rho1,L,z,dp_100m[1])
                df_result['Weymouth'] = Summary_calculations(q,D,G,mu,0.95,p1,P2[2],t,m_wt,k,rho1,L,z,dp_100m[2])
                st.dataframe(df_result)
                st.download_button("Click to download your calculations table!", convert_data(df_result.reset_index()),"pressure_drop_calculations.csv","text/csv", key = "download3")
                st.warning('Neqsim simulator couldnt converge')
    elif s1 == 'Liquid pressure drop/NPSHa':
        st.write('## Estimation of Equivalent Length') 
        st.write("""When the piping layout is not available, the equivalent length (Le) of the piping will be estimated based on the straight length (Ls) as follows:\n 1. Process area: 3.0 times Ls\n 2. Common area: 1.5 times Ls\n 3. Offsite area: 1.3 times Ls""")
        st.write("""Note that Ls is the sum of XYZ coordinate length. \n For large size or high pressure piping, it is recommended to estimate the number of elbows tees and valves and evaluate the equivalent length, assuming piping layout.""")
        
        def Darcy_equation_liq(Q,L,D,rho_liq,mu):
            
                Q = Q /3600
                D = fluids.nearest_pipe(NPS=find_nearest(D))[1]
                mu = mu*0.001
                A = np.pi * (D**2) * 0.25
                
                v_liq = Q/A 
                f=0.02
                e = 0.00005
                error = 10
                Re = (rho_liq*v_liq*D)/mu
                if Re < 2040:
                            f1 = 64/Re 
                            error = 0
                else:
                    while error > 0.0001:
                        f1 = (1/(-2*np.log10((e/(3.7*D))+(2.51/(Re*(f**0.5))))))**2 
                        error = abs(f - f1)
                        f = f1
                dp = (f*L*rho_liq*(v_liq**2))/(D*2)*0.001*0.0101972

                return dp,v_liq,Re,f,e
        def Nelson_equation(Q,L,D,rho_liq,mu):
            Q = Q /3600
            mu = mu
            D1 = fluids.nearest_pipe(NPS=find_nearest(D))[1]
            A = np.pi * (D1**2) * 0.25
            
            v_liq = (Q/A)
            
            Re = (rho_liq*v_liq*D1)/(mu*0.001)
            D2 = fluids.nearest_pipe(NPS=find_nearest(D))[0]
            
            L = L*3.28084
            v_liq = (Q/A )*3.28084
            S_gr = rho_liq/(16.01846336974*62.4)
            if Re <2040:
                 f = 15.017*(Re**(-0.995))
            else : 
                 f=0.048*(Re**(-0.183))

            e = np.nan
            
            dp = (0.323*f*(v_liq**2)*L*S_gr)/D2

            v_liq = v_liq/3.28084
            dp =dp*0.070307
            return dp,v_liq,Re,f,e
        df_liq['input'] = 0.00
        s2 = st.selectbox('Calculate NPSHa?',('No','Yes'), key = 'NPSHa')
        if s2 == 'No':
            edited_df = st.experimental_data_editor(df_liq.iloc[:7,:])
            p1,t,Q,rho_liq, mu,L,D = edited_df.iloc[0,0],edited_df.iloc[1,0],edited_df.iloc[2,0],edited_df.iloc[3,0],edited_df.iloc[4,0],edited_df.iloc[5,0],edited_df.iloc[6,0]
            
            
            
            if st.button("Reveal Calculations", key = 'calculations_tableLiq'):
                Nelson_equation(Q,L,D,rho_liq,mu)
                dp,v,Re,f,e = [0,0],[0,0],[0,0],[0,0],[0,0]
                p2 = [0,0]
                dp[0],v[0],Re[0],f[0],e[0]=Darcy_equation_liq(Q,L,D,rho_liq,mu)
                p2[0]= p1-dp[0]
                dp[1],v[1],Re[1],f[1],e[1]=Nelson_equation(Q,L,D,rho_liq,mu)
                p2[1]= p1-dp[1]
                df_liq.iloc[[0,1,2,3,4,5,6,9,10,11,12,13,14],0] = [p1,t,Q,rho_liq, mu,L,D,p2[0],dp[0],v[0],Re[0],f[0],e[0]]
                df_liq['Nelson (fannings Equation)'] = [p1,t,Q,rho_liq, mu,L,D,np.nan,np.nan,p2[1],dp[1],v[1],Re[1],f[1],e[1],np.nan]
                df_liq.rename(columns={'input': 'Darcy Equation'}, inplace=True)
                st.dataframe(df_liq.iloc[[0,1,2,3,4,5,6,9,10,11,12,13,14],:])
                st.download_button("Click to download your calculations table!", convert_data(df_liq.iloc[[0,1,2,3,4,5,6,9,10,11,12,13,14],:].reset_index()),"pressure_drop_calculations.csv","text/csv", key = "download4")
        else: 
            edited_df = st.experimental_data_editor(df_liq.iloc[:9,:])
            p1,t,Q,rho_liq, mu,L,D,H,Vp = edited_df.iloc[0,0],edited_df.iloc[1,0],edited_df.iloc[2,0],edited_df.iloc[3,0],edited_df.iloc[4,0],edited_df.iloc[5,0],edited_df.iloc[6,0],edited_df.iloc[7,0],edited_df.iloc[8,0]
            
            if st.button("Reveal Calculations", key = 'calculations_tableLiq'):
                dp,v,Re,f,e,NPSHa = [0,0],[0,0],[0,0],[0,0],[0,0],[0,0]
                p2 = [0,0]
                dp[0],v[0],Re[0],f[0],e[0]=Darcy_equation_liq(Q,L,D,rho_liq,mu)
                p2[0]= p1-dp[0]
                NPSHa[0] = (p2[0] - Vp + 1.03323)*10/(rho_liq*0.001)+H
                dp[1],v[1],Re[1],f[1],e[1]=Nelson_equation(Q,L,D,rho_liq,mu)
                p2[1]= p1-dp[1]
                NPSHa[1] = (p2[1] - Vp + 1.03323)*10/(rho_liq*0.001)+H
                df_liq.rename(columns={'input': 'Darcy Equation'}, inplace=True)
                df_liq['Darcy Equation'] = [p1,t,Q,rho_liq, mu,L,D,H,Vp,p2[0],dp[0],v[0],Re[0],f[0],e[0],NPSHa[0]]
                df_liq['Nelson (fannings Equation)'] = [p1,t,Q,rho_liq, mu,L,D,H,Vp,p2[1],dp[1],v[1],Re[1],f[1],e[1],NPSHa[1]]
                st.dataframe(df_liq)
                st.download_button("Click to download your calculations table!", convert_data(df_liq.reset_index()),"pressure_drop_calculations.csv","text/csv", key = "download4")
    elif s1 == 'Estimate Equivalent Length':
            url4 = 'http://raw.githubusercontent.com/Ahmedhassan676/pressure_drop/main/equ_length.csv'
            df_Le = pd.read_csv(url4, index_col=[0])
            Length = st.number_input('Insert pipe straight length (m)', key = 'pipe')
            if "df" not in st.session_state:
                st.session_state.df = pd.DataFrame(columns=[0,1,2,3])
            
            rw = -1
            fittings_list=[]
            with st.form("my_form"):
            
               
                
                
                selected_columns = st.selectbox('Select fitting', options=df_Le.columns)
                selected_indices = st.selectbox('Select diameter/d', options=df_Le.index)
                count = st.slider('Count of fittings', 1, 10, 1)
                fittings_list.append([selected_indices,selected_columns,df_Le.loc[selected_indices,selected_columns],count])
                
                
                
                # create an "Add" button to add selected indices to the list
                add_button = st.form_submit_button('Add')
                if add_button:
                    
                    rw = st.session_state.df.shape[0] 
                    st.session_state.df.loc[rw] = fittings_list[0]
            df_fitting = pd.DataFrame(st.session_state.df).rename(columns={0:'D', 1:'fitting', 2:'Le (ft)', 3: 'count'})
            df_fitting['Le (m)'] = df_fitting['Le (ft)']*0.3048
            st.dataframe(df_fitting.loc[:,['D','fitting','Le (m)','count']])
            
            st.write('Sum of All fittings Length = {} m and Equivalent Length is {} m'.format(round(sum(df_fitting['Le (m)']*df_fitting['count']),2),round(sum(df_fitting['Le (m)']*df_fitting['count']),2)+Length))
            
           
    elif s1=='Check Standards for Line Sizing':
       
        from pandas.api.types import (is_categorical_dtype,is_datetime64_any_dtype,is_numeric_dtype,is_object_dtype,)



        def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
            """
            Adds a UI on top of a dataframe to let viewers filter columns

            Args:
                df (pd.DataFrame): Original dataframe

            Returns:
                pd.DataFrame: Filtered dataframe
            """
            modify = st.checkbox("Add filters")

            if not modify:
                return df

            df = df.copy()

            # Try to convert datetimes into a standard format (datetime, no timezone)
            for col in df.columns:
                if is_object_dtype(df[col]):
                    try:
                        df[col] = pd.to_datetime(df[col])
                    except Exception:
                        pass

                if is_datetime64_any_dtype(df[col]):
                    df[col] = df[col].dt.tz_localize(None)

            modification_container = st.container()

            with modification_container:
                to_filter_columns = st.multiselect("Filter dataframe on", df.columns)
                for column in to_filter_columns:
                    left, right = st.columns((1, 20))
                    # Treat columns with < 10 unique values as categorical
                    if is_categorical_dtype(df[column]) or df[column].nunique() < 10:
                        user_cat_input = right.multiselect(
                            f"Values for {column}",
                            df[column].unique(),
                            default=list(df[column].unique()),
                        )
                        df = df[df[column].isin(user_cat_input)]
                    elif is_numeric_dtype(df[column]):
                        _min = float(df[column].min())
                        _max = float(df[column].max())
                        step = (_max - _min) / 100
                        user_num_input = right.slider(
                            f"Values for {column}",
                            min_value=_min,
                            max_value=_max,
                            value=(_min, _max),
                            step=step,
                        )
                        df = df[df[column].between(*user_num_input)]
                    else:
                        user_text_input = right.text_input(
                            f"Search in {column}",
                        )
                        if user_text_input:
                            df = df[df[column].astype(str).str.contains(user_text_input.lower())]

            return df
        df = pd.read_csv('http://raw.githubusercontent.com/Ahmedhassan676/pressure_drop/main/criteria.csv', index_col=[0])
        st.dataframe(filter_dataframe(df))
        st.write('Based On an Excel sheet Compiled by: Ajay S. Satpute')
if __name__ == '__main__':
    main()