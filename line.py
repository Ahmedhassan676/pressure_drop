import pandas as pd
import numpy as np
import streamlit as st
import pandas as pd
import neqsim
from neqsim.thermo.thermoTools import *
from neqsim.process.processTools import *
from neqsim import methods
from neqsim.thermo import fluid, TPflash, createfluid2
from neqsim.process import pipe, pipeline, clearProcess, stream, runProcess


url ='https://raw.githubusercontent.com/Ahmedhassan676/pressure_drop/main/table.csv'
df_gas = pd.read_csv(url, index_col=[0])
url_1 = 'http://raw.githubusercontent.com/Ahmedhassan676/pressure_drop/main/comp.csv'
df_comp_table = pd.read_csv(url_1)
url_2 ='http://raw.githubusercontent.com/Ahmedhassan676/pressure_drop/main/summary.csv'
df_summary = pd.read_csv(url_2, index_col=[0])
url_3 ='https://raw.githubusercontent.com/Ahmedhassan676/pressure_drop/main/table_liq.csv'
df_liq = pd.read_csv(url_3, index_col=[0])

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
def Summary_calculations(Q_std,D,G,mu,f_E,p1,p2,t,m_wt,k,rho2,L,z):
        tb = 273.15+15.55556
        Pb = 101.325
        D_mm = 25.4 * D
        A = np.pi*0.25*((D_mm/1000)**2)
        Q_normal = Q_std*(273.15/(273.15+15))
        Q_actual = Q_std*(1.033023/p2)*((t+273.15)/tb)
        v = Q_actual/(A*3600) #Velocity (m/s)
        sonic_velocity =((9.81*m_wt*847.9*(t+273.15))/k)**0.5
        Re = 0.5134*(Pb/tb)*((G*Q_std*24)/(mu*0.01*D_mm))
        mach = v/sonic_velocity
        dp_percent = ((p1-p2)/p1)*100
        rho_v_2 = rho2*(v**2)
        summary_list = [p1,p2,t,L,D,Q_std,Q_normal,Q_actual,dp_percent,f_E,Re,v,sonic_velocity,mach,rho_v_2,m_wt,z,k]
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
        
        return Z, m_wt,Pr,Tr
def choose_composition():
            
            df = pd.DataFrame({'Composition/property':df_comp_table.columns[1:],'mol%':np.zeros(len(df_comp_table.columns)-1), 'm.wt':df_comp_table.iloc[0,1:],'Pc':df_comp_table.iloc[1,1:],'Tc':df_comp_table.iloc[2,1:]})
            try:
                sum_of_comp = 0 
                c1,c2,c3,c15,c4,c5,c6,c7,c8,c9,c16,c10,c11,c13,c14,nh3,h2o = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
                options_list = [df_comp_table.columns[i] for i in [22,1,4,3,6,11,10,13,15,20,21,24,25,23,18,17,19]]
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
                    if df_comp_table.columns[3] in options:
                        c15 = st.number_input('ethylene%', key = 'c15')
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
                    if df_comp_table.columns[18] in options:
                        c14 = st.number_input('air%', key = 'c14')
                    if df_comp_table.columns[17] in options:
                        nh3 = st.number_input('Ammonia%', key = 'nh3')
                    if df_comp_table.columns[19] in options:
                        h2o = st.number_input('Water vapor%', key = 'h2o')
                    if c1 or c2 or c3 or c15 or c4 or c5 or c6 or c7 or c8 or c9 or c16 or c10 or c11 or c13 or c14 or nh3 or h2o:
                        c = []
                        for i in (c1,c2,c3,c15,c4,c5,c6,c7,c8,c9,c16,c10,c11,c13,c14,nh3,h2o):
                            c.append(i)
                        
                        for (i,j) in zip(options_list,c):
                            if j != None:
                                    df.loc[i,'mol%'] = j
                        
                        sum_of_comp = np.sum(df['mol%'])
                        
                st.success('Composition in Mol. percent completed!', icon="✅")
                
                return df[df['mol%'] != 0]

            except (ValueError, st.errors.DuplicateWidgetID): pass
            except (TypeError, KeyError, ZeroDivisionError):st.write('Please Check your data')
def NeqSim_calculations(q,D,df_comp,t,p1,p2,L):
        p1 = p1*98.066
        q = q*24/(1000000)
        D = D*25.4*0.001
        #Creating inlet fluid using SRK-EoS
            
        names =  list(df_comp.index.str.lower())
        molefractions = list(df_comp['mol%']*0.01)
        fluid1 = createfluid2(names, molefractions)
        fluid1.setMixingRule('classic')
        fluid1.setTemperature(t, "C")
        fluid1.setPressure(p1, "bara")
        fluid1.setTotalFlowRate(q, "MSm3/day")
        
        TPflash(fluid1)

        stream1 = stream(fluid1)

        deltaElevation = 0.0
        pipeLength = L*1000
        roughness= 5.0e-5
        diameter = D
        error = 10
        while error > 0.001:
            
            clearProcess()
            pipeSimple = pipe(stream1, pipeLength, deltaElevation, diameter, roughness)
            runProcess()
            
            error = p2 - pipeSimple.getOutStream().getFluid().getPressure('bara')
            if error >0:
                 q = q - q*0.01
            elif error < 0:
                 q = q + q*0.01
            else: error = error

        return q*1000000/24
def detailed_NeqSim(q,D,df_comp,t,p1,L):
    p1 = p1*0.98066
    q = q*24/(1000000)
    D = D*25.4*0.001
     #Creating inlet fluid using SRK-EoS
    names = list(df_comp.index.str.lower())
    molefractions = list(df_comp['mol%']*0.01)
    fluid1 = createfluid2(names, molefractions)
    
    #Creating stream and pipeline
    clearProcess()
    stream1 = stream(fluid1)
    stream1.setFlowRate(q, "MSm3/day")
    stream1.setTemperature(t, "C")
    stream1.setPressure(p1, "bara")

    diameter = [D,D] #meter
    roughnes = [5.0e-5,5.0e-5] #meter
    position = [0,L] #meter
    height = [0.0, 0.0] #meter
    outtemperatures =[288.15,288.15] #Kelvin
    outHeatU = [0, 0] #W/m2K
    wallHeatU = [0,0] #W/m2K

    pipe1 = pipeline(stream1, position, diameter, height, outtemperatures, roughnes,outHeatU,wallHeatU)
    pipe1.setNumberOfNodesInLeg(100)
    runProcess()
    numberofnodes = pipe1.getPipe().getTotalNumberOfNodes()
    pres = []
    temp = []
    length = []
    height = []
    calcdensity = []
    gasvelocity = []

    for node in range (0,pipe1.getPipe().getTotalNumberOfNodes()-1):
        pres.append(pipe1.getPipe().getNode(node).getBulkSystem().getPressure('bara'))
        temp.append(pipe1.getPipe().getNode(node).getBulkSystem().getTemperature('C'))
        height.append(pipe1.getPipe().getNode(node).getVerticalPositionOfNode())
        length.append(pipe1.getPipe().getNode(node).getDistanceToCenterOfNode())
        calcdensity.append(pipe1.getPipe().getNode(node).getBulkSystem().getDensity('kg/m3'))
        gasvelocity.append(pipe1.getPipe().getNode(node).getVelocity()) 
    fig,axs = plt.subplots(5,1)
    fig.set_figheight(12)
    fig.set_figwidth(12)
    axs[0].plot(length, pres, '-')
    axs[0].set_ylabel('Pressure [bara]')
    axs[0].set_xlabel('Length [meter]')

    
    axs[1].plot(length, temp)
    axs[1].set_xlabel('Length [meter]')
    axs[1].set_ylabel('Temperature[C]')

   
    axs[2].plot(length, calcdensity, '-')
    axs[2].set_ylabel('Density [kg/m3]')
    axs[2].set_xlabel('Length [meter]')

    
    axs[3].plot(length, gasvelocity, '-')
    axs[3].set_ylabel('gasvelocity [m/sec]')
    axs[3].set_xlabel('Length [meter]')

   
    axs[4].plot(length, height, '-')
    axs[4].set_ylabel('height [meter]')
    axs[4].set_xlabel('position [meter]')

    return st.pyplot(fig)
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

def general_gas_equation(p1,p2,D,G,z,L,t,mu):
    f = 0.02
    e = 0.00015*12 #Roughness
    tb = 273+15.55556
    Pb = 101.325
    calc_port =G*(t+273.15)*L*z
    dp_2 = ((p1*98.066)**2)-((p2*98.066)**2)
    D_mm = 25.4 * D
    error =10
    while error > 0.0001:
        Q_std = ((11.4946*(10**-4))*(1/(f**0.5))*(tb/Pb)*((dp_2/calc_port)**0.5)*((D*25.4)**2.5))/24 #*0.947942947917463/24
        Re = 0.5134*(Pb/tb)*((G*Q_std*24)/(mu*0.01*D_mm))
        f1 = (1/(-2*np.log10((e/(3.7*D))+(2.51/(Re*(f**0.5))))))**2 
        error = abs(f - f1)
        f = f1
    
    return Q_std,f
def panhandleA_equation(p1,p2,D,G,z,L,t):
    E = 0.95
    tb = 273+15.55556
    Pb = 101.325
    calc_port =(G**0.8539)*(t+273.15)*L*z
    dp_2 = ((p1*98.066)**2)-((p2*98.066)**2)
    Q_std = ((4.5965*10**-3)*(E)*((tb/Pb)**1.0788)*((dp_2/calc_port)**0.5394)*((D*25.4)**2.6182))/24 #*(273.15/(273.15+15))/24
    return Q_std,E
def panhandleB_equation(p1,p2,D,G,z,L,t):
    E = 0.95
    tb = 273+15.55556
    Pb = 101.325
    calc_port =(G**0.961)*(t+273.15)*L*z
    dp_2 = ((p1*98.066)**2)-((p2*98.066)**2)
    Q_std = ((1.002*(10**-2))*(E)*((tb/Pb)**1.02)*((dp_2/calc_port)**0.51)*((D*25.4)**2.53))/24 #*(273.15/(273.15+15))/24
    return Q_std,E
def weymouth_equation(p1,p2,D,G,z,L,t):
    E = 0.95
    tb = 273+15.55556
    Pb = 101.325
    calc_port =(G)*(t+273.15)*L*z
    dp_2 = ((p1*98.066)**2)-((p2*98.066)**2)
    Q_std = ((3.7435*(10**-3))*(E)*((tb/Pb))*((dp_2/calc_port)**0.5)*((D*25.4)**(8/3)))/24 #*(273.15/(273.15+15))/24
    return Q_std,E

def main():
    html_temp="""
    <div style="background-color:lightblue;padding:16px">
    <h2 style="color:black"; text-align:center> Line sizing tool </h2>
    </div>
    
        """
    st.markdown(html_temp, unsafe_allow_html=True)
    s1 = st.selectbox('Gas or liquid line sizing?',('Gas','Liquid'), key = 'type')
    if s1 == 'Gas':
        df_gas['input'] = 0.00
        edited_df = st.experimental_data_editor(df_gas)
        
        
        
        p1,p2,t,L,D = edited_df.iloc[0,0],edited_df.iloc[1,0],edited_df.iloc[2,0],edited_df.iloc[3,0]*0.001,edited_df.iloc[4,0]
        Q_std = [0,0,0,0,0]
        try:
                df_comp = choose_composition()
                
                z1, m_wt,Pr1,Tr1 = Z_calculations(df_comp,t,p1)
                z2, m_wt,Pr2,Tr2 = Z_calculations(df_comp,t,p2)
                G = m_wt/29
                z = (z1+z2)*0.5
                mu,rho1 = get_viscosity(df_comp,p1,t)
                mu,rho2 = get_viscosity(df_comp,p2,t)
                k = k_calculations(df_comp,df_comp_table,t,t)
                
        except (ValueError,TypeError, KeyError, ZeroDivisionError):st.write('your total mol. percent should add up to 100')
        except UnboundLocalError: pass

        if st.button("Reveal Calculations", key = 'calculations_table22'):
            
            Q_std[0],f = general_gas_equation(p1,p2,D,G,z,L,t,mu)
            Q_std[1],E_w = weymouth_equation(p1,p2,D,G,z,L,t)
            Q_std[2],E_a = panhandleA_equation(p1,p2,D,G,z,L,t)
            Q_std[3],E_b = panhandleB_equation(p1,p2,D,G,z,L,t)
            Q_std[4] = NeqSim_calculations(Q_std[0],D,df_comp,t,p1,p2,L)
            df_result = pd.DataFrame(df_summary)
            df_result['General Gas'] = Summary_calculations(Q_std[0],D,G,mu,f,p1,p2,t,m_wt,k,rho2,L,z)
            df_result['Weyouth'] = Summary_calculations(Q_std[1],D,G,mu,E_w,p1,p2,t,m_wt,k,rho2,L,z)
            df_result['Panhandle_A'] = Summary_calculations(Q_std[2],D,G,mu,E_a,p1,p2,t,m_wt,k,rho2,L,z)
            df_result['Panhandle_B'] = Summary_calculations(Q_std[3],D,G,mu,E_b,p1,p2,t,m_wt,k,rho2,L,z)
            df_result['NeqSim Simulator'] = Summary_calculations(Q_std[4],D,G,mu,0,p1,p2,t,m_wt,k,rho2,L,z)
            st.dataframe(df_result)
            #detailed_NeqSim(Q_std[4],D,df_comp,t,p1,L)
    else:
        def Darcy_equation(Q,L,D,rho_liq,mu):
                A = np.pi * (D**2) * 0.25
                v_liq = Q/A 
                f=0.02
                e = 0
                error = 10
                while error > 0.0001:
                    Re = (rho_liq*v_liq*D)/mu
                    f1 = (1/(-2*np.log10((e/(3.7*D))+(2.51/(Re*(f**0.5))))))**2 
                    error = abs(f - f1)
                    f = f1
                dp = (f*L*rho_liq*(v_liq**2))/(D*2)*0.001*0.0101972
        
                return dp,v_liq,Re,f,e
        df_liq['input'] = 0.00
        edited_df = st.experimental_data_editor(df_liq.iloc[:7,:])
        p1,t,Q,rho_liq, mu,L,D = edited_df.iloc[0,0],edited_df.iloc[1,0],edited_df.iloc[2,0],edited_df.iloc[3,0],edited_df.iloc[4,0],edited_df.iloc[5,0],edited_df.iloc[6,0]
        Q = Q /3600
        D = D*25.4*0.001
        mu = mu*0.001
        if st.button("Reveal Calculations", key = 'calculations_tableLiq'):
            dp,v,Re,f,e = Darcy_equation(Q,L,D,rho_liq,mu)
            p2 = p1-dp
            df_liq['input'] = [p1,t,Q,rho_liq, mu,L,D,p2,dp,v,Re,f,e]
            st.dataframe(df_liq)
if __name__ == '__main__':
    main()
