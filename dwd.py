'''DWD Data Sampler
collection of tools to easily get DWD station data
(cc) conrad.jackisch@tbt.tu-freiberg.de 2021
'''

#todo bugs: if no hourly values exists it raises an error nobody will understand
# let people choose from stations in case of ambiguity
# let people choose also about aggregation?

import numpy as np
import pandas as pd
import requests
import zipfile
import os
import numba


def dwd_stations():
    http_fi = 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/KL_Tageswerte_Beschreibung_Stationen.txt'
    r = requests.get(http_fi, allow_redirects=True)
    open('stat.txt', 'wb').write(r.content)

    dummy = pd.read_fwf('stat.txt',skiprows=2,header=None,encoding='latin1')
    dummy.columns = ['station_id','von_datum','bis_datum','Stationshoehe','geoBreite','geoLaenge','Stationsname','Bundesland']
    dummy.von_datum = pd.to_datetime(dummy.von_datum,format='%Y%m%d')
    dummy.bis_datum = pd.to_datetime(dummy.bis_datum,format='%Y%m%d')
    dummy.index = dummy.station_id
    return dummy.iloc[:,1:]


def dwd_stationno(st_str):
    dwdstat = dwd_stations()
    s_id = dwdstat.index[dwdstat.Stationsname.str.contains(st_str).values]
    if len(s_id)>1:
        if any(dwdstat.loc[s_id,'Stationsname']==st_str):
            s_id = (dwdstat.loc[s_id,'Stationsname']==st_str).index[(dwdstat.loc[s_id,'Stationsname']==st_str).values]
        else:
            dummyc = (pd.Timestamp.now() - dwdstat.loc[s_id,'bis_datum'])
            s_id = dummyc.iloc[(dummyc == dummyc.min()).values].index
            if len(s_id)>1:
                dummyc = dwdstat.loc[s_id,'bis_datum']-dwdstat.loc[s_id,'von_datum']
                s_id = dummyc.idxmax()
    return s_id.values[0]
            
def get_dwdfi(http_fo=None,station_no=3815):
    if http_fo==None:
        http_fo = 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/air_temperature/historical/'
    r = requests.get(http_fo, allow_redirects=True)
    open('dummy.html', 'wb').write(r.content)

    try:
        dummyf = pd.read_table('dummy.html',skiprows=9,sep='"',header=None).iloc[:-2,1]
        os.remove('dummy.html')
        return dummyf[dummyf.str.contains(str(station_no).zfill(5)).values].values[0]
    except:
        return 'no such file'


def dwd_recent(stationname='Emden',freq='daily',param='precip',ti='recent'):
    '''Read recent DWD data for stations'''
    
    sid = dwd_stationno(stationname)
    
    if freq == 'daily':
        http_fi = 'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/'
        if ti == 'recent':
            if param == 'precip':
                http_fi = http_fi + 'more_precip/recent/tageswerte_RR_'
            elif param == 'kl':
                http_fi = http_fi + 'kl/recent/tageswerte_KL_'
            else:
                print('check '+http_fi+' and implement request in this function')
                return 
            http_fi = http_fi+str(sid).zfill(5)+'_akt.zip'
        elif ti == 'historical':
            if param == 'precip':
                http_fi = http_fi + 'more_precip/historical'
            elif param == 'kl':
                http_fi = http_fi + 'kl/historical'
            elif param == 'solar':
                http_fi = http_fi + 'solar'
            else:
                print('check '+http_fi+' and implement request in this function')
                return             
            fidummy = get_dwdfi(http_fi,sid)
            if fidummy == 'no such file':
                return 'File not found!'
            else:
                http_fi = http_fi + '/' + fidummy
                
            
        
    elif freq == 'hourly':
        http_fi='https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/hourly/'
        if ti == 'recent':
            if param == 'precip':
                http_fi = http_fi + 'precipitation/recent/stundenwerte_RR_'
                http_fi = http_fi+str(sid).zfill(5)+'_akt.zip'
            elif param == 'temp':
                http_fi = http_fi + 'air_temperature/recent/stundenwerte_TU_'
                http_fi = http_fi+str(sid).zfill(5)+'_akt.zip'
            elif param == 'wind':
                http_fi = http_fi + 'wind/recent/stundenwerte_FF_'
                http_fi = http_fi+str(sid).zfill(5)+'_akt.zip'
            elif param == 'pressure':
                http_fi = http_fi + 'pressure/recent/stundenwerte_P0_'
                http_fi = http_fi+str(sid).zfill(5)+'_akt.zip'
            elif param == 'solar':
                http_fi = http_fi + 'solar/stundenwerte_ST_'
                http_fi = http_fi+str(sid).zfill(5)+'_row.zip'
            else:
                print('check '+http_fi+' and implement request in this function')
                return             
        elif ti == 'historical':
            if param == 'precip':
                http_fi = http_fi + 'precipitation/historical'
            elif param == 'temp':
                http_fi = http_fi + 'air_temperature/historical'
            elif param == 'wind':
                http_fi = http_fi + 'wind/historical'
            elif param == 'pressure':
                http_fi = http_fi + 'pressure/historical'
            elif param == 'solar':
                http_fi = http_fi + 'solar'
            elif param == 'kl':
                http_fi = http_fi + 'kl/historical/'
            else:
                print('check '+http_fi+' and implement request in this function')
                return             
            fidummy = get_dwdfi(http_fi,sid)
            if fidummy == 'no such file':
                return 'File not found!'
            else:
                http_fi = http_fi + '/' + fidummy            
        
    #get station file from server    
    r = requests.get(http_fi, allow_redirects=True)
    open('dummy.zip', 'wb').write(r.content)

    if os.stat('dummy.zip').st_size<500:
        print('Variable '+param+' not given at Station '+str(sid))
        dummy = pd.DataFrame([])
    else:
        myzipfile = zipfile.ZipFile('dummy.zip')
    
        dummy = pd.read_csv(myzipfile.open(myzipfile.namelist()[-1]),sep=';',na_values=-999)
        if freq == 'daily':
            dummy.index = pd.to_datetime(dummy.MESS_DATUM.values,format='%Y%m%d')
        elif freq == 'hourly':
            if param == 'solar':
                dummy.index = pd.to_datetime(dummy.MESS_DATUM.values,format='%Y%m%d%H:%M')
            else:
                dummy.index = pd.to_datetime(dummy.MESS_DATUM.values,format='%Y%m%d%H')
        
        os.remove('dummy.zip')
    
    return dummy
    #return http_fi

@numba.jit
def Te_opt(T_e,gammax, vabar):
    maxdeltaT_e = 1.
    maxit = 9999999
    itc = 0
    while maxdeltaT_e < 0.01:
        v_e = 0.6108 * np.exp(17.27 * T_e / (T_e + 237.3))  # saturated vapour pressure at T_e (S2.5)
        T_enew = gammax * (v_e - vabar)  # rearranged from S8.8
        deltaT_e = T_enew - T_e
        T_e = T_enew
        maxdeltaT_e = np.abs(np.max(deltaT_e))
        if itc > maxit:
            break
        itc += 1
    return T_e

def ET_SzilagyiJozsa(data, Elev=120., lat=53.5, windfunction_ver='1948', alpha=0.23, zerocorr=True):
    # Taken from R package Evapotranspiration >> Danlu Guo <danlu.guo@adelaide.edu.au>
    # Daily Actual Evapotranspiration after Szilagyi, J. 2007, doi:10.1029/2006GL028708
    # data is assumed to be a pandas data frame with at least:
    # T or Tmin/max - daily temperature in degree Celcius,
    # RH or RHmin/max - daily  relative humidity in percentage,
    # u2 - daily wind speed in meters per second
    # Rs - daily solar radiation in Megajoule per square meter.
    # Result in [mm/day] for EToSJ, EToPM, and EToPT reference ET

    # update alphaPT according to Szilagyi and Jozsa (2008)
    alphaPT = 1.31
    lat_rad = lat * np.pi / 180.
    alphaPT = 1.31
    sigma = 4.903e-09
    Gsc = 0.082
    lambdax = 2.45

    # add julian Days
    J = data.index.dayofyear

    # Calculating mean temperature
    if ('Ta' in data.columns):
        Ta = data['Ta']
    elif ('T' in data.columns):
        Ta = data['T']
    elif ('Temp' in data.columns):
        Ta = data['Temp']
    else:
        Ta = (data['Tmax'] + data[
            'Tmin']) / 2  # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998.

    # Saturated vapour pressure
    vs_Tmax = 0.6108 * np.exp(17.27 * data['Tmax'] / (data['Tmax'] + 237.3))  # Equation S2.5
    vs_Tmin = 0.6108 * np.exp(17.27 * data['Tmin'] / (data['Tmin'] + 237.3))  # Equation S2.5
    vas = (vs_Tmax + vs_Tmin) / 2.  # Equation S2.6

    # Vapour pressure
    if 'RHmax' in data.columns:
        vabar = (vs_Tmin * data['RHmax'] / 100. + vs_Tmax * data['RHmin'] / 100.) / 2.  # Equation S2.7
        # mean relative humidity
        RHmean = (data['RHmax'] + data['RHmin']) / 2

    else:
        vabar = (vs_Tmin + vs_Tmax) / 2 * data['RH'] / 100.

    if 'Rs' in data.columns:
        R_s = data['Rs']
    elif ('Rad' in data.columns):
        R_s = data['Rad']*0.01
    else:
        print('Radiation data missing')
        return

    # Calculations from data and constants for Penman

    P = 101.3 * ((293. - 0.0065 * Elev) / 293.) ** 5.26  # atmospheric pressure (S2.10)
    delta = 4098 * (0.6108 * np.exp((17.27 * Ta) / (Ta + 237.3))) / (
                (Ta + 237.3) ** 2)  # slope of vapour pressure curve (S2.4)
    gamma = 0.00163 * P / lambdax  # psychrometric constant (S2.9)
    d_r2 = 1 + 0.033 * np.cos(2 * np.pi / 365 * J)  # dr is the inverse relative distance Earth-Sun (S3.6)
    delta2 = 0.409 * np.sin(2 * np.pi / 365 * J - 1.39)  # solar dedication (S3.7)
    w_s = np.arccos(-1. * np.tan(lat_rad) * np.tan(delta2))  # sunset hour angle (S3.8)
    N = 24 / np.pi * w_s  # calculating daily values
    R_a = (1440 / np.pi) * d_r2 * Gsc * (
                w_s * np.sin(lat_rad) * np.sin(delta2) + np.cos(lat_rad) * np.cos(delta2) * np.sin(
            w_s))  # extraterristrial radiation (S3.5)
    R_so = (0.75 + (2 * 10 ** - 5) * Elev) * R_a  # clear sky radiation (S3.4)

    R_nl = sigma * (0.34 - 0.14 * np.sqrt(vabar)) * ((data['Tmax'] + 273.2) ** 4 + (data['Tmin'] + 273.2) ** 4) / 2 * (
                1.35 * R_s / R_so - 0.35)  # estimated net outgoing longwave radiation (S3.3)
    # For vegetated surface
    R_nsg = (1 - alpha) * R_s  # net incoming shortwave radiation (S3.2)
    R_ng = R_nsg - R_nl  # net radiation (S3.1)

    if 'u2' in data.columns:
        u2 = data['u2']
        if windfunction_ver == "1948":
            f_u = 2.626 + 1.381 * u2  # wind function Penman 1948 (S4.11)
        elif windfunction_ver == "1956":
            f_u = 1.313 + 1.381 * u2  # wind function Penman 1956 (S4.3)

        Ea = f_u * (vas - vabar)  # (S4.2)

        Epenman_Daily = delta / (delta + gamma) * (R_ng / lambdax) + gamma / (
                    delta + gamma) * Ea  # Penman open-water evaporation (S4.1)

    else:

        Epenman_Daily = 0.047 * R_s * np.sqrt(Ta + 9.5) - 2.4 * (R_s / R_a) ** 2 + 0.09 * (Ta + 20) * (
                    1 - RHmean / 100)  # Penman open-water evaporation without wind data by Valiantzas (2006) (S4.12)

    # Iteration for equilibrium temperature T_e
    T_e = Ta
    gammax = Ta - 1 / gamma * (1 - R_ng / (lambdax * Epenman_Daily))

    T_e = Te_opt(T_e.values, gammax.values, vabar.values)
    T_e = pd.Series(T_e, index=Ta.index)

    deltaT_e = 4098 * (0.6108 * np.exp((17.27 * T_e) / (T_e + 237.3))) / (
                (T_e + 237.3) ** 2)  # slope of vapour pressure curve (S2.4)
    E_PT_T_e = alphaPT * (deltaT_e / (deltaT_e + gamma) * R_ng / lambdax)  # Priestley-Taylor evapotranspiration at T_e
    E_SJ_Act_Daily = 2 * E_PT_T_e - Epenman_Daily  # actual evapotranspiration by Szilagyi and Jozsa (2008) (S8.7)

    if zerocorr:
        ET_Daily = np.fmax(E_SJ_Act_Daily, 0.)
    else:
        ET_Daily = E_SJ_Act_Daily

    return ET_Daily,Epenman_Daily,E_PT_T_e


def resample_DWD(stat='Norderney',calc_et=False, histo=False, ETclean0=True, rss='1D'):
    if (rss=='1H'):
        get_freq='hourly'
    else:
        get_freq='daily'

    if get_freq == 'hourly':
        dummy3 = dwd_recent(stat,'hourly','temp')
        dummy1 = dwd_recent(stat,'hourly','precip')
        try:
            dummy4 = dwd_recent(stat,'hourly','solar')
        except:
            dummy4 = dwd_recent('Norderney','hourly','solar')
            print('Used Norderney Solar radiation because no radiation data at station!')
        dummy5 = dwd_recent(stat,'hourly','wind')
        dummy6 = dwd_recent(stat,'hourly','pressure')
    
        try:
            dummyx = pd.concat([dummy3.TT_TU.resample(rss).mean(), # > AirT_2m °C
                        dummy3.TT_TU.resample(rss).min(),
                        dummy3.TT_TU.resample(rss).max(),
                        dummy1['  R1'].resample(rss).sum(), # > precip mm
                        dummy4.FG_LBERG.resample(rss).sum().loc[dummy3.index[0]:], # > hourly sum of solar incoming radiation J/cm^2 
                        dummy3.RF_TU.resample(rss).min(), # > RH_2m
                        dummy3.RF_TU.resample(rss).max(),
                        dummy5['   F'].resample(rss).mean(), # > mean windspeed m/s
                        dummy6['   P'].resample(rss).mean()],axis=1) # > mean air pressure hPa
            dummyx.columns = ['T','Tmin','Tmax','Prec','Rad','RHmin','RHmax','u2','aP']
        except:
            dummyx = pd.concat([dummy3.TT_TU.resample(rss).mean(), # > AirT_2m °C
                        dummy3.TT_TU.resample(rss).min(),
                        dummy3.TT_TU.resample(rss).max(),
                        dummy1['  R1'].resample(rss).sum(), # > precip mm
                        dummy3.RF_TU.resample(rss).min(), # > RH_2m
                        dummy3.RF_TU.resample(rss).max(),
                        dummy5['   F'].resample(rss).mean(), # > mean windspeed m/s
                        dummy6['   P'].resample(rss).mean()],axis=1) # > mean air pressure hPa
            dummyx.columns = ['T','Tmin','Tmax','Prec','RHmin','RHmax','u2','aP']


    if histo == True:
        if get_freq == 'hourly':
            try:
                dummy4h = dwd_recent(stat, get_freq, 'solar', 'historical')
                rad = dummy4h.FG_LBERG.resample(rss).sum()
            except:
                dummy4h = dwd_recent('Norderney', get_freq, 'solar','historical')
                print('Used Norderney Solar radiation because no radiation data at station!')
                rad = dummy4h.FG_LBERG.resample(rss).sum()
            dummy1h = dwd_recent(stat,get_freq,'precip','historical')
            dummy3h = dwd_recent(stat,get_freq,'temp','historical')
            dummy5h = dwd_recent(stat,get_freq,'wind','historical')
            dummy6h = dwd_recent(stat,get_freq,'pressure','historical')

            try:
                dummyxh = pd.concat([dummy3h.TT_TU.resample(rss).mean(), # > AirT_2m °C
                             dummy3h.TT_TU.resample(rss).min(),
                             dummy3h.TT_TU.resample(rss).max(),
                             dummy1h['  R1'].resample(rss).sum(), # > precip mm
                             rad.loc[dummy3h.index[0]:], # > hourly sum of solar incoming radiation J/cm^2
                             dummy3h.RF_TU.resample(rss).min(), # > RH_2m
                             dummy3h.RF_TU.resample(rss).max(),
                             dummy5h['   F'].resample(rss).mean(), # > mean windspeed m/s
                             dummy6h['   P'].resample(rss).mean()],axis=1) # # > mean air pressure hPa

                dummyxh.columns = ['T','Tmin','Tmax','Prec','Rad','RHmin','RHmax','u2','aP']
            except:
                dummyxh = pd.concat([dummy3h.TT_TU.resample(rss).mean(), # > AirT_2m °C
                             dummy3h.TT_TU.resample(rss).min(),
                             dummy3h.TT_TU.resample(rss).max(),
                             dummy1h['  R1'].resample(rss).sum(), # > precip mm
                             dummy3h.RF_TU.resample(rss).min(), # > RH_2m
                             dummy3h.RF_TU.resample(rss).max(),
                             dummy5h['   F'].resample(rss).mean(), # > mean windspeed m/s
                             dummy6h['   P'].resample(rss).mean()],axis=1) # # > mean air pressure hPa

                dummyxh.columns = ['T','Tmin','Tmax','Prec','RHmin','RHmax','u2','aP']

            dummyx = pd.concat([dummyxh.loc[:dummyx.index[0]],dummyx])
            dummyx = dummyx[~pd.Index.duplicated(dummyx.index,keep='first')]

        elif get_freq == 'daily':
            dummy3h = dwd_recent(stat,get_freq,'kl','historical')
            dummy4h = dwd_recent(stat,get_freq,'solar','historical')
            
            try:
                dsx = dummy4h.FG_STRAHL
            except:
                dsx = pd.Series(np.ones(len(dummy3h))*np.nan,index=dummy3h.index)

            dummyxh = pd.concat([dummy3h[' TMK'], # > AirT_2m °C
                             dummy3h[' TNK'],
                             dummy3h[' TXK'],
                             dummy3h[' RSK'], # > precip mm
                             dsx, # > daily sum of solar incoming radiation J/cm^2 
                             dsx*0.01,
                             dummy3h[' UPM'], # > RH_2m
                             dummy3h['  FM'], # > mean windspeed m/s
                             dummy3h['  FX'], # > mean windspeed m/s
                             dummy3h[' VPM'], # > mean vapor pressure hPa
                             dummy3h['  PM']],axis=1) # # > mean air pressure hPa
        

            dummyxh.columns = ['T','Tmin','Tmax','Prec','Rad','Rs','RH','u2','u2mx','vap','aP']
        
            dummyx = dummyxh

    if ((calc_et==True) & (rss=='1D')):
        import pyeto as pt
        if dummyxh.columns[6] == 'RH':
            EToPM = pt.fao56_penman_monteith(dummyx.Rs.values,dummyx['T'].values+273.15,dummyx.u2.values,pt.svp_from_t(dummyx['T'].values),dummyx.vap*0.1,pt.delta_svp(dummyx['T'].values),pt.psy_const(dummyx.aP.values*0.1))
        else:    
            EToPM = pt.fao56_penman_monteith(dummyx.Rs.values,dummyx['T'].values+273.15,dummyx.u2.values,pt.svp_from_t(dummyx['T'].values),pt.avp_from_rhmin_rhmax(pt.svp_from_t(dummyx.Tmin.values), pt.svp_from_t(dummyx.Tmax.values), dummyx.RHmin.values, dummyx.RHmax.values),pt.delta_svp(dummyx['T'].values),pt.psy_const(dummyx.aP.values*0.1))
        EToPM = pd.Series(EToPM)
        EToPM.index = dummyx.index
        EToHG = pt.hargreaves(dummyx.Tmin.values,dummyx.Tmax.values,dummyx['T'].values,pt.et_rad(52. * np.pi/180.,pt.sol_dec(dummyx.index.dayofyear.values),pt.sunset_hour_angle(52. * np.pi/180.,pt.sol_dec(dummyx.index.dayofyear.values)),pt.inv_rel_dist_earth_sun(dummyx.index.dayofyear.values)))
        EToHG = pd.Series(EToHG)
        EToHG.index = dummyx.index
        
        EToSJ,EToPM2,EToPT = ET_SzilagyiJozsa(dummyx , 3. , 53.4, zerocorr=ETclean0) 

        #if dummyxh.columns[5] == 'RH':
        #    ETact,vap = ET_SJ(dummyx.Tmax,dummyx.Tmin,dummyx.RH,None,dummyx.Rad*0.01,dummyx.u2*3.6,dummyx.index.dayofyear,53.,1.,0.23,0.02)
        #else:
        #    ETact,vap = ET_SJ(dummyx.Tmax,dummyx.Tmin,dummyx.RHmax,dummyx.RHmin,dummyx.Rad*0.01,dummyx.u2*3.6,dummyx.index.dayofyear,53.,1.,0.23,0.02)
        #if ETclean0:
        #    ETact[ETact>0.] = 0.
        #    vap[vap>10.] = np.nan
        dummyx = pd.concat([dummyx,EToPM,EToHG,EToSJ,EToPM2,EToPT],axis=1)
        if (get_freq == 'daily'):
            dummyx.columns = ['T','Tmin','Tmax','Prec','Rad','Rs','RH','u2','u2mx','vap','aP','EToPM','EToHG','EToSJ','EToPM2','EToPT']
        else:
            dummyx.columns = ['T','Tmin','Tmax','Prec','Rad','RHmin','RHmax','u2','aP','EToPM','EToHG','EToSJ','EToPM2','EToPT']
        
    return dummyx