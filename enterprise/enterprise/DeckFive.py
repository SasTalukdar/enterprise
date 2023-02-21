#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 01:39:01 2021
@author: picard

This module contains data input output related functions
"""

import re
import numpy as np
import datetime as dt
import pandas as pd

######################################################################################################
def afloat(x):
    try :
        x=float(x)
    except :
        x=np.nan
    return x
######################################################################################################

class rsrw :
    # functions to handle university of wyoming rsrw data
    def extract_indices(file):
        with open(file,'r') as f:
            data=f.readlines()
        mon_name=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov', 'Dec']
        variables=['station','date','utc','lat','lon','elev','show','li','lift','sweat','ki','cti',
            'vti','tti','cape','cape_vir','cin','cin_vir','eqlb_level','eql_vir','lfc','lfc_vir',
            'brn','brn_cape','lcl_t','lcl_p','m_pot_t','m_mx_rat','prec']
        df = pd.DataFrame(columns=variables)
        start=False
        for i, line in enumerate(data):
            if 'Observations at' in line :
                if start :
                    df = pd.DataFrame(np.insert(df.values, len(df.index), values = [station,dt.date(year,mon,day),
                                     time,lat,lon,elev,show,li,lift,sweat,ki,cti,vti,tti,cape,cape_vir,cin,cin_vir,
                                     eqlb_level,eql_vir,lfc,lfc_vir,brn,brn_cape,lcl_t,lcl_p,m_pot_t,m_mx_rat,prec], axis=0))
                    df.columns = variables
                t_ind=line.index('Observations at')+len('Observations at')+1
                obs_time=line[t_ind:t_ind+15]
                time=int(obs_time[:2])
                day=int(obs_time[4:6])
                mon=int(mon_name.index(obs_time[7:10])+1)
                year=int(obs_time[-4:])
                [lat,lon,elev,show,li,lift,sweat,ki,cti,vti,tti,cape,cape_vir,
                  cin,cin_vir,eqlb_level,eql_vir,lfc,lfc_vir,brn,brn_cape,lcl_t,
                  lcl_p,m_pot_t,m_mx_rat,prec] = np.nan*np.ones(26)
                start = True
            elif 'Station number' in line and '****' not in line:
                station=float(line[line.index('number:')+len('number:'):-1])
            elif 'Station latitude' in line and '****' not in line:
                lat=float(line[line.index('latitude:')+len('latitude:'):-1])
            elif 'Station longitude' in line and '****' not in line:
                lon=float(line[line.index('longitude:')+len('longitude:'):-1])
            elif 'Station elevation' in line and '****' not in line:
                elev=float(line[line.index('elevation:')+len('elevation:'):-1])
            elif 'Showalter index' in line and '****' not in line:
                show=float(line[line.index('index:')+len('index:'):-1])
            elif 'Lifted index' in line and '****' not in line:
                li=float(line[line.index('index:')+len('index:'):-1])
            elif 'LIFT computed using virtual temperature:' in line and '****' not in line:
                lift=float(line[line.index('temperature:')+len('temperature:'):-1])
            elif 'SWEAT index' in line and '****' not in line:
                sweat=float(line[line.index('index:')+len('index:'):-1])
            elif 'K index' in line and '****' not in line:
                ki=float(line[line.index('index:')+len('index:'):-1])
            elif 'Cross totals index:' in line and '****' not in line:
                cti=float(line[line.index('index:')+len('index:'):-1])
            elif 'Vertical totals index:' in line and '****' not in line:
                vti=float(line[line.index('index:')+len('index:'):-1])
            elif 'Totals totals index:' in line and '****' not in line:
                tti=float(line[line.index('index:')+len('index:'):-1])
            elif 'Convective Available Potential Energy:' in line and '****' not in line:
                cape=float(line[line.index('Energy:')+len('Energy:'):-1])
            elif 'CAPE using virtual temperature:' in line and '****' not in line:
                cape_vir=float(line[line.index('temperature:')+len('temperature:'):-1])
            elif 'Convective Inhibition:' in line and '****' not in line:
                cin=float(line[line.index('Inhibition:')+len('Inhibition:'):-1])
            elif 'CINS using virtual temperature:' in line and '****' not in line:
                cin_vir=float(line[line.index('temperature:')+len('temperature:'):-1])    
            elif 'Equilibrum Level:' in line and '****' not in line:
                eqlb_level=float(line[line.index('Level:')+len('Level:'):-1])
            elif 'Equilibrum Level using virtual temperature:' in line and '****' not in line:
                eql_vir=float(line[line.index('temperature:')+len('temperature:'):-1])
            elif 'Level of Free Convection:' in line and '****' not in line:
                lfc=float(line[line.index('Convection:')+len('Convection:'):-1])
            elif 'LFCT using virtual temperature:' in line and '****' not in line:
                lfc_vir=float(line[line.index('temperature:')+len('temperature:'):-1])
            elif 'Bulk Richardson Number:' in line and '****' not in line:
                brn=float(line[line.index('Number:')+len('Number:'):-1])
            elif 'Bulk Richardson Number using CAPV:' in line  and '****' not in line:
                brn_cape=float(line[line.index('CAPV:')+len('CAPV:'):-1])
            elif 'Temp [K] of the Lifted Condensation Level:' in line and '****' not in line:
                lcl_t=float(line[line.index('Level:')+len('Level:'):-1])
            elif 'Pres [hPa] of the Lifted Condensation Level:' in line and '****' not in line:
                lcl_p=float(line[line.index('Level:')+len('Level:'):-1])
            elif 'Mean mixed layer potential temperature:' in line and '****' not in line:
                m_pot_t=float(line[line.index('temperature:')+len('temperature:'):-1])
            elif 'Mean mixed layer mixing ratio:' in line and '****' not in line:
                m_mx_rat=float(line[line.index('ratio:')+len('ratio:'):-1])
            elif 'Precipitable water [mm] for entire sounding:' in line and '****' not in line:
                prec=float(line[line.index('sounding:')+len('sounding:'):-1])
        return df
    
    def extract_prof(file):
        with open(file,'r') as f:
            data=f.readlines()
        mon_name=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov', 'Dec']
        variables=['station','date','utc','PRES','HGHT','TEMP','DWPT','RELH',
                   'MIXR','DRCT','SKNT','THTA','THTE','THTV']
        df = pd.DataFrame(columns=variables)
        for i, line in enumerate(data):
            if 'Observations at' in line :
                station = line[:5]
                t_ind=line.index('Observations at')+len('Observations at')+1
                obs_time=line[t_ind:t_ind+15]
                time=int(obs_time[:2])
                day=int(obs_time[4:6])
                mon=int(mon_name.index(obs_time[7:10])+1)
                year=int(obs_time[-4:])
                [lat,lon,elev,PRES,HGHT,TEMP,DWPT,RELH,MIXR,DRCT,SKNT,THTA,THTE,THTV] = np.nan*np.ones(14)
            elif not re.findall('[a-zA-Z]', line) and len(re.findall('[0-9]', line)) > 0 :
                PRES=afloat(line[0:7])
                HGHT=afloat(line[7:14])
                TEMP=afloat(line[14:21])
                DWPT=afloat(line[21:28])
                RELH=afloat(line[28:35])
                MIXR=afloat(line[35:42])
                DRCT=afloat(line[42:49])
                SKNT=afloat(line[49:56])
                THTA=afloat(line[56:63])
                THTE=afloat(line[63:70])
                THTV=afloat(line[70:77])
                df = pd.DataFrame(np.insert(df.values, len(df.index), values = [station,dt.date(year,mon,day),time,
                                                                                    PRES,HGHT,TEMP,DWPT,RELH,
                                                                                    MIXR,DRCT,SKNT,THTA,THTE,THTV], axis=0))
        df.columns = variables     
        return df
