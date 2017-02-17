# -*- coding: utf-8 -*-
"""
Auteur : Kamel Chikhar

Objet : Controle de qualite des données de bouées IABP
        (site de téléchargement : ftp://iabp.apl.washington.edu/pub/IABP/C/)

Input : Parametres de controle contenus dans le fichier PARAM
"""

import numpy as np
from datetime import datetime, timedelta
#import Fonctions as rod
from os import listdir
from os.path import isfile, join
import pandas as pd
import scipy.interpolate as interpolate
import rpnpy.librmn.all as rmn
import os, sys
from collections import defaultdict
import params_bouees2 as params
from timeit import default_timer as timer
import cPickle as pickle
import math
from math import radians

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, -1)

# sample usage
#save_object(company1, 'company1.pkl')


def haversine(lat1, lon1, lat2, lon2):
    '''
    Calcule la distance entre deux points sur la sphere (terre)
    :param lat1: latitude du premier point
    :param lon1: longitude du premier point
    :param lat2: latitude du deuxieme point
    :param lon2: longitude du deuxieme point
    :return: distance en metres
    '''
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371*1000  # Radius of earth in kilometers. Use 3956 for miles
    return c * r


def grille(date):
    # Function that reads a standard file and builds an interpolator (for ice concentration variable) given a file date

    filein = '/fs/cetus/fs3/mrb/armn/armnkch/maestro_suites/creg025/sorties/cice5_avec_formdrag/cice_output/'+ str(date) + '_000d'  #2004091800_000d
    fileName = os.path.join(filein)
    try:
        fileId = rmn.fstopenall(fileName)
    except:
        sys.stderr.write("Problem opening the file: %s\n" % fileName)
    try:
        keylist = rmn.fstinl(fileId, nomvar='>>') # Get the list of record keys matching nomvar='>>'
        for k in keylist:
            # Get the record meta data
            r = rmn.fstluk(k, dtype=np.float32)
            grille_lon = r['d']

        keylist = rmn.fstinl(fileId, nomvar='^^') # Get the list of record keys matching nomvar='^^'
        for k in keylist:
            # Get the record meta data
            r = rmn.fstluk(k, dtype=np.float32)
            grille_lat = r['d']

        keylist = rmn.fstinl(fileId, nomvar='aice') # Get the list of record keys matching nomvar='aice'
        k = keylist[0]
        r = rmn.fstluk(k, dtype=np.float32)
        concentration = r['d']

    except:
        pass
    finally:
        # Close file even if an error occured above
        rmn.fstcloseall(fileId)
    gl_interpolator = interpolate.LinearNDInterpolator(np.vstack((grille_lon.flatten(), grille_lat.flatten())).T,
                                                       concentration.flatten())
    return gl_interpolator

def interpolation(date, lat, lon):
    # Function that interpolates value given lat/lon

    interpolator = dd.values()[dates.index(date)]
    if lon < 0:
        data_int = interpolator(lon + 360, lat)
    else:
        data_int = interpolator(lon, lat)

    print 'Interpolation pour la date : ', int(date)
    print 'Coordonnees Lat, Lon : ', lat, lon, data_int

    return data_int


start = timer()

startDate = params.startDate
endDate = params.endDate
raw_data = params.path_bouee
latMin = params.Latmin
latMax = params.Latmax
lonMin = params.Lonmin
lonMax = params.Lonmax
seuil_concentration = params.seuil_de_concentration


data = pd.read_csv(raw_data, sep="\s+", header = None, names=['annee', 'mois', 'jour','heure', 'id', 'sLat','sLon'])


data['temps'] = data.apply(lambda arow: int(datetime(int(arow['annee']), int(arow['mois']), int(arow['jour']), int(arow['heure'])).strftime("%Y%m%d%H")), axis=1)


data = data.drop(data.columns[[0, 1, 2, 3]], axis=1)

data = data[['id', 'temps','sLat','sLon']]

# Filtrer les dates

data[8] = data['temps'].apply(lambda x:  datetime.strptime(str(x),'%Y%m%d%H'))
data[9] = data[8].map(lambda x: x.hour)
data = data[(data[8] >= datetime.strptime(startDate, '%Y%m%d')) & (data[8] <= datetime.strptime(endDate, '%Y%m%d')) & (data[9] == 0)]
print ' Nombre d''observations apres filtrage des dates : ', len(data)


#**************************************************************************************
# Filtre 1 : on elimine toutes les bouées situées à l'exterieur de la region considérée
#        latMin <=  lat <= latMax
#        lonMin <=  lon <= lonMax
#*************************************************************************************

indf = data[(data['sLat'] >= latMin) & (data['sLat'] <= latMax) & (data['sLon'] >= lonMin) & (data['sLon'] <= lonMax)]
print ' Nombre d''observations apres filtrage des lat/lon : ', len(data['id'])

indf = indf.drop(indf.columns[[4, 5]], axis=1)


day = timedelta(days=1)
date_format = "%Y%m%d%H"
outdf = indf.copy()
outdf['temps+24'] = indf['temps'].map(lambda d: int((datetime.strptime(str(d), date_format) + day).strftime(date_format)))
st_date_to_lat_lon = {}
for st, data in outdf.groupby('id'):
    st_date_to_lat_lon[st] = defaultdict(lambda: (np.nan, np.nan))
    for row_i, row in data.iterrows():
        st_date_to_lat_lon[st][row['temps']] = (row['sLat'], row['sLon'])

outdf['eLat'] = outdf.apply(lambda arow: st_date_to_lat_lon[arow['id']][arow['temps+24']][0], axis=1)
outdf['eLon'] = outdf.apply(lambda arow: st_date_to_lat_lon[arow['id']][arow['temps+24']][1], axis=1)
# Some dates do not have the corresponding data a day later (so drop nans),

outdf = outdf.dropna()
print 'Nombre de bouees apres filtrages des bouees avec H+24: ', len(outdf)

#Filtrage des doublons
for i in range(len(outdf)):
    if (outdf['sLat'].values[i] == outdf['eLat'].values[i]) and (outdf['sLon'].values[i] == outdf['eLon'].values[i]):
        outdf.drop(outdf.index[[i]])

print 'apres filtrage des doublons: ', len(outdf)

# Filtrage des vitesses
sim_length_S = 24 * 3600
outdf['vitesse'] = outdf.apply(lambda row: haversine(row['sLat'], row['sLon'], row['eLat'], row['eLon']) / sim_length_S, axis=1)

outdf2 = outdf[(outdf['vitesse'] >= 1)]
outdf = outdf[(outdf['vitesse'] < 1)]


global dd, dates
dates = list(np.unique(outdf["temps"]))

dd = dict((h, None) for h in dates)
for i in range(len(dates)):
    print dates[i]
    dd[dates[i]] = grille(dates[i])


outdf['aice'] = outdf.apply(lambda row: interpolation(row['temps'], row['sLat'], row['sLon']), axis=1)




# On filtre les bouees avec concentration < 0.5
outdf = outdf[outdf['aice'] >= seuil_concentration]
print ' Nombre de bouees apres filtrage des concentrations: ', len(outdf)



format="%10.0f %10.0f %6.3f %8.3f %10.0f %6.3f %8.3f %5.3f"
with open('/users/dor/armn/kch/scripts/controled_buoys_avec_2.csv', 'a') as f:
    np.savetxt(f,zip(outdf['id'].values,
                     outdf['temps'].values,
                     outdf['sLat'].values,
                     outdf['sLon'].values,
                     outdf['temps+24'].values,
                     outdf['eLat'].values,
                     outdf['eLon'].values,
                     outdf['aice'].values),
               fmt=format)
            
if len(outdf2) > 0:
    with open('/users/dor/armn/kch/scripts/controled_buoys_avec_rejected.csv', 'a') as f2:
         np.savetxt(f2,zip(outdf2['id'].values,
                     outdf2['temps'].values,
                     outdf2['sLat'].values,
                     outdf2['sLon'].values,
                     outdf2['temps+24'].values,
                     outdf2['eLat'].values,
                     outdf2['eLon'].values,
                     outdf2['vitesse'].values),
                 fmt=format)
               
              

dates_f = list(np.unique(outdf["temps"]))
nb_bouees=[]
for d in dates_f :
    nb_bouees.append(len(outdf[outdf['temps'] == d]['id']))
for i in range(len(dates_f)):
    print 'Date: ', dates[i],'Nobre de bouees considerees: ', nb_bouees[i]

#df_2004_11 = outdf

#save_object(df_2004_10_12, 'QualityControl.pkl')
#with open('QualityControl.pkl', 'ab') as f:
#    pickle.dump(df_2004_11, f)

end = timer()
print ' Finished in : ', (end - start), ' s'
#..fin..


