# -*- coding: utf-8 -*-
"""
Utilité du code :
    Ce code vise à lire un fichier d'observation de bouées afin d'en ressortir
    les positions de départ et d'arrivée de différents points d'observation. Le fichier d'observation
    doit d'abord avoir été traité à l'aide de QualityControlBuoy1.py pour qu'il respecte une certaine structure.
    Il prend également quelques paramètres en entrée qui peuvent être modifiés dans le fichier d'entrée PARAM.
    
    Une fois que les vitesses des observations ont été calculées, l'écriture d'un fichier d'entrée
    pour MLDP commence. Chaque position de départ du fichier d'observation est écrite à l'intérieur
    du fichier d'entrée. MLDP est ensuite lancé. Il est également possible de "bypasser" cette étape. 
    (bypassMLDP=True). Si c'est le cas, les fichiers de sorties de MLDP DOIVENT déjà exister, 
    sinon il y aura une erreur.
    
    Finalement, la lecture des fichiers de sorties de MLDP est exécuté et des statistiques de vitesse
    et d'angle sont creées. L'écriture d'un fichier contenant les vitesses et positions des observations
    et du modèle est faite afin de pouvoir exécuter BasemapBouee.py Un graphique de comparaison de vitesse
    et d'angle est généré
    
Mode d'emploi:
    S'assurer que le fichier texte PARAM est bien rempli selon le format nécessaire. Choisir ensuite
    les options de bypass MLDP . 
    Pour obtenir les résultats pour toutes les dates,  mettre la variable "intervalle_temps"=24. 
    Pour obtenir les résultats qui commencent les mercredis,  mettre la variable "intervalle_temps"=168.


    Important : on doit executer au prealable le script QualityControlBuoy_v1.py
"""    
from timeit import default_timer as timer

import os, sys
import numpy as np
import datetime
import shutil
import math
import rpnpy.librmn.all as rmn
import rpnpy.rpndate as rpndate
#import Fonctions as rod
from os import listdir
from os.path import isfile, join
import pandas as pd
import subprocess
#import matplotlib.pyplot as plt
from shutil import copyfile
import params_bouees as params
from math import radians
import cPickle as pickle
start = timer()
sys.stdout = open('/fs/cetus/fs3/mrb/armn/armnkch/VerifGlace/stdout.txt', 'w')

def replace_input_MLDP(file, date, id, lat, lon, sim_length):
    """
    Prepare input file for MLDP
    :param file: filename
    :param date: date
    :param lat: starting latitude
    :param lon: starting longitude
    :param id: particule name
    :param sim_length: simulation length
    """
    d = datetime.datetime.strptime(date, "%Y%m%d%H")
    delta = datetime.timedelta(hours=3)
    prev_dates = [d + i * delta for i in range(sim_length / 3 + 1)]

    with open(file) as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i].strip()[0:12] == 'SRC_NAME   =':
                lines[i] = lines[i].strip()[0:12] + ' ' + id + '                 # Buoy identification\n'
            if lines[i].strip()[0:12] == 'SRC_TIME   =':
                lines[i] = lines[i].strip()[0:12] + ' ' + date + '00' + lines[i][30:]
            if lines[i].strip()[0:12] == 'SRC_COORD  =':
                lines[i] = lines[i].strip()[0:12] + ' ' + str(lat) + ' ' + str(lon) + '     #  Starting coordinates \n'
            if lines[i].strip()[0:11] == 'MET_FILES =':
                for k in range(sim_length / 3 + 1):
                    str1 = datetime.datetime.strftime(prev_dates[k], "%Y-%m-%d")
                    h = datetime.datetime.strftime(prev_dates[k], "%H").zfill(3)
                    str2 = datetime.datetime.strftime(prev_dates[k], "%Y%m%d") + '00'
                    if (prev_dates[k] - d).days ==1:
                        k1 = prev_dates[k] - datetime.timedelta(days=1)
                        str1 = datetime.datetime.strftime(k1, "%Y-%m-%d")
                        h = '024'
                        str2 = datetime.datetime.strftime(k1, "%Y%m%d")+'00'
                    lines[i+k+1] = '/fs/cetus/fs3/mrb/armn/armnkch/VerifGlace/MLDP_champs/cice5_avec_formdrag/%s/%s_%s_pn5km\n' %(str1, str2, h)
    f.close()
    if os.path.isfile('tmpfile'):
        os.remove('tmpfile')

    with open('tmpfile','w') as f:
        f.writelines(lines)
    f.close()
    shutil.move('tmpfile', file)
    return()

def rad(x):
    y =  x * math.pi / 180
    return y

def deg(x):
    y = x * 180 / math.pi
    return y

def angle(a,b):
    '''
    Calcule l'angle a partir du Nord
    :param a: coordonnee longitudinale
    :param b: coordonnee zonale
    :return: angle en degres
    '''
    global alpha
    case0 = (a > 0 and b > 0)
    case1 = (a < 0 and b > 0)
    case2 = (a < 0 and b < 0)
    case3 = (a > 0 and b < 0)
    case4 = (a == 0 and b > 0)
    case5 = (a == 0 and b < 0)
    case6 = (b == 0 and a > 0)
    case7 = (b == 0 and a < 0)

    if case0:
        alpha = deg(np.arctan(np.abs(b)/np.abs(a)))
    if case1:
        alpha = 180 - deg(np.arctan(np.abs(b)/np.abs(a)))
    if case2:
        alpha = 180 + deg(np.arctan(np.abs(b)/np.abs(a)))
    if case3:
        alpha = 360 - deg(np.arctan(np.abs(b)/np.abs(a)))
    if case4:
        alpha = 90
    if case5:
        alpha = 270
    if case6:
        alpha = 360
    if case7:
        alpha = 180

    return alpha


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


def replace_lanceur_MLDP(date):
    """
    Prepare le lancement MLDP
    :param date: date  initiale pour le lanceur
    :param id : indcatif de la bouee
    """
    with open('lancement_MLDPn_buoys') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i].strip()[0:6] == 'rm -rf':
                lines[i] = lines[i].strip()[0:6] + ' ' + '/fs/cetus/fs3/mrb/armn/armnkch/VerifGlace/MLDP_Output/SortiesBouees/' + date + ' \n'
            if lines[i].strip()[0:8] == 'mkdir -p':
                lines[i] = lines[i].strip()[0:8] + ' ' + '/fs/cetus/fs3/mrb/armn/armnkch/VerifGlace/MLDP_Output/SortiesBouees/' + date + ' \n'
            if lines[i].strip()[0:8] == 'fileOut=':
                lines[i] = lines[i].strip()[0:8] + '/fs/cetus/fs3/mrb/armn/armnkch/VerifGlace/MLDP_Output/SortiesBouees/' + date + ' \n'
        f.close()
    if os.path.isfile('lanceur_MLDPn_buoys'):
        os.remove('lanceur_MLDPn_buoys')
    with open('lanceur_MLDPn_buoys', 'w') as f:
        f.writelines(lines)
    f.close()
    os.chmod('lanceur_MLDPn_buoys', 0775)
    return()

bypassMLDP = False
sim_length = 24
sim_length_S = 24 * 3600  # en secondes
path_B='/users/dor/armn/kch/scripts/controled_buoys_avec_2.csv'
MLDP_input_template = '/users/dor/armn/kch/scripts/sources_python/validation_cice/inputMLDP_buoys'

Rterre=6371.0*1000

with open(MLDP_input_template) as f:
    gen_lines= f.readlines()
for i in range(len(gen_lines)):
    print i , gen_lines[i]


if not os.path.isfile(path_B):
    print "***** Erreur : Veuillez utiliser le QualityControlBuoy_v1.py avant de lancer ce script!!! *****"

data = pd.read_csv(path_B, sep="\s+", header = None, names=['id', 'temps', 'sLat', 'sLon', 'temps+24', 'eLat', 'eLon', 'aice'])

# Selectionner les bouees ayant des vitesses superieures a 5 cm/s

#data['vitesse'] = data.apply(lambda row: haversine(row['sLat'], row['sLon'], row['eLat'], row['eLon']) / sim_length_S, axis=1)

#data = data[(data['vitesse'] >= 0.05)]

s0 = params.startDate + '00'
s1 = params.endDate + '00'
d0 = datetime.datetime.strptime(s0, "%Y%m%d%H")
d1 = datetime.datetime.strptime(s1, "%Y%m%d%H")
oneday = datetime.timedelta(days=1)
datesStr = [d0 + oneday * i for i in range((d1 - d0).days + 1)]
work_dates = [datetime.datetime.strftime(d, "%Y%m%d%H") for d in datesStr] 
dates = np.unique(data['temps'].values)
#print datesStr 
#sys.exit('test')
delta = datetime.timedelta(hours=3)

Nombre_bouees = np.zeros(len(dates))

# Nombre de bouees par date:
for i in range(len(dates)):
    Nombre_bouees[i] = len(data[data['temps'] == dates[i]])
    print 'Date: ', dates[i], '   Nombre de bouees: ', int(Nombre_bouees[i])

# Initialisation:
vitesse_B_moy = []
vitesse_B_std = []
vitesse_B_min = []
vitesse_B_max = []
vitesse_B2_moy = []

angle_B_moy = []
angle_B_std = []
angle_B_min = []
angle_B_max = []

vitesse_MLDP_moy = []
vitesse_MLDP_std = []
vitesse_MLDP_min = []
vitesse_MLDP_max = []

angle_MLDP_moy = []
angle_MLDP_std = []
angle_MLDP_min = []
angle_MLDP_max = []

distance_moy = []
distance_std = []
distance_min = []
distance_max = []

angle_diff_moy = []
angle_diff_std = []
angle_diff_min = []
angle_diff_max = []

nb_bouees = []

vitesse_bias_moy = []
vitesse_bias_std = []
angle_bias_moy = []
angle_bias_std = []


# Extract working dates
x = dates[dates >= int(s0)]
x = x[x <= int(s1)]
dates = x
print dates
print type(dates)
#sys.exit('test')
for date in dates:
    wf = data[data['temps'] == date]
    print 'On travaille sur la date:  ', date
    d = datetime.datetime.strptime(str(date), "%Y%m%d%H")
    prev_dates = [d + i * delta for i in range(sim_length / 3 + 1)]
    tt = str(date) + '00'
    if bypassMLDP == False:   # Dans ce cas on utilise MLDP pour calculer les trajectoies des particules
        if os.path.isfile('/fs/cetus/fs3/mrb/armn/armnkch/VerifGlace/MLDP_tmp_buoys'):
            os.remove('/fs/cetus/fs3/mrb/armn/armnkch/VerifGlace/MLDP_tmp_buoys')
        with open('/fs/cetus/fs3/mrb/armn/armnkch/VerifGlace/MLDP_tmp_buoys', 'w') as InputFile:
            InputFile.writelines(gen_lines[:14])
            for indicatif in wf['id'].values:
                startLat = wf[wf['id'] == indicatif]['sLat'].values
                startLon = wf[wf['id'] == indicatif]['sLon'].values
                InputFile.writelines(gen_lines[14])
                InputFile.write('SRC_NAME   =    %s     # Buoy identification  \n' % str(indicatif))
                InputFile.write('SRC_TIME   =    %s     # Emission date  \n'  % tt)
                InputFile.writelines(gen_lines[17])
                InputFile.write('SRC_COORD  =    %s   %s  \n' %(float(startLat), float(startLon)))
                InputFile.writelines(gen_lines[19:23])
            InputFile.writelines(gen_lines[23:40])                            #2004122800_021_pn5km
            for k in range(sim_length / 3 + 1):
                str1 = datetime.datetime.strftime(prev_dates[k], "%Y-%m-%d")
                h = datetime.datetime.strftime(prev_dates[k], "%H").zfill(3)
                str2 = datetime.datetime.strftime(prev_dates[k], "%Y%m%d") + '00'
                if (prev_dates[k] - d).days == 1:
                    k1 = prev_dates[k] - datetime.timedelta(days=1)
                    str1 = datetime.datetime.strftime(k1, "%Y-%m-%d")
                    h = '024'
                    str2 = datetime.datetime.strftime(k1, "%Y%m%d") + '00'
                InputFile.write('/fs/cetus/fs3/mrb/armn/armnkch/VerifGlace/MLDP_champs/cice5_avec_formdrag/%s/%s_%s_pn5km \n' % (
                str1, str2, h))

        InputFile.close()
        replace_lanceur_MLDP(str(date))
        x=subprocess.call('/users/dor/armn/kch/scripts/sources_python/validation_cice/lanceur_MLDPn_buoys', shell=True)

    # On recupere les lat/lon simulees par MLDP
    pathMLDP = '/fs/cetus/fs3/mrb/armn/armnkch/VerifGlace/MLDP_Output/SortiesBouees/' + str(date)
    filesMLDP = [f for f in listdir(pathMLDP) if isfile(join(pathMLDP, f))]
    for files in filesMLDP:
        if files.endswith('024'):
            fileName=pathMLDP + '/' + files
    try:
        fileId = rmn.fstopenall(fileName)
    except:
        sys.stderr.write("Problem opening the file: %s\n" % fileName)
    try:
        # Get the list of record keys matching nomvar='LA', 'LO'
        keylist1 = rmn.fstinl(fileId, nomvar='LA')
        keylist2 = rmn.fstinl(fileId, nomvar='LO')
        MLDP_lat = np.zeros(len(keylist1))
        MLDP_lon = np.zeros(len(keylist1))
        MLDP_etiket = []
        for k in keylist1:
            r = rmn.fstluk(k)
            MLDP_lat[keylist1.index(k)] = r['d']
            MLDP_etiket.append(r['etiket'].strip())
        for k in keylist2:
            r = rmn.fstluk(k)
            MLDP_lon[keylist2.index(k)] = r['d'] - 360
    except:
        pass
    finally:
        # Close file even if an error occured above
        rmn.fstcloseall(fileId)

    # Calcul des vitesses (Haversine formula), angles et statistiques (bouees et sorties MLDP)
    # Bouees
    ll = len(wf)
    vitesse_B = np.zeros(ll)
    angle_B = np.zeros(ll)
    vitesse_B2 = np.zeros(ll)
    vitesse_MLDP = np.zeros(ll)
    angle_MLDP = np.zeros(ll)
    angle_diff = np.zeros(ll)
    distance = np.zeros(ll)


    for k in range(ll):
        vitesse_B[k] = haversine(wf['sLat'].values[k], wf['sLon'].values[k], wf['eLat'].values[k], wf['eLon'].values[k]) / sim_length_S
    vitesse_lat_B = rad(wf['eLat'].values - wf['sLat'].values) * Rterre / sim_length_S
#    vitesse_B = 2 * Rterre * np.arcsin(np.sqrt( (np.sin(rad(wf['eLat'].values - wf['sLat'].values)/2))**2 \
#                  + np.cos(rad(wf['eLat'].values)) * np.cos(rad(wf['sLat'].values)) * (np.sin(rad(wf['eLon'].values - wf['sLon'].values)/2))**2)) / sim_length_S
    vitesse_lon_B = (vitesse_B**2-vitesse_lat_B**2)**0.5
    vitesse_lon_B[wf['eLon'].values < wf['sLon'].values] = -1 * vitesse_lon_B[wf['eLon'].values < wf['sLon'].values] # changer le signe pour eLon<sLon
    for b in range(ll):
        angle_B[b] = angle(vitesse_lat_B[b], vitesse_lon_B[b])
        vitesse_B2[b] = haversine(wf['sLat'].values[b], wf['sLon'].values[b], wf['eLat'].values[b], wf['eLon'].values[b]) / sim_length_S
    # MLDP
    vitesse_lat_MLDP = rad(MLDP_lat - wf['sLat'].values) * Rterre / sim_length_S

    for k in range(ll):
        vitesse_MLDP[k] = haversine(wf['sLat'].values[k], wf['sLon'].values[k], MLDP_lat[k], MLDP_lon[k]) / sim_length_S
#    vitesse_MLDP = 2 * Rterre * np.arcsin(np.sqrt((np.sin(rad(MLDP_lat - wf['sLat'].values) / 2)) ** 2 \
#                    + np.cos(rad(MLDP_lat)) * np.cos(rad(wf['sLat'].values)) * (np.sin(rad(MLDP_lon - wf['sLon'].values) / 2)) ** 2)) / sim_length_S
    vitesse_lon_MLDP = (vitesse_MLDP ** 2 - vitesse_lat_MLDP ** 2) ** 0.5
    vitesse_lon_MLDP[MLDP_lon < wf['sLon'].values] = -1 * vitesse_lon_MLDP[MLDP_lon < wf['sLon'].values]  # changer le signe pour eLon<sLon


    for b in range(len(wf)):
        angle_MLDP[b] = angle(vitesse_lat_MLDP[b], vitesse_lon_MLDP[b])
        angle_diff[b] = angle_B[b] - angle_MLDP[b]
        distance[b] = haversine(wf['eLat'].values[b], wf['eLon'].values[b], MLDP_lat[b], MLDP_lon[b])
#    distance = 2 * Rterre * np.arcsin(np.sqrt((np.sin((rad(MLDP_lat - wf['eLat'].values)) / 2)) ** 2 \
#                   + np.cos(rad(wf['eLat'].values)) * np.cos(rad(MLDP_lat)) * (np.sin(rad(MLDP_lon - wf['eLon'].values) / 2)) ** 2))

    # Stats
    vitesse_B_moy.append(np.mean(vitesse_B))
    vitesse_B_std.append(np.std(vitesse_B))
    vitesse_B_min.append(np.min(vitesse_B))
    vitesse_B_max.append(np.max(vitesse_B))

    vitesse_MLDP_moy.append(np.mean(vitesse_MLDP))
    vitesse_MLDP_std.append(np.std(vitesse_MLDP))
    vitesse_MLDP_min.append(np.min(vitesse_MLDP))
    vitesse_MLDP_max.append(np.max(vitesse_MLDP))

    vitesse_bias = vitesse_MLDP - vitesse_B
    vitesse_bias_moy.append(np.mean(vitesse_bias))
    vitesse_bias_std.append(np.std(vitesse_bias))
 
    angle_B_moy.append(np.mean(angle_B))
    angle_B_std.append(np.std(angle_B))
    angle_B_min.append(np.min(angle_B))
    angle_B_max.append(np.max(angle_B))
  
    angle_MLDP_moy.append(np.mean(angle_MLDP))
    angle_MLDP_std.append(np.std(angle_MLDP))
    angle_MLDP_min.append(np.min(angle_MLDP))
    angle_MLDP_max.append(np.max(angle_MLDP))
    
    angle_bias = np.abs(angle_MLDP - angle_B)
    for i in range(len(angle_bias)):
         if angle_bias[i] > 180:
	   angle_bias[i] = 360 - angle_bias[i]

    angle_bias_moy.append(np.mean(angle_bias))
    angle_bias_std.append(np.std(angle_bias))

    distance_moy.append(np.mean(distance / 1000))  # en km
    distance_std.append(np.std(distance / 1000))  # en km
    distance_min.append(np.min(distance / 1000))  # en km
    distance_max.append(np.max(distance / 1000))  # en km

    nb_bouees.append(ll)
 
    vitesse_B2_moy.append(np.mean(vitesse_B2))
    format="%8.3f %8.3f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f"
    with open('/fs/cetus/fs3/mrb/armn/armnkch/VerifGlace/basemap/data_for_basemap_avec_' + str(date) + '.csv', 'w') as f:
        np.savetxt(f,zip(wf['sLon'].values,
                         wf['sLat'].values,
                         vitesse_lon_B,
                         vitesse_lat_B,
                         vitesse_B,
                         vitesse_lon_MLDP,
                         vitesse_lat_MLDP,
                         vitesse_MLDP),
                         fmt=format)
#creation de la liste de dates pour Pandas
liste_dates=[]
for date in dates:
    liste_dates.append(datetime.datetime.strptime(str(date), "%Y%m%d%H"))

#création des dictionnaires Pandas
OBS1_0 = pd.DataFrame(vitesse_B_moy, index=liste_dates, columns= ['vitesseB_moy'])
OBS1_1 = pd.DataFrame(vitesse_B_std, index=liste_dates, columns= ['vitesseB_std'])
OBS1_2 = pd.DataFrame(vitesse_B_min, index=liste_dates, columns= ['vitesseB_min'])
OBS1_3 = pd.DataFrame(vitesse_B_max, index=liste_dates, columns= ['vitesseB_max'])
OBS1_4 = pd.DataFrame(vitesse_bias_moy, index=liste_dates, columns= ['vitesse_bias_moy'])
OBS1_5 = pd.DataFrame(vitesse_bias_std, index=liste_dates, columns= ['vitesse_bias_std'])

OBS2_0 = pd.DataFrame(angle_B_moy, index=liste_dates, columns= ['angleB_moy'])
OBS2_1 = pd.DataFrame(angle_B_std, index=liste_dates, columns= ['angleB_std'])
OBS2_2 = pd.DataFrame(angle_B_min, index=liste_dates, columns= ['angleB_min'])
OBS2_3 = pd.DataFrame(angle_B_max, index=liste_dates, columns= ['angleB_max'])

OBS3 = pd.DataFrame(np.array(nb_bouees), index=liste_dates, columns=['nombreB'])
OBS4 = pd.DataFrame(vitesse_B2_moy, index=liste_dates, columns= ['vitesse_B2'])

OBS=pd.concat([OBS1_0, OBS1_1, OBS1_2, OBS1_3, OBS1_4, OBS1_5, OBS2_0, OBS2_1, OBS2_3, OBS3, OBS4], axis=1)

MLDP1_0 = pd.DataFrame(vitesse_MLDP_moy, index=liste_dates, columns= ['vitesseM_moy'])
MLDP1_1 = pd.DataFrame(vitesse_MLDP_std, index=liste_dates, columns= ['vitesseM_std'])
MLDP1_2 = pd.DataFrame(vitesse_MLDP_min, index=liste_dates, columns= ['vitesseM_min'])
MLDP1_3 = pd.DataFrame(vitesse_MLDP_max, index=liste_dates, columns= ['vitesseM_max'])

MLDP2_0 = pd.DataFrame(angle_MLDP_moy, index=liste_dates, columns= ['angleM_moy'])
MLDP2_1 = pd.DataFrame(angle_MLDP_std, index=liste_dates, columns= ['angleM_std'])
MLDP2_2 = pd.DataFrame(angle_MLDP_min, index=liste_dates, columns= ['angleM_min'])
MLDP2_3 = pd.DataFrame(angle_MLDP_max, index=liste_dates, columns= ['angleM_max'])

MLDP3_0 = pd.DataFrame(angle_bias_moy, index=liste_dates, columns= ['angle_bias_moy'])
MLDP3_1 = pd.DataFrame(angle_bias_std, index=liste_dates, columns= ['angle_bias_std'])

MLDP=pd.concat([MLDP1_0, MLDP1_1, MLDP1_2, MLDP1_3, MLDP2_0, MLDP2_1, MLDP2_2, MLDP2_3, MLDP3_0, MLDP3_1], axis=1)

distance_pd_moy=pd.DataFrame(distance_moy, index=liste_dates, columns= ['distance_moy'])
distance_pd_std=pd.DataFrame(distance_std, index=liste_dates, columns= ['distance_std'])
distance_pd_min=pd.DataFrame(distance_min, index=liste_dates, columns= ['distance_min'])
distance_pd_max=pd.DataFrame(distance_max, index=liste_dates, columns= ['distance_max'])

df_global = pd.concat([OBS, MLDP, distance_pd_moy, distance_pd_std, distance_pd_min, distance_pd_max], axis = 1)

with open('analyse_bouee_data_avec_v0_1.pkl', 'wb') as output:
  pickle.dump(df_global, output, -1)

df_global.to_csv('analyse_bouee_data_avec_v0_1.csv')

end = timer()
print 'Fin execution : ', (end - start), ' s'

