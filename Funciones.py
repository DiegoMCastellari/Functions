#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.ops import nearest_points

import matplotlib.pyplot as plt

# import warnings
# warnings.filterwarnings('ignore')


# In[3]:


# the function calculta the distance between a point and other points
# it return the ID and distance of a specific number (cantidad) of chosen nearest points 

def punto_cercano(DF1, DF2, cantidad):
# DF1: starting point
# DF2: destiny points to calculate distance
# cantidad: number of nearest points to return 

    DF11= DF1.copy()
    DF22= DF2.copy()

    # transform gdf in geoseries to calculate distnaces:    
    GS1 = DF1['geometry']
    GS2 = DF2['geometry']
    
    # Calculate the distances and save them in an array
    point_count = np.empty((len(GS1),len(GS2))) # Create an empty array
    for i, point in enumerate(GS1):
        DF22[i] = ([point.distance(points) for points in GS2]) # Calculate the distances and save them in an array
    
    # Orden df by distance y keep de fisrts 'cantidad' requered
    DF22.rename(columns = {DF22.columns[-1]: 'cercania'}, inplace=True)
    cantidad_index = cantidad - 1
    DF22_sorted = DF22.sort_values('cercania').iloc[cantidad_index]

    return DF22_sorted[['ID_ZONA', 'cercania']]


# In[97]:


# the function merge small polygons into bigger ones
# takes all polygon with area less than Sup_Limit (or Sup_Limit[0]) and iterates over them. It takes all the neightbours with 
# the same values in the columns in filtro_con, calculates the nearest neighbour by centroid, calculate the area of both polygons 
# and if this sum is less than Sup_Limit (or Sup_Limit[1]), they are merge.
# this repeats until the sum is greater than Sup_Limit
# polygons that already have an area greater than Sup_Limit stay the same
# the function returns a gdf with all the polygons merged into new zones, with the sum of their atributes, and an ID with the 
# same name as the first polygon

def union_poligonos(df, ID_Col, Sup_Col, Sup_Limit, filtro_con=[]):
    # df: gdf proyected[m]
    # ID_Col (str): id column name                ** pasa a ser 'ID_ZONA'
    # Sup_Col (str): area column name             ** pasa a ser 'AREA'
    # Sup_Limit (int o int list): area limit value/s - can be a unique value o a tuple (acept: 3000, [3000], [1000, 3000])
    #                 [area máx that a pol must have to be analyzed, area máx of the generated zone]
    # filtro_con (list of str): filter columns - the list can have any lenght
    
    if (df.crs == {'init' :'epsg:4326'}):
        print('ERRRROORRRR: cambiar proyeccción')
    
    df=df.rename(columns = {ID_Col:'ID_ZONA', Sup_Col: 'AREA'})
    ID_Col_original = ID_Col
    ID_Col = 'ID_ZONA'
    
    if type(Sup_Limit)==int:    # if Sup_Limit is int creates a tuple of values
        Sup_Limit = [Sup_Limit, Sup_Limit]
    elif len(Sup_Limit)==1:     # if Sup_Limit is list with one value creates one with two
        Sup_Limit = [Sup_Limit[0], Sup_Limit[0]]
    else:
        Sup_Limit = [Sup_Limit[0], Sup_Limit[1]]
  
    df['ZONA_NUEVA'] = '' # new field to gather the new zone name for the polygons to merge 
    df = df.sort_values('AREA')
        
    df.reset_index(inplace=True)
    df.rename(columns = {df.columns[0]: 'index_drop'}, inplace=True)
    df.drop(columns={'index_drop'}, inplace=True)
    
    print("INICIO -  Limite: ", Sup_Limit, " - i totales: ", len(df.loc[df.AREA <= Sup_Limit[0]]))
    contador = 0
    for i in list(df.loc[df.AREA <= Sup_Limit[0]].index): # i only take values of the indexes with area less than Sup_Limit
        
        contador = contador + 1
        if (contador%1000 == 0):
            print('Contador= ', contador)
                  
        if (df.loc[i, 'ZONA_NUEVA'] == '') & (df.loc[i, 'AREA']<= Sup_Limit[0]):  
            
            # aplying the filter:
            # generates a df 'DF_PoligonosSinZona' with all polygons not asigned to a zone yet 
            # and which filtro_con columns values are iqual to i
            mask =['(df.ZONA_NUEVA == "")']    # mask to save conditions 
            for c in range(len(filtro_con)):
                mask.append('(df["' + filtro_con[c] + '"] == df.loc[' + str(i) + ', "' + filtro_con[c] + '"])') # se crea una condición para columna del filtro
            mask = ' & '.join(mask) # creates a string with the whole condition
            DF_PoligonosSinZona = df.loc[eval(mask)] # filters the df with the mask

            List_Poligonos_Zonas = [DF_PoligonosSinZona.loc[i, ID_Col]] # polygons selected
            Poligono_Selec = DF_PoligonosSinZona.loc[DF_PoligonosSinZona.ID_ZONA.isin(List_Poligonos_Zonas)] # polygons in the zone

            cantidad_vecinos = 1 # define a this different than 0 to start the 'while', then it´s modify
            
            # start itarating adding polygons until reaching the limit area (Sup_Limit)
            while (Poligono_Selec['AREA'].sum() <= Sup_Limit[1]) & (cantidad_vecinos !=0):

                # DISSOLVE POLYGONS SELECTED (zone)
                Pol_SEL = Poligono_Selec
                Pol_SEL['ID_UNION'] = 'unionselected'
                Pol_SEL = Pol_SEL.dissolve(by='ID_UNION')

                # CREATE A BUFFER TO RUN OVERLAY
                Pol_SEL_buffer = pd.DataFrame(Pol_SEL.buffer(5))
                Pol_SEL_buffer.rename(columns = {Pol_SEL_buffer.columns[0]: 'geometry'}, inplace=True)
                Pol_SEL_buffer = gpd.GeoDataFrame(Pol_SEL_buffer, geometry='geometry')

                # NEIGHBOURS
                Inter_Vecinos = gpd.overlay(Pol_SEL_buffer, DF_PoligonosSinZona, how='intersection') # All polygons: neighbours and zone
                Vecinos = DF_PoligonosSinZona.loc[DF_PoligonosSinZona.ID_ZONA.isin(list(Inter_Vecinos.ID_ZONA))] # All polygons: neighbours and zone
                Vecinos = Vecinos.loc[~Vecinos.ID_ZONA.isin(List_Poligonos_Zonas)] # only neighbours
                cantidad_vecinos = len(Vecinos) # number of neighbours, if 0 then don´t look for neighbour and the 'while' will break

                if (cantidad_vecinos !=0): # IF THERE´S AT LEAST ONE NEIGHTBOUR 

                    # CENTROIDS
                    # calculate the centroid of the dissolved zone
                    points1 = Pol_SEL
                    points1['geometry'] = points1['geometry'].centroid
                    # Calculate the centroid of each neighbour
                    points2 = Vecinos
                    points2['geometry'] = points2['geometry'].centroid

                    # CALCULATE TEH DISTANCE AND KEEP THE NEAREST 
                    df_func_cercania = punto_cercano(points1, points2, 1) # returns the ID of the nearest zone
                    df_cercania = df_func_cercania['ID_ZONA'] # keep the nearest neighbour
                    df_distancia_min = df_func_cercania['cercania']

                    List_Poligonos_Zonas_Eval = List_Poligonos_Zonas.copy() # create a copy of the zone polygons
                    List_Poligonos_Zonas_Eval.append(df_cercania)    # add the nearest polygon to the list
                    Sup_Resultante = DF_PoligonosSinZona.loc[DF_PoligonosSinZona.ID_ZONA.isin(List_Poligonos_Zonas_Eval), 'AREA'].sum() # calculate area

                    if (Sup_Resultante <= Sup_Limit[1]): # If area is less than the limit (Sup_Limit) 
                        List_Poligonos_Zonas.append(df_cercania) # add the neighbour to the list of plygons of the zone
                        Poligono_Selec = DF_PoligonosSinZona.loc[DF_PoligonosSinZona.ID_ZONA.isin(List_Poligonos_Zonas)] # Polygons in the zone
                    else: # if the area is grater than the limit 
                        cantidad_vecinos=0    # modify the variable and 'while' will break
            
            # to all the polygons in the new zone, create a new ID
            df.loc[df.ID_ZONA.isin(List_Poligonos_Zonas), 'ZONA_NUEVA'] = df.loc[i, ID_Col] 
    
    # if poylons area is alerady greater than the limit, complete the new zone field with it´s ID
    df.loc[df.AREA > Sup_Limit[0], 'ZONA_NUEVA'] = df.loc[df.AREA > Sup_Limit[0], ID_Col] 
    
    # create a list with all the not numerial columns to keep
    if len(filtro_con)!=0:
        lista_merge = filtro_con
        lista_merge += ['ID_ZONA']
    else:
        lista_merge = ['ID_ZONA']
    
    # dissolve the polygons by the new zone ID field, adding numerical columns
    dfzonas = df.dissolve(by=['ZONA_NUEVA'], aggfunc='sum').reset_index() 
    # merge with the original df not numerical columns to recover them
    dfzonasfinal = pd.merge(dfzonas, df[lista_merge], left_on='ZONA_NUEVA', right_on='ID_ZONA', how='left') 
          
    # fix the gdf to be returned and able to itarate the function
    dfzonasfinal = gpd.GeoDataFrame(dfzonasfinal, geometry='geometry')
    dfzonasfinal.rename(columns = {dfzonasfinal.columns[0]: 'ID_ZONA', dfzonasfinal.columns[-1]: 'ID_ZONA_dup'}, inplace=True)
    dfzonasfinal.drop(columns={'ID_ZONA_dup'}, inplace=True)
    dfzonasfinal.rename(columns = {'ID_ZONA':ID_Col_original, 'AREA':Sup_Col}, inplace=True)
    dfzonasfinal.crs = df.crs

    print("FIN -  Limite: ", Sup_Limit, " - Polígonos finales: ", len(dfzonas), " / ", len(df))
    return dfzonasfinal


# In[ ]:


# the function create a dataframe with a sample of pixels inside a window with a specific class name
# it´s also posible to pass a distance(d) to add the neighbours data to the central pixel

def muestreo_imagen_vecinos (imagen, xini, xfin, yini, yfin, d, clase):  
    # imagen: image from where sample will be taken
    # xini, xfin, yini, yfin: define image rectangle of the sample
    # clase: class asigned to pixels
    # d: number of pixels around the central pixel ("radio")
    
    img = imagen               
    pixel_list = []            # list where gather pixels and their neighbours data

    for yy in range(yini, yfin,1):                  # itarate over row and columns    
            for xx in range(xini, xfin, 1):
                pixel_fila = []                           # list where will be gather the data of every pixel in the window search (around the central)
                pixel_pos = img.getpixel(xy=(xx,yy))      # index central pixel
                pixel_fila += pixel_pos                   # put data into list

                for d1 in range(-d, (d+1), 1):                   # iterate over rows and columns around central pixel 
                    for d2 in range(-d, (d+1), 1):
                        if (d1 != 0) or (d2 != 0):                # exclude central pixel

                            pixel_around = img.getpixel(xy=(xx+d1,yy+d2))       # index pixels values
                            pixel_fila += pixel_around                          # add to parcial list
                            
                pixel_list.append(pixel_fila)        # add row to list 


    df = pd.DataFrame(pixel_list)       # generate dataframe
    
    posind = ["R", "G", "B"]            # generate a list with column names
    for ind in range(1, int(len(df.columns)/3)):        # iterate over columns to complete their names
        for color in ["R", "G", "B"]:
            posind.append(color + str(ind))
    df.columns = posind                                 # replace columns names
    
    df['clase'] = clase             # create a 'class' field with the class name 
    
    return df


# In[ ]:


# the function, for every pixel in the image, calculate the R, G and B mean and std of all the neighbour pixels

def imagen_info_vecinos (d_img_compl):
    # d_img_compl: number of pixels around the central pixel ("radio")
    
    pixel_list = []   # list where info will be gather

    for yy in range(0, img_compl.size[1],1):                  # iterate over rows and columns         
            for xx in range(0, img_compl.size[0], 1):
                pixel_fila2 = []                           # list where will be gather the data of every pixel in the window search (around the central)
                pixel_pos = img_compl.getpixel(xy=(xx,yy))      # index central pixel
                pixel_fila2 += pixel_pos                   # put data into list

                for d1 in range(-d_img_compl, (d_img_compl+1), 1):                   # iterate over rows and columns around central pixel
                    for d2 in range(-d_img_compl, (d_img_compl+1), 1):
                        if (d1 != 0) or (d2 != 0):                                          # exclude central pixel
                            if (0 <xx+d1< img_compl.size[0]) and (0 <yy+d2< img_compl.size[1]):   # limit image range (don´t get outside in corner and edge)           
                                pixel_around = img_compl.getpixel(xy=(xx+d1,yy+d2))       # index pixels values
                                pixel_fila2 += pixel_around                               # add to parcial list

                a = list(range(0,len(pixel_fila2),3))
                b = list(np.array(a) + 1)
                c = list(np.array(a) + 2)

                # Separate list by color and transform to array 
                pixel_fila2_R = np.array(pixel_fila2)[a]
                pixel_fila2_G = np.array(pixel_fila2)[b]
                pixel_fila2_B = np.array(pixel_fila2)[c]

                listaRGB = []
                # add each color the mean and std to the list 
                listaRGB += (pixel_fila2_R.mean(), pixel_fila2_R.std())        
                listaRGB += (pixel_fila2_G.mean(), pixel_fila2_G.std())   
                listaRGB += (pixel_fila2_B.mean(), pixel_fila2_B.std())  

                pixel_list.append(listaRGB)
                
    return pixel_list

