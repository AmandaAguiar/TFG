#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 16:15:29 2020
@author: amagalaguiar
"""
from astropy.io import fits #Se importan los módulos necesarios
import numpy as np
import matplotlib.pyplot as plt
import os

plt.close('all') #Al ejecutar el programa, borra las figuras anteriores 

#Datos de entrada
toper0_1 = 60000 #Tiempos durante el cual el detector "espera" en ms
toper0_2 = 300000
toper0_3 = 600000
topers0 = np.array([toper0_1,toper0_2,toper0_3]) #Array con los tres toper
nfrsec0 = 200 #Imágenes por cada exposición del detector
nimgobbl0 = 1 #Número de secuencias 
readmode = 'Multi-FOWLER' #Modo de lectura

#Se crean elementos para almacenar los datos recogidos
data = np.zeros((2048, 2048)) #Cuentas recogidas por cada píxel
medias = np.zeros((nfrsec0,32)) #Medias con las que se trabajará
medias_v = np.zeros(32) #Varianzas con las que se trabajará
ts = np.empty(nfrsec0) #Tiempos asociados a cada imagen      
DC = np.zeros(32) #DC de cada canal
RON = np.zeros(32) #RON de cada canal
datmin = 1000 #Valores mín y máx para píxeles buenos en primera lectura     
datmax = 10000                           
good = np.empty((4,8),dtype='object') #Índice de píxeles válidos
file_id0 = 1225372 #Número asociado a primer filename (toper=60000)
fecha = '20170518' #Fecha asociada a las primeras imágenes
file_idf = 1226975 #Número asociado al último filename 
file_id1 = 1226112 #Número del filename a partir del cual cambia la fecha

#Se define el número de píxeles que distan del píxel central en cada canal 
#dependiendo de su orientación
xv = 25 ; yv = 50 #Canales verticales       
xh = 50 ; yh = 25 #Canales horizontales  
                                                       
path = '/media/amagalaguiar/Elements/TFG/gan_lin_ron/' #Directorio

nfig = 1 #Índice para las figuras
i = 0 #Índice del bucle

#Se crean varios txt para recoger los valores de interés 
file = open('RON.txt', 'w')
file.write(u'Toper Canal RON')
file.write(u'\n')

file2 = open('DC_RON.txt', 'w')
file2.write(u'Toper Canal Grupo DC')
file2.write(u'\n')

file3 = open('RON_medio.txt','w')
file3.write(u'Toper RON Error')
file3.write(u'\n')

#Se emplea un bucle para los cálculos y representaciones
while (i+file_id0)<file_idf: #Se toman todos los archivos considerados
    frsec=0
    
    while frsec<nfrsec0: #Se ejecuta esta parte para cada imagen        
        if (file_id0+i)>file_id1:
            fecha = '20170519' #Cambio de fecha en el nombre del archivo
            
        #Se escribe el nombre de cada archivo fit             
        filename = path+'000'+str('%d'%(file_id0+i))+\
        '-' + fecha +'-EMIR-TEST0.fits'
        
        i+=1 #Se incrementa el índice en una unidad 
       
        if not os.path.exists(filename): #Si no existe un archivo con ese 
            continue                     #nombre, que siga pasando por el resto

        frsec = fits.getval(filename,'frsec') #Contador imágenes de exposición 
        
        if frsec==1: #Se recogen los valores de los headers para la primera
                     #imagen de la exposición          
            nfrsec = fits.getval(filename, 'nfrsec') #Número de imágenes por
                                                     #exposición                                                
            toper = fits.getval(filename, 'toper') #Tiempo "de espera" en 
                                                   #misma secuencia            
            nrdil = fits.getval(filename, 'nrdil') #Número de lecturas 
                                                   #productivas            
            nrdil2 = nrdil/2.            
            nrdilm = nrdil - nrdil2 #Se toman las imágenes de la segunda parte
                                    #de cada serie para calcular la DC
                                                                                                                            
        ts[frsec-1] = fits.getval(filename, 'tsutc2') #Se recoge el tiempo    
        data = fits.getdata(filename).astype(float) #Se recogen las cuentas 
                 
        for j in range(8): #Número de canales por cuadrante
            
            if frsec==1: #Se seleccionan los píxeles válidos en la primera
                         #imagen de la exposición              
                dat1 = data[512-xv:512+xv, 64+j*128-yv:64+128*j+yv].flatten()
                npun1 = len(dat1)
                ind1 = [m for m in range(npun1) if datmin<dat1[m]<datmax]
                good[0,j] = ind1
                
                dat2 = data[64+j*128-xh:64+j*128+xh, 512+1024-yh:512+1024+\
                            yh].flatten()
                npun2 = len(dat2)
                ind2 = [m for m in range(npun2) if datmin<dat2[m]<datmax]                
                good[1,j] = ind2
                
                dat3 = data[512+1024-xv:512+1024+xv, 1024+64+j*128-yh:1024+64+\
                            j*128+yh].flatten()
                npun3 = len(dat3)
                ind3 = [m for m in range(npun3) if datmin<dat3[m]<datmax]                
                good[2,j] = ind3
                
                dat4 = data[1024+64+128*j-xh:1024+64+128*j+xh, 512-yh:\
                            512+yh].flatten()
                npun4 = len(dat4)
                ind4 = [m for m in range(npun4) if datmin<dat4[m]<datmax]                
                good[3,j] = ind4
            
            #Se recogen los datos de cada canal y se promedia la señal de los
            #píxeles tomados como válidos
            dat1 = data[512-xv:512+xv, 64+j*128-yv:64+128*j+yv].flatten()
            dat2 = data[64+j*128-xh:64+j*128+xh, 512+1024-yh:512+1024+\
                        yh].flatten()
            dat3 = data[512+1024-xv:512+1024+xv, 1024+64+j*128-yh:1024+64+\
                        j*128+yh].flatten()
            dat4 = data[1024+64+128*j-xh:1024+64+128*j+xh, 512-yh:512\
                        +yh].flatten()
            
            medias[frsec-1, j] = np.mean(dat1[good[0,j]])
            medias[frsec-1, j + 8] = np.mean(dat2[good[1,j]])
            medias[frsec-1, j + 16] = np.mean(dat3[good[2,j]])
            medias[frsec-1, j + 24] = np.mean(dat4[good[3,j]])
            
            #Se halla la varianza de las diferencia entre las dos últimas
            #imágenes para los píxeles válidos
            if frsec==nfrsec-1: 
                
                data2 = data
            
            elif frsec == nfrsec:
                
                dat2_1 = data2[512-xv:512+xv, 64+j*128-yv:64+128*j+\
                               yv].flatten()
                dat2_2 = data2[64+j*128-xh:64+j*128+xh, 512+1024-yh:512+1024\
                               +yh].flatten()
                dat2_3 = data2[512+1024-xv:512+1024+xv, 1024+64+j*128-yh:\
                               1024+64+j*128+yh].flatten()
                dat2_4 = data2[1024+64+128*j-xh:1024+64+128*j+xh, 512-yh:\
                               512+yh].flatten()
                
                medias_v[j] = np.var(dat1[good[0,j]]-dat2_1[good[0,j]])
                medias_v[j + 8] = np.var(dat2[good[1,j]]-dat2_2[good[1,j]])
                medias_v[j + 16] = np.var(dat3[good[2,j]]-dat2_3[good[2,j]])
                medias_v[j + 24] = np.var(dat4[good[3,j]]-dat2_4[good[3,j]])  
               
    ts-=ts[0] #Tiempo asociado a cada imagen
                
    fig=plt.figure(figsize=(8.3,11.7)) #Se define el tamaño de la imagen
    u=np.where(topers0==toper)[0][0] #Índice escogido en función del toper
    n_g = int(nfrsec/nrdil) #Número de grupos
    
    for k in range(32): #Número de canales del detector

        plt.subplot(8,4,k+1)        
        RON[k] = np.sqrt(medias_v[k]/2.) #Cálculo del RON
        
        if k==0 or k==4 or k==8 or k==12 or k==16 or k==20 or k==24 or k==28:
            plt.ylabel(u'ADU', fontsize = 8)                    
        
        if k > 27:     
            plt.xlabel(u'T[s]', fontsize = 8)        
            plt.xticks(np.linspace(int(ts[0]),int(ts[-1]),5),fontsize = 8)
        
        else:         
            plt.tick_params(axis='x',bottom=False, labelbottom=False)
            
        if k<8:
            c = 'tab:pink'
        elif k>=8 and k<=15:
            c = 'tab:cyan'
        elif k>15 and k<=23:
            c = 'darkgoldenrod'
        else:
            c = 'forestgreen'
                
        plt.plot(ts, medias[:,k],'*',color=c)
        plt.title(u'RON = %.2f [ADU]' % (RON[k]), fontsize=8)
        plt.yticks(np.linspace(int(medias[0,k]),int(medias[-1,k]),2),\
                   fontsize = 8,rotation=90)
        plt.xticks(np.linspace(int(ts[0]),int(ts[-1]),5), fontsize = 8)   
        
        file.write('%s %s %s\n' % (toper,k+1,RON[k]))
        
        for w in range(n_g): #Se halla la DC a partir de la segunda mitad de  
                             #cada serie de medidas          
            w1 = nrdil2 + nrdil*w #Índices para seleccionar estos valores
            w2 = w1 + nrdilm
            p = np.polyfit(ts[int(w1):int(w2)],medias[int(w1):int(w2),k],1)       
            DC[k] = p[0]    
            
            file2.write('%s %s %s %s\n' % (toper,k+1,w,DC[k]))          
            
    exptime = fits.getval(filename,'exptime') #Tiempo de exposición     
    
    file3.write('%d %.5f %.5f\n' % (toper,np.mean(RON[:]),np.std(RON[:])))       
            
    plt.suptitle(r'$T_{exp}=%.1f\:s,\:NFRSEC=%d,\:TOPER=%d\:ms,\:NIMGOBBL=%d$'\
                 %  (exptime,nfrsec0,toper,nimgobbl0) +\
                 '$\; RON = %.2f \pm %.3f\:[ADU]$' % (np.mean(RON[:]),  
                 np.std(RON[:])),x=0.5,y=0.99)     
    
    plt.tight_layout(pad=1, w_pad=0.8, h_pad=0.8)
    plt.subplots_adjust(top=0.95)      

    fig.savefig('RON_'+str('%d'%nfig)+'.png',dpi=200)
    plt.close(fig)
    nfig+=1
    
file.close()
file2.close()
file3.close()
