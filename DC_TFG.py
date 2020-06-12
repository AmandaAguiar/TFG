#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 20:00:07 2020
@author: amagalaguiar
Este es un programa para calcular la DC tanto en modo rampa como en CDS 
"""
from astropy.io import fits #Se importan los módulos necesarios
import numpy as np
import matplotlib.pyplot as plt
import os

plt.close('all') #Al ejecutar el programa, borra las figuras anteriores 

fig=plt.figure(figsize=(8.3,11.7)) #Se define el tamaño de la imagen

#Se definen parámetros para hacer el programa más sencillo de modificar
#Ahora mismo se distinguen dos casos de modo que se debe comentar con 
#comillas el que no se va a utilizar
'''
#Rampa
toper0 = 750 #Tiempo durante el cual el detector "espera" en ms
nexp0 = 10 #Número de exposiciones del detector
nfrsec0 = 28 #Imágenes por cada exposición del detector
nimgobbl0 = 5 #Número de secuencias de imágenes 
primer_filename = 2345204 #Parte del nombre del archivo que varía para cada fit
readmode = 'RAMP' #Modo de lectura
'''

#CDS
toper0 = 443 #Tiempo durante el cual el detector "espera" en ms
nexp0 = 10 #Número de exposiciones del detector
nfrsec0 = 2 #Imágenes por cada exposición del detector
nimgobbl0 = 5 #Número de secuencias de imágenes                                           
primer_filename = 2344704 #Parte del nombre del archivo que varía para cada fit
readmode = 'CDS' #Modo de lectura


#Se crean elementos para almacenar los datos recogidos
data = np.zeros((2048, 2048)) #Cuentas recogidas por cada píxel
medias = np.zeros((nfrsec0,nexp0,32)) #Medias con las que se trabajará                 
ts = np.empty(nfrsec0*nexp0) #Tiempos asociados a cada imagen                                                                  
DC = np.zeros((nexp0,32)) #DC de cada exposición en cada canal
medias_canal = np.zeros((32)) #DC de cada canal
desvest_canal = np.zeros((32)) #Error en DC de cada canal

#Se define el número de píxeles que distan del píxel central en cada canal 
#dependiendo de su orientación
xv = 25 ; yv = 50 #Canales verticales       
xh = 50 ; yh = 25 #Canales horizontales  
                                                       
path = '/media/amagalaguiar/Elements/TFG/calib_31oct/' #Directorio de archivos

files = os.listdir(path) #Lista con los nombres de los archivos

#Se emplean varios bucles para llevar a cabo los procesos necesarios
for i in range(len(files)): #Se establece un rango muy superior al necesario 
                            #para asegurar que pase por todos los archivos
                            #En GTC se almacenan por orden y podría haber
                            #haber archivos de otros instrumentos
   
   #Se escribe el nombre de cada archivo fit 
   filename = path+'000'+str('%d'%(primer_filename+i))+\
   '-20191031-EMIR-TEST0.fits'
     
   if not os.path.exists(filename): #Si no existe un archivo con ese nombre
       continue                     #que siga pasando por los demás

   #Se recogen algunos parámetros que caracterizan a cada archivo
   frsec = fits.getval(filename,'frsec') #Contador imágenes de una exposición
   
   if frsec==1:
       nfrsec = fits.getval(filename, 'nfrsec') #Imágenes por exposición                                                     
       exp = fits.getval(filename, 'exp') #Contador de exposiciones
       nexp = fits.getval(filename, 'nexp') #Número de exposiciones
       imgobbl = fits.getval(filename, 'imgobbl') #Contador de secuencias
       nimgobbl = fits.getval(filename, 'nimgobbl') #Número de secuencias 
       toper = fits.getval(filename, 'toper') #Tiempo "de espera"   

   if nfrsec != nfrsec0 or toper != toper0:
       continue #Si el archivo no tiene las características indicadas, que siga
                #ejecutándose  
                
   if exp > nexp0 or imgobbl > nimgobbl0: 
       break #Si el archivo se asocia a un contador de exposición o secuencia 
             #mayor al máximo establecido, esta parte deja de ejecutarse       
 
   data = fits.getdata(filename) #Se recogen las cuentas de cada fit 
   
   if imgobbl == 1: #Se toma el tiempo en la primera secuencia, pues se quiere 
                    #tomar el tiempo asociado a cada imagen de las exposiciones
                    #totales              
       m = (exp-1)*nfrsec + frsec #Índice para almacenar tiempos                           
       ts[m-1] = fits.getval(filename, 'tsutc2') #Tiempo tras leer el detector
       
       if exp == 1 and frsec == 1: 
           ts0 = ts[0]  #Tiempo asociado a la primera imagen de la primera         
                        #exposición de la primera secuencia 
                        
   for j in range(8): #Número de canales por cuadrante 
       
       #Se distinguen las operaciones necesarias en cada cuadrante
       cuadrante1 = np.mean(data[512-xv:512+xv, 64+j*128-yv:64+128*j+yv])
       cuadrante2 = np.mean(data[64+j*128-xh:64+j*128+xh, 512+1024-yh:\
                                 512+1024+yh]) 
       cuadrante3 = np.mean(data[512+1024-xv:512+1024+xv, 1024+64+j*128-yh:\
                                 1024+64+j*128+yh])
       cuadrante4 = np.mean(data[1024+64+128*j-xh:1024+64+128*j+xh, \
                                 512-yh:512+yh]) 
       
       #Se hace un sumatorio para todas las secuencias
       medias[frsec-1, exp-1, j]+=  cuadrante1
       medias[frsec-1, exp-1, j + 8]+= cuadrante2
       medias[frsec-1, exp-1, j + 16]+= cuadrante3
       medias[frsec-1, exp-1, j + 24]+= cuadrante4 
       
#Se divide entre el número de secuencias para obtener el valor medio asociado        
medias/= nimgobbl0

#Se hace un cambio en estos valores obtenidos para visualizar más claramente
#los resultados
medias[:,:,:]+=1000.-medias[0,:,:] 

#Se resta tiempo inicial para obtener el tiempo transcurrido desde el inicio
ts-=ts0 
        
exptime = fits.getval(path +'/000'+str('%d'%(primer_filename))+ \
        '-20191031-EMIR-TEST0.fits','exptime') #Tiempo de exposición 
                                               
#Se define el número de puntos que se van a descartar para los cálculos
if readmode == 'CDS':
    eliminados = 0
    
elif readmode == 'RAMP': 
    eliminados = 5

#Por último, se emplea un bucle para hacer las representaciones y los cálculos   
for k in range(32): #Número de canales del detector 
    
    for r in range(nexp0): #Número de exposiciones
        
        #Ajuste lineal hecho eliminando cierto número de puntos
        p = np.polyfit(ts[nfrsec0*r+eliminados:nfrsec0*(r+1)],\
                          medias[eliminados:,r,k],1)                                                                              
        
        DC[r,k] = p[0] #Se recoge la pendiente, que se corresponde con la DC
        
        #Se hace una gráfica para cada uno de los 32 canales
        ax = fig.add_subplot(8,4,k+1)       
        plt.yticks(np.linspace(int(medias[0,0,k]),int(medias[-1,0,k]),2),\
                   rotation = 90, fontsize = 8)
        
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
                    
        plt.plot(ts[nfrsec0*r:nfrsec0*(r+1)], medias[:,r,k], '*', color= c, \
                    markersize = 3) 
        
    medias_canal[k] = np.mean(DC[:,k]) #Medias para cada canal
    desvest_canal[k]=np.std(DC[:,k]) #Desviación estándar
    plt.title('DC = %.2f' % (medias_canal[k]) + str('$\pm$') + '%.3f [ADU]' % (desvest_canal[k]), loc='center',fontsize = 8)
                        
plt.tight_layout(pad=1, w_pad=0.8, h_pad=0.8)
plt.subplots_adjust(top=0.95) 

#Se hallan los promedios 
media_total = np.average(medias_canal,weights=1./desvest_canal**2)
desvest_total = np.sqrt(1./(np.sum(1./(desvest_canal**2))))

plt.suptitle(r'$ %s:$'%(readmode)+ 
             r'$T_{exp}=%.1f\:s,\:NFRSEC=%d,\:TOPER=%d\:ms,\:NIMGOBBL=%d,$'
             %  (exptime,nfrsec0,toper0,nimgobbl0)  + \
             '$\; DC = %.2f$' % (media_total) + str('$\pm$') +\
             '$%.3f\:[ADU]$'%(desvest_total),x=0.5,y=0.99)

plt.savefig('DC_%s_'%(readmode)+str('%.1f'%exptime)+'.png',dpi=200)