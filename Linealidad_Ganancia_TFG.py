#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 20:39:33 2020
@author: amagalaguiar
"""
from astropy.io import fits #Se importan los módulos necesarios
import numpy as np
import matplotlib.pyplot as plt
import os

plt.close('all') #Al ejecutar el programa, borra las figuras anteriores 

#Datos de entrada de la rampa
toper0 = 668 #Tiempo durante el cual el detector "espera" en ms
nfrsec0 = 30 #Imágenes por cada exposición del detector
ptos0 = 6 #Número inicial de puntos para los cálculos
n0 = 1 #Número de puntos no empleados para el cálculo de la pendiente
nimgobbl0 = 1

#Se escribe el número asociado al primer archivo de la primera rampa y al 
#último de la segunda
primer_filename1 = 1225310 
ultimo_filename2 = primer_filename1 + 29 + 1 + 30

#Se crean elementos para almacenar los datos recogidos
data1 = np.zeros((2048, 2048)) #Cuentas recogidas por cada píxel en rampa 1
data2 = np.zeros((2048, 2048)) #Cuentas recogidas por cada píxel en rampa 2
medias = np.zeros((nfrsec0,32)) #Medias con las que se trabajará  
medias_v = np.zeros((nfrsec0,32)) #Varianzas con las que se trabajará
ts = np.zeros((30)) #Tiempos asociados a cada imagen
                                                                      
tol1 = 0.01 #1% de tolerancia
tol2 = 0.02 #2% de tolerancia

#Se define el número de píxeles que distan del píxel central en cada canal 
#dependiendo de su orientación
xv = 25 ; yv = 50 #Canales verticales       
xh = 50 ; yh = 25 #Canales horizontales  

#Se crean arrays para determinar en qué frsec deja de cumplirse la condición
#y para almacenar las cuentas asociadas a la saturación en cada canal
frsec_utilizados1 = np.zeros((32)) #Asociados a tol1
cuentas1 = np.zeros((32))
frsec_utilizados2 = np.zeros((32)) #Asociados a tol2
cuentas2 = np.zeros((32))

g = np.zeros((32)) #Array para almacenar los factores de conversión   
error_g = np.zeros((32)) #Array para almacenar error de g                                                 

path = '/media/amagalaguiar/Elements/TFG/gan_lin_ron/' #Directorio

#Para saber cuántas veces se debe ejecutar el bucle, se crea una lista con
#los nombres de los archivos
files = os.listdir(path)

#Se emplean varios bucles para llevar a cabo los procesos necesarios
for i in range(len(files)): #Se establece un rango muy superior al necesario 
                            #para asegurar que pase por todos los archivos
                            #En GTC se almacenan por orden y podría haber
                            #haber archivos de otros instrumentos
   #Se escribe el nombre de cada archivo fit asociado a cada rampa 
   filename1 = path+'000'+str('%d'%(primer_filename1+i))+\
   '-20170518-EMIR-TEST0_raw.fits' #Rampa 1
   filename2 = path+'000'+str('%d'%(primer_filename1+nfrsec0+1+i))+\
   '-20170518-EMIR-TEST0_raw.fits' #Rampa2
  
   if not os.path.exists(filename1) or not os.path.exists(filename2): 
       continue #Si no existe el archivo con el nombre indicado, que se siga
                #ejecutando
   if (primer_filename1+nfrsec0+1+i)>ultimo_filename2:
       break #Se detiene cuando se terminan de tomar los datos de las rampas

   #Se recogen datos de interés de los archivos para cada rampa 
   frsec1 = fits.getval(filename1,'frsec') #Contador imágenes de una exposición
   frsec2 = fits.getval(filename2,'frsec') 
   
   if frsec1==1 and frsec2==2:
       nfrsec1 = fits.getval(filename1, 'nfrsec') #Imágenes por exposición  
       toper1 = fits.getval(filename1, 'toper') #Tiempo "de espera"
       nfrsec2 = fits.getval(filename2, 'nfrsec') 
       toper2 = fits.getval(filename2, 'toper')
       
   if nfrsec1 != nfrsec0 or toper1 != toper0:
       continue #Si el archivo no tiene las características indicadas que siga
                #ejecutándose
   if nfrsec2 != nfrsec0 or toper2 != toper0:
       continue 
   
   #Se recogen las cuentas de cada píxel para cada rampa      
   data1 = fits.getdata(filename1).astype(float)
   data2 = fits.getdata(filename2).astype(float)
   #Se recoge el tiempo
   ts[frsec1-1] = fits.getval(filename1, 'tsutc2')  
                                                    
   for j in range(8): #Número de canales por cuadrante del detector
       #Se hace el promedio de las dos rampas
       #Se distinguen las operaciones necesarias en cada cuadrante
       cuadrante1 = np.mean(data1[512-xv:512+xv, 64+j*128-yv:64+128*j+yv]+\
                            data2[512-xv:512+xv, 64+j*128-yv:64+128*j+yv])
       cuadrante2 = np.mean(data1[64+j*128-xh:64+j*128+xh, 512+1024-yh:\
                                  512+1024+yh]+data2[64+j*128-xh:64+j*128+xh,\
                                             512+1024-yh:512+1024+yh]) 
       cuadrante3 = np.mean(data1[512+1024-xv:512+1024+xv, 1024+64+j*128-yh:\
                                  1024+64+j*128+yh]+data2[512+1024-xv:\
                                  512+1024+xv, 1024+64+j*128-yh:1024+64+j*128+\
                                  yh])
       cuadrante4 = np.mean(data1[1024+64+128*j-xh:1024+64+128*j+xh, 512-yh:\
                                  512+yh]+data2[1024+64+128*j-xh:1024+64+128*j\
                                        +xh, 512-yh:512+yh])
                            
       medias[frsec1-1, j]+=  cuadrante1/2. 
       medias[frsec1-1, j + 8]+= cuadrante2/2.
       medias[frsec1-1, j + 16]+= cuadrante3/2.
       medias[frsec1-1, j + 24]+= cuadrante4/2. 
       
       #Se halla la varianza asociada a la resta de las dos rampas
       cuadrante1_g = np.var(data1[512-xv:512+xv, 64+j*128-yv:64+128*j+yv]-\
                            data2[512-xv:512+xv, 64+j*128-yv:64+128*j+yv])
       cuadrante2_g = np.var(data1[64+j*128-xh:64+j*128+xh, 512+1024-yh:\
                                   512+1024+yh]-data2[64+j*128-xh:64+j*128+xh,\
                                              512+1024-yh:512+1024+yh]) 
       cuadrante3_g = np.var(data1[512+1024-xv:512+1024+xv, 1024+64+j*128-yh:\
                                   1024+64+j*128+yh]-data2[512+1024-xv:\
                                   512+1024+xv,1024+64+j*128-yh:1024+64+\
                                   j*128+yh])
       cuadrante4_g = np.var(data1[1024+64+128*j-xh:1024+64+128*j+xh, 512-yh:\
                                   512+yh]-data2[1024+64+128*j-xh:\
                                         1024+64+128*j+xh, 512-yh:512+yh])

       medias_v[frsec1-1, j]+=  cuadrante1_g/2.
       medias_v[frsec1-1, j + 8]+= cuadrante2_g/2.
       medias_v[frsec1-1, j + 16]+= cuadrante3_g/2.
       medias_v[frsec1-1, j + 24]+= cuadrante4_g/2. 

exptime = fits.getval(path +'/000'+str('%d'%(primer_filename1))+ \
        '-20170518-EMIR-TEST0_raw.fits','exptime') #Tiempo de exposición 

#Se resta tiempo inicial para obtener el tiempo transcurrido desde el inicio           
ts[:]-=ts[0]    

#LINEALIDAD      

#Se define el tamaño de la imagen y el bucle para los cálculos y representar
plt.figure(figsize=(8.3,11.7))

for k in range(32): #Número de canales del detector     
    
    tol1_bool=True
    plt.subplot(8,4,k+1)
    for n in range(nfrsec0-ptos0): #Para ir añadiendo puntos       
        #Ajuste lineal empezando por los 5 primeros puntos
        #Posteriormente, se van añadiendo
        p = np.polyfit(np.arange(n0,ptos0+n),medias[n0:ptos0+n,k],1)
        #Se recoge el valor predicho para el siguiente punto
        punto_recta = np.polyval(p,ptos0+n)
        
        #Se comprueba si cumplen la condición impuesta
        if tol1_bool==True and (abs(punto_recta-medias[ptos0+n,k]) > \
                                tol1*medias[ptos0+n,k]): 
            plt.plot(np.arange(ptos0+n),np.polyval(p,np.arange(ptos0+n)),'-',\
                     linewidth=7, color='silver') #Se representa punto a punto
            
            frsec_utilizados1[k] = n+ptos0 
            cuentas1[k]=(medias[ptos0+n,k]+medias[ptos0+n-1,k])/2.
            
            tol1_bool=False #Cuando deja de cumplirse la condición, deja de
                            #ejecutarse esta parte del bucle
        if abs(punto_recta-medias[ptos0+n,k]) > tol2*medias[ptos0+n,k]:
            break #Deja de ejecutarse el bucle cuando deja de cumplirse la 
                  #condición
    #Se representa la recta cuya pendiente se obtiene con todos los puntos 
    #que se dan por válidos
    frsec_utilizados2[k]= n+ptos0
    cuentas2[k] = (medias[ptos0+n,k]+medias[ptos0+n-1,k])/2.

    plt.plot(np.arange(ptos0+n),np.polyval(p,np.arange(ptos0+n)),'-',\
             linewidth=3,color='black')   
    plt.title('S1 = %.2f [ADU]' % (cuentas1[k])+'\n'+'S2 = %.2f [ADU]'%\
              (cuentas2[k]),fontsize=8, loc='center')   
    plt.xticks([])        
    plt.yticks(np.linspace(int(medias[0,k]),int(medias[-1,k]),2),rotation =90,\
               fontsize = 8)
    
    if k==0 or k==4 or k==8 or k==12 or k==16 or k==20 or k==24 or k==28:           
        plt.ylabel(u'ADU', fontsize = 8) 
        plt.ticklabel_format(axis="y") 
                            
    if k > 27:         
        plt.xlabel(u'T [s]', fontsize = 8)
        plt.xticks(np.linspace(int(ts[0]),int(ts[-1]),5),fontsize = 8) 
    if k<8:
        c = 'tab:pink'
    elif k>=8 and k<=15:
        c = 'tab:cyan'
    elif k>15 and k<=23:
        c = 'darkgoldenrod'
    else:
        c = 'forestgreen'   
        
    plt.plot(np.arange(nfrsec0),medias[:,k],'*', markersize=2, color=c)    
        

plt.tight_layout(pad=1, w_pad=0.8, h_pad=0.8)

plt.subplots_adjust(top=0.9)    
    
  
plt.suptitle(r'$T_{exp}=%.1f\:s,\:NFRSEC=%d,\:TOPER=%d\:ms,\:NIMGOBBL=%d$'
             %  (exptime,nfrsec0,toper0,nimgobbl0)+ 
             '\n'+'$S_{tol1}=%.2f \pm %.2f\:[ADU],\:$' \
             % (np.mean(cuentas1),np.std(cuentas1))+
             '$S_{tol2}=%.2f \pm %.2f\:[ADU]$' %\
             (np.mean(cuentas2),np.std(cuentas2)))
    
plt.savefig('Linealidad_2tolerancias',dpi=200)  

#GANANCIA

#Se define el tamaño de la imagen
plt.figure(figsize=(8.3,11.7),dpi=200)

#Se establece un límite en el número de cuentas para el cálculo de la g para
#asegurar que, para todos los canales, se opera en el rango lineal
lim_cuentas = 25000.

for k in range(32): #Número de canales del detector 
    plt.subplot(8,4,k+1)  
    #Se halla el índice donde se supera el número de cuentas para cada canal
    lim = np.argmax(medias[medias[:,k]<lim_cuentas,k])+1
    
    #Se hace un ajuste lineal de la varianza frente a la media hasta el límite
    p,r = np.polyfit(medias[n0:lim,k],medias_v[n0:lim,k],1,cov=True) 
    
    #Se recoge la inversa de la pendiente, correspondiente al g
    g[k] = 1./p[0]
    e_pendiente = np.sqrt((np.diag(r)))
    error_g[k] = e_pendiente[0]/p[0]**2
    
    plt.title('g = %.2f' % (g[k])+ str('$\pm$') + '%.2f [e-/ADU]' % \
              (error_g[k]), fontsize=8, loc='center')
    plt.xticks(np.linspace(int(medias[n0,k]),int(medias[lim,k]),3),\
               fontsize = 8)       
    plt.yticks(np.linspace(int(medias_v[n0,k]),int(medias_v[lim,k]),3),\
               rotation = 90, fontsize = 8)
        
    if k==0 or k==4 or k==8 or k==12 or k==16 or k==20 or k==24 or k==28:
        plt.ylabel(u'Varianza [ADU]', fontsize = 8) 
        plt.ticklabel_format(axis="y") 
                            
    if k > 27:            
        plt.xlabel(u'Medias [ADU]', fontsize = 8) 
        
    if k<8:
        c = 'tab:pink'
    elif k>=8 and k<=15:
        c = 'tab:cyan'
    elif k>15 and k<=23:
        c = 'darkgoldenrod'
    else:
        c = 'forestgreen'   
        
    plt.plot(medias[n0:lim,k],medias_v[n0:lim,k],'*',color=c)    
        
plt.tight_layout(pad=1, w_pad=0.8, h_pad=0.8)
plt.subplots_adjust(top=0.93)    

#Se hallan los promedios
media_total = np.average(g, weights=1./error_g**2)
desvest_total = np.sqrt(1./(np.sum(1./error_g)**2))    
  
plt.suptitle(r'$T_{exp}=%.1f\:s,\:NFRSEC=%d,\:TOPER=%d\:ms,\:NIMGOBBL=%d,$'
             %  (exptime,nfrsec0,toper0,nimgobbl0) + \
             '$\; g = %.3f$' % (media_total) + str('$\pm$') +\
             '$%.3f\:[e-/ADU]$'%(desvest_total))
    
plt.savefig('Ganancia',dpi=200)

