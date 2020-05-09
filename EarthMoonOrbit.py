from datetime import datetime
import time
import numpy as np

fmt = '%d-%m-%Y %H:%M:%S'

def DateDelta(Date1,Month1,Year1,Hrs1,Min1,Sec1,Date2,Month2,Year2,Hrs2,Min2,Sec2):
    fmt = '%d-%m-%Y %H:%M:%S'
    d1 = str(Date1)+"-"+str(Month1)+"-"+str(Year1)+" "+str(Hrs1)+":"+str(Min1)+":"+str(Sec1)
    d2 = str(Date2)+"-"+str(Month2)+"-"+str(Year2)+" "+str(Hrs2)+":"+str(Min2)+":"+str(Sec2)

    date1 = datetime.strptime(d1,fmt)
    date2 = datetime.strptime(d2,fmt)
    
    a = time.mktime(date1.timetuple())
    b = time.mktime(date2.timetuple())
    
    seconds = int(b-a)     #seconds
    minutes = seconds/60
    hrs = minutes/60
    days = hrs/24
    
    return days

def Delta_D(Date1,Date2):
    d1 = str(Date1)
    d2 = str(Date2)

    fmt = '%d-%m-%Y %H:%M:%S'
    
    date1 = datetime.strptime(d1,fmt)
    date2 = datetime.strptime(d2,fmt)
    
    a = time.mktime(date1.timetuple())
    b = time.mktime(date2.timetuple())
    
    seconds = int(b-a)     #seconds
    minutes = seconds/60
    hrs = minutes/60
    days = hrs/24
    
    return days

def Solve_M(val):
    E = 2*val
    i = 1
    while i < 10:

        der = .5*(1 - e*np.cos(E))
        Fun = .5*(E - e*np.sin(E)) - val

        h = -1*Fun/der
        E = E + h

        i = i + 1
    return E    

def Angle_Dates_Perigee_degree(Date):
    PerigeeMoon = '2-1-2018 11:05:00'
    days = Delta_D(PerigeeEarth,Date)
    AreaSwept = Area_Time_Ratio*days
    while AreaSwept > np.pi:
        AreaSwept = AreaSwept - np.pi

    E = Solve_M(AreaSwept)
    theta = 2*np.arctan(np.tan(E/2)/factor)
    thetad = theta*(180/np.pi)

    if thetad < 0:
        thetad = 360 + thetad

    return thetad

def Degree_Minutes(angle):
    degree = np.floor(angle)
    rashi = np.floor(degree/30) + 1
    degree_rashi = degree - 30*(rashi-1)
    minutes = np.floor((angle - degree)*60)
    seconds = np.floor(((angle - degree)*60 - minutes)*60)
    D_M = str(rashi)+ "  " + str(degree_rashi) + "  " + str(minutes) + "  " + str(seconds)
    return D_M
e = 0.0167

factor = np.sqrt((1-e)/(1+e))

PerigeeEarth = '3-1-2018 11:05:00'
ApogeeEarth = '6-7-2018 22:17:00'

VernalEquinox = '20-3-2018 21:45:00' 
AutumnalEquinox = '23-9-2018 7:24:00'

SummerSolstice = '21-6-2018 15:37:00'
WinterSolstice = '22-12-2018 3:53:00'

Sankranti = '14-1-2018 13:47:00'
Aries = '14-4-2018 8:14:00'

SiderialDay = 0.99726956    # Solar Days (24 hr)
SiderialYear = 365.2563629  # Solar Days (24 hr)

AreaTotal = np.pi 

Area_Time_Ratio = AreaTotal/SiderialYear
Time_Area_Ratio = SiderialYear/AreaTotal

TestDate = '12-10-2018 14:25:00'

ImpDates = [Sankranti,Aries,VernalEquinox,AutumnalEquinox,SummerSolstice,WinterSolstice,ApogeeEarth,TestDate]

for i in ImpDates:
    A = Angle_Dates_Perigee_degree(i)
    D = Delta_D(PerigeeEarth,i)
    
    Angle_Aries = A - Angle_Dates_Perigee_degree(Aries)
    if Angle_Aries < 0:
        Angle_Aries = 360 + Angle_Aries
        
    k = Degree_Minutes(Angle_Aries)

    Nakshtra = np.floor(Angle_Aries*(27/360)) + 1
    #print(A,"     ",Angle_Aries,"     " ,k,"   ",i,"    ",D," Days    ",Nakshtra," Nakshtra")
    #print(A,"     ",i)
    print(A,"     ",Angle_Aries,"     " ,k,"   ",i)

#print(Degree_Minutes(174.86821901602067))s

