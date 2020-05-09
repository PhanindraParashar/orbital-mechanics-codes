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
    print(d1)
    seconds = int(b-a)     #seconds
    minutes = seconds/60
    hrs = minutes/60
    days = hrs/24
    
    return days

perigeeMoon = '6-10-2018 3:56:00'
AriesMoon = '27-9-2018 1:54:00'
perigeeMoonAngle = 125.35166  #degree
perigeeMoonAngleR = perigeeMoonAngle*(np.pi/180)

emoon = 0.0549
factor = np.sqrt((1 - emoon)/(1 + emoon))

theta = np.pi/6

E = 2*np.arctan(np.tan(theta/2)*factor)

area = (1/2)*(E - emoon*np.sin(E))

Tarea = np.pi
OrbitPeriod = 27.321661

const = OrbitPeriod/Tarea

timeDelta = const*area          #days


d = DateDelta(6,10,2018,3,56,00,8,10,2018,5,46,00)

error = (timeDelta - d)*24  #Hrs
#print(error)

theta2 = perigeeMoonAngleR + np.pi/6
E2 = 2*np.arctan(np.tan(theta2/2)*factor)
area2 = (1/2)*(E2 - emoon*np.sin(E2))

print(area2)

theta1 = perigeeMoonAngleR
E1 = 2*np.arctan(np.tan(theta1/2)*factor)
area1 = (1/2)*(E1 - emoon*np.sin(E1))

print(area1)

area = area2 - area1

print(area)

timeDelta = const*area

print(timeDelta)
d = DateDelta(27,9,2018,1,54,00,29,9,2018,8,26,00)
error = (timeDelta - d)*24
print(error)



