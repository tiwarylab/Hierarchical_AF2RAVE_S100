#!/usr/bin/env python3

center=(-6.88,-0.805,51.142)
size=(19.755,13.30,23.85)
dels=(size[0]/2,size[1]/2,size[2]/2)
print("HEADER    CORNERS OF BOX")
print("REMARK    CENTER (X Y Z)   "+str(center[0])+" "+str(center[1])+" "+str(center[2]))
print("REMARK    DIMENSIONS (X Y Z) "+str(size[0])+" "+str(size[1])+" "+str(size[2]))
x,y,z=center[0]-dels[0],center[1]-dels[1],center[2]-dels[2]
def fill_space(n,spaces=7):
    s=str(round(n,3))
    if len(s)>spaces:
        return s
    else:
        return " "*(spaces-len(s))+s
def getStr(x,y,z):
    return fill_space(x)+" "+fill_space(y)+" "+fill_space(z)


print("ATOM      1  DUA BOX     1     "+getStr(x,y,z))
x+=2*dels[0]
print("ATOM      2  DUA BOX     1     "+getStr(x,y,z))
y+=2*dels[1]
print("ATOM      3  DUA BOX     1     "+getStr(x,y,z))
x-=2*dels[0]
print("ATOM      4  DUA BOX     1     "+getStr(x,y,z))

y-=2*dels[1]
z+=2*dels[2]

print("ATOM      5  DUA BOX     1     "+getStr(x,y,z))
x+=2*dels[0]
print("ATOM      6  DUA BOX     1     "+getStr(x,y,z))
y+=2*dels[1]
print("ATOM      7  DUA BOX     1     "+getStr(x,y,z))
x-=2*dels[0]
print("ATOM      8  DUA BOX     1     "+getStr(x,y,z))

print("CONECT    1    2    4    5")
print("CONECT    2    1    3    6")
print("CONECT    3    2    4    7")
print("CONECT    4    1    3    8")
print("CONECT    5    1    6    8")
print("CONECT    6    2    5    7")
print("CONECT    7    3    6    8")
print("CONECT    8    4    5    7")
