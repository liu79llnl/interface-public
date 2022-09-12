import math

from main.geoms.geoms import getDistance

# Input: stencil = 3x3 array of area fractions
# Output: Young's normal as unit vector
def getYoungsNormal(stencil):
    youngs_alpha = 2
    f_e = 1/(2+youngs_alpha)*(stencil[2][0] + youngs_alpha*stencil[2][1] + stencil[2][2])
    f_w = 1/(2+youngs_alpha)*(stencil[0][0] + youngs_alpha*stencil[0][1] + stencil[0][2])
    f_n = 1/(2+youngs_alpha)*(stencil[0][2] + youngs_alpha*stencil[1][2] + stencil[2][2])
    f_s = 1/(2+youngs_alpha)*(stencil[0][0] + youngs_alpha*stencil[1][0] + stencil[2][0])
    normal = [(f_e - f_w)/2, (f_n - f_s)/2]
    normal_magnitude = getDistance([0,0], normal)
    normal = [normal[0]/normal_magnitude, normal[1]/normal_magnitude]

    return normal

# Input: stencil = 3x3 array of area fractions
# Output: Young's normal as unit vector
# TODO finish this
def getElviraNormal(stencil, norm="l-2"):
    s = [0 for _ in range(6)]
    n = [None for _ in range(24)]
    lines = [None for _ in range(24)]
    l2s = [0 for _ in range(24)]

    for i in range(-1, 2):
        s[0] += stencil[1][1+i] - stencil[0][1+i]
        s[1] += 0.5*(stencil[2][1+i] - stencil[0][1+i])
        s[2] += stencil[2][1+i] - stencil[1][1+i]
        s[3] += stencil[1+i][1] - stencil[1+i][0]
        s[4] += 0.5*(stencil[1+i][2] - stencil[1+i][0])
        s[5] += stencil[1+i][2] - stencil[1+i][1]

    n[0] = [s[0]/math.sqrt(1+s[0]**2), -1/math.sqrt(1+s[0]**2)]
    n[1] = [s[1]/math.sqrt(1+s[1]**2), -1/math.sqrt(1+s[1]**2)]
    n[2] = [s[2]/math.sqrt(1+s[2]**2), -1/math.sqrt(1+s[2]**2)]
    n[3] = [-1/math.sqrt(1+s[3]**2), s[3]/math.sqrt(1+s[3]**2)]
    n[4] = [-1/math.sqrt(1+s[4]**2), s[4]/math.sqrt(1+s[4]**2)]
    n[5] = [-1/math.sqrt(1+s[5]**2), s[5]/math.sqrt(1+s[5]**2)]

    for i in range(6):
        n[i+6] = [-n[i][0], -n[i][1]]
        n[i+12] = [-n[i][0], n[i][1]]
        n[i+18] = [n[i][0], -n[i][1]]
    
    curmin = 0
    if norm == "l-2": #use l2 norm
        for option in range(24):
            #print(option)
            l1, l2 = getLinearFacetFromNormal(opolys[x][y], areas[x][y], n[option], threshold)
            for i in range(-1, 2):
                for j in range(-1, 2):
                    if not(x < 0 or y < 0 or x >= len(opolys) or y >= len(opolys[0])) and l2s[option] <= l2s[curmin]:
                        #this is a valid neighbor square
                        #print(getPolyLineArea(opolys[x+i][y+j], l1, l2)/getArea(opolys[x+i][y+j]))
                        l2s[option] += (areas[x+i][y+j] - getPolyLineArea(opolys[x+i][y+j], l1, l2)/getArea(opolys[x+i][y+j]))**2

            if l2s[option] < l2s[curmin]:
                curmin = option
                
    elif norm == "l-inf": #use l-inf norm
        for option in range(12):
            l1, l2 = getLinearFacetFromNormal(opolys[x][y], areas[x][y], n[option], threshold)
            for i in range(-1, 2):
                for j in range(-1, 2):
                    if not(x < 0 or y < 0 or x >= len(opolys) or y >= len(opolys[0])):
                        linf = (areas[x+i][y+j] - getPolyLineArea(opolys[x+i][y+j], l1, l2)/getArea(opolys[x+i][y+j]))
                        if linf > l2s[curmin]:
                            break
                        else:
                            l2s[option] = max(l2s[option], linf)

            if l2s[option] < l2s[curmin]:
                curmin = option

    else:
        raise ValueError(f"Unknown norm {norm} in getElviraNormal")

    return n[curmin]