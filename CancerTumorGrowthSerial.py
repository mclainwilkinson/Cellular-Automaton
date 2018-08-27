import numpy as np
import random
import math

# grid shape and grid origin
ROWS = 101
COLS = 101	
origin = (COLS//2, ROWS//2)

# constants
k1, k2, k3, k4 = (0.7, 0.2, 0.3, 0.3)
phi = ROWS * COLS
rho = 3.85 

# set generations
generations = 50 

# cell types: N, C, E, D (normal, cancerous, effector, dead)

# initialize grid of normal cells
M = np.empty([ROWS, COLS],dtype=object)
M.fill('N')

# initialize cancer cell(s) in middle of grid
M[COLS//2, ROWS//2] = 'C'
M[COLS//2 + 1, ROWS//2] = 'C'
M[COLS//2 - 1, ROWS//2] = 'C'
M[COLS//2, ROWS//2 - 1] = 'C'
M[COLS//2, ROWS//2 + 1] = 'C'

# functions
def sumCellType(cell_type):
    # number of type of cell in grid
    return sum(sum(M==cell_type))

def mitosisProb(k, n, p):
    # probability of cell mitosis
    return (k * (1 - n / p))

def originDistance(r, c):
    # calculate grid distance from origin
    return math.sqrt((r - origin[0])**2 + (c - origin[1])**2)

def densityDevelopment(ROWS, COLS):
    # calculate density development of tumor
    c = sumCellType('C')
    e = sumCellType('E')
    d = sumCellType('D')
    nPrime = c + e + d
    R = 0
    for i in range(ROWS):
        for j in range(COLS):
	    if M[i, j] == 'C':
                R += originDistance(i, j)
    R = R / nPrime
    return nPrime / R**2

def getQuadrant(r, c, origin):
    # get quadrant of coordinates relative to origin
    if r <= origin[0] and c > origin[1]:
	return 'I'
    elif r <= origin[0] and c <= origin[1]:
	return 'II'
    elif r > origin[0] and c <= origin[1]:
	return 'III'
    else:
	return 'IV'
    
def mitosis(choice, dense, quadrant):
    # model cell division w/ density development
    up = (r-1, c)
    rt = (r, c+1)
    dn = (r+1, c)
    lt = (r, c-1) 
    denseMap = {'I':[up,rt],'II':[lt,up],'III':[dn,lt],'IV':[rt,dn]}
    notDenseMap = {'I':[dn,lt],'II':[rt,dn],'III':[up,rt],'IV':[lt,up]} 
    notNormal = ('E','D')
    if dense:
	if choice < 0.5 and newM[denseMap[quadrant][0]] not in notNormal:
	    newM[denseMap[quadrant][0]] = 'C'
	elif choice >= 0.5 and newM[denseMap[quadrant][1]] not in notNormal:
	    newM[denseMap[quadrant][1]] = 'C'
    else:
	if choice < 0.5 and newM[notDenseMap[quadrant][0]] not in notNormal:
	    newM[notDenseMap[quadrant][0]] = 'C' 
	elif choice >= 0.5 and newM[notDenseMap[quadrant][1]] not in notNormal: 
	    newM[notDenseMap[quadrant][1]] = 'C'

# show starting matrix
print(M)
print('')

# create copy of matrix to apply changes to
newM = np.copy(M)

# model expansion of tumor for g generations
for g in range(generations):
    dense = (densityDevelopment(ROWS, COLS) > rho) 
    for r in range(1, ROWS):
        for c in range(1, COLS):
	    quadrant = getQuadrant(r, c, origin)
            if M[r, c] == 'C':
                if random.random() < mitosisProb(k1, sumCellType('C'), phi):
                    choice = random.random()
		    mitosis(choice, dense, quadrant) 
                elif random.random() < k2:
                    newM[r, c] = 'E'
            elif M[r, c] == 'E':
                if random.random() < k3:
                    newM[r, c] = 'D'
            elif M[r, c] == 'D':
                if random.random() < k4:
                    newM[r, c] = 'N'
    M = np.copy(newM)

# show final matrix and show cell type totals
print(M)
print('Normal Cells: ',sumCellType('N'),'Cancer cells: ',sumCellType('C'),'Effector Cells: ',sumCellType('E'),'Dead Cells: ',sumCellType('D'))
